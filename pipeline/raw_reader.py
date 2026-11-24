"""
Fast reading of raw voltage from a DigitalRF archive on spinning disks.

DigitalRF stores a channel as one HDF5 file per second.  In the Millstone Hill
recordings each file holds its 1,000,000 complex int16 samples as a single
uncompressed, contiguous block at the same byte offset, and ext4 lays one
channel's consecutive files next to each other on the disk.  Reading a pulse
through DigitalRFReader.read_vector_1d costs several small reads (the HDF5
superblock, object headers, the sample index, then the samples), each a disk
seek, so a pipeline that reads pulse by pulse is limited by seeks, not by the
disk's transfer rate.

RawReader is a drop-in for read_vector_1d.  It learns the block's offset once
per channel, checks every file against it by its size, and reads a whole second
with one pread straight into a numpy buffer, keeping the last few seconds so
that consecutive pulses cost nothing.  With whole_seconds=False it instead
reads exactly the requested samples, still with one read and no HDF5, which is
better for pulses scattered far apart.  A file that does not match the layout
(a partial second at the start or end of a recording, a gap) is read through
digital_rf, so the result is always what DigitalRFReader would return.

    r = RawReader("/path/to/recording/rf_data/")
    z = r.read_vector_1d(key, 10000, "zenith-l")      # complex64, as digital_rf
"""

import collections
import datetime
import os

import h5py
import numpy as n


class RawReader:
    def __init__(self, rf_dir, cache_seconds=3, whole_seconds=True, prefetch=0):
        """
        rf_dir         the rf_data directory, as for DigitalRFReader
        cache_seconds  seconds kept in memory per channel (4 MB each)
        whole_seconds  read a whole second when a sample in it is asked for;
                       False reads only the requested samples
        prefetch       ask the kernel to start reading this many following
                       seconds in the background (posix_fadvise WILLNEED)
        """
        self.rf_dir = rf_dir
        self.cache_seconds = cache_seconds
        self.whole_seconds = whole_seconds
        self.prefetch = prefetch
        self._ch = {}
        self._cache = {}
        self._drf = None
        self.n_fast = 0          # files read directly
        self.n_fallback = 0      # requests handed to digital_rf

    # ---- layout -------------------------------------------------------------
    def _channel(self, ch):
        if ch in self._ch:
            return self._ch[ch]
        with h5py.File(os.path.join(self.rf_dir, ch, "drf_properties.h5"), "r") as h:
            a = h.attrs
            sr_num, sr_den = int(a["sample_rate_numerator"]), int(a["sample_rate_denominator"])
            fc_ms = int(a["file_cadence_millisecs"])
            sub_s = int(a["subdir_cadence_secs"])
            ok = sr_den == 1 and (sr_num * fc_ms) % 1000 == 0 and int(a["num_subchannels"]) == 1 \
                and int(a["is_complex"]) == 1
        c = dict(spf=sr_num * fc_ms // 1000, fc_ms=fc_ms, sub_s=sub_s, ok=ok, offset=None, size=None, dtype=None)
        self._ch[ch] = c
        return c

    def _path(self, ch, c, fi):
        t_ms = fi * c["fc_ms"]
        s = t_ms // 1000
        sub = datetime.datetime.fromtimestamp(s - s % c["sub_s"], datetime.timezone.utc)
        return os.path.join(self.rf_dir, ch, sub.strftime("%Y-%m-%dT%H-%M-%S"),
                            "rf@%d.%03d.h5" % (s, t_ms % 1000))

    def _learn(self, c, p):
        """The block's offset, dtype and the file size, from one complete file."""
        with h5py.File(p, "r") as h:
            d = h["rf_data"]
            idx = h["rf_data_index"][()]
            off = d.id.get_offset()
            if off is None or d.shape != (c["spf"], 1) or d.chunks is not None or d.compression is not None \
                    or idx.shape != (1, 2) or d.dtype != n.dtype([("r", "<i2"), ("i", "<i2")]):
                return False
            c.update(offset=off, dtype=d.dtype, size=os.path.getsize(p))
        return True

    # ---- reading ------------------------------------------------------------
    def _fast_fd(self, ch, c, fi):
        """An open file descriptor of file fi if it has the known layout, else None."""
        if not c["ok"]:
            return None
        p = self._path(ch, c, fi)
        try:
            fd = os.open(p, os.O_RDONLY)
        except FileNotFoundError:
            return None
        if c["offset"] is None and not self._learn(c, p):
            os.close(fd)
            c["ok"] = False
            return None
        if os.fstat(fd).st_size != c["size"]:
            os.close(fd)
            return None
        return fd

    def _pread(self, fd, c, i0, i1):
        buf = n.empty((i1 - i0, 2), dtype=n.int16)
        nb = os.preadv(fd, [buf], c["offset"] + 4 * i0)
        if nb != buf.nbytes:
            raise IOError("short read")
        return buf

    def _advise(self, ch, c, fi):
        for a in range(1, self.prefetch + 1):
            try:
                fd = os.open(self._path(ch, c, fi + a), os.O_RDONLY)
            except FileNotFoundError:
                return
            os.posix_fadvise(fd, c["offset"], 4 * c["spf"], os.POSIX_FADV_WILLNEED)
            os.close(fd)

    def _second(self, ch, c, fi):
        """int16 (spf, 2) samples of file fi, from the cache or one read; None if not fast."""
        cache = self._cache.setdefault(ch, collections.OrderedDict())
        if fi in cache:
            cache.move_to_end(fi)
            return cache[fi]
        fd = self._fast_fd(ch, c, fi)
        if fd is None:
            return None
        try:
            buf = self._pread(fd, c, 0, c["spf"])
        finally:
            os.close(fd)
        self.n_fast += 1
        if self.prefetch:
            self._advise(ch, c, fi)
        cache[fi] = buf
        while len(cache) > self.cache_seconds:
            cache.popitem(last=False)
        return buf

    def _piece(self, ch, c, fi, i0, i1):
        if self.whole_seconds or fi in self._cache.get(ch, ()):
            b = self._second(ch, c, fi)
            return None if b is None else b[i0:i1]
        fd = self._fast_fd(ch, c, fi)
        if fd is None:
            return None
        try:
            self.n_fast += 1
            return self._pread(fd, c, i0, i1)
        finally:
            os.close(fd)

    def _fallback(self, start_sample, vector_length, ch):
        if self._drf is None:
            from digital_rf import DigitalRFReader
            self._drf = DigitalRFReader(self.rf_dir)
        self.n_fallback += 1
        return self._drf.read_vector_1d(start_sample, vector_length, ch)

    def read_vector_1d(self, start_sample, vector_length, channel_name, sub_channel=0):
        """vector_length complex64 samples from start_sample, as DigitalRFReader.read_vector_1d."""
        c = self._channel(channel_name)
        start_sample, vector_length = int(start_sample), int(vector_length)
        spf = c["spf"]
        pieces = []
        s = start_sample
        end = start_sample + vector_length
        while s < end:
            fi, i0 = divmod(s, spf)
            i1 = min(spf, i0 + end - s)
            p = self._piece(channel_name, c, fi, i0, i1)
            if p is None:
                return self._fallback(start_sample, vector_length, channel_name)
            pieces.append(p)
            s += i1 - i0
        raw = pieces[0] if len(pieces) == 1 else n.concatenate(pieces)
        return raw.astype(n.float32).view(n.complex64).ravel()


# ---- the pulse index ------------------------------------------------------------
# The id_metadata of a recording is one small HDF5 file per second, 256 kB and
# some 450 objects each (four scalars per pulse).  Reading it pulse by pulse
# through DigitalMetadataReader is CPU bound (30 ms a file) and, with many
# processes at once, turns the disk's access pattern random.  pulse_index()
# reads it once in two stages, the files whole and in time order by a few
# processes (disk bound), then decoded from the page cache by all CPUs, and
# keeps the result as one small file for later runs.

INDEX_FIELDS = ("sweepid", "modeid", "sweepnum")
CACHE_DIR = os.path.expanduser("~/.cache/isr_analysis")


def _hour_dirs(md_dir):
    return sorted(e.path for e in os.scandir(md_dir) if e.is_dir())


def _warm(dirs):
    nb = 0
    for d in dirs:
        for f in sorted(os.listdir(d)):
            with open(os.path.join(d, f), "rb", buffering=0) as fh:
                nb += len(fh.read())
    return nb


def _decode(args):
    md_dir, t0, t1 = args
    from digital_rf import DigitalMetadataReader
    rec = DigitalMetadataReader(md_dir).read(t0, t1)
    if not rec:
        return n.zeros((0, 1 + len(INDEX_FIELDS)), dtype=n.int64)
    keys = sorted(rec)
    return n.array([[k] + [int(rec[k][f]) for f in INDEX_FIELDS] for k in keys], dtype=n.int64)


def pulse_index(dirname, cache_dir=CACHE_DIR, n_readers=16, n_proc=None, rebuild=False):
    """
    Every pulse of a recording: {"key": sample index of the pulse (int64,
    sorted), "sweepid", "modeid", "sweepnum"}, from dirname/metadata/id_metadata.

    Built on the first call, which takes a few minutes for two days of data,
    and kept in cache_dir/<recording>/id_index.npz; later calls load it.  The
    cache is rebuilt when the metadata's bounds or number of hour directories
    change, or with rebuild=True.
    """
    import multiprocessing
    from digital_rf import DigitalMetadataReader
    md_dir = os.path.join(dirname, "metadata", "id_metadata")
    bounds = DigitalMetadataReader(md_dir).get_bounds()
    dirs = _hour_dirs(md_dir)
    path = os.path.join(cache_dir, os.path.basename(os.path.normpath(dirname)), "id_index.npz")
    if not rebuild and os.path.exists(path):
        z = n.load(path)
        if tuple(z["bounds"]) == tuple(bounds) and int(z["n_dirs"]) == len(dirs):
            return {k: z[k] for k in ("key",) + INDEX_FIELDS}
    ctx = multiprocessing.get_context("fork")
    with ctx.Pool(n_readers) as pool:
        pool.map(_warm, [[d] for d in dirs], chunksize=1)
    step = 600 * 10**6
    spans = [(md_dir, t, min(t + step - 1, bounds[1])) for t in range(bounds[0], bounds[1] + 1, step)]   # read() includes its end
    with ctx.Pool(n_proc or os.cpu_count()) as pool:
        a = n.concatenate(pool.map(_decode, spans, chunksize=1))
    a = a[n.argsort(a[:, 0], kind="stable")]
    out = {"key": a[:, 0]}
    out.update({f: a[:, 1 + i] for i, f in enumerate(INDEX_FIELDS)})
    os.makedirs(os.path.dirname(path), exist_ok=True)
    n.savez(path + ".tmp.npz", bounds=n.array(bounds), n_dirs=len(dirs), **out)
    os.replace(path + ".tmp.npz", path)
    return out
