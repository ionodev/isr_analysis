"""
Where the analysis scripts find their data, from environment variables:

    ISR_DATA_DIR               the recording (the folder with rf_data/ and metadata/)
    ISR_SATELLITE_CATALOGUE    the satellite catalogue at P_FA 1e-10 (the detector's default)
    ISR_QUIET_CATALOGUE        the catalogue at P_FA 1e-6 of the quiet reference cycle

e.g. export ISR_DATA_DIR=/path/to/usrp-rx0-r_20240407T100000_20240409T110000.
A script that needs one that is not set stops with a message naming it.
"""
import os


class _Unset:
    """Stands in for a path whose environment variable is not set; using it
    as a path or a string raises an error that names the variable."""

    def __init__(self, var, what):
        self.var, self.what = var, what

    def _fail(self, *args, **kwargs):
        raise RuntimeError("set the environment variable %s to %s" % (self.var, self.what))

    __str__ = __fspath__ = __format__ = _fail

    def __repr__(self):
        return "<%s not set>" % self.var


def _path(var, what):
    v = os.environ.get(var)
    return v if v else _Unset(var, what)


DATA_DIR = _path("ISR_DATA_DIR", "the recording's folder (with rf_data/ and metadata/)")
SATELLITE_CATALOGUE = _path("ISR_SATELLITE_CATALOGUE", "the satellite catalogue at P_FA 1e-10")
QUIET_CATALOGUE = _path("ISR_QUIET_CATALOGUE", "the P_FA 1e-6 satellite catalogue of the quiet reference cycle")
