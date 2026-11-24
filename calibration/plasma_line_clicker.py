import numpy as n
import glob
import h5py
import matplotlib.pyplot as plt
from digital_rf import DigitalRFReader, DigitalMetadataReader, DigitalMetadataWriter
import sys

# range window, in km, over which each record's plasma-line spectrum is
# maximised.  It used to be the gate indices 50:150, which in the 288-gate
# records of the 2023 and 2024 recordings (3.745 km gates from 0 km) are
# 187.4 to 558.4 km; this window selects exactly those gates.
PLASMA_LINE_RANGES_KM=(185.0,560.0)
# records with fewer gates than this inside the window are skipped (the
# 25-gate long-range MISA records, which have 2)
MIN_GATES=50

def click_spec(misa_t,misa_s,freqs,dirname,calname="cal.h5"):
    # pip install mpl_point_clicker
    from mpl_point_clicker import clicker

    fig, ax = plt.subplots(constrained_layout=True)

    c=ax.pcolormesh(misa_t,freqs,misa_s.T,vmin=0,vmax=0.1,cmap="plasma")
#    fig.colorbar(c,ax=ax)
    ax.set_title(calname)
    

    #zoom_factory(ax)
    #ph = panhandler(fig, button=2)
    klicker = clicker(
        ax,
       ["plasma_frequency"],
       markers=["x"]
    )

    from typing import Tuple
    def point_added_cb(position: Tuple[float, float], klass: str):
        x, y = position
        print(f"New point of class {klass} added at {x=}, {y=}")

        kp=klicker.get_positions()
        print(kp)
        print("writing %s/%s.h5"%(dirname,calname))
        ho=h5py.File("%s/%s.h5"%(dirname,calname),"w")
        ho["plasma_frequency"]=kp["plasma_frequency"]
        ho.close()

    def point_removed_cb(position: Tuple[float, float], klass: str, idx):
        x, y = position

        suffix = {'1': 'st', '2': 'nd', '3': 'rd'}.get(str(idx)[-1], 'th')
        print(
            f"The {idx}{suffix} point of class {klass} with position {x=:.2f}, {y=:.2f}  was removed"
        )

    klicker.on_point_added(point_added_cb)
    klicker.on_point_removed(point_removed_cb)
    plt.show()


def gate_window(ranges,ranges_km=PLASMA_LINE_RANGES_KM):
    """Indices of the range gates inside ranges_km, as a slice (they are contiguous)."""
    r=n.asarray(ranges)
    idx=n.where((r>=ranges_km[0])&(r<ranges_km[1]))[0]
    if len(idx)==0:
        return slice(0,0)
    return slice(int(idx[0]),int(idx[-1])+1)


def build_spectra(dirname,ranges_km=PLASMA_LINE_RANGES_KM,min_gates=MIN_GATES):
    """
    Per minute of the recording, the mean over the records of each antenna of
    the plasma-line spectrum's maximum over the gates in ranges_km (after
    removing each gate's median), minus its median.  Returns
    (zenith_t, zenith_s, misa_t, misa_s, freqs).
    """
    plmd=DigitalMetadataReader("%s/metadata/integrated_plasma_line_metadata_hires/"%(dirname))
    b=plmd.get_bounds()

    dt=plmd.get_file_cadence_secs()

    n_files=int((b[1]-b[0])/dt)

    freqs=n.fft.fftshift(n.fft.fftfreq(12000,d=1/25e6))#d["freqs"][0]

    misa_t=[]
    zenith_t=[]
    misa_s=[]
    zenith_s=[]
    for fi in range(n_files):
        print("%d/%d"%(fi,n_files))
        try:
            d = plmd.read_flatdict(b[0]+fi*60,b[0]+fi*60+60)
        except:
            continue
        if "spec" in d.keys():
            n_misa=0
            n_zenith=0
            z_t=0
            m_t=0
            z_s=n.zeros(len(freqs))
            m_s=n.zeros(len(freqs))

            for si in range(len(d["spec"])):
                g=gate_window(d["ranges"][si],ranges_km)
                use=(g.stop-g.start)>=min_gates
                if d["antenna"][si] == "MISA":
                    print(d["spec"][si].shape)
                    if use:
                        m_s+=n.nanmax(d["spec"][si][:,g]-n.nanmedian(d["spec"][si][:,g],axis=0),axis=1)
                        m_t+=d["t1"][si]
                        n_misa+=1
                    else:
                        print("ignoring")
                    print("misa")
                else:
                    if use:
                        z_s+=n.nanmax(d["spec"][si][:,g]-n.nanmedian(d["spec"][si][:,g],axis=0),axis=1)
                        z_t+=d["t1"][si]
                        n_zenith+=1.0
                    else:
                        print("ignoring")
                    print("zenith")

            if n_misa>0:
                m_s=m_s/n_misa
                m_s=m_s-n.nanmedian(m_s)
                m_t=m_t/n_misa
                misa_s.append(m_s)
                misa_t.append(m_t)
            if n_zenith>0:
                print(n_zenith)
                z_s=z_s/n_zenith
                z_s=z_s-n.nanmedian(z_s)
                z_t=z_t/n_zenith
                print(z_t)
                zenith_s.append(z_s)
                zenith_t.append(z_t)

        else:
            print("no spec key")
    return zenith_t,zenith_s,misa_t,misa_s,freqs


if __name__ == "__main__":
    dirname=sys.argv[1]#"/media/j/fee7388b-a51d-4e10-86e3-5cabb0e1bc13/isr/2023-09-05/usrp-rx0-r_20230905T214448_20230906T040054/"
    zenith_t,zenith_s,misa_t,misa_s,freqs=build_spectra(dirname)
    zenith_s=n.vstack(zenith_s)
    click_spec(zenith_t,zenith_s,freqs,dirname,calname="zenithcal")
    if len(misa_s) > 0:
        misa_s=n.vstack(misa_s)
        click_spec(misa_t,misa_s,freqs,dirname,calname="misacal")
