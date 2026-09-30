import numpy as n
import matplotlib.pyplot as plt
from digital_rf import DigitalRFReader, DigitalMetadataReader, DigitalMetadataWriter
import scipy.interpolate as sint

radar_lat=42.61932878636544
radar_lon=-71.49124624803031
radar_hgt=146.0

# The zenith antenna does not point at the geodetic zenith: its beam is at
# elevation 88.16 deg, azimuth 172.9 deg (memo 7), from satellite echoes and
# confirmed by the transit of Cygnus A. The metadata report 90 deg. The angular
# uncertainty of this fit has not been estimated (memo 23 depends on it).
ZENITH_BEAM_EL_DEG=88.16
ZENITH_BEAM_AZ_DEG=172.9

def get_tx_power_model(dirn,plot=False):
    print("Reading transmit power meter metadata. Might take a few seconds")
    dmd=DigitalMetadataReader(dirn)
    b=dmd.get_bounds()

    zenith_t=[]
    zenith_pwr=[]
    misa_t=[]
    misa_pwr=[]
    sid = dmd.read(b[0],b[1],"zenith_power")    
    for keyi,key in enumerate(sid.keys()):
        zenith_t.append(key)
        zenith_pwr.append(sid[key])
#        print("%1.2f %1.2f"%(key,sid[key]))
        
    sid = dmd.read(b[0],b[1],"misa_power")    
    for keyi,key in enumerate(sid.keys()):
        misa_t.append(key)
        misa_pwr.append(sid[key])
#        print("%1.2f %1.2f"%(key,sid[key]))

    if plot:
        plt.plot(zenith_t,zenith_pwr,label="Zenith")
        plt.plot(misa_t,misa_pwr,label="MISA")
        plt.legend()
        plt.show()

    zenith_t[0]=zenith_t[0]-3600
    zenith_t[-1]=zenith_t[-1]+3600
    misa_t[0]=misa_t[0]-3600
    misa_t[-1]=misa_t[-1]+3600
    zenith_pwrf=sint.interp1d(zenith_t,zenith_pwr)
    misa_pwrf=sint.interp1d(misa_t,misa_pwr)
    
    return(zenith_pwrf,misa_pwrf)


def get_antenna_select(dirn,plot=False):
    """
    Transmit and receive antenna selection as functions of time.

    1 = misa
    -1 = zenith

    Returns (tx_sel, rx_sel). Both are step functions of time in microseconds,
    holding the most recent switch, and holding the first and last value
    outside the recorded range.
    """
    print("Reading transmit power meter metadata. Might take a few seconds")
    dmd=DigitalMetadataReader(dirn)
    b=dmd.get_bounds()

    tx_t=[]
    tx_v=[]
    rx_t=[]
    rx_v=[]

    sid = dmd.read(b[0],b[1],"rx_antenna")    
    for keyi,key in enumerate(sid.keys()):

        rx_t.append(key)                    
        if sid[key]==b'MISA':
            rx_v.append(1.0)
        else:
            rx_v.append(-1.0)

    sid = dmd.read(b[0],b[1],"tx_antenna")    
    for keyi,key in enumerate(sid.keys()):

        tx_t.append(key)
  #      print("%s %d"%(sid[key],key))
        if sid[key]==b'MISA':
            tx_v.append(1.0)
        else:
            tx_v.append(-1.0)
            
    rx_t=n.array(rx_t)
    tx_t=n.array(tx_t)
    rx_t0=n.copy(rx_t)
    tx_t0=n.copy(tx_t)    
    
    # the antenna selection is a step function: it holds a value until the next
    # switch. interpolating it linearly makes the value sweep through zero at
    # every switch, so tests like tx_ant(t)<=-0.99 fail for a slice of each
    # transition and those pulses are discarded.
    #
    # the previous version also moved the first sample 24 hours earlier and the
    # last 24 hours later to extend the range. with linear interpolation that
    # ramps from the first value towards the second across the whole 24 hours,
    # so times before the first event came back close to the *second* event's
    # value, which is the opposite answer whenever the antenna switched. holding
    # the first and last value is what extending the range should mean.
    rx_sel=sint.interp1d(rx_t,rx_v,kind="previous",
                         bounds_error=False,fill_value=(rx_v[0],rx_v[-1]))
    tx_sel=sint.interp1d(tx_t,tx_v,kind="previous",
                         bounds_error=False,fill_value=(tx_v[0],tx_v[-1]))

    if plot:
        t=n.linspace(b[0],b[1],num=100000)
        plt.plot(t/1e6,rx_sel(t),label="rx")
 #       for rxt in rx_t0:
#            plt.axvline(rxt/1e6)
        plt.plot(t/1e6,tx_sel(t),label="tx")
        plt.ylim([-1.1,1.1])
        plt.legend()
        plt.show()
    # tx first: every caller unpacks this as tx_ant,rx_ant. the previous order
    # was the other way round, which was harmless only because both are tested
    # for the same value.
    return(tx_sel,rx_sel)



def get_misa_pointing(dirn):
    """
    MISA's commanded pointing at any sample index, read as the experiment
    schedule the antenna control metadata are: each record cycle starts with
    an event that gives MISA's azimuth and elevation, and holds until the next
    event.  Returns f(keys) -> (az, el, cycle, stationary) for sample indices
    (microseconds): azimuth (0-360) and elevation in degrees of the latest
    event at or before each key, the cycle name of the latest named event, and
    whether the latest event opens one of MISA's own record cycles (its name
    starts with "misa"), during which MISA holds still.  In the zenith cycles
    and the gaps between cycles MISA may be slewing to its next position.

    get_misa_az_el_model interpolates between events, which puts a period near
    a cycle change between two positions; this is a step function.
    """
    acmd=DigitalMetadataReader(dirn)
    b=acmd.get_bounds()
    ev=acmd.read(b[0],b[1])
    keys=n.array(sorted(ev.keys()),dtype=n.int64)
    az=n.array([float(ev[k]["misa_azimuth"])%360.0 for k in keys])
    el=n.array([float(ev[k]["misa_elevation"]) for k in keys])
    name=[ev[k].get("cycle_name","") for k in keys]
    name=n.array([x.decode() if isinstance(x,bytes) else str(x) for x in name])
    named=n.array(name,dtype=object)
    last=""
    for i in range(len(named)):
        if named[i]:
            last=named[i]
        named[i]=last
    stat=n.array([x.startswith("misa") for x in name])
    def f(k):
        k=n.atleast_1d(n.asarray(k,dtype=n.int64))
        i=n.clip(n.searchsorted(keys,k,side="right")-1,0,len(keys)-1)
        return az[i],el[i],named[i],stat[i]
    return f


def get_misa_az_el_model(dirn):
    acmd=DigitalMetadataReader(dirn)
    b=acmd.get_bounds()
    azk = acmd.read(b[0],b[1],"misa_azimuth")
    az=[]
    azt=[]
    el=[]
    elt=[]
    for keyi,key in enumerate(azk.keys()):
        az.append(azk[key])
        azt.append(key/1e6)
        
    elk = acmd.read(b[0],b[1],"misa_elevation")        
    for keyi,key in enumerate(elk.keys()):
        el.append(elk[key])
        elt.append(key/1e6)

    azt=n.array(azt)
    elt=n.array(elt)
    az=n.array(az)
    el=n.array(el)
    azt[0]=azt[0]-24*3600.0
    azt[-1]=azt[-1]+24*3600.0    
    elt[0]=elt[0]-24*3600.0
    elt[-1]=elt[-1]+24*3600.0

    azf=sint.interp1d(azt,az)
    elf=sint.interp1d(elt,el)
    return(azf,elf,[b[0]/1e6,b[1]/1e6])
    
        
        
        
    

if __name__ == "__main__":
    import sys


    #    get_tx_power_model(sys.argv[1],plot=True)    
    get_tx_power_model("%s/metadata/powermeter"%(sys.argv[1]),plot=True)
    az,el,b=get_misa_az_el_model("%s/metadata/antenna_control_metadata"%(sys.argv[1]))
    t=n.linspace(b[0],b[1],num=10000)
    plt.plot(t,az(t))
    plt.plot(t,el(t))
    plt.show()

    get_antenna_select("%s/metadata/antenna_control_metadata"%(sys.argv[1]),plot=True)

    
  #  import jcoord
   # llh=jcoord.az_el_r2geodetic(radar_lat,radar_lon,146,az(b[0]+3600),el(b[0]+3600),400e3)
    #    print(llh)
    #    
