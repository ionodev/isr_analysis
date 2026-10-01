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


# seconds before a cycle-end event from which the antenna is treated as
# unknown when the next cycle uses the other antenna (memo 30)
ANTENNA_SWITCH_GUARD_S=2.5

def get_antenna_select(dirn,plot=False,switch_guard_s=ANTENNA_SWITCH_GUARD_S):
    """
    Transmit and receive antenna selection as functions of time.

    1 = misa
    -1 = zenith
    0 = unknown: around a change of antenna, see below

    Returns (tx_sel, rx_sel). Both are step functions of time in microseconds,
    holding the most recent switch, and holding the first and last value
    outside the recorded range.

    The metadata change the antenna at the event that opens the next record
    cycle, but the radar changes it at the event that closes the previous
    one: in the eclipse2024 recording the first pulse on the new antenna comes
    2.3-20.8 s (median 8.6 s) before the metadata say so, after a pause of
    about 1.25 s without pulses, and the first new pulse comes up to 2 s
    before the closing event (memo 30).  So from switch_guard_s seconds
    before the closing event until the opening event, where the next cycle
    uses the other antenna, both functions return 0, and tests like
    tx_ant(t)<=-0.99 reject those pulses.  The guard was chosen on this
    recording.  Before the first event the functions return the first
    recorded value, which is the antenna of the cycle in progress when the
    recording starts, also when that first event closes a cycle.
    switch_guard_s=None gives the metadata as recorded, and is also used,
    with a message, when the metadata have no cycle_name for every event.
    """
    print("Reading transmit power meter metadata. Might take a few seconds")
    dmd=DigitalMetadataReader(dirn)
    b=dmd.get_bounds()

    tx_t=[]
    tx_v=[]
    rx_t=[]
    rx_v=[]

    if switch_guard_s is not None:
        try:
            names=dmd.read(b[0],b[1],"cycle_name")
        except Exception as e:
            print("antenna select: no cycle_name metadata (%s), "
                  "using the antenna metadata as recorded"%e)
            switch_guard_s=None
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
    if switch_guard_s is not None and not (set(tx_t)|set(rx_t))<=set(names.keys()):
        print("antenna select: cycle_name is missing for some antenna events, "
              "using the antenna metadata as recorded")
        switch_guard_s=None
    # held before the first event: the antenna of the cycle in progress
    tx_v0=tx_v[0]
    rx_v0=rx_v[0]
    if switch_guard_s is not None:
        tx_t,tx_v=_unknown_at_switches(tx_t,tx_v,names,switch_guard_s)
        rx_t,rx_v=_unknown_at_switches(rx_t,rx_v,names,switch_guard_s)
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
    rx_sel=_step(rx_t,rx_v,rx_v0)
    tx_sel=_step(tx_t,tx_v,tx_v0)

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



def _step(t,v,v_before):
    """
    Step function of time holding the most recent value of v, v_before before
    the first time and the last value after the last.
    """
    return sint.interp1d(t,v,kind="previous",
                         bounds_error=False,fill_value=(v_before,v[-1]))


def _unknown_at_switches(t,v,names,guard_s):
    """
    Value 0 from guard_s before each cycle-end event (an event with an empty
    cycle name) until the next event, where that next event changes the
    antenna.  t (microseconds) and v are one field's events, sorted; names
    maps event times to cycle names.  Returns the new (t, v).
    """
    t=n.asarray(t,dtype=n.int64)
    v=n.asarray(v,dtype=float)
    def unnamed(k):
        nm=names.get(k,b"")
        return (nm.decode() if isinstance(nm,bytes) else str(nm))==""
    sw=[i for i in range(1,len(t)) if v[i]!=v[i-1] and unnamed(t[i-1])]
    n_change=int(n.sum(v[1:]!=v[:-1]))
    if len(sw)<n_change:
        print("antenna select: %d of %d antenna changes do not follow a "
              "cycle-end event and keep the metadata's time"%(n_change-len(sw),n_change))
    v=n.copy(v)
    v[[i-1 for i in sw]]=0.0
    # the guard starts no earlier than the event before the closing one
    g=[max(t[i-1]-int(guard_s*1e6),t[i-2]+1 if i>=2 else t[i-1]-int(guard_s*1e6)) for i in sw]
    tt=n.concatenate([t,n.array(g,dtype=n.int64)])
    vv=n.concatenate([v,n.zeros(len(g))])
    o=n.argsort(tt,kind="stable")
    return tt[o],list(vv[o])


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
