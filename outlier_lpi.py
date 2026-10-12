import os
# mpi does paralelization, both multithread matrix operations
# just one per process to avoid cache trashing. otherwise each mpi process will
# try to use all cpus for the linear algebra, which slows things to a halt.
#os.system("export OMP_NUM_THREADS=1")
os.environ["OMP_NUM_THREADS"] = "1"

import numpy as n
import matplotlib.pyplot as plt
import pyfftw
import stuffr
import scipy.signal as s
from digital_rf import DigitalRFReader, DigitalMetadataReader, DigitalMetadataWriter

import h5py
from scipy.ndimage import median_filter
from scipy import sparse
from mpi4py import MPI
import scipy.signal as ss
import scipy.constants as c
import traceback
import time

import millstone_radar_state as mrs
# transmit and receive gate timing, shared with avg_range_doppler_spec.py and
# tx_delay.py rather than copied into each
from radar_timing import TMM as tmm, T_INJECTION as T_injection
# design matrix columns for detected satellite echoes
import satellite_columns as satcol
from raw_reader import RawReader
from impulse_blanking import BlankingReader

comm=MPI.COMM_WORLD
size=comm.Get_size()
rank=comm.Get_rank()

fft=pyfftw.interfaces.numpy_fft.fft
ifft=pyfftw.interfaces.numpy_fft.ifft

# we really need to be narrow band to avoid interference. there is plenty of it.
# but we can't go too narrow, because then we lose the ion-line.
# 1000 m/s is 3 kHz and the code itself is 66.6 kHz
# 66.6/2 + 3 = 70 kHz, and thus the maximum frequency offset will be 35 kHz.
# that is tight, but hopefully enough to filter out the interference.
#pass_band=0.1e6
def ideal_lpf(z,sr=1e6,f0=1.2*0.1e6,L=200):
    m=n.arange(-L,L)+1e-6
    om0=n.pi*f0/(0.5*sr)
    h=s.windows.hann(len(m))*n.sin(om0*m)/(n.pi*m)

    Z=n.fft.fft(z)
    H=n.fft.fft(h,len(Z))
    z_filtered=n.roll(n.fft.ifft(Z*H),-L)

    return(z_filtered)

class simple_decimator:
    def __init__(self,L=10000,dec=10):
        self.L=L
        self.dec=dec
        decL=int(n.floor(L/dec))
        self.idxm=n.zeros([decL,dec],dtype=int)
        for ti in range(decL):
            self.idxm[ti,:]=n.arange(dec,dtype=int) + ti*dec
            
    def decimate(self,z):
        """ 
        this decimate has to be a sum to ensure all range resolutions have the same 
        magic constant.
        """
        decL=int(n.floor(len(z)/self.dec))
        return(n.sum(z[self.idxm[0:decL,:]],axis=1))

class fft_lpf:
    def __init__(self,z_len=10000,sr=1e6,f0=1.2*0.1e6,L=20):
        m=n.arange(-L,L)+1e-6
        om0=n.pi*f0/(0.5*sr)
        h=s.windows.hann(len(m))*n.sin(om0*m)/(n.pi*m)
        # normalize to impulse response to unity.
        #h=n.array(h/n.sum(n.abs(h)**2.0),dtype=n.complex64)
        # unity gain at DC
        h=n.array(h/n.sum(h),dtype=n.complex64)
        self.h=h
        #pyfftw.interfaces.numpy_fft.fft()
        self.H=fft(h,z_len)
        self.L=L
        
    def lpf(self,z):
        return(n.roll(ifft(self.H*fft(z)),-self.L))
        

def ideal_lpf_h(sr=1e6,f0=1.2*0.1e6,L=200):
    m=n.arange(-L,L)+1e-6
    om0=n.pi*f0/(0.5*sr)
    h=s.windows.hann(len(m))*n.sin(om0*m)/(n.pi*m)
    return(h)

def estimate_dc(d_il,tmm,sid,channel):
    # estimate dc offset first
    z_dc=n.zeros(10000,dtype=n.complex64)
    n_dc=0.0
    for keyi,key in enumerate(sid.keys()):
        if sid[key] not in tmm.keys():
            print("pulse code %d is not in the timing table, ignoring"%(sid[key]))
            continue
        # fftw "allocated vector"
        z_echo = d_il.read_vector_1d(key, 10000, channel).astype("c8", casting="unsafe", copy=False)
        last_echo=tmm[sid[key]]["last_echo"]        
        gc=tmm[sid[key]]["gc"]
        
        z_echo[0:(gc+4000)]=n.nan
        z_echo[last_echo:10000]=n.nan
        z_dc+=z_echo
        n_dc+=1.0
    z_dc=z_dc/n_dc

    if False:
        # plot the dc offset estimate
        plt.plot(z_dc.real)
        plt.plot(z_dc.imag)
        plt.axhline(n.nanmedian(z_dc.real))
        plt.axhline(n.nanmedian(z_dc.imag))
        print(n.nanmedian(z_dc.real))
        print(n.nanmedian(z_dc.imag))        
        plt.show()
    
    z_dc=n.complex64(n.nanmedian(z_dc.real)+n.nanmedian(z_dc.imag)*1j)
    return(z_dc)
    

def invert_normal(ATA):
    """
    Error covariance (the inverse of the normal matrix ATA), and the indices
    of the unknowns no measurement constrains.

    Such an unknown (a zero column of the design matrix: at some lags the
    lowest range gate of the first periods of a recording) makes ATA
    singular, which used to lose the whole lag. Invert for the others and
    give it NaN on the diagonal only, so that the NaN reaches the estimate
    and variance of that unknown and nothing else. Only an exactly zero
    diagonal counts as unconstrained: a NaN in ATA takes the plain inverse,
    as before. With every unknown constrained this is the plain inverse.
    """
    d=n.real(n.diag(ATA))
    unconstrained=n.where(d==0)[0]
    if len(unconstrained)==0:
        return n.linalg.inv(ATA),unconstrained
    fi=n.where(d!=0)[0]
    Sigma=n.zeros(ATA.shape,dtype=ATA.dtype)
    Sigma[n.ix_(fi,fi)]=n.linalg.inv(ATA[n.ix_(fi,fi)])
    Sigma[unconstrained,unconstrained]=n.nan
    return Sigma,unconstrained


def convolution_matrix(envelope, rmin=0, rmax=100):
    """
    we imply that the number of measurements is equal to the number of elements
    in code

    Use the index matrix (idxm) to efficiently grab the numbers from a 1d array to build a matrix
    A = z_tx_envelope[idxm]
    """
    L = len(envelope)
    ridx = n.arange(rmin, rmax,dtype=int)
    A = n.zeros([L, rmax - rmin], dtype=n.complex64)
    idxm = n.zeros([L, rmax - rmin], dtype=int)    
    for i in n.arange(L):
        A[i, :] = envelope[(i - ridx) % L]
        idxm[i,:] = n.array(n.mod( i-ridx, L ),dtype=int)
    result = {}
    result["A"] = A
    result["ridx"] = ridx
    result["idxm"] = idxm
    return(result)

#
# tbd: add range gates of different sizes
#
def lpi_files(dirname="/media/j/fee7388b-a51d-4e10-86e3-5cabb0e1bc13/isr/2023-09-05/usrp-rx0-r_20230905T214448_20230906T040054",
              avg_dur=10,  # n seconds to average
              channel="zenith-l",
              rg=60,       # how many microseconds is one range gate
              min_tx_frac=0.5,  # how much of the pulse can be missing due to ground clutter clipping, defines the minimum range gate
              reanalyze=False,
              pass_band=0.1e6,
              filter_len=20,
              use_long_pulse=True,
              maximum_range_delay=7000,    # microseconds. defines the highest range to analyze
              save_acf_images=True,
              min_tx_pwr=400e3,
              fft_len=1024,                 # store diagnostic spectrum for RFI identification
              lags=n.arange(1,46,dtype=int)*10,
              lag_avg=1,
              output_base=None,
              max_time_s=None,
              # detected satellite echoes, modelled instead of rejected. An
              # object whose for_period(i0, i1) returns {pulse key: [(delay in
              # samples relative to the transmit template, Doppler in Hz),
              # ...]}, e.g. satellite_columns.CatalogueDetections. Each echo
              # gets its own design matrix column with a free complex amplitude
              # per lag. None leaves the inversion as it was.
              satellite_detections=None,
              # transmit template the satellite columns are built from:
              # "average" the leakage of its pulse code averaged over the period,
              # whose lower noise lets brighter echoes be modelled; "pulse" each
              # pulse's own leakage, as the plasma columns use
              satellite_template="average",
              # compute the per row noise weights in double precision; see the
              # comment where it is used for why a bright echo needs it. False
              # reproduces the output of the code before this option, which
              # differs only at roundoff level when no bright echo is present.
              precise_weights=True,
              # how the DC offset is removed before the noise injection
              # calibration: "period" estimates it once from all the pulses of
              # the integration period; "window" subtracts each pulse's own
              # background window mean, which biased T_sys and alpha low by a
              # signal dependent factor (memo 7) and reproduces the earlier
              # output exactly.
              noise_dc="period",
              # reject lagged products by the ratio test and the local power
              # cut below. Off when satellite echoes are to be modelled, as
              # those tests reject the very products the columns explain.
              outlier_rejection=True,
              # indices of the integration periods to analyse, None for all
              periods=None,
              # read the raw voltage with raw_reader.RawReader, one read per
              # second of data, instead of digital_rf's several small reads per
              # pulse. The samples are identical; False uses DigitalRFReader.
              fast_read=True,
              # set impulsive interference (power-line sparks, memo 21) to zero
              # in the raw voltage before the inversion: impulse_blanking.
              # BlankingReader. False leaves the samples as recorded.
              blank_impulses=False,
              ):
    if output_base is None:
        output_base = dirname
    print("mkdir -p %s/lpi_%d/%s"%(output_base,rg,channel))
    os.system("mkdir -p %s/lpi_%d/%s"%(output_base,rg,channel))
    
        
    id_read = DigitalMetadataReader("%s/metadata/id_metadata"%(dirname))
    d_il = RawReader("%s/rf_data/"%(dirname)) if fast_read else DigitalRFReader("%s/rf_data/"%(dirname))
    # the reader as recorded, for the calibration windows (impulse_blanking)
    d_raw = d_il
    if blank_impulses:
        d_il = BlankingReader(d_il, dirname)

    zpm,mpm=mrs.get_tx_power_model("%s/metadata/powermeter"%(dirname))
    tx_ant,rx_ant=mrs.get_antenna_select("%s/metadata/antenna_control_metadata"%(dirname))    

    idb=id_read.get_bounds()
    # sample rate for metadata
    idsr=1000000
    # sample rate for ion line channel
    sr=1000000

    plot_voltage=False
    use_ideal_filter=True
    debug_gc_rem=False

    # how many integration cycles do we have
    n_times = int(n.floor((idb[1]-idb[0])/idsr/avg_dur))
    if max_time_s is not None:
        n_times = min(n_times, int(max_time_s / avg_dur))

    # which lags to calculate
    
    # how many lags do we average together?
    

    # calculate the average lag value
    n_lags=len(lags)-lag_avg+1
    mean_lags=n.zeros(n_lags)
    for i in range(n_lags):
        mean_lags[i]=n.mean(lags[i:(i+lag_avg)])

    # maximum number of microseconds of delay, which we analyze
    # this is experiment specific; lpi_files takes it as an argument, and
    # run_analysis.py from the configuration (max_range_delay_us)
    n_rg=int(n.floor(maximum_range_delay/rg))
    rgs=n.arange(n_rg)*rg
    rmax=n_rg

    # round trip speed of light in vacuum propagation, one microsecond
    rg_1us=c.c/1e6/2.0/1e3

    # range gates
    rgs_km=rgs*rg_1us

    # first entry in tx pulse metadata
    i0=idb[0]

    lpf=fft_lpf(10000,f0=1.2*pass_band,L=filter_len)

    decim=simple_decimator(L=10000,dec=rg)

    pwr_spec=n.zeros(fft_len,dtype=n.float32)
    n_pwr_spec=0.0
    spec_window=ss.windows.hann(fft_len)

    

    period_list = list(range(n_times)) if periods is None else [p for p in periods if p < n_times]

    # go through one integration window at a time
    for ai in period_list[rank::size]:

        
        i0 = ai*int(avg_dur*idsr) + idb[0]

        if os.path.exists("%s/lpi_%d/%s/lpi-%d.png"%(output_base,rg,channel,int(i0/1e6))) and reanalyze==False:
            print("already analyzed %d"%(i0/1e6))
            continue

        # get info on all the pulses transmitted during this averaging interval
        # get some extra for gc
        sid = id_read.read(i0,i0+int(avg_dur*idsr)+40000,"sweepid")

        n_pulses=len(sid.keys())

        # satellite echoes detected in this period, and one entry per satellite
        # column: (pulse key, delay in samples, Doppler in Hz)
        sat_det={}
        if satellite_detections is not None:
            sat_det=satellite_detections.for_period(i0,i0+int(avg_dur*idsr)+40000)
        sat_cols=[]

        # USRP DC offset bug due to truncation instead of rounding.
        # Ryan Volz has a fix for firmware in USRPs.
        # note that this appears to change as a function of time
        # we can probably only estimate this from the estimated autocorrelation functions
        z_dc=n.complex64(-0.212-0.221j)
        # usrp n200 is fixed
        if channel == "zenith-l2":
            z_dc=0.0

        # the per code average transmit template for the satellite columns
        sat_avg_tx={}
        if len(sat_det) > 0 and satellite_template == "average":
            def keep_pulse(k):
                if channel in ("zenith-l","zenith-l2"):
                    return (tx_ant(k) <= -0.99) and (rx_ant(k) <= -0.99) and (zpm(k/1e6) >= min_tx_pwr)
                return (tx_ant(k) >= 0.99) and (rx_ant(k) >= 0.99) and (mpm(k/1e6) >= min_tx_pwr)
            skeys=list(sid.keys())
            sat_avg_tx=satcol.average_templates(d_il,skeys[3:len(skeys)-3],sid,channel,
                                                z_dc,tmm,keep_pulse)
            for code in sorted(sat_avg_tx.keys()):
                print("satellite template, code %d: %d pulses, scatter %.1f dB"%(
                    code,sat_avg_tx[code][1],10*n.log10(sat_avg_tx[code][2]+1e-30)))
        
        bg_samples=[]
        bg_plus_inj_samples=[]
        z_dc_samples=[]

        sidkeys=list(sid.keys())

        A=[]
        # satellite column vectors, one list per lag holding for every entry of
        # A[li] the [(column index, column), ...] of that pulse
        S=[]
        mgs=[]
        mes=[]
#        sigmas=[]
        idxms=[]
        rmins=[]

        sample0=800
        sample1=8200
        rdec=rg
        m0=int(n.round(sample0/rdec))
        m1=int(n.round(sample1/rdec))

        n_meas=m1-m0

        # count the number of good measurements encountered as a function of delay
        # in lagged products
        ok_count = n.zeros(n_meas,dtype=int)
        meas_count = n.zeros(n_meas,dtype=int)
        meas_delays_us = n.arange(m0,m1)*rdec
        
        pwr_spec[:]=0.0
        n_pwr_spec=0.0


        for li in range(n_lags):
            # determine what is the lowest range that can be estimated
            # rg0=(gc - txstart - 0.6*pulse_length + lag)/range_decimation
            rmin=int(n.round((sample0-111-480*min_tx_frac+lags[li])/rdec))
            cm=convolution_matrix(n.zeros(m1),rmin=rmin,rmax=rmax)
            rmins.append(rmin)
            idxms.append(cm["idxm"])
            A.append([])
            S.append([])
            mgs.append([])
            mes.append([])
        n_good_estimates=0

        avg_pwr=0.0
        avg_pwr_n=0

        # pulse codes seen in this integration period that the timing table
        # does not cover, and how many pulses each cost us
        unknown_codes={}
        long_codes={}

        # start at 3, because we may need to look back for GC
        for keyi in range(3,n_pulses-3):
            
            t0=time.time()
            key=sidkeys[keyi]

            zenith_pwr=zpm(key/1e6)
            misa_pwr=mpm(key/1e6)

            if (channel == "zenith-l") or (channel=="zenith-l2"):
                if (tx_ant(key) > -0.99) or (rx_ant(key) > -0.99) or (zenith_pwr < min_tx_pwr):
                    print("no zenith data. P_tx %1.2f (MW) skipping"%(zenith_pwr/1e6))
                    continue
                else:
                    avg_pwr+=zenith_pwr
                    avg_pwr_n+=1            

            if channel == "misa-l":
                if (tx_ant(key) < 0.99) or (rx_ant(key) < 0.99) or (misa_pwr < min_tx_pwr):
                    print("no misa data. skipping")
                    continue
                else:
                    avg_pwr+=misa_pwr
                    avg_pwr_n+=1
                    

            if sid[key] not in tmm.keys():
                # counted rather than printed here: an experiment with a mode
                # the timing table does not know would otherwise print one line
                # per pulse, hundreds per integration period
                unknown_codes[sid[key]] = unknown_codes.get(sid[key], 0) + 1
                continue
            if tmm[sid[key]]["read_length"] > 10000:
                # the inversion reads 10000 samples per pulse; a mode whose echo
                # window reaches further (MISA's low elevation mode 800, to 30000
                # us) belongs to avg_range_doppler_spec.py, and used to crash
                # the period on an empty slice (memo 20)
                long_codes[sid[key]] = long_codes.get(sid[key], 0) + 1
                continue

            z_echo=None
            zd=None
            next_key=None

            try:
                z_echo = d_il.read_vector_1d(key, 10000, channel).astype("c8", casting="unsafe", copy=False) - z_dc
            except:
                traceback.print_exc()
                print("couldn't read echo")
                continue

            # no filtering of tx to get better ambiguity function.
            # note that the transmit envelope is the pulse leaking into this same
            # echo channel, not the tx-h channel, so the tx-h to echo channel
            # delay (see tx_delay.py) cancels here and needs no correction.
            z_tx=n.copy(z_echo)

            if sid[key] == 300:
                if use_long_pulse == False:
                    # ignore long pulse
                    continue
                # if long pulse, then take the next long pulse
                next_key = sidkeys[keyi+3]
                try:
                    z_echo1 = d_il.read_vector_1d(next_key, 10000, channel).astype("c8", casting="unsafe", copy=False) - z_dc
                except:
                    traceback.print_exc()
                    print("couldn't read echo")
                    continue
                    
                
            elif sid[key] == sid[sidkeys[keyi+1]]:
                # if first AC, subtract next one
                next_key = sidkeys[keyi+1]
                try:                
                    z_echo1 = d_il.read_vector_1d(next_key, 10000, channel).astype("c8", casting="unsafe", copy=False) - z_dc
                except:
                    traceback.print_exc()
                    continue
                    

            elif sid[key] == sid[sidkeys[keyi-1]]:
                # if second AC, subtract previous one.
                next_key = sidkeys[keyi-1]
                try:
                    z_echo1 = d_il.read_vector_1d(next_key, 10000, channel).astype("c8", casting="unsafe", copy=False) - z_dc
                except:
                    traceback.print_exc()
                    continue

                if debug_gc_rem:
                    plt.plot(zd.real+2000)
                    plt.plot(zd.imag+2000)
                    plt.plot(z_echo.real)
                    plt.plot(z_echo.imag)            
                    plt.title(sid[key])            
                    plt.show()


            noise0=tmm[sid[key]]["noise0"]
            noise1=tmm[sid[key]]["noise1"]
            last_echo=tmm[sid[key]]["last_echo"]
            tx0=tmm[sid[key]]["tx0"]
            tx1=tmm[sid[key]]["tx1"]
            gc=tmm[sid[key]]["gc"]
            e_gc=tmm[sid[key]]["e_gc"]



            # filter noise injection.
            if blank_impulses:
                # the calibration takes its windows as recorded: blanking
                # cannot treat the background and the brighter injection
                # window alike, and impulses left in both cancel in alpha
                z_noise=d_raw.read_vector_1d(key, 10000, channel).astype("c8", casting="unsafe", copy=False) - z_dc
            else:
                z_noise=n.copy(z_echo)
            z_noise=lpf.lpf(z_noise)

            # the dc offset changes. The mean of this pulse's background window
            # is kept, together with the mean power and mean of each window, so
            # the offset can be removed after the loop using all the pulses.
            z_dc_noise=n.mean(z_noise[(last_echo-500):last_echo])
            z_dc_samples.append(z_dc_noise)
            if noise_dc == "window":
                bg_samples.append( n.mean(n.abs(z_noise[(last_echo-500):last_echo]-z_dc_noise)**2.0) )
                bg_plus_inj_samples.append( n.mean(n.abs(z_noise[(noise0):noise1]-z_dc_noise)**2.0) )
            else:
                bg_samples.append( (n.mean(n.abs(z_noise[(last_echo-500):last_echo])**2.0), z_dc_noise) )
                bg_plus_inj_samples.append( (n.mean(n.abs(z_noise[(noise0):noise1])**2.0), n.mean(z_noise[(noise0):noise1])) )

            z_tx[0:tx0]=0.0
            z_tx[tx1:10000]=0.0

            # normalize tx pwr
            z_tx=z_tx/n.sqrt(n.sum(n.real(z_tx*n.conj(z_tx))))
            z_echo[last_echo:10000]=0.0
            z_echo1[last_echo:10000]=0.0

            z_echo[0:gc]=0.0
            z_echo1[0:gc]=0.0


            if False:
                plt.subplot(121)
                plt.plot(z_tx.real)
                plt.plot(z_tx.imag)
                plt.subplot(122)
                plt.plot(z_echo.real)
                plt.plot(z_echo.imag)
                plt.show()

            if False:
                # testing notching of frequencies.
                ZE=fft(z_echo)
                ZE1=fft(z_echo1)
                z_fftfreq=n.fft.fftfreq(len(z_echo),d=1/sr)

                if False:
                    plt.plot(n.fft.fftshift(z_fftfreq),n.fft.fftshift(10.0*n.log10(n.abs(ZE))**2.0))
                    plt.show()

                for freq_range in notch_freq_range:
                    fridx0=n.argmin(n.abs(z_fftfreq-freq_range[0]))
                    fridx1=n.argmin(n.abs(z_fftfreq-freq_range[1]))
                    noise_std=n.sqrt(0.25*(n.mean(n.abs(ZE[(fridx0-200):(fridx0-100)])**2.0)+n.mean(n.abs(ZE[(fridx1+100):(fridx1+200)])**2.0)+n.mean(n.abs(ZE1[(fridx0-200):(fridx0-100)])**2.0)+n.mean(n.abs(ZE1[(fridx1+100):(fridx1+200)])**2.0)))
                    nrand=fridx1-fridx0

                    ZE[fridx0:fridx1]=noise_std*(n.random.randn(nrand)+n.random.randn(nrand)*1j)/n.sqrt(2.0)
                    ZE1[fridx0:fridx1]=noise_std*(n.random.randn(nrand)+n.random.randn(nrand)*1j)/n.sqrt(2.0)

                if False:
                    plt.plot(n.fft.fftshift(z_fftfreq),n.fft.fftshift(10.0*n.log10(n.abs(ZE))**2.0))                
                    plt.show()

                z_echo=ifft(ZE)
                z_echo1=ifft(ZE1)

                

            # calculate power spectrum after notch
            Z=n.fft.fftshift(fft(spec_window*z_echo[(last_echo-fft_len):(last_echo)]))
            pwr_spec+=n.real(Z*n.conj(Z))
            n_pwr_spec+=1.0

            z_echo=lpf.lpf(z_echo)
            z_echo1=lpf.lpf(z_echo1)

            zd=z_echo-z_echo1

            zd[0:gc]=n.nan
            z_echo[0:gc]=n.nan
            z_echo[last_echo:10000]=n.nan
            zd[last_echo:10000]=n.nan

            # satellite templates for this pulse. The plain and the ground
            # clutter subtracted lagged products share one design matrix, and
            # the latter carry the echoes of the pulse this one is differenced
            # against, so an echo detected only there gets a column here too;
            # its amplitude in the plain products then comes out near zero.
            sat_here=sat_det.get(key,[])
            if len(sat_here)==0 and next_key is not None:
                sat_here=sat_det.get(next_key,[])
            sat_templates=[]
            sat_tx=sat_avg_tx[sid[key]][0] if sid[key] in sat_avg_tx else z_tx
            for sat_delay,sat_doppler in sat_here:
                sat_templates.append(satcol.satellite_template(sat_tx,sat_delay,sat_doppler,lpf,gc,last_echo,sr=sr))
                sat_cols.append((key,sat_delay,sat_doppler))
            sat_col0=len(sat_cols)-len(sat_templates)
            t1=time.time()
            read_time=t1-t0
            t0=time.time()
            for li in range(n_lags):
                for ai in range(lag_avg):
                    amb=decim.decimate(z_tx[0:(len(z_tx)-lags[li+ai])]*n.conj(z_tx[lags[li+ai]:len(z_tx)]))

                    # gc removal by the T. Turunen subtraction of two pulses with the same code, transmitted in
                    # close proximity to one another.
                    measg=decim.decimate(zd[0:(len(z_echo)-lags[li+ai])]*n.conj(zd[lags[li+ai]:len(z_echo)]))

                    # no gc removal
                    mease=decim.decimate(z_echo[0:(len(z_echo)-lags[li+ai])]*n.conj(z_echo[lags[li+ai]:len(z_echo)]))

                    #TM=amb[idxms[li]]
                    # add a column of ones to allow an additional noise process that is independent of range
                    O=n.ones(m1,dtype=n.complex64)
                    O.shape=(m1,1)
                    TM=n.hstack([amb[idxms[li]],O])
                    TM=sparse.csc_matrix(TM[m0:m1,:])

                    mgs[li].append(measg[m0:m1])
                    mes[li].append(mease[m0:m1])
                    A[li].append(TM)
                    S[li].append([(sat_col0+si,satcol.template_lag_column(st,lags[li+ai],decim,m0,m1)) for si,st in enumerate(sat_templates)])
            t1=time.time()
            ambiguity_time=t1-t0
            print("prep %d/%d ambiguity time %1.2f read time %1.2f (s)"%(keyi,n_pulses,ambiguity_time,read_time))            

        for code in sorted(unknown_codes.keys()):
            print("pulse code %d is not in the timing table, skipped %d pulses"%(code,unknown_codes[code]))
        for code in sorted(long_codes.keys()):
            print("pulse code %d reads beyond the inversion's 10000 samples, skipped %d pulses"%(code,long_codes[code]))

        # append the satellite columns, now that their number is known. Each
        # column is nonzero only in the rows of its own pulse.
        n_sat=len(sat_cols)
        if n_sat > 0:
            print("%d satellite columns"%(n_sat))
            for li in range(n_lags):
                for pi in range(len(A[li])):
                    cols=S[li][pi]
                    if len(cols) > 0:
                        sdata=n.concatenate([col for _,col in cols])
                        srow=n.tile(n.arange(n_meas,dtype=int),len(cols))
                        scol=n.repeat(n.array([ci for ci,_ in cols],dtype=int),n_meas)
                        SM=sparse.csc_matrix((sdata,(srow,scol)),shape=(n_meas,n_sat),dtype=n.complex64)
                    else:
                        SM=sparse.csc_matrix((n_meas,n_sat),dtype=n.complex64)
                    A[li][pi]=sparse.hstack([A[li][pi],SM],format="csc")
        sat_amp_e=n.full([n_sat,n_lags],n.nan,dtype=n.complex64)
        sat_amp_g=n.full([n_sat,n_lags],n.nan,dtype=n.complex64)
        sat_amp_var=n.full([n_sat,n_lags],n.nan,dtype=n.float32)

        acfs_g=n.zeros([rmax,n_lags],dtype=n.complex64)
        acfs_e=n.zeros([rmax,n_lags],dtype=n.complex64)
        
        # store noise autocorrelation function
        noise_e=n.zeros(n_lags,dtype=n.complex64)
        noise_g=n.zeros(n_lags,dtype=n.complex64)
        
        acfs_g[:,:]=n.nan
        acfs_e[:,:]=n.nan    

        acfs_var=n.zeros([rmax,n_lags],dtype=n.float32)
        acfs_var[:,:]=n.nan

        if noise_dc != "window" and len(z_dc_samples) > 0:
            # one DC offset for the period, the median over all its pulses of the
            # background window mean. A single window's mean is an estimate from
            # only N = 2*1.2*pass_band*500 us independent samples, 21.6 at 18 kHz,
            # and subtracting it per pulse biased T_sys low (memo 7): it removed
            # 1/N of the background power and added P_bg/N to the injection
            # window. Estimated over the period the error is negligible. Each
            # window's power about it follows from its mean power and mean,
            # mean|x-d|^2 = mean|x|^2 - 2 Re(d* mean x) + |d|^2.
            z_dc_period=n.median(n.real(z_dc_samples))+1j*n.median(n.imag(z_dc_samples))
            about_dc=lambda s: s[0]-2.0*n.real(n.conj(z_dc_period)*s[1])+n.abs(z_dc_period)**2.0
            bg_samples=[about_dc(s) for s in bg_samples]
            bg_plus_inj_samples=[about_dc(s) for s in bg_plus_inj_samples]
        noise=n.median(bg_samples)
        alpha=(n.median(bg_plus_inj_samples)-n.median(bg_samples))/T_injection
        T_sys=noise/alpha

        for li in range(n_lags):
            print(li)
            if len(A[li]) < 16:
                print("not enough measurements. skipping")
                continue
            else:
                n_good_estimates+=1
            
            
            AA=sparse.vstack(A[li])
            #print(AA.shape)
            mm_g=n.concatenate(mgs[li])
            mm_e=n.concatenate(mes[li])
            sigma_lp_est=n.zeros(len(mm_g))
            sigma_lp_est[:]=1.0
            
            n_ipp=0
            # remove outliers and estimate standard deviation 
            if True:
                print("ratio test")
                # the fourth moment term below was measured on 20 consecutive
                # files, 8602 range-lag cells: pooled and aligned to twice the
                # ACF phase, the pseudo variance is +0.024 +- 0.003 of the
                # variance, 8.4 sigma, and it grows with signal exactly as
                # R(tau)^2 requires: +0.004 at weak cells, +0.157 at the
                # strongest tenth. So the doubt below is correct, and the effect
                # is confined to the bright gates, where it makes the error
                # ellipse eccentric by about 17 per cent. Treating it properly
                # needs a widely linear least squares here and a non-circular
                # weighting in fit_acf; see memo 6. Not done: measured, bounded,
                # and judged not worth the rewrite for now.
                #
                # tbd: estimate the fourth moments for lagged products
                # 
                # <(m_t m_{t+\tau}^*) (m_t^* m_{t+\tau})>
                # but also for this one:
                # <(m_t m_{t+\tau}^*) (m_t m_{t+\tau}^*)>
                # as it might not be zero when snr is high!!!
                # this would require doing the least-squares with
                # a slightly different method
                # 
                mm_gm=n.copy(mm_g)
                mm_em=n.copy(mm_e)
                n_ipp=int(len(mm_gm)/n_meas)
                mm_gm.shape=(n_ipp,n_meas)
                mm_em.shape=(n_ipp,n_meas)

                
                sigma_lp_est=n.sqrt(n.percentile(n.abs(mm_em[:,:])**2.0,34,axis=0)*2.0)
                sigma_lp_est_g=n.sqrt(n.percentile(n.abs(mm_gm[:,:])**2.0,34,axis=0)*2.0)            

                ratio_test=n.abs(mm_em)/sigma_lp_est
                ratio_test_g=n.abs(mm_gm)/sigma_lp_est_g

                localized_sigma=n.abs(n.copy(mm_em))**2.0
                if precise_weights:
                    # the lagged products are complex64, so this running mean
                    # would otherwise be a single precision FFT, whose roundoff
                    # (~1e-7 of the largest value) swamps every pulse's weight
                    # in the rows a bright echo reaches, clean pulses included
                    localized_sigma=localized_sigma.astype(n.float64)
                wf=n.repeat(1/10,10)
                WF=fft(wf,localized_sigma.shape[0])
                for ri in range(mm_em.shape[1]):
                    # we need to wrap around, to avoid too low values.
                    localized_sigma[:,ri]=n.roll(n.sqrt(ifft(WF*fft(localized_sigma[:,ri])).real),-5)

                # make sure we don't have a division by zero
                msig=n.nanmedian(localized_sigma)
                if msig<0:
                    msig=1.0
                localized_sigma[localized_sigma<msig]=msig


                if False:
                    plt.pcolormesh(localized_sigma.T)
                    plt.colorbar()
                    plt.show()

                    plt.pcolormesh(ratio_test.T)
                    plt.colorbar()
                    plt.show()
                    plt.pcolormesh(ratio_test_g.T)
                    plt.colorbar()
                    plt.show()

                debug_outlier_test=False
                if debug_outlier_test:
                    plt.pcolormesh(mm_em.real.T)
                    plt.colorbar()
                    plt.show()


                if outlier_rejection:
                    # is this threshold too high?
                    # maybe 6-7 might still be possible.
                    mm_em[ratio_test > 10]=n.nan
                    mm_gm[ratio_test_g > 10]=n.nan

                    # these will be shit no matter what
                    mm_em[localized_sigma > 100*msig]=n.nan
                    mm_gm[localized_sigma > 100*msig]=n.nan

                ok_count+=n.sum((n.isnan(mm_em)!=True)*(n.isnan(mm_gm)!=True),axis=0)
                meas_count+=n_ipp
                
                if debug_outlier_test:            
                    plt.pcolormesh(mm_em.real.T)
                    plt.colorbar()
                    plt.show()

                    plt.pcolormesh(localized_sigma.T)
                    plt.colorbar()
                    plt.show()

                sigma_lp_est=localized_sigma
                sigma_lp_est.shape=(len(mm_g),)

                mm_gm.shape=(len(mm_g),)
                mm_em.shape=(len(mm_e),)
                mm_g=mm_gm
                mm_e=mm_em

            mm_g=mm_g/sigma_lp_est
            mm_e=mm_e/sigma_lp_est

            gidx = n.where( (n.isnan(mm_e)==False) & (n.isnan(mm_g)==False) & (n.isnan(sigma_lp_est) == False) )[0]
            print("%d/%d measurements good"%(len(gidx),len(mm_g)))

            # take outliers and bad measurements
            AA=AA[gidx,:]
            mm_g=mm_g[gidx]
            mm_e=mm_e[gidx]
            
            # at this point, we could add regularization to reduce range resolution on the top-side
            #
            # acf(rg[i])**rg[i]**2.0 = acf(rg[i+1])**rg[i+1]**2.0
            #
            # Something like this:
            # acf(rg[i]) - acf(rg[i+1])*(rg[i+1]**2.0/rg[i]**2.0) = 0
            #
            # n_rgs_this_lag = rmax-rmins[li]
            #

            
            srow=n.arange(len(gidx),dtype=int)
            scol=n.arange(len(gidx),dtype=int)
            sdata=1/sigma_lp_est[gidx]

            Sinv = sparse.csc_matrix( (sdata, (srow,scol)) ,shape=(len(gidx),len(gidx)))


            if len(gidx) < n_rg:
                print("not enough measurements. skipping")
                continue
                
            

            try:
                t0=time.time()
                # scale the design matrix by 1/sigma once, rather than applying
                # Sinv on each side of every product below. With B = Sinv A,
                # B^H B is the same Fisher information matrix A^H Sinv^2 A and
                # B^H mm the same A^H Sinv^2 m, because mm already carries one
                # factor of 1/sigma and Sinv is real and diagonal.
                B=Sinv.dot(AA)
                if n_sat > 0:
                    # a bright satellite's rows are weighted down by orders of
                    # magnitude, which leaves the columns of B on very different
                    # scales and the normal matrix too ill conditioned to invert
                    # in double precision. Scale every column to unit norm,
                    # solve, and scale back: the same estimate, computed stably.
                    col_norm=n.sqrt(n.asarray(B.multiply(B.conj()).real.sum(axis=0))).ravel()
                    col_norm[(col_norm==0) | ~n.isfinite(col_norm)]=1.0
                    if li in (0,n_lags//2):
                        BTu=n.conj(B.T)
                        print("lag %d normal matrix condition number: %.2e as it was, %.2e with the columns scaled"%(
                            li,n.linalg.cond(BTu.dot(B).toarray()),
                            n.linalg.cond((BTu.dot(B).toarray())/n.outer(col_norm,col_norm))))
                    B=B.dot(sparse.diags(1.0/col_norm))
                # (Sinv A)^H
                BT=n.conj(B.T)
                # A^H S^{-1} S^{-1} A (Fisher information matrix)
                ATA=BT.dot(B).toarray()

                # A^H \Sigma^{-1} m_g with ground clutter mitigation
                # note that 1/sigma is taken earlier when forming mm_g and mm_e
                # here we add a 1/sigma to get 1/sigma^2 on the diagonal of \Sigma^{-1}
                ATm_g=BT.dot(mm_g)
                # A^H \Sigma^{-1} m_e no ground clutter mitigation
                # note that 1/sigma is taken earlier when forming mm_g and mm_e
                ATm_e=BT.dot(mm_e)

                # error covariance, NaN for unknowns no measurement constrains
                Sigma,unconstrained=invert_normal(ATA)
                if len(unconstrained)>0:
                    print("lag %d: %d unknowns unconstrained (%s), left NaN"%(li,len(unconstrained),unconstrained[:8].tolist()))

                # ML estimate for ACF lag without ground clutter mitigation measures in place
                xhat_e=n.dot(Sigma,ATm_e)

                # ML estimate for ACF lag with ground clutter mitigation measures            
                xhat_g=n.dot(Sigma,ATm_g)

                if n_sat > 0:
                    # back from the scaled unknowns
                    xhat_e=xhat_e/col_norm
                    xhat_g=xhat_g/col_norm
                    Sigma=Sigma/n.outer(col_norm,col_norm)

                t1=time.time()
                t_simple=t1-t0        
                print("simple %1.2f"%(t_simple))
                # unknowns: the plasma gates, the background, then the satellites
                n_pl=rmax-rmins[li]
                acfs_e[ rmins[li]:rmax, li ]=xhat_e[0:n_pl]
                noise_e[li]=xhat_e[n_pl]
                acfs_g[ rmins[li]:rmax, li ]=xhat_g[0:n_pl]
                noise_g[li]=xhat_g[n_pl]

                acfs_var[ rmins[li]:rmax, li ] = n.diag(Sigma.real)[0:n_pl]
                if n_sat > 0:
                    sat_amp_e[:,li]=xhat_e[(n_pl+1):]
                    sat_amp_g[:,li]=xhat_g[(n_pl+1):]
                    sat_amp_var[:,li]=n.diag(Sigma.real)[(n_pl+1):]
            except:
                traceback.print_exc()
                print("something went wrong.")

        if n_good_estimates > 0:
            print("saving")
            if save_acf_images:
                # plot real part of acf
                acf_std=1.77*n.nanmedian(n.abs(acfs_e.real))
                plt.pcolormesh(mean_lags,rgs_km[0:rmax],acfs_e.real,vmin=-acf_std,vmax=2*acf_std)
                plt.xlabel(r"Lag ($\mu$s)")
                plt.ylabel("Range (km)")
                plt.colorbar()
                plt.title("%s T_sys=%1.0f K"%(stuffr.unix2datestr(i0/sr),T_sys))
                plt.tight_layout()
                plt.savefig("%s/lpi_%d/%s/lpi-%d.png"%(output_base,rg,channel,i0/sr))
                plt.close()
                plt.clf()

            #
            # tbd: determine if this could be done better with digital_metadata
            #
            ho=h5py.File("%s/lpi_%d/%s/lpi-%d.h5"%(output_base,rg,channel,i0/sr),"w")
            ho["acfs_g"]=acfs_g       # pulse to pulse ground clutter removal
            ho["acfs_e"]=acfs_e       # no ground clutter removal
            ho["noise_e"]=noise_e     # store estimated noise ACF
            ho["noise_g"]=noise_g     # store estimated noise ACF   
            ho["acfs_var"]=acfs_var   # variance of the acf estimate
            ho["rgs_km"]=rgs_km[0:rmax]
            ho["channel"]=channel
            ho["P_tx"]=avg_pwr/avg_pwr_n
            ho["lags"]=mean_lags/sr
            # the span this integration period covers. i0 is kept because
            # existing files and readers use it; t0 and t1 say where the period
            # ends, which a reader could otherwise only guess at from the start
            # of the following file.
            ho["i0"]=i0/sr
            ho["t0"]=i0/sr
            ho["t1"]=(i0 + avg_dur*idsr)/sr
            ho["T_sys"]=T_sys     # T_sys = alpha*noise_power
            ho["alpha"]=alpha     # This can scale power to T_sys (e.g., noise_power = T_sys/alpha)   T_sys * power/noise_pwr = T_pwr
            #
            ho["z_dc"]=n.median(z_dc_samples)
            if noise_dc != "window":
                # how the DC offset was removed for T_sys and alpha
                ho["noise_dc"]=noise_dc
            ho["pass_band"]=pass_band        # sort of important to store this, as this defines the low pass filter  
            ho["filter_len"]=filter_len      #
            # keep track of how many lagged products are rejected as bad as a function of time delay
            ho["retained_measurement_fraction"]=n.array(ok_count/meas_count,dtype=n.float32)
            ho["meas_delays_us"]=meas_delays_us
            ho["diagnostic_pwr_spec"]=pwr_spec/n_pwr_spec
            # written only when not at their defaults, so a default run writes
            # the same file as before
            if satellite_detections is not None or not outlier_rejection:
                ho["outlier_rejection"]=outlier_rejection
            if satellite_detections is not None:
                # one entry per satellite column: the pulse, and the delay
                # (samples, relative to the transmit template) and Doppler (Hz)
                # the column was built with
                ho["sat_keys"]=n.array([sc[0] for sc in sat_cols],dtype=n.int64)
                ho["sat_delay_samples"]=n.array([sc[1] for sc in sat_cols],dtype=n.float64)
                ho["sat_doppler_hz"]=n.array([sc[2] for sc in sat_cols],dtype=n.float64)
                # fitted complex amplitude of each column per lag, for the plain
                # and the ground clutter subtracted lagged products, and its
                # variance. The magnitude is the echo energy of the pulse and
                # should not change with lag.
                ho["sat_amp_e"]=sat_amp_e
                ho["sat_amp_g"]=sat_amp_g
                ho["sat_amp_var"]=sat_amp_var
                ho["satellite_template"]=satellite_template
                if len(sat_avg_tx) > 0:
                    # per code: pulses averaged, and their scatter about the average
                    ho["sat_template_codes"]=n.array(sorted(sat_avg_tx.keys()),dtype=int)
                    ho["sat_template_pulses"]=n.array([sat_avg_tx[c][1] for c in sorted(sat_avg_tx.keys())],dtype=int)
                    ho["sat_template_scatter"]=n.array([sat_avg_tx[c][2] for c in sorted(sat_avg_tx.keys())])
            ho.close()
        else:
            print("no estimates in this integration period")

if __name__ == "__main__":

    
    if True:
        datadir="/media/j/fee7388b-a51d-4e10-86e3-5cabb0e1bc13/isr/2021-12-01/usrp-rx0-r_20211201T230000_20211202T160100"
        lpi_files(dirname=datadir,
                  avg_dur=10,  # n seconds to average
                  channel="misa-l",
                  rg=30,       # how many microseconds is one range gate
                  min_tx_frac=0.2, # of the pulse can be missing
                  pass_band=0.018e6, # +/- 50 kHz 
                  filter_len=100,    # short filter, less problems with correlated noise, more problems with RFI
                  maximum_range_delay=7200,
                  save_acf_images=True,
                  lag_avg=1,
                  reanalyze=True)
        lpi_files(dirname=datadir,
                  avg_dur=10,  # n seconds to average
                  channel="zenith-l",
                  rg=30,       # how many microseconds is one range gate
                  min_tx_frac=0.2, # of the pulse can be missing
                  pass_band=0.018e6, # +/- 50 kHz 
                  filter_len=100,    # short filter, less problems with correlated noise, more problems with RFI
                  maximum_range_delay=7200,
                  save_acf_images=True,
                  lag_avg=1,
                  reanalyze=True)
        exit(0)

    
    if False:
        datadir="/media/j/4df2b77b-d2db-4dfa-8b39-7a6bece677ca/eclipse2024/usrp-rx0-r_20240407T100000_20240409T110000"
        lpi_files(dirname=datadir,
                  avg_dur=10,  # n seconds to average
                  channel="zenith-l",
                  rg=30,       # how many microseconds is one range gate
                  min_tx_frac=0.2, # of the pulse can be missing
                  pass_band=0.018e6, # +/- 50 kHz 
                  filter_len=100,    # short filter, less problems with correlated noise, more problems with RFI
                  maximum_range_delay=7200,
                  save_acf_images=True,
                  lag_avg=1,
                  reanalyze=False)
        exit(0)

        lpi_files(dirname=datadir,
                  avg_dur=10,  # n seconds to average
                  channel="misa-l",
                  rg=30,       # how many microseconds is one range gate
                  min_tx_frac=0.2, # of the pulse can be missing
                  pass_band=0.018e6, # +/- 50 kHz 
                  filter_len=100,    # short filter, less problems with correlated noise, more problems with RFI
                  maximum_range_delay=7200,
                  save_acf_images=True,
                  lag_avg=1,
                  reanalyze=True)
        
        
        datadir="/media/j/fee7388b-a51d-4e10-86e3-5cabb0e1bc13/isr/2023-09-24/usrp-rx0-r_20230924T200050_20230925T041059/"
        # bottom-side
        lpi_files(dirname=datadir,
                  avg_dur=10,  # n seconds to average
                  channel="zenith-l",
                  rg=30,       # how many microseconds is one range gate
                  output_prefix="%s/lpi_30"%(datadir),
                  min_tx_frac=0.5, # of the pulse can be missing
                  pass_band=0.018e6, # +/- 50 kHz 
                  filter_len=100,    # short filter, less problems with correlated noise, more problems with RFI
                  maximum_range_delay=6000,
                  save_acf_images=True,
                  lag_avg=1,
                  reanalyze=True)
        exit(0)
        
    if False:        
        datadir="/media/j/fee7388b-a51d-4e10-86e3-5cabb0e1bc13/isr/2023-09-05/usrp-rx0-r_20230905T214448_20230906T040054"    
        # Top-side
        lpi_files(dirname=datadir,
                  avg_dur=10,  # n seconds to average
                  channel="zenith-l",
                  rg=240,       # how many microseconds is one range gate
                  output_prefix="%s/lpi_240"%(datadir),
                  min_tx_frac=0.0, # of the pulse can be missing
                  pass_band=0.05e6, # +/- 50 kHz 
                  filter_len=10,    # short filter, less problems with correlated noise, more problems with RFI
                  maximum_range_delay=7200,
                  save_acf_images=True,
                  lag_avg=1,
                  reanalyze=False)
        exit(0)
    if True:
        dirname="/media/j/fee7388b-a51d-4e10-86e3-5cabb0e1bc13/isr/2023-09-28/usrp-rx0-r_20230928T211929_20230929T040533"
        lpi_files(dirname=dirname,
                  avg_dur=10,  # n seconds to average
                  channel="zenith-l",
                  rg=240,       # how many microseconds is one range gate
                  output_prefix="%s/lpi_240"%(dirname),
                  min_tx_frac=0.0, # how much of the pulse can be missing
                  filter_len=20,
                  pass_band=0.05e6,
                  maximum_range_delay=7200,
                  save_acf_images=False,
                  reanalyze=False,
                  lag_avg=1
                  )
        lpi_files(dirname=dirname,
                  avg_dur=10,  # n seconds to average
                  channel="zenith-l",
                  rg=120,       # how many microseconds is one range gate
                  output_prefix="%s/lpi_120"%(dirname),
                  min_tx_frac=0.0, # how much of the pulse can be missing
                  filter_len=20,
                  pass_band=0.05e6,
                  maximum_range_delay=7200,
                  save_acf_images=False,
                  reanalyze=False,
                  lag_avg=1
                  )
        lpi_files(dirname=dirname,
                  avg_dur=10,  # n seconds to average
                  channel="zenith-l",
                  rg=60,       # how many microseconds is one range gate
                  output_prefix="%s/lpi_60"%(dirname),
                  min_tx_frac=0.0, # how much of the pulse can be missing
                  filter_len=20,
                  pass_band=0.05e6,
                  maximum_range_delay=7200,
                  save_acf_images=False,
                  reanalyze=False,
                  lag_avg=1
                  )
        lpi_files(dirname=dirname,
                  avg_dur=10,  # n seconds to average
                  channel="zenith-l",
                  rg=30,       # how many microseconds is one range gate
                  output_prefix="%s/lpi_30"%(dirname),
                  min_tx_frac=0.5, # how much of the pulse can be missing
                  filter_len=20,
                  pass_band=0.05e6,
                  maximum_range_delay=7200,
                  save_acf_images=False,
                  reanalyze=False,
                  lag_avg=1
                  )

        
        exit(0)



    if False:
        dirname="/media/j/fee7388b-a51d-4e10-86e3-5cabb0e1bc13/isr/2023-09-24/usrp-rx0-r_20230924T200050_20230925T041059"
        lpi_files(dirname=dirname,
                  avg_dur=10,  # n seconds to average
                  channel="zenith-l",
                  rg=240,       # how many microseconds is one range gate
                  output_prefix="%s/lpi_240"%(dirname),
                  min_tx_frac=0.0, # how much of the pulse can be missing
                  filter_len=20,
                  pass_band=0.05e6,
                  maximum_range_delay=7200,
                  save_acf_images=True,
                  reanalyze=False,
                  lag_avg=1
                  )
        exit(0)        
        # E-region analysis
        # newly acquired fact: spectrally wide ambiguity functions mix more with out of band interference.
        # lower range resolutions works better in the presence of noise.
        # 30 microsecond gating is worse than 60 us or 120 us gating!
        lpi_files(dirname=dirname,
                  avg_dur=10,  # n seconds to average
                  channel="zenith-l",
                  rg=30,       # how many microseconds is one range gate
                  output_prefix="lpi_2023-09-24_30",
                  min_tx_frac=0.5, # how much of the pulse can be missing
                  filter_len=10,
                  pass_band=0.1e6,
                  maximum_range_delay=5000,
                  save_acf_images=False,
                  reanalyze=False,
                  lag_avg=3
                  )
        exit(0)
    # F1-region analysis
    lpi_files(dirname=dirname,
              avg_dur=10,  # n seconds to average
              channel="zenith-l",
              rg=60,       # how many microseconds is one range gate
              output_prefix="lpi_2023-09-24_60",
              min_tx_frac=0.3, # how much of the pulse can be missing
              filter_len=10,
              pass_band=0.1e6,
              maximum_range_delay=7200,
              save_acf_images=False,
              reanalyze=True
              )
    exit(0)
    
    if True:
        # F-region analysis
        lpi_files(dirname=datadir,
                  avg_dur=10,  # n seconds to average
                  channel="zenith-l",
                  rg=120,       # how many microseconds is one range gate
                  output_prefix="lpi_120",
                  min_tx_frac=0.0, # of the pulse can be missing
                  pass_band=0.1e6, # +/- 100 kHz 
                  filter_len=10,    # short filter, less problems with correlated noise, more problems with RFI
                  maximum_range_delay=7200,
                  save_acf_images=False,                  
                  reanalyze=True)
        
    
    if True:
        # F1-region analysis
        lpi_files(dirname=datadir,
                  avg_dur=10,  # n seconds to average
                  channel="zenith-l",
                  rg=60,       # how many microseconds is one range gate
                  output_prefix="lpi_60",
                  min_tx_frac=0.3, # how much of the pulse can be missing
                  filter_len=10,
                  pass_band=0.1e6,
                  maximum_range_delay=7200,
                  save_acf_images=False,
                  reanalyze=True
                  )

    if True:
        # E-region analysis
        # newly acquired fact: spectrally wide ambiguity functions mix more with out of band interference.
        # lower range resolutions works better in the presence of noise.
        # 30 microsecond gating is worse than 60 us or 120 us gating!
        lpi_files(dirname=datadir,
                  avg_dur=10,  # n seconds to average
                  channel="zenith-l",
                  rg=30,       # how many microseconds is one range gate
                  output_prefix="lpi_30",
                  min_tx_frac=0.5, # how much of the pulse can be missing
                  filter_len=10,
                  pass_band=0.1e6,
                  maximum_range_delay=4000,
                  save_acf_images=False,
                  reanalyze=True
                  )



