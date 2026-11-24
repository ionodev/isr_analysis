import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "pipeline"))  # the pipeline modules
import outlier_lpi as olpi
import avg_range_doppler_spec as ards
import traceback

# the recordings to analyse, on the command line:
#     python3 analysis/analyze_dirs.py <data dir> [<data dir> ...]
dirs=sys.argv[1:]
for d in dirs:
    try:
        olpi.lpi_files(dirname=d,
                       avg_dur=10,  # n seconds to average
                       channel="zenith-l2",
                       rg=30,       # how many microseconds is one range gate
                       min_tx_frac=0.5, # of the pulse can be missing
                       pass_band=0.018e6, # +/- 50 kHz 
                       filter_len=100,    # short filter, less problems with correlated noise, more problems with RFI
                       maximum_range_delay=7200,
                       save_acf_images=True,
                       lag_avg=1,
                       reanalyze=False)
    except:
        traceback.print_exc()
    try:
        olpi.lpi_files(dirname=d,
                       avg_dur=10,  # n seconds to average
                       channel="zenith-l",
                       rg=30,       # how many microseconds is one range gate
                       min_tx_frac=0.5, # of the pulse can be missing
                       pass_band=0.018e6, # +/- 50 kHz 
                       filter_len=100,    # short filter, less problems with correlated noise, more problems with RFI
                       maximum_range_delay=7200,
                       save_acf_images=True,
                       lag_avg=1,
                       reanalyze=False)
    except:
        traceback.print_exc()
        
    try:
        olpi.lpi_files(dirname=d,
                       avg_dur=10,  # n seconds to average
                       channel="misa-l",
                       rg=30,       # how many microseconds is one range gate
                       min_tx_frac=0.5, # of the pulse can be missing
                       pass_band=0.018e6, # +/- 50 kHz 
                       filter_len=100,    # short filter, less problems with correlated noise, more problems with RFI
                       maximum_range_delay=7200,
                       save_acf_images=True,
                       lag_avg=1,
                       reanalyze=False)
    except:
        traceback.print_exc()
    
    try:
        ards.avg_range_doppler_spectra(dirname=d,
                                       channel="zenith-l",
                                       save_png=True,
                                       avg_dur=10,
                                       reanalyze=False
                                       )
    except:
        traceback.print_exc()
    try:
        ards.avg_range_doppler_spectra(dirname=d,
                                       channel="misa-l",
                                       save_png=True,
                                       avg_dur=10,
                                       reanalyze=False
                                       )
    except:
        traceback.print_exc()
    
    
