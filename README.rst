ADCPyProc
#########

Version 0.0.1: *Under development. Locally working code. Not extensively 
tested or QCed.* 

Overview
--------

Post-processing and basic analysis of ocean ADCP data using python.

(Currently only has functionality for RDI ADCP data - starting from .mat 
files output by WinADCP).

(At the moment: primarily functionality for processing, not much for
analysis/inspection. Will add more of this later. )

Creates a ``RdiObj`` object containing the data. This object can be easily 
manipulated with associated methods for cleaning, flagging, chopping and other 
useful post-processing operations. The processed data can then be exported to 
a python dictionary, a matlab *.mat* file, or a netCDF file (the latter is not 
implemented yet).

A key aspect of the algorithm is that it keeps track of which processing steps
have been applied, in what order, and with what effect on the dataset. This is 
useful in order to keep track of the processing applied to a dataset, and it can
be a useful reference when i.e. writing a Methods section of a scientific 
manuscript or report.

The code manipulates data in the form of *masked numpy arrays*, which use an 
auxilliary array of masks to keep track of which values are considered valid. Upon
export to matlab or netcdf formats, invalid values are converted to NaNs.
The "masking" and "flagging" are used somewhat interchangeably in the documentation. 

Dependencies
-------------

*Adcpyproc* is a Python package, and requires Python 3 (will not work on 2. 
and has currently only been tested on 3.8).

Beyond standard Python libraries such as *numpy*, *scipy*, *matplotlib*, etc, 
*adcpyproc* depends on:

- `pickle <https://docs.python.org/3/library/pickle.html>`_ (for saving python dictionaries).
- `geomag <https://pypi.org/project/geomag/>`_ (for computing time-and location-dependent magnetic declination based on the World Magnetic Model forcompass corrections).  
 
  - **Note**: *This dependency should be removed. Clunky application, and not ideal to rely on this small and poorly maintained package. In the future: Let the user supply magnetic declination* (easily obtained from WMM online calulators etc.)  

- `netCDF4 <https://unidata.github.io/netcdf4-python/>`_ for exporting to netCDF.
 
To process raw *.000* files output from RDI ADCPs, you need to use RDI's WinADCP
software.


- **Note**: It would probably be a good idea to change the base of the sorftware from 
  numpy arrays to xarray Datasets, as in the *SigPyProc* module.

Installing
----------

To install, you should be able to do something like this:

1. Obtain a copy of *adcpyproc*. Either: 
    - Clone the git repository into a suitable location.
    - Download the zip file and unpack into a suitable folder.

2. Navigate to the top folder of the pproject (the one containing ``setup.py`` ).
3. If you haven't already: Install *pip*. For conda systems: ``conda install pip``.
4. Install *adcpyproc* using ``pip install -e .``.
    - The ``-e`` flag should be used while the code is still in the initial
      development phase.
5. You should now be able to load the module using, e.g.:
    - ``import adcpyproc``
    - ``from adcpyproc import rdi_adcp``

Example usage
-------------

Minimal example
+++++++++++++++

An example of processing steps which may be considered a minimum for ADCP
processing.


::

    # Import adcpyproc
    from adcpyproc import rdi_adcp

    # Set the path to the .mat file output from WinADCP
    fn_mat = 'path/to/matfile/file_from_winadcp.mat'

    # Load data from a WinADCP-produced .mat-file into an RdiObj object
    d = rdi_adcp.RdiObj(fn_mat)           

    # Chop away deployment and recovery
    d.remove_ship_time()   

    # Reject bins where the median depth is above the surface.
    d.reject_surface_bins()     

    # Set lat and lon and correct current direction for magnetic declination.
    d.set_latlon(80, 30)                   
    d.correct_magdec()   

    # MASKING/FLAGGING BASED ON THRESHOLDS 
    d.mask_umax(thr_uvamp=100)        # Mask entries where u or v exceed 100 cm/s amplitude.
    d.mask_pg(thr_pg=75)              # Mask entries where percent good (PG1 + PG4) is below 75%.
    d.mask_amp_a(thr_amp=64)          # Mask entries where beam mean amplitude is less than 64 db.
    d.mask_errvel(thr_errvel=50)      # Mask entries where error velocity is greater than 50 cm/s.
    d.mask_ca(thr_cor=30)             # Mask entries where the mean beam correlation in two or  
                                      # more beams is below 45 counts.
    d.mask_tilt(thr_tilt=20)          # Mask entries where tile exceeds 20 degrees.

    # Print a string showing a summary of the dataset (time/depth ranges, mean velocities, etc).
    d.print_summary()                      

    # Print a processing history listing the individual steps applied to the dataset.   
    d.print_proc()                                              

    # Save as matfile with the typically most important parameters (t, depth, u, v, ..)
    d.to_matfile('test_fn.mat', sparse = True)       


Extended example
+++++++++++++++++


::

    # Import adcpyproc
    from adcpyproc import rdi_adcp

    # Set the path to the .mat file output from WinADCP
    fn_mat = 'path/to/matfile/file_from_winadcp.mat'

    # Load data from a WinADCP-produced .mat-file into an RdiObj object
    d = rdi_adcp.RdiObj(fn_mat)           

    # -------------------------------------------------------------------------------------

    # Print some system parameters (instrument configuration)                 
    d.print_system_info()                  

    # Chop away deployment and recovery
    d.remove_ship_time()    

    # Adjust transducer and bin depths 3.2 m *downwards*               
    d.apply_depth_offset(3.2)              

    # Reject bins where the median depth is above the surface.
    d.reject_surface_bins()                

    # Reject the two bins nearest to transducer(rows 0 and 1)). Will prompt y/n.
    d.reject_rows([0, 1])                     

    # Set lat and lon and correct current direction for magnetic declination.
    d.set_latlon(80, 30)                   
    d.correct_magdec()                     

    # -------------------------------------------------------------------------------------

    ## Masking based on criteria (apply the relevant ones and modify the criteria) ##
    ## (Masks will end up as NaNs when exporting to matlab)

    d.mask_umax(thr_uvamp=100)           # Mask entries where u or v exceed 100 cm/s amplitude.
    # MASK_SURF_SIDELOBE NOT APPLIED YET!
    #d.mask_surf_sidelobe() (X)          # Mask entries falling within the estimated range of 
                                         # sidelobe interference of the surface.
    d.mask_pg(thr_pg=75)                 # Mask entries where percent good (PG1 + PG4) is below 75%.
    d.mask_amp_a(thr_amp=64)             # Mask entries where beam mean amplitude is less than 64 db.
    d.mask_errvel(thr_errvel=50)         # Mask entries where error velocity is greater than 50 cm/s.
    d.mask_ca(thr_cor=30)                # Mask entries where the mean beam correlation in two or  
                                         # more beams is below 45 counts.
    d.mask_ca_mean(thr_cor=30)           # Mask entries where the mean beam correlation is below
                                         # 45 counts.
    d.mask_w(thr_w=30)                   # Mask entries where the mean vertical is below 30 cm/s.
    d.mask_amp_jump(max_amp_increase=30) # Masking entries where the beam amplitude of any beam has 
                                         # a jump of 30 db or more (masking *after* the jump).
    d.mask_amp_jump(max_amp_increase=30, # Same, but also masks all entries *above* such jumps.
                    mask_above=True)
    d.mask_tilt(thr_tilt=20)             # Mask entries where tile exceeds 20 degrees.
    
    # -------------------------------------------------------------------------------------
    
    # calculate vertical shear (s2, shu, shv)
    d.calculate_shear()                   

    # Reject rows with less than 50% valid (unmasked) entries.
    d.reject_rows(masked_max=50)      

    # Print a string showing a summary of the dataset (time/depth ranges, mean velocities, etc).
    d.print_summary()                      

    # Print a processing history listing the individual steps applied to the dataset.   
    d.print_proc()                          

    # Export the dataset to a python Bunch.
    b = d.to_dict()                       

    # Save as matfile with the typically most important parameters (t, depth, u, v, ..)
    d.to_matfile('test_fn.mat', sparse=True)       

    # Save as pickled python dictionary (all parameters).
    d.to_matfile('test_fn_full.mat')      

     # Save as pickled python dictionary (all parameters).
    d.to_pickle('test_fn.p')              
    
    # Save as netcdf file (all parameters). TO_NETCDF4 NOT APPLIED YET!
    #d.to_netcdf('test_fn.nc')            