'''
##################
### _RDI_DEFS  ###
##################

Defining some quantities used in the rdi_adcp.py module.
(stuff that is too bulky to put in the main module)

_VARDICT_RDI
-----------

Defining the variables which will be read from the WinADCP matfiles.

- nm_out: The name that will be used once loaded into python (e.g.
            we will use "u" for "SerEmmpersec").
- fact_in_to_out: Scale factor used when transforming from WinADCP to
                    python (e.g. 1e-1 to convert from mm/s to cm/s).
- unit_out: The unit of used in the python processing (e.g. "cm/s").
- desc: A quick description of the field (e.g. 'Eastward velocity').
- dim: The dimantions of the quantity (1 for time, 2 for depth-time).

_PROC_STR_DICT
-------------

Strings to output when reading the procesing history.

_OUTVARS_RDI
-----------

Lists of variables to exports when using the RdiObj export functions.

_EXAMPLE_STR
------------

Strings to output to show processing workflow examples using the 
rdi_adcp.example() function-

_EXCEP_STR
----------

Various Exception error message strings.

_WARN_STR
---------

Various Warning message strings.  

_NETCDF_DICT
-------------

Contains conventions used when translating data to netCDF.

'''

# _VARDICT_RDI
# Dictionary of names of variables output by the RDI software.
# Also contains the expected dimensions of this variable, 
# a human readable description, the units we use for this variable,
# and the scale factor applied to convert to this unit.  

_vardict_rdi = {
'SerEmmpersec' : {'nm_out':'u', 'unit_out':'cm/s', 'fact_in_to_out':1e-1,
                    'desc':'Eastward velocity', 'dim':2},
'SerNmmpersec' : {'nm_out':'v', 'unit_out':'cm/s', 'fact_in_to_out':1e-1,
                    'desc':'Northward velocity', 'dim':2},
'SerVmmpersec' : {'nm_out':'w', 'unit_out':'cm/s', 'fact_in_to_out':1e-1,
                    'desc':'Upward velocity', 'dim':2},
'SerErmmpersec' : {'nm_out':'errvel', 'unit_out':'cm/s', 'fact_in_to_out':1e-1,
                    'desc':'Error velocity', 'dim':2},

'AnDepthmm' : {'nm_out':'tdepth', 'unit_out':'m', 'fact_in_to_out':1e-3,
                'desc':'Transducer depth', 'dim':1},
'AnT100thDeg' : {'nm_out':'T', 'unit_out':'deg C', 'fact_in_to_out':1e-2,
                'desc':'Temperature', 'dim':1},
'AnR100thDeg' : {'nm_out':'roll', 'unit_out':'degrees', 'fact_in_to_out':1e-2,
                'desc':'Roll', 'dim':1},
'AnP100thDeg' : {'nm_out':'pitch', 'unit_out':'degrees', 'fact_in_to_out':1e-2,
                'desc':'Pitch', 'dim':1},
'AnH100thDeg' : {'nm_out':'heading', 'unit_out':'degrees', 'fact_in_to_out':1e-2,
                'desc':'Heading', 'dim':1},
'SerDir10thDeg' : {'nm_out':'direction', 'unit_out':'degrees', 'fact_in_to_out':1e-1,
                    'desc':'Velocity direction', 'dim':2},
'AnOrienUP' : {'nm_out':'ori_ud', 'unit_out':'1/0 (y/n) binary', 'fact_in_to_out':1, 
                'desc':'Up/down orientation', 'dim':1},
'AnBatt' : {'nm_out':'batt', 'unit_out':'V?', 'fact_in_to_out':1, 
            'desc':'Battery level', 'dim':1},

'SerEA1cnt' : {'nm_out':'amp1', 'unit_out':'db', 'fact_in_to_out':1, 
                'desc':'Return amplitude, beam 1', 'dim':2},
'SerEA2cnt': {'nm_out':'amp2', 'unit_out':'db', 'fact_in_to_out':1, 
                'desc':'Return amplitude, beam 2', 'dim':2},
'SerEA3cnt' : {'nm_out':'amp3', 'unit_out':'db', 'fact_in_to_out':1, 
                'desc':'Return amplitude, beam 3', 'dim':2},
'SerEA4cnt' : {'nm_out':'amp4', 'unit_out':'db', 'fact_in_to_out':1, 
                'desc':'Return amplitude, beam 4', 'dim':2},
'SerEAAcnt' : {'nm_out':'amp_a', 'unit_out':'db', 'fact_in_to_out':1, 
                'desc':'Return amplitude, averaged over all functional beams', 
                'dim':2},

'SerPG1' : {'nm_out':'pg1', 'unit_out':'perc', 'fact_in_to_out':1, 
            'desc':'Percent good, beam 1', 'dim':2},
'SerPG2' : {'nm_out':'pg2', 'unit_out':'perc', 'fact_in_to_out':1, 
            'desc':'Percent good, beam 2', 'dim':2},
'SerPG3' : {'nm_out':'pg3', 'unit_out':'perc', 'fact_in_to_out':1, 
            'desc':'Percent good, beam 3', 'dim':2},
'SerPG4' : {'nm_out':'pg4', 'unit_out':'perc', 'fact_in_to_out':1,
            'desc':'Percent good, beam 4', 'dim':2},

'SerC1cnt' : {'nm_out':'ca1', 'unit_out':'counts, 0-255', 
            'fact_in_to_out':1, 'desc':'Correlation, beam 1', 'dim':2},
'SerC2cnt' : {'nm_out':'ca2', 'unit_out':'counts, 0-255', 
            'fact_in_to_out':1, 'desc':'Correlation, beam 2', 'dim':2},
'SerC3cnt' : {'nm_out':'ca3', 'unit_out':'counts, 0-255', 
            'fact_in_to_out':1, 'desc':'Correlation, beam 3', 'dim':2},
'SerC4cnt' : {'nm_out':'ca4', 'unit_out':'counts, 0-255', 
            'fact_in_to_out':1, 'desc':'Correlation, beam 4', 'dim':2},
'SerCAcnt' : {'nm_out':'ca_a', 'unit_out':'counts, 0-255', 'fact_in_to_out':1, 
            'desc':'Correlation, averaged over all functional beams', 'dim':2},
}

# PROC_STR_DICT_RDI
# Dictionary of strings associated with specific operations.
# Used in the documentation side of rdi_adcp, 
# can be desplayed using RdiObj.print_proc_steps() 

_proc_str_dict_rdi = {
    'shape_init' : ('Initial dataset shape: (Bins = %i, Profiles '
        '= %i) -> %i data points.'),
    'mask_init' : 'Initially masked entries: %i (%.2f %% of total).\n',
    'uvrot' : ' - Rotated currents by %.1f degrees CCW.',
    'magdec' : (' - Applied magnetic declination correction to horizontal currents:'
                '\n   Rotated current vectors CW by a mean angle %.2f deg.'),
    'shear' : ('Computed shear (shu, shv, s2) from first differences of'
                ' every 2nd bin.%s'),
    'rowrem' : (' - Removed %i rows (indices: %s, depths: %s).'
                '\n   %i -> %i rows (%.1f %% reduction).'),
    'rowrem_mask' : (' - Removed %i rows (indices: %s, depths: %s) based on mask'
                ' criterion\n  (rejecting rows with >%.1f %% masked entries).   '
                '\n   %i -> %i rows (%.1f %% reduction).'),
    'rowrem_surf' : (' - Removed %i rows (indices: %s, depths: %s) with median depth'
                     ' above the surface.'
                '\n   %i -> %i rows (%.1f %% reduction).'),
    'colrem' : (' - Removed %i columns.'
                '\n   %i -> %i rows (%.1f %% reduction).'),
    'colrem_wrongway' : (' - Removed %i columns where the ori_ud field indicates '
                  'that the instrument was facing %s.'
                '\n   %i -> %i rows (%.1f %% reduction).'),
    'depoffs' : (' - Applied %.1f m depth offset. \n   Median transducer depth '
                'changed from %.1f m to %.1f m.'),
    'shiptrem' : (' - Removed deployment (%i columns) and recovery (%i '
                  'columns) based on pressure record.'
                '\n   %i -> %i columns (%.1f %% reduction).'),
    'mask_umax' : (' - Flagged entries where u²+v² > (%.1f cm/s)².'
                '\n   -> Masked %i of %i entries (%.2f%% of previously unmasked).'),
    'mask_pg' : (' - Flagged entries where the "percent good" criterion (for '
                 '3- and/or 4-beam\n   solutions) is less than %.1f %%.'
                '\n   -> Masked %i of %i entries (%.2f%% of previously unmasked).'),
    'mask_amp_a' : (' - Flagged entries where the average backscatter amplitude'
                 ' is less than\n   %.1f db.'
                '\n   -> Masked %i of %i entries (%.2f%% of previously unmasked).'),
    'mask_errvel' : (' - Flagged entries where the error velocity is greater '
                 'than %.1f cm/s\n   (4-beam solutions only).'
                '\n   -> Masked %i of %i entries (%.2f%% of previously unmasked).'), 
    'mask_w' : (' - Flagged entries where the vertical velocity is greater '
                 'than %.1f cm/s.'
                '\n   -> Masked %i of %i entries (%.2f%% of previously unmasked).'), 
    'mask_ca' : (' - Flagged entries where the pulse-to-pulse correlation '
                 'amplitude in two or more beams is more than %.1f (0-255).'
                '\n   -> Masked %i of %i entries (%.2f%% of previously unmasked).'), 
    'mask_ca_mean' : (' - Flagged entries where the beam average pulse-to-pulse  '
                 'correlation amplitude is more than %.1f (0-255).'
                '\n   -> Masked %i of %i entries (%.2f%% of previously unmasked).'), 
    'mask_amp_jump' : (' - Flagged entries directly above amplitude jumps '
                 '(in any beam) exceeding\n   %.1f db.'
                '\n   -> Masked %i of %i entries (%.2f%% of previously unmasked).'), 
    'mask_amp_jump_abv' : (' - Flagged ALL entries above amplitude jumps '
                 '(in any beam) exceeding\n   %.1f db.'
                '\n   -> Masked %i of %i entries (%.2f%% of previously unmasked).'),
    'mask_tilt' : (' - Flagged profiles where tilt '
                 'exceeds %.1f degrees.'
                '\n   -> Masked %i profiles - %i of %i entries'
                ' (%.2f%% of previously unmasked).'),
    'heading_corr' : (' - Applied heading-dependent correction to compass heading and '
                'current vector angle.\n   Sinusiodal correction parameters:\n   '
                '[amplitude  phase  (offset)]  )   %s.\n   STD of sinusoidal fit: '
                '%.1f  degrees.\n   Average heading correction: %.1f degrees.'),
    'shape_final' : ('\nFinal dataset shape: (Bins = %i, Profiles '
        '= %i) -> %i data points.'),
    'mask_final' : 'Final masked entries: %i (%.2f %% of total).',

}

# _OUTVARS_RDI_FULL
# Dictionary containing lists of variables which we want to include 
# when exporting RdiObj objects after post-processing is done.
# (full: extensive list of parameters)
# (sparse: sparse list of parameters)

_outvars_rdi =  {
    'full': [
    'meta', 'binsize',  'u', 'v', 'w', 'errvel', 'tdepth', 't_datetime', 
    'T', 'roll', 'pitch', 'heading', 'direction', 'ori_ud', 'batt', 'amp1', 'amp2', 
    'amp3', 'amp4', 'amp_a', 'pg1', 'pg2', 'pg3', 'pg4', 'ca1', 'ca2', 'ca3', 'ca4', 
    'ca_a', 't_mpl', 'dep', 'dep_med', 'tilt', 'Ndep', 'Nt', 'shu', 'shv',
    'SN', 'lon', 'lat',
    ], 

    'sparse' : [
    'meta', 'binsize',  'u', 'v', 'tdepth', 'amp_a', 't_mpl', 't_datetime', 
    'dep', 'shu', 'shv', 'SN', 'lon', 'lat',
    ],
    }


# _EXAMPLE STR
# Dictionary containing strings with example usages of the module

_example_str =  {
    'extended': '''\
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
    ''',

    'minimal' : '''\

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
    '''}



# _EXCEP_STR
# Dictionary containing exception strings

_excep_str = {
    'magdec_lonlat': (
        'Correcting for magnetic declination (*correct_magdec()* method): \n'
        'Lat/lon need to be specified using the set_latlon() method.'),

    'shear_bins' : (
        'Cannot compute shear after we have removed individual rows (bins). ' 
        'Either compute shear *before* removing bins, or compute the shear '
        'yourself taking varying depth interval into account.'),

    'rrej_nodef' : (
        'Either *rows* or *masked_max* must be specified to '
         'RdiObj.reject_rows()'),

    'rrej_bothdef' : (
        'Both *rows* and *masked_max* specified to RdiObj.reject_rows(). Use '
        'only one.'),
}


# _WARNING_STR
# Dictionary containing warning strings

_warn_str = {
    'magdec_2024' : (
        'WMM2020 is only valid through 2024. Proceed with caution..'),

    'magdec_2010' : (
        'WMM2010 is only valid starting in 2010. Proceed with caution..'),                         
}
