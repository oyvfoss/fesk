'''

VARDICT_RDI

Defining the variables which will be read from the WinADCP matfiles.

- nm_out: The name that will be used once loaded into python (e.g.
            we will use "u" for "SerEmmpersec").
- fact_in_to_out: Scale factor used when transforming from WinADCP to 
                    python (e.g. 1e-1 to convert from mm/s to cm/s for u).
- unit_out: The unit of used in the python processing (e.g. "cm/s"). 
- desc: A quick description of the field (e.g. 'Eastward velocity').
- dim: The dimantions of the quantity (1 for time, 2 for depth-time).

PROC_STR_DICT

Strings to output when reading the procesing history.
'''

# VARDICT_RDI
# Dictionary of names of variables output by the RDI software.
# Also contains the expected dimensions of this variable, 
# a human readable description, the units we use for this variable,
# and the scale factor applied to convert to this unit.  

vardict_rdi = {
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
                    'desc':'Velocity direction', 'dim':1},
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

proc_str_dict_rdi = {
    'shape_init' : ('Initial dataset shape: (Bins = %i, Profiles '
        '= %i) -> %i data points.'),
    'mask_init' : 'Initially masked entries: %i (%.2f %% of total).\n',
    'uvrot' : ' - Rotated currents by %.1f degrees CCW.',
    'magdec' : (' - Applied magnetic declination correction to horizontal currents:'
                '\n   Rotated current vectors CW by a mean angle %.2f deg.'),
    'rowrem' : (' - Removed %i rows (indices: %s, depths: %s).'
                '\n   %i -> %i bins (%.1f %% reduction).'),
    'colrem' : (' - Removed %i columns.'
                '\n   %i -> %i bins (%.1f %% reduction).'),
    'depoffs' : (' - Applied %.1f m depth offset. \n   Median transducer depth '
                'changed from %.1f m to %.1f m.'),
    'shiptrem' : (' - Removed deployment (%i columns) and recovery (%i '
                  'columns) based on pressure record.'
                '\n   %i -> %i columns (%.1f %% reduction).'),
    'mask_umax' : (' - Flagged entries where u²+v² > (%.1f cm/s)².'
                '\n   -> Masked %i of %i entries (%.2f%% of previously unmasked).'),
    'mask_pg' : (' - Flagged entries where the "percent good" criterion (for '
                 '3- and/or 4-beam solutions) is less than %.1f %%.'
                '\n   -> Masked %i of %i entries (%.2f%% of previously unmasked).'),
    'mask_amp_a' : (' - Flagged entries where the average backscatter amplitude'
                 ' is less than %.1f db.'
                '\n   -> Masked %i of %i entries (%.2f%% of previously unmasked).'),
    'mask_errvel' : (' - Flagged entries where the error velocity is greater '
                 'than %.1f cm/s (4-beam solutions only).'
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
                 '(in any beam) exceeding %.1f db.'
                '\n   -> Masked %i of %i entries (%.2f%% of previously unmasked).'), 
    'mask_amp_jump_abv' : (' - Flagged ALL entries above amplitude jumps '
                 '(in any beam) exceeding %.1f db.'
                '\n   -> Masked %i of %i entries (%.2f%% of previously unmasked).'),
    'mask_tilt' : (' - Flagged profiles where tilt '
                 'exceeds %.1f degrees.'
                '\n   -> Masked %i profiles - %i of %i entries'
                ' (%.2f%% of previously unmasked).'),
}

# OUTVARS_RDI_FULL
# List of variables which we want to include when exporting RdiObj objects
# after post-processing is done.
# (_FULL: extensive list of parameters)
# (_SPARSE: sparse list of parameters)

outvars_rdi_full = [
    'meta', 'binsize',  'u', 'v', 'w', 'errvel', 'tdepth', 't_datetime', 
    'T', 'roll', 'pitch', 'heading', 'direction', 'ori_ud', 'batt', 'amp1', 'amp2', 
    'amp3', 'amp4', 'amp_a', 'pg1', 'pg2', 'pg3', 'pg4', 'ca1', 'ca2', 'ca3', 'ca4', 
    'ca_a', 't_mpl', 't_datetime', 'dep', 'dep_med', 'tilt', 'Ndep', 'Nt', 'shu', 'shv'
]

outvars_rdi_sparse = [
    'meta', 'binsize',  'u', 'v', 'tdepth', 'amp_a', 't_datetime', 'dep', 'shu', 'shv'
]
