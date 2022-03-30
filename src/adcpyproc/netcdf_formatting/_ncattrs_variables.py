# Dictionaries containing conventions for translating to netCDF



# _NETCDF_DICT_RDI
# For RDI ADCPs

_ncattrs_variables = {
    'u' :{'ncname': 'UCUR', 
          'attrs' : {
          'long_name': 'Eastward current velocity',
          'standard_name': 'eastward_sea_water_velocity',
          'coverage_content_type':  'physicalMeasurement',
          'processing_level': (
                'Known bad data has been replaced with null values'),
          'grid_mapping' : 'crs',
          'source' : 'mooring',
          },},

    'v' :{'ncname': 'UCUR',
          'attrs' : {
          'long_name': 'Northward current velocity',
          'standard_name': 'northward_sea_water_velocity',
          'coverage_content_type':  'physicalMeasurement', 
          'processing_level': (
                'Known bad data has been replaced with null values'),
          'grid_mapping' : 'crs',
          'source' : 'mooring',
          },},
          
    'dep': {'ncname': 'DEPTH', 
    'attrs' : {
            'long_name' : 'Time-varying depth of each bin',
            'standard_name' : 'depth',
            'ancillary_variables':  'CRS INSTRUMENT PLATFORM',
            'axis' : 'Z', 'units' : 'm', 'positive' : 'down',
            }}
}



