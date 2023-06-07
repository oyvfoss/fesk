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
          'source' : 'mooring',
          },},

      'v' :{'ncname': 'VCUR',
          'attrs' : {
          'long_name': 'Northward current velocity',
          'standard_name': 'northward_sea_water_velocity',
          'coverage_content_type':  'physicalMeasurement', 
          'processing_level': (
                'Known bad data has been replaced with null values'),
          'source' : 'mooring',
          },},
          
      'dep': {'ncname': 'DEPTH', 
            'attrs' : {
            'long_name' : 'Time-varying depth of each bin',
            'standard_name' : 'depth',
            'coverage_content_type':  'physicalMeasurement', 
            'axis' : 'Z', 'units' : 'm', 'positive' : 'down',
            }},
      'lon': {'ncname': 'LONGITUDE',
            'attrs' : {'long_name':'longitude', 'standard_name':'longitude',
              'units':'degrees_east'}},
      'lat': {'ncname': 'LATITUDE',
            'attrs' : {'long_name':'latitude', 'standard_name':'latitude',
              'units':'degrees_north'}},
}



