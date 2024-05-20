### EXAMPLE NETCDF ATTRIBUTE DICTIONARY ###
'''
When exporting your dataset to netCDF, modify 

Some mandatory keywords (lat/lon ranges, time period, etc are created
automatically from the data).
'''



ncattr_vars = {
    # These must be filled out for each individual dataset. 
    'title' : 'Ocean currents from Fram Expedition',
    'id': '?',
    'naming authority': 'Norwegian Polar Institute', # 
    'summary' : 'Currents ', # Summary of the dataset - analogous to a paper abstract

    'creator_type' : 
    'creator_institution' :
    'creator_name' : 
    'creator_email' : 
    'creator_url' :
    'institution' : 
    'publisher_name' : 
    'publisher_email' : 
    'publisher_url' : 
    'project' : 
    'platform' : 
    'activity_type' : 
    'operational_status' : 
    'acknowledgements' : 'Funded by XX agency. Jane Doe collected the data.'
    'ancillary_keywords' : 'In-Situ,Mooring,Current, ADCP',

    
    'area' : 'Arctic Ocean',#
    'platform': 'In Situ Ocean-based Platforms>MOORINGS>MOORINGS'


    # These are standard parameters that you typically don't have to change.

    'featureType' : 'timeSeriesProfile', 
    'keywords' : 'EARTH SCIENCE>OCEANS>OCEAN CIRCULATION>OCEAN CURRENTS,EARTH SCIENCE>OCEANS>OCEAN PRESSURE>WATER PRESSURE',
    'keywords_vocabulary' : 'NASA/GCMD Science Keywords 9.1.5',
    'standard_name_vocabulary' : 'CF Standard Name Table v77'


    'conventions': 'NPI-0.2, CF-1.7, NMDC-0.6, ACDD-1.3',
    'QC_indicator' : 'good_data',
    'history': '',
    'source'; '',
    'processing_level': '', 
    'keywords' : '',
    'license': 'CC-BY 4.0 (https://creativecommons.org/licenses/by/4.0/)',
    'platform_volabulary' : 'NASA/GCMD Platform Keywords 9.1.5',
}

