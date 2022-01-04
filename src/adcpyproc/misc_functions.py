#######################################################################
# MISC_FUNCTIONS
#######################################################################
'''
Various functions that are not specific to the classes defined in
adcpyproc:

- Find magnetic declination from lon/lat/time using geomag.
- Convert mpl time to matlab time
'''
import numpy as np
import datetime as dt

try:
    import geomag
except ModuleNotFoundError:
    raise Warning("Couldn't load module *geomag* -> won't be "
          "able to apply magnetic declination correction.")

from adcpyproc import __file__ as _initloc

# Get the path of the folder containing the .COF file (a little hacky..)
cofdir = _initloc.replace('__init__.py', 'WMM_COF_files/')

#######################################################################

def declination(dlat, dlon, time, h=0, model=None, explicit=True):
    """
    Wrapper for function in the *geomag* module by Christopher
    Weiss (https://pypi.org/project/geomag/).

    Used to calculate magnetic declination based on World
    Magnetic Model coefficients allowing different WMM models.

    Allows interation over several times (and lat/lons).
    (implentation is a bit clunky but works fine for limited datasets).

    Calculate magnetic declination in degrees
    dlat = Latitude in degrees (float or array)
    dlon = Longitude in degrees (float or array)
    h = Altitude in feet, default=0
    time = Date for computing declination (float or array)
    model = Magnetic model (WMM). Default behaviour is to select the
            model based on the time. Possible inputs are 'WMM2010',
            'WMM2015', 'WMM2020'.
    
    Returns: Magnetic declination in degrees (float or array depending 
             on the input)

    NOTE 1: Each WMM model is valid for five years. The three models
    included here are valid between 1-1-2010 and 31-12-2024. Caution
    is advised outside of this time range.
    NOTE 2: Discontinuities can occur when switching between models
    (e.g. at 1-1-2015 00:00). Consider using a single WMM model for
    a single dataset.  

    The World Magnetic Model is a representation of the Earth's magntic
    field developed by the US National Geophysical Data Center and the
    British Geological Survey. .COF files with coefficients were
    downloaded from www.ncei.noaa.gov.
    """

    if not model:
        if time.year > 2019:
            model = 'WMM2020'
            if time.year > 2024:
                raise Warning('WMM2020 is only valid through 2024.'
                              ' Proceed with caution..')
        elif time.year < 2015:
            model = 'WMM2010'
            if time.year < 2010:
                raise Warning('WMM2010 is only valid starting in 2010.'
                              ' Proceed with caution..')
        else:
            model = 'WMM2015'

    if explicit:
        print('Using %s to compute declination..' % model)
        
    wmm_filename='%s%s.COF'% (cofdir, model)

    __singleton__ = geomag.geomag.GeoMag(wmm_filename=wmm_filename)

    # For a single time entry 
    if not hasattr(time, '__iter__'): 
        time_date = dt.date(*time.timetuple()[:-6])
        magdec = __singleton__.GeoMag(dlat, dlon, h, time_date).dec
        
    # For many time entries
    else:
        magdec = np.array([])
        time_dates = [dt.date(*t_.timetuple()[:-6]) for t_ in time]

        N = len(time)

        # For a single lat/lon entry:
        if not hasattr(dlon, '__iter__'): 
            for nn, time_ in enumerate(time_dates):
                magdec_ = __singleton__.GeoMag(dlat, dlon, h, time_).dec
                magdec = np.append(magdec, magdec_)
                if explicit:
                    print('%.2f %%..\r'%(100*nn/N), end = '')

        # For multiple lat/lon entries:
        else:
            if time.shape != dlat.shape:
                raise Exception('magdec calculation: lat/lon must be constant' 
                                'or have the same shape as t!')

            for time_, dlon_, dlat_ in zip(time_dates, dlon, dlat):
                magdec_ = __singleton__.GeoMag(dlat_, dlon_, h, time_).dec
                magdec = np.append(magdec, magdec_)
                if explicit:
                    print('%.2f %%..\r'%(100*nn/N), end = '')
        
        print('...  done.')

    return magdec


#######################################################################


def t_python_to_mat(pytime):
    '''
    Convert matlab datenum (days) to Matplotlib dates (days).

    MATLAB base: 00-Jan-0000
    Matplotlib base: 01-Jan-1970
    '''

    mattime = pytime + 719529.0

    return mattime