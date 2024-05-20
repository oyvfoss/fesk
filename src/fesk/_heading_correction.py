#######################################################################
# HEADING_CORRECTIONS
#######################################################################
'''
Occasionally, we do not trust the instrument compass calibration and therefore
create a "lookup table" where we orient the instrument in various directions and
compare the instrument's reported heading with the known "true" heading (taken
from a compass external to the instrument).   
'''

from scipy.optimize import curve_fit
import numpy as np

def fit_lookup_table(heading_true, heading_obs, allow_offset = True):
    '''
    Least squares best fit of the observed difference between trua and observed
    heading to the function f = AMP * sin(deg_obs/180*pi + PHASE) + OFFSET.

    If allow_offset is *False*, we assume a zero offset and only fit two
    parameters.
    
    All parameters (including phase) are in degrees.

    Returns f, std_degs, where:

    f is the function

        f(heading)

    which is the best fit to approximate

        heading_true = heading_obs + f(heading_obs)

    std_degs is the standard deviation between 
    
        f(heading_obs) and heading_true - heading_obs
    '''

    degdiff = heading_true - heading_obs
    degdiff[degdiff>180] -= 360
    degdiff[degdiff<-180] += 360

    if allow_offset:
        def sinfunc(ddiff, amp, ph, offs):
            F = amp * np.sin(ddiff/180*np.pi + ph/180*np.pi) + offs
            return F
    else:
        def sinfunc(ddiff, amp, ph):
            F = amp * np.sin(ddiff/180*np.pi + ph/180*np.pi)
            return F

    best_fit = curve_fit(sinfunc, heading_obs, degdiff)

    f = lambda deg: sinfunc(deg, *best_fit[0])

    std_degs = np.std(degdiff - f(heading_obs))

    params = best_fit[0]

    return f, std_degs, params

