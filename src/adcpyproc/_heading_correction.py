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


def _apply_heading_correction(self, heading_true, heading_obs, 
            allow_offset = True, plot_fit = True, auto_accept = False):
    '''
    '''

    f, std_degs, params = fit_lookup_table(heading_true, heading_obs, 
            allow_offset = allow_offset)
    
    if plot_fit:

        degdiff = heading_true - heading_obs
        degdiff[degdiff>180] -= 360
        degdiff[degdiff<-180] += 360

        fig, ax = plt.subplots(figsize = (6, 3))
        ax.plot(heading_obs, degdiff, 'o', label = 'Table')
        xdeg = np.arange(0, 361)
        ax.plot(xdeg, f(xdeg), '-', label = 'Fit')
        plt.legend()
        plt.show()

    if auto_accept:
        accept = 'y'
    if not auto_accept:
        print('Fit parameters [amplitude  phase  (offset)]: ', params)
        print('Std of fit: %.1f degrees.'%std_degs)
        accept = input('Accept and apply heading correction? (y/n): ')
    
    if accept == 'y':
        heading_orig = self.heading
        heading_correction = f(heading_orig)
        heading_new = heading_orig + heading_correction
        heading_new[heading_new>360] -= 360
        heading_new[heading_new<0] += 360
        self.heading = heading_new

        direction_new =  self.direction + heading_correction
        direction_new[direction_new>360] -= 360
        direction_new[direction_new<0] += 360
        self.direction = direction_new

        self._rotate_uv(-heading_correction, write_to_proc_dict=False)
        self._record_proc_step('heading_corr', 
            (params, std_degs, heading_correction.mean()))


