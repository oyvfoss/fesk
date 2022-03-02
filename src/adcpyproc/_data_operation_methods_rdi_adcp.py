# Internal methods to manipulate/add data for the rdi_adcp.RdiObj class

import numpy as np
import scipy.io
import datetime as dt
from collections import OrderedDict
from adcpyproc import _rdi_defs

#######################################################################
# SETUP THE OBJECT FROM A .MAT FILE
#######################################################################


def _setup(self, matfn):
    """
    Load matfile (typically processed from binary file using
    WinADCP) and transfer variables to RdiObject with conventions
    (see _rdi_defs.py).

    Convention: Depth decreases with increasing 0th index.
    """

    print('Initializing postprocessing by creading RdiObj object from'
            ' this file:\n %s..\n' % matfn)
    m = scipy.io.matlab.loadmat(matfn)

    print('Writing some metadata and creating some useful structures..')
    self.meta = (
        'Created on %s using *pyrdi_adcp* from the source file %s.' % (
            dt.datetime.now().strftime('%d.%m.%y'), matfn))
    self.binsize = m['RDIBinSize']
    self.var_desc, self.units = {}, {}

    # Documentation
    self._hist_dict = {'procdate': dt.datetime.now().strftime('%d.%m.%Y'),
                        'source_matfile': matfn,
                        'source_rawfile': m['RDIFileName'][0]}
    self._proc_dict = OrderedDict()
    self._extract_config(m)

    print('Loading data to numpy arrays and converting to '
            'preferred names/units..')
    self._load_data_from_matfile(m, _rdi_defs._vardict_rdi,
                                 explicit=False)

    # Store dimensions
    self.Ndep, self.Nt = self.u.shape

    # Record keeping: Initial data details
    self._record_proc_step('shape_init', (
        self.Nt, self.Ndep, self.Nt * self.Ndep))
    self._record_proc_step('mask_init', (
        self.u.mask.sum(), 100 * self.u.mask.sum()
        / (self.Nt * self.Ndep)))


    # Calculate time
    self._calculate_time(m)

    # Check orientation and remove any profiles where the instrument
    # faced the wrong way
    self._get_orientation()

    print('Calculating bin depths..')
    self._calculate_bin_depths(m)

    self.dep_med = np.median((self.dep), axis=1)


    # Flipping depth dimension if orientation is downward
    self._flip_to_ori()

    print('Calculating tilt from pitch and roll..')
    self._calculate_tilt()

    # Create summary of the deployment
    self._update_summary()

    print('\nInitialization complete. Now for the processing.'
            '\nRun rdi_adcp.example() for an example.')


#######################################################################
# FUNCTIONS FOR MANIPULATING DATA
#######################################################################

def _apply_mask(self, mask_arr, include_errvel=True, return_ns=False):
    """
    Mask the True indices of the *mask_arr* array.
    Velocity components (u, v, w) are always masked.

    include_errvel: Specify whether the error velocity should also
    be masked.

    if return_ns = True:
        Returns n, n_umasked_before, n_umasked_after, where:

    n: Total number of data points (masked and unmasked) at the
        time of the operation.
    n_umasked_before: Number of unmasked points after the
                        operation.
    n_umasked_after: Number of unmasked points after the
                        operation.
    """
    n = np.sum(np.ones(self.u.shape))
    n_umasked_before = np.sum(~self.u.mask)
    if include_errvel:
        add_field = [
            'errvel']
    else:
        add_field = []
    if hasattr(self, 'shu'):
        add_field += ['shu', 'shv', 's2']
        raise Warning('Mask applied to shear. In general it is better' 
                        'to complete the masking of the velocity components'
                        ' *before* computing shear!')
    for nm in [
        'u', 'v', 'w'] + add_field:
        v0 = getattr(self, nm)
        v0.mask = v0.mask + mask_arr
        setattr(self, nm, v0)
    else:
        n_umasked_after = np.sum(~self.u.mask)
        return (
            n, n_umasked_before, n_umasked_after)


#######################################################################

def _rotate_uv(self, ang_deg, write_to_proc_dict=True):
    """
    Rotates the horizontal currents clockwise.

    ang_deg: Angle (degrees) - can be 1D iterable with dimension 
                                (time).

    Rotates u and v (and shu and shv if applicable), and applies 
    an offset to *direction*.
    """
    uvc_orig = self.u + 1j * self.v
    uvc_rot = uvc_orig * np.exp(1j * ang_deg*np.pi/180)
    self.u, self.v = uvc_rot.real, uvc_rot.imag
    self.direction = self.direction + ang_deg
    if hasattr(self, 'shu'):
        shuvc_orig = self.shu + 1j * self.shv
        shuvc_rot = shuvc_orig * np.exp(1j * ang_deg)
        self.shu, self.shv = shuvc_rot.real, shuvc_rot.imag
    if write_to_proc_dict:
        self._record_proc_step('uvrot', ang_deg)

#######################################################################

def _remove_rows(self, cut_ind, write_to_proc_dict=True, 
                    return_val_tuple=False):
    """
    Removes the selected rows (bins) from the dataset altogether.
    Applies to all variables with a depth dimension.

    rowinds: Row(s) to remove. Can be one of the following:
        - An integer (4)
        - A slice (slice(12, 15))
        - An iterable ( [1, 3, 4], np.array([1, 3, 4]) )
    """
    row_inds = np.arange(self.Ndep)
    rows_to_cut = row_inds[cut_ind]
    rowdeps_to_cut = np.round(self.dep_med[cut_ind], 1)
    Ndep_0 = self.Ndep
    if hasattr(rows_to_cut, '__iter__'):
        n_rows_to_cut = len(rows_to_cut)
    else:
        n_rows_to_cut = 1

    remain_bool = np.bool_(np.ones(self.Ndep))
    remain_bool[cut_ind] = False
    
    for attrname in dir(self):
        if isinstance(getattr(self, attrname), np.ndarray):
            if getattr(self, attrname).shape[0] == self.Ndep:
                setattr(self, attrname, 
                        getattr(self, attrname)[remain_bool])

    self.Ndep = len(self.dep_med)

    val_tuple = (
            n_rows_to_cut, rows_to_cut, rowdeps_to_cut, Ndep_0,
            self.Ndep, -100 * (1 - Ndep_0 / self.Ndep))

    if write_to_proc_dict:
        self._record_proc_step('rowrem', val_tuple)
    if return_val_tuple:
        return val_tuple

#######################################################################

def _remove_cols(self, cut_ind, write_to_proc_dict=True, 
                    return_val_tuple=False):
    """
    Removes the selected columns (profiles) from the dataset 
    altogether.
    
    Applies to all variables with a time dimension.

    cut_ind: Columns(s) to remove. Can be one of the following:
        - An integer (4)
        - A slice (slice(12, 15))
        - An iterable ( [1, 3, 4], np.array([1, 3, 4]) )
    """
    col_inds = np.arange(self.Nt)
    cols_to_cut = col_inds[cut_ind]
    if hasattr(cols_to_cut, '__iter__'):
        n_cols_to_cut = len(cols_to_cut)
    else:
        n_cols_to_cut = 1
    remain_bool = np.bool_(np.ones(self.Nt))
    remain_bool[cut_ind] = False

    for attrname in dir(self):
        if isinstance(getattr(self, attrname), np.ndarray):
            if getattr(self, attrname).shape[(-1)] == self.Nt:
                index = tuple([slice(None)] 
                    * (getattr(self, attrname).ndim - 1) 
                    + [remain_bool])
                setattr(self, attrname, getattr(self, attrname)[index])

    Nt_0 = self.Nt
    self.Nt = len(self.t_mpl)
    val_tuple = (
        n_cols_to_cut, Nt_0,
        self.Nt, -100 * (1 - Nt_0 / self.Nt))
    if write_to_proc_dict:
        self._record_proc_step('colrem', val_tuple)
    if return_val_tuple:
        return val_tuple


#######################################################################

def _flip_to_ori(self):
    '''
    Flips all variables with a depth dimensions if orientation is 
    downward.
    '''

    if self._ori == 'down':
        for attrname in dir(self):
            if isinstance(getattr(self, attrname), np.ndarray):
                if getattr(self, attrname).shape[0] == self.Ndep:
                    setattr(self, attrname, 
                            getattr(self, attrname)[::-1])
    else: pass


#######################################################################

def _calculate_tilt(self):
    """
    Calculate tilt from pitch/roll components.
    See Mantovanelli et al 2014 and Woodgate et al 2011.
    """
    P_rad = self.pitch / 180 * np.pi
    R_rad = self.roll / 180 * np.pi
    T_rad = np.arccos(np.sqrt(1 - np.sin(P_rad) ** 2
                        - np.sin(R_rad) ** 2))
    self.tilt = T_rad * 180 / np.pi
    self.units['tilt'] = 'degrees'
    self.var_desc['tilt'] = 'Tilt (deg)'


#######################################################################
# Functions for record keeping / documentation
#######################################################################

def _update_summary(self):
    '''
    Update summary dictionary containing some basic statistics of the
    dataset.
    '''
    sd_ = {}

    timefmt = '%d-%b-%Y %H:%M'
    time_span = self.t_mpl[-1] - self.t_mpl[0]
    sd_['time_range'] = ('Time range: %s to %s (%.1f days)' 
        % (self.t_datetime[0].strftime(timefmt), 
            self.t_datetime[-1].strftime(timefmt),
            time_span))

    sd_['depth_range'] = ('Median bin depth: %.1f m to %.1f m'
        % (self.dep_med.min(), self.dep_med.max()))
   
    d_t_min = np.median(np.diff(self.t_mpl))*24*60
    d_t_sec = d_t_min*60
    sd_['time_step'] = 'Median time step: %.1f min / %.0f sec'%(
                        d_t_min, d_t_sec)

    sd_['tilt_mean'] = 'Mean tilt: %.1f deg'%(self.tilt.mean())
    sd_['tdep_mean'] = 'Mean transducer depth: %.1f m'%(self.tdepth.mean())

    sd_['u_mean'] = 'MEAN U: %.2f cm/s'%(self.u.mean())
    sd_['v_mean'] = 'MEAN V: %.2f cm/s'%(self.v.mean())
    sd_['sp_mean'] = 'Mean speed: %.2f cm/s'%(
            np.ma.sqrt((self.u**2 + self.v**2)).mean())

    uvdir = np.arctan2(self.v.mean(), self.u.mean())*180/np.pi
    sd_['uv_dir'] = '%.1f deg CCW of east'%(uvdir)

    self._sum_dict = sd_


#######################################################################

def _record_proc_step(self, operation_name, value):
    """
    Records processing steps to the proc_dict dictionary.

    Keeps track of multiple interations of the same operation by
    appending a number (_2, _3, ..) to the dictionary key (e.g.
    "uvrot2").

    Assigns a tuple to the dictionary entry. The first index of the
    tuple is the operation type (e.g. "uvrot"), the second is
    *value*, which may be a number or an interable, depending on the 
    operation.

    The function *print_proc_steps* knows how to read these entries
    and produce a human-readable processing history string. 
    """
    opnm = operation_name

    # If we have already performed this step, we distinguish
    # subsequent steps of the same operation with _1, _2, etc.
    if opnm in self._proc_dict.keys():
        nn = 1
        while opnm in self._proc_dict.keys():
            nn = nn + 1
            opnm = '%s_%i' % (operation_name, nn)

    self._proc_dict[opnm] = (operation_name, value)