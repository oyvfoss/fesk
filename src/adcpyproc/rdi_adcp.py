
#######################################################################
# rdi_adcp
#######################################################################
"""
- Loads a matfile created using WinADCP and generates an RdiObj object
  with various specialized functions to perform post-processing of the
  dataset.
- The result can be exported to matfiles, python dictionaries, and
  (TBD) netCDF.
- Run rdi_adcp.exaple() to print an example of processing using
  these tools.
"""

#######################################################################
# IMPORTS
#######################################################################

# Standard libraries
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import xarray as xr
import datetime as dt
import pickle
from matplotlib.dates import num2date, date2num
import scipy.io
from collections import OrderedDict
import time
import pickle
import warnings
import os

# Internal libraries
from adcpyproc.misc_functions import declination, t_python_to_mat
# Dictionaries, long strings, warming messages etc
from adcpyproc import _rdi_defs

# Imported internal methods for the RdiObj class
# (We put most internal class methods in these external modules and reserve
# the present module for user relevant functions)
from adcpyproc import _read_rdifile_methods_rdi_adcp
from adcpyproc import _data_operation_methods_rdi_adcp
from adcpyproc._netcdf_export import _create_netcdf
from adcpyproc import __file__ as _initloc
from adcpyproc import rdi_toolbox


# Enable interactive mode for pyplot
# (necessary for correct rendering of figures in loops, etc)
plt.ion()


#######################################################################
# CLASS DEFINITION AND CONSTRUCTOR METHOD
# (Loading an ADCP file into this processing framework amounts to
# generating an RdiObj, e.g.: d = rdi_adcp.RdiObj('adcpdata.mat')
#######################################################################

class RdiObj:

    def __init__(self, matfile, orientation = 'auto'):
        '''
        Load matfile (typically processed from binary file using
        WinADCP) and transfer variables to RdiObject with conventions
        (see _rdi_defs.py).

        If 'orientation' is 'auto', we assume an instrument orientation 
        as specified in the file. Can be overridden by setting orientation
        as 'up' or 'down'

        Convention: Depth decreases with increasing 0th index.

        Usage example:

        d = RdiObj('sourcefile_from_winadcp.mat')
        '''

        # Adding methods from the internal modules imported above
        for method_module in (
            _read_rdifile_methods_rdi_adcp,
            _data_operation_methods_rdi_adcp,
        ):

            [setattr(RdiObj, f.__name__, f) for _, f in
             method_module.__dict__.items() if callable(f)]

        # Setup the object, populating various fields
        self._setup(matfile, orientation=orientation)


    #######################################################################
    ####### PROCESSING FUNCTIONS
    #######################################################################

    def add_serial_number(self, SN):
        """
        Add serial number information.
        """
        self.SN = SN
        print('Set serial number: %i' % SN)

    #######################################################################

    def set_latlon(self, lat=None, lon=None):
        """
        Add lat and lon (decimal degrees). Can be single floats or arrays of
        floats.
        """

        if not lat:
            lat = np.float(input('Input latitude (decimal degrees): '))
        if not lon:
            lon = np.float(input('Input longitude (decimal degrees): '))
        self.lon, self.lat = lon, lat

        self.latlon_single = not hasattr(self.lon, 'len')

        print('Set lat and lon.')

    #######################################################################

    def correct_magdec(self, model=None):
        """
        Correction for time-dependent magnetic declination using
        WMM2010, WMM2015 or WMM2020 using functionality from the
        *geomag* module.

        Using one single model for each deployment (choosign based on
        average time stamp of the data).

        Rotates *u*, *v*, and shear components *shu*, *shv* if they are
        available. Also adds the offset to the *direction* field.

        (see documentation for *adcpyproc.misc_functions.declination*)
        """

        if hasattr(self, 'lat'):
            lat, lon = self.lat, self.lon
        else:
            raise Exception(_rdi_defs._excep_str['magdec_lonlat'])

        if not model:
            t_mid = num2date(self.t_mpl.mean())
            if t_mid.year > 2019:
                model = 'WMM2020'
            if t_mid.year > 2024:
                warnings.warn(_rdi_defs._warn_str['magdec_2024'])
            elif t_mid.year < 2015:
                model = 'WMM2010'
            if t_mid.year < 2010:
                warnings.warn(_rdi_defs._warn_str['magdec_2010'])
        else:
            model = 'WMM2015'

        self.magdec = declination(lat, lon, (self.t_datetime), model=model,
                                  explicit=True)
        self._rotate_uv((-self.magdec), write_to_proc_dict=False)
        magdec_mean = np.angle(np.mean(np.exp(1j * self.magdec
                               * np.pi / 180))) * 180 / np.pi
        self._record_proc_step('magdec', magdec_mean)

    #######################################################################

    def apply_depth_offset(self, offset_m):
        """
        Apply an offset to depths (transducer and bin depths)

        positive *offset_m* -> shift downwards
        """
        med_tdepth_bef = np.median(self.tdepth)
        self.tdepth = self.tdepth + offset_m
        self.dep += offset_m
        self.dep_med = np.median((self.dep), axis=1)
        med_tdepth_aft = np.median(self.tdepth)
        depoffs_values = (
            offset_m, med_tdepth_bef, med_tdepth_aft)
        self._record_proc_step('depoffs', depoffs_values)

    #######################################################################

    def remove_ship_time(self, thr_m=3, autoaccept=False):
        """
        Reads the transducer pressure (depth), suggests
        deployment/recovery range, and, if accepted by the user,
        removes the edges.

        The suggested range is from the first to the last index where
        transducer depth is within *thr_m* meters of the deployment
        median depth.

        autoaccept = True: Skip user prompt - just accept the chopping
        suggested by the algorithm.
        """
        med_tdep = np.ma.median(self.tdepth)
        ind_ = np.ma.where(med_tdep - self.tdepth < thr_m)[0]
        indices = np.arange(self.Nt)

        if not autoaccept:
            fig, ax = plt.subplots(figsize=(6, 2))
            ax.plot(indices, self.tdepth, '.-k')
            ax.set_xlabel('Index')
            ax.set_ylabel('Depth [m]')
            plt.tight_layout()
            plt.show()

            fig_handle = display(fig, display_id=True )

        accept_range = False
        first = True

        while accept_range == False:
            if first:
                ind0, ind1 = ind_[0], ind_[(-1)]
            else:
                ind0 = int(input('Previous start index was %i.' % ind0
                        + ' Input new: '))
                ind1 = int(input('Previous end index was %i.' % ind1
                        + ' Input new: '))

            sl = slice(ind0, ind1 + 1)
            if not autoaccept:
                tx_range = ax.plot((indices[sl]), (self.tdepth[sl]),
                                   '.r', ms=2)
#                plt.pause(0.1)
                fig_handle.update(fig)

                accept_range_yn = input(
                    'Current range (%i, %i).' % (ind0, ind1)
                    + ' Accept? (y/n): ')
                print(accept_range_yn)
            else:
                accept_range_yn = 'y'
            if accept_range_yn == 'y':
                accept_range = True
            else:
                tx_range[0].remove()
                first = False
        print('Range accepted -> shortening time series to range '
              '(%i, %i). ' % (ind0, ind1))

        Nt_0 = self.Nt
        self._remove_cols((slice(ind1 + 1, None)), write_to_proc_dict=False)
        Nt_1 = self.Nt
        self._remove_cols((slice(None, ind0)), write_to_proc_dict=False)
        Nt_2 = self.Nt
        chop_ship_values = (
            Nt_1 - Nt_2, Nt_0 - Nt_1, Nt_0, Nt_2, 100 * (1 - Nt_2 / Nt_0))
        self._record_proc_step('shiptrem', chop_ship_values)


    #######################################################################

    def reject_surface_bins(self):
        '''
        Reject entire rows where the median depth is above the surface.
        '''

        rrows = np.where(self.dep_med<0)
        val_tuple = self._remove_rows(rrows, return_val_tuple = True,
                                      write_to_proc_dict=False)
        self._record_proc_step('rowrem_surf', val_tuple)

    #######################################################################

    def reject_rows(self, rows = False, masked_max = None, autoaccept = False):
        '''
        Reject entire rows.
        
        Can be specified, or determined based on where the percent masked
        exceeds *masked_max* (best to apply *after* masking/flagging).

        rows: 
        autoaccept = False: Prompts the user before removing.
        '''

        row_inds = np.arange(self.Ndep)
        mask_prof = self.u.mask.mean(axis=1)*100

        if rows:
            rej_flag = row_inds[rows]

        else:
            if not masked_max:
                raise Exception(_rdi_defs._excep_str['rrej_nodef'])
        if masked_max:
            rej_flag = mask_prof > masked_max

            if rows:
                raise Exception(_rdi_defs._excep_str['rrej_bothdef'])
            else:
                rows = np.where(rej_flag)[0]

        rows_to_cut = row_inds[rows]
        nrej = len(rows_to_cut)

        if not autoaccept:
            if masked_max:
                print('Found %.0f rows where percent masked exceeds %.1f pct.'
                %(nrej, masked_max))   
                return_vals = True  
            
            line_top = ''#'-'*80 + '\n'
            line1 = '='*80 + '\nIndex: '
            line2 = 'Depth: '
            line3 = 'Masked:'
            line4 = 'Reject:'

            for ind_ in row_inds:
                line1 += ('%.0f'%ind_).rjust(6) 
                line2 += ('%.0f'%self.dep_med[ind_]).rjust(6)
                if mask_prof[ind_]>1:
                    line3 += ('%.1f'%mask_prof[ind_]).rjust(6)
                else:
                    line3 += '    <1'
                if ind_ in rows_to_cut:
                    l4str = '     Y'
                else:
                    l4str = '     -'
                line4 += l4str

                if len(line4)>75:
                    line_top += '%s\n%s\n%s\n%s\n'%(line1, line2, line3, line4)
                    line1 = '-'*len(line4) + '\nIndex: '
                    line2 = 'Depth: '
                    line3 = 'Masked:'
                    line4 = 'Reject:'

            print('%s%s\n%s\n%s\n%s'%(line_top, line1, line2, line3, line4))
            print('='*80+'\n')

            print('Rejecting rows: %s'%rows_to_cut) 
            print('Median depths:  %s'%np.round(self.dep_med[rows], 1))
            reject = input('\nREMOVE THESE ROWS ALTOGETHER? (y/n): ')

        else: 
            reject = 'y'

        if reject == 'y':
            if not masked_max:
                self._remove_rows(rows, write_to_proc_dict=True)
            else:
                val_tuple = self._remove_rows(rows, return_val_tuple = True,
                            write_to_proc_dict=False)
                self._record_proc_step('rowrem_mask', 
                        (val_tuple[0], val_tuple[1],  val_tuple[2], masked_max,
                        val_tuple[3], val_tuple[4], val_tuple[5]))
            print('\n-> %i rows rejected'%nrej)

        else:
            print('\n-> No rows rejected.')


    #######################################################################
    # FLAGGING FUNCTIONS
    # Masking entries based on certain criteria
    #######################################################################

    def mask_umax(self, thr_uvamp = 100):
        '''
        Mask measurements where horizontal current speed
        is greater than *thr_uvamp* (cm/s). 
        '''

        ns = self._apply_mask(np.ma.sqrt(self.u**2 + self.v**2)>thr_uvamp, 
                include_errvel=True, return_ns=True)

        perc_masked = (1 - ns[2]/ns[1])*100 # Pct reduction in unmasked entries
        entries_masked = ns[1]-ns[2] 
        self._record_proc_step('mask_umax', (thr_uvamp, entries_masked, ns[1], 
                                perc_masked))

    def mask_surf_sidelobe(self, beam_angle = 20, thr_perc = None):
        '''
        TBD
        Mask entries collected within a certain distance of the surface
        (given as a percentage of the transducer depth).

        Calculated based on beam angle (default: 20 degrees), 
        or specified as *thr_perc* percent (the latter overrides the former).
        '''
        warnings.warn(
            '"mask_surf_sidelobe()" not implemented yet - > did nothing')
        pass

    #######################################################################

    def mask_pg(self, thr_pg = 75):
        '''
        Mask where the total number of pings which come from good beams 
        is less than *thr_pg* (percent).
        '''

        ns = self._apply_mask((self.pg1+self.pg4)<thr_pg,
                include_errvel=True, return_ns=True)

        perc_masked = (1 - ns[2]/ns[1])*100 # Pct reduction in unmasked entries
        entries_masked = ns[1]-ns[2] 
        self._record_proc_step('mask_pg', (thr_pg, entries_masked, ns[1], 
                                perc_masked))


    #######################################################################

    def mask_amp_a(self, thr_amp = 64):
        '''
        Mask where the average backscatter amplitude is less than *thr_amp* 
        (db).
        '''

        ns = self._apply_mask((self.amp_a)<thr_amp,
                include_errvel=True, return_ns=True)

        perc_masked = (1 - ns[2]/ns[1])*100 # Pct reduction in unmasked entries
        entries_masked = ns[1]-ns[2] 
        self._record_proc_step('mask_amp_a', (thr_amp, entries_masked, ns[1], 
                                perc_masked))

    #######################################################################

    def mask_errvel(self, thr_errvel = 50):
        '''
        Mask where the error velocity is greater than *thr_amp* (cm/s).
        Where errvel is invalid (masked or non-existent due to 3-beam solution),
        no masking is applied.        
        '''

        errvel_above_thr = np.abs(self.errvel)>thr_errvel
        errvel_above_thr[self.errvel.mask] = False

        # When no error velocity is available (e.g. for a 3-beam solution), 
        # RDI software sets this to -3276.8 or similar. Don't want to reject 
        # these entries as the 3-beam solution could still be valid. 
        errvel_above_thr[self.errvel<-3200] = False

        ns = self._apply_mask(errvel_above_thr,
                include_errvel=True, return_ns=True)

        perc_masked = (1 - ns[2]/ns[1])*100 # Pct reduction in unmasked entries
        entries_masked = ns[1]-ns[2] 
        self._record_proc_step('mask_errvel', (thr_errvel, entries_masked, ns[1], 
                                perc_masked))


    #######################################################################

    def mask_w(self, thr_w = 30):
        '''
        Mask where the vertical velocity is greater than *thr_w* (cm/s).
        '''
        ns = self._apply_mask(np.abs(self.w)>thr_w,
                include_errvel=True, return_ns=True)

        perc_masked = (1 - ns[2]/ns[1])*100 # Pct reduction in unmasked entries
        entries_masked = ns[1]-ns[2] 
        self._record_proc_step('mask_w', (thr_w, entries_masked, ns[1], 
                                perc_masked))


    #######################################################################

    def mask_ca(self, thr_cor = 30):
        '''
        Mask where the pulse-to-pulse correlation amplitude in two or more 
        beams is less than *thr_cor* (counts, 0-255).
        '''
        cf1, cf2 = np.int_(self.ca1<thr_cor), np.int_(self.ca2<thr_cor) 
        cf3, cf4 = np.int_(self.ca3<thr_cor), np.int_(self.ca4<thr_cor)

        ns = self._apply_mask((cf1+cf2+cf3+cf4)>1,
                include_errvel=True, return_ns=True)

        perc_masked = (1 - ns[2]/ns[1])*100 # Pct reduction in unmasked entries
        entries_masked = ns[1]-ns[2] 
        self._record_proc_step('mask_ca', (thr_cor, entries_masked, ns[1], 
                                perc_masked))


    #######################################################################

    def mask_ca_mean(self, thr_cor = 30):
        '''
        Mask where the beam average pulse-to-pulse correlation
        is less than *thr_amp* (counts, 0-255).        
        '''
        
        ns = self._apply_mask(self.ca_a<thr_cor,
                include_errvel=True, return_ns=True)

        perc_masked = (1 - ns[2]/ns[1])*100 # Pct reduction in unmasked entries
        entries_masked = ns[1]-ns[2] 
        self._record_proc_step('mask_ca_mean', (thr_cor, entries_masked, ns[1],
                                perc_masked))

    #######################################################################

    def mask_amp_jump(self, max_amp_increase = 30, mask_above = False):
        '''
        Masking bins where the beam amplitude of any beam has a jump exceeding 
        *max_amp_increase* (dB). 
        
        mask_above = False -> Masks only the exact bins above the jump.
        mask_above = True -> Also masks all bins *above* the jump.
        '''
        
        amp_stack = np.float_(np.dstack([self.amp1, self.amp2, self.amp3, 
                                         self.amp4]))
        d_amp_stack_bool = (np.diff(amp_stack, axis = 0).max(axis = -1) 
                            > max_amp_increase)

        # Array which is 1 where there is an amplitude jump (on the index 
        # *after* the jump)
        jump_arr = np.zeros(self.u.shape)
        jump_arr[1:] = d_amp_stack_bool

        if not mask_above:
            mask_arr = np.bool_(jump_arr)
            jump_str = 'mask_amp_jump' 
        else:
            after_jump_arr = np.cumsum(jump_arr, axis = 0)
            mask_arr = after_jump_arr > 0
            jump_str = 'mask_amp_jump_abv' 

        ns = self._apply_mask(mask_arr,
                include_errvel=True, return_ns=True)

        perc_masked = (1 - ns[2]/ns[1])*100 # Pct reduction in unmasked entries
        entries_masked = ns[1]-ns[2] 
        self._record_proc_step(jump_str, (max_amp_increase, 
                                entries_masked, ns[1],
                                perc_masked))

    #######################################################################

    def mask_tilt(self, thr_tilt = 20):
        '''
        Mask profiles where tilt is greater than *thr_tilt* (degrees). 
        '''

        ns = self._apply_mask(self.tilt>thr_tilt,
                include_errvel=True, return_ns=True)

        nprofs = np.ma.sum(self.tilt>thr_tilt)
        perc_masked = (1 - ns[2]/ns[1])*100 # Pct reduction in unmasked entries
        entries_masked = ns[1]-ns[2] 
        self._record_proc_step('mask_tilt', (thr_tilt, nprofs, 
                                entries_masked, ns[1],
                                perc_masked))



    #######################################################################
    # CALCULATE SHEAR 
    #######################################################################

    def calculate_shear(self):
        """
        Calculate vertical shear.
        """
        for shkey in ('shu', 'shv'):
            empty_arr = self.u.copy() * 0
            empty_arr.mask = True
            setattr(self, shkey, empty_arr)

        # Check that we haven't chopped out any individual bins..
        if (-np.round(np.diff(self.dep_med)) != self.binsize).any():
            raise Exception(_rdi_defs._excep_str['shear_bins'])
        
        self.shu[1:-1] = (0.01 * (self.u[2:] - self.u[:-2]) 
                        / (2 * self.binsize))
        self.shv[1:-1] = (0.01 * (self.v[2:] - self.v[:-2]) 
                        / (2 * self.binsize))
        self.s2 = self.shu ** 2 + self.shv ** 2
        self.units['shu'], self.units['shv'], self.units['s2'] = (
                '1/s', '1/s', '1/s2')

        self._record_proc_step('shear', '')

    #######################################################################
    # PRINTING/HISTORY/DOC/PLOTTING FUNCTIONS
    #######################################################################

    def add_note(self, note_str=None):
        """
        Add a note to be stored with the data. Could be any useful detail
        (processing info, version control, purpose, whatever).

        Can be supplied as a string (*note_str* - otherwise gives an input
        prompt).
        """
        if not note_str:
            note_str = input('Input processing note: ..\n ')
        if hasattr(self, 'proc_notes'):
            self.proc_notes = self.proc_notes + '\n\n%s' % note_str
        else:
            self.proc_notes = note_str

    #######################################################################

    def print_proc(self, return_str=False):
        """
        Print a processing history showing the steps applied to the dataset.
        """
        proc_str = ''
        for key in self._proc_dict.keys():
            type_str = self._proc_dict[key][0]
            desc_str = (_rdi_defs._proc_str_dict_rdi[type_str] 
                        % self._proc_dict[key][1])
            proc_str += desc_str + '\n'
        
        if return_str:
            return proc_str
        else:
            print(proc_str)

    #######################################################################

    def print_notes(self):
        """
        Print any processing notes (processig notes can be added using
        the *add_note()* method).
        """
        if hasattr(self, 'proc_notes'):
            print('PROCESSING NOTES:\n\n%s' % self.proc_notes)
        else:
            print('(no processing notes).')

    #######################################################################

    def print_system_info(self, return_str=False):
        """
        Print details about the source file, instrument, and
        configuration extracted from the source mat file.
        """

        sys_str = ''
        for key in self._conf_dict:
           sys_str += '%s:\n  %s\n'%(key, self._conf_dict[key])
        
        if return_str:
            return sys_str
        else:
            print(sys_str)

    #######################################################################

    def print_summary(self):
        """
        Print a basic summary of the dataset and some average values.
        """

        self._update_summary()

        sum_str = 'SUMMARY OF DATASET AFTER CURRENT PROCESSING STEPS:\n\n'
        for key in ['time_range', 'depth_range', 'time_step', 'tdep_mean']:
           sum_str += '%s\n'%(self._sum_dict[key])

        sum_str += '\n%s\n'%(self._sum_dict['u_mean'])
        sum_str += ('%s'%self._sum_dict['v_mean'] + ' '*5 
                    + '(%s)\n'%self._sum_dict['uv_dir'])

        sum_str += '%s\n'%(self._sum_dict['sp_mean'])

        print(sum_str)

    # Plotting functions from rdi_toolbox
    #######################################################################

    def hist(self, varnm, nbins = 50, line = None, return_ax = False):
        '''
        Show histogram(s) of a variable.
        
        For 1D parameters (tdepth, heading, tilt..): 
            Showing a simple histogram
        For 2D parameters (u, amp_a, ca1..)
            Simple histogram + histograms by depth
        
        varnm: Short variable name ('u', 'amp1', 'tilt'..)
        nbins: Number of bins to use for histogram
        line: Plot a vertical line here. Useful for looking at 
            thresholds.
        return_ax: return the axis object
        '''

        ax = rdi_toolbox.hist(self, varnm, nbins = nbins, line = line,
                        return_ax = True)
        plt.show()

        if return_ax:
            return ax

    #######################################################################

    def plot_ellipse(self, dep0 = None, dep1 = None, lp_days = 5, 
                    return_ax = False):
        '''
        Plot of u and v components averaged between *dep0* and *dep1* and
        low pass filtered with a running mean of *lp_days*.

        Showing the mean current vector, the low-pass-filtered 
        and subsampled currents, and the semi-major and -minor axes
        of the variance ellipse. 

        return_ax: return the axis object

        '''
        ax = rdi_toolbox.plot_ellipse(self, dep0 = dep0, dep1 = dep1, 
                lp_days = lp_days)
        plt.show()
        if return_ax:
            return ax

    #######################################################################

    def uv_panels(self, wlen = None, subsamp = None, cmap_saturate = 99, 
                save_file = None, return_ax = False):
        '''
        Quick panels of U and V.

        wlen: Window length for running mean (days). Default ~1 day
        subsamp: Subsampling period (for reducing the array sizes) (days).
                Default is around 1/4 of wlen.s
        cmap_saturate: Saturating the colormap at this percentile. 
        save_file: Enter a valid file path to save the figure.
        return_ax: return the axis object
        '''

        ax = rdi_toolbox.uv_panels(self, wlen_days = wlen, 
            subsamp_days = subsamp,
            cmap_saturate = cmap_saturate, save_file = save_file,
            return_ax = True)
        plt.show()

        if return_ax:
            return ax


    #######################################################################
    # EXPORTING FUNCTIONS
    #######################################################################

    def to_dict(self, sparse = True, params_list = None):
        '''
        Create a python dictionary with the relevant data.
        For export to other formats.

        sparse = True skips a bunch of parameters and only exports
                 (t_mpl, dep, u, v, amp_a, tdepth, (shu, shv)) as 
                 well as some documentation parameters. 

        params_list: Iterable of variable names if these should be
                     specified exactly (overrides *sparse*) 

        Returns the dictionary.
        '''

        # Obtain the list names of the variables we want to include 
        if params_list:
            outvars = params_list
        else:
            if sparse:
                outvars = _rdi_defs._outvars_rdi['sparse']
            else:
                outvars = _rdi_defs._outvars_rdi['full']

        # Loop through and transfer variables and units/descriptions
        D = {}
        D['var_desc'] = {}
        D['units'] = {}

        for var_ in outvars:
            if hasattr(self, var_):
                if var_ in self.var_desc.keys():
                    D['var_desc'][var_] = self.var_desc[var_]
                if var_ in self.units.keys():
                    D['units'][var_] = self.units[var_]
                D[var_] = getattr(self, var_)
        
        # Record keeping: Final data details
        self._record_proc_step('shape_final', (
            self.Nt, self.Ndep, self.Nt * self.Ndep))
        self._record_proc_step('mask_final', (
            self.u.mask.sum(), 100 * self.u.mask.sum() 
            / (self.Nt * self.Ndep)))

        D['processing_history'] = self.print_proc(return_str = True)
        D['system_info'] = self.print_system_info(return_str = True)

        # Removing the final data details in case we want to continue 
        # processing on the RdiObj
        self._proc_dict.pop('shape_final')
        self._proc_dict.pop('mask_final')

        return D

    #######################################################################

    def to_matfile(self, file_name, sparse = True, matlab_format = '5'):
        '''
        Export to a matfile (replacing matplotlib time with a time stamp
        string parseable by matlab).
        '''

        if sparse:
            outvars = _rdi_defs._outvars_rdi['full'].copy()
        else:
            outvars = _rdi_defs._outvars_rdi['sparse'].copy()

        outvars.remove('t_mpl')
        outvars.remove('t_datetime')
        
        D = self.to_dict(params_list = outvars)

        # Time string
        t_fmt = "%d-%b-%Y %H:%M:%S"
        D['t_str'] = np.array([t_.strftime(t_fmt) for t_ in self.t_datetime])
        D['var_desc']['t_str'] = 'Time string'

        # Matlab time
        D['t_mat'] = t_python_to_mat(self.t_mpl)
        D['var_desc']['t_mat'] = 'Time, matlab serial date number'
        D['units']['t_mat'] = 'Days since 00-Jan-0000'

        # Convert masks to NaNs
        for key in D.keys():
            if isinstance(D[key], np.ma.core.MaskedArray):
                var_ = np.array(D[key])
                var_[D[key].mask] = np.nan
                D[key] = var_
        
        # Save
        scipy.io.matlab.savemat(file_name, D, format = matlab_format)
        print('SAVED TO : %s'%file_name)

    #######################################################################

    def to_pickle(self, file_name, sparse = True):
        '''
        Save to a pickled dictionary.
        '''

        D = self.to_dict(sparse = sparse)

        with open(file_name, 'wb') as f:
            pickle.dump(D, f)

        print('SAVED TO : %s'%file_name)

    #######################################################################

    def prepare_for_netcdf(self, netcdf_dir,
                    attr_file='netcdf_global_attributes.py'):
        '''
        Create a file with netCDF global attributes to be filled out
        before creatign the netcdf
        '''

        # Get the path of the default file with netCDF attributes
        ncattrfile_default = _initloc.replace('__init__.py', 
            'netcdf_formatting/netcdf_global_attributes.txt')

        # Make a copy of this file in the destination'
        cp_cmd = 'cp %s %snetcdf_global_attributes.py'%(
                ncattrfile_default, netcdf_dir)
        try:
            os.popen(cp_cmd) 
        except: 
            raise Exception(
                'Could not execute the following command: %s'%cp_cmd)

        self.prepared_for_nc = True


    #######################################################################


    def to_netcdf(self, netcdf_dir, netcdf_file, 
                    attr_file='netcdf_global_attributes.py'):
        '''
        Save to a NetCDF file (TBD).
        '''

        if hasattr(self, 'prepared_for_nc'):
            _create_netcdf(self, netcdf_dir, netcdf_file, 
                    attr_file='netcdf_global_attributes.py')
        else:
            raise Exception('Before you can make a netCDF file, ' 
            'you need to run *prepare_for_netcdf()* and edit the file '
            '*netcdf_global_attributes.py* with details of your dataset.')


    #######################################################################




#######################################################################

def example(extended = False):
    '''
    Print a quick example of processing steps.

    extended = True/False: Print an extended/minimal example.
    '''
    if extended:
        example_str = _rdi_defs._example_str['extended']
    else:
        example_str = _rdi_defs._example_str['minimal']

    print(example_str)