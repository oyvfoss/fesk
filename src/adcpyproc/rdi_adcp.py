
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
import datetime as dt, pickle
from matplotlib.dates import num2date, date2num
import scipy.io
from collections import OrderedDict
import time
import pickle

# Internal libraries
from adcpyproc.misc_functions import declination, t_python_to_mat
from adcpyproc import variables_dictionary 


#######################################################################
# CLASS DEFINITION AND CONSTRUCTOR METHODS
# (Loading an ADCP file into this processing framework amounts to
# generating an RdiObj, e.g.: d = rdi_adcp.RdiObj('adcpdata.mat')
#######################################################################

class RdiObj:

    def __init__(self, matfn):
        """
    Load matfile (typically processed from binary file using
    WinADCP) and transfer variables to RdiObject with conventions
    (see variables_dictionary.py).

    Convention: Depth decreases with increasing 0th index.
    """
        print('Initializing postprocessing by creading RdiObj object from'
              ' this file:\n %s..\n' % matfn)
        m = scipy.io.matlab.loadmat(matfn)

        print('Writing some metadata and creating some useful structures..')
        self.meta = ('Created on %s using *pyrdi_adcp* from the source file %s.' % (
            dt.datetime.now().strftime('%d.%m.%y'), matfn))
        self.binsize = m['RDIBinSize']
        self.var_desc, self.units = {}, {}

        # Documentation    
        self._hist_dict = {'procdate':dt.datetime.now().strftime('%d.%m.%Y'), 
            'source_matfile':matfn, 
            'source_rawfile':m['RDIFileName'][0]}
        self._proc_dict = OrderedDict()
        self._extract_config(m)

        print('Loading data to numpy arrays and converting to '
              'preferred names/units..')
        self._load_data_from_matfile(m, variables_dictionary.vardict_rdi, 
                                    explicit = False)

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
    #######################################################################
    ####### PROCESSING FUNCTIONS 
    #######################################################################
    #######################################################################


    def add_serial_number(self, SN):
        """
        Add serial number information.
        """
        self.SN = SN
        print('Set serial number: %i' % SN)

    #######################################################################


    def set_latlon(self, lat, lon):
        """
        Add lat and lon (decimal degrees).
        """
        self.lon, self.lat = lon, lat
        print('Set lat and lon.')

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
            raise Exception(
            'Correcting for magnetic declination (*correct_magdec()*'
            'method): \n'
            'Lat/lon need to be specified using the  ' 
            'set_latlon() method.')

        if not model:
            t_mid = num2date(self.t_mpl.mean())
            if t_mid.year > 2019:
                model = 'WMM2020'
            if t_mid.year > 2024:
                raise Warning('WMM2020 is only valid through 2024.' 
                    ' Proceed with caution..')
            elif t_mid.year < 2015:
                model = 'WMM2010'
            if t_mid.year < 2010:
                raise Warning('WMM2010 is only valid starting in 2010.' 
                              'Proceed with caution..')
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

        accept_range = False
        first = True

        while accept_range == False:
            if first:
                ind0, ind1 = ind_[0], ind_[(-1)]
            else:
                ind0 = int(input('Previous start index was %i.'%ind0+
                        ' Input new: '))
                ind1 = int(input('Previous end index was %i.'%ind1+
                        ' Input new: '))

            sl = slice(ind0, ind1 + 1)
            if not autoaccept:
                tx_range = ax.plot((indices[sl]), (self.tdepth[sl]), 
                                    '.r', ms=2)
                plt.pause(0.1)
                plt.show(block=False)
                accept_range_yn = input(
                    'Current range (%i, %i).'%(ind0, ind1) 
                    + ' Accept? (y/n): ')
                print(accept_range_yn)

            if accept_range_yn=='y':
                accept_range = True

            if accept_range == False:
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
                raise Exception('Either *rows* or *masked_max* must be' 
                            'specified to RdiObj.reject_rows()')
        if masked_max:
            rej_flag = mask_prof > masked_max

            if rows:
                raise Exception('Both *rows* and *masked_max* specified'
                    ' to RdiObj.reject_rows(). Use only one.')
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
    # FUNCTIONS FOR READING INFORMATION FROM MATFILE
    #######################################################################

    def _load_data_from_matfile(self, mfile, vardict, explicit = False):
        """
        Loads data from WinADCP-derived matfile and adds it as 
        attributes (numpy arrays / masked arrays) to the RdiObj-
        """
        n_1dvars, n_2dvars = (0, 0)
        for nm in vardict.keys():
            if nm in mfile.keys():
                vd_ = vardict[nm]
                outnm_ = vd_['nm_out']
                if explicit:
                    print(('Reading data arrays: (field %i: %s..)%s\r' % (
                        n_2dvars + n_1dvars + 1, outnm_, ' '*10)), end = '')

                try:
                    if vd_['dim'] == 2:
                        outvar_ = mfile[nm].T * vd_['fact_in_to_out']
                        n_2dvars += 1
                    elif vd_['dim'] == 1:
                        outvar_ = mfile[nm] * vd_['fact_in_to_out']
                        n_1dvars += 1
                except Exception:
                    Warning('Issue with variable "%s"/"%s"' % (
                        outnm_, nm) +
                        '(skipping this variable).')

                setattr(self, outnm_, np.ma.squeeze(outvar_))
                vstr = '%s (%s). From field "%s".' % (vd_['desc'],
                    vd_['unit_out'], nm)
                self.var_desc[vd_['nm_out']] = vstr
                self.units[vd_['nm_out']] = vd_['unit_out']

        print('Reading data arrays: -> done.         '
                '(%i 2D variables)' % (n_2dvars + n_1dvars))
        for nm in ('u', 'v', 'w', 'direction', 'errvel'):
            try:
                outvar_masked_ = np.ma.masked_invalid(getattr(self, nm))
                outvar_masked_.mask[outvar_masked_ == -3276.8] = True
                setattr(self, nm, outvar_masked_)
            except AttributeError:
                pass

    #######################################################################

    def _get_orientation(self):
        '''
        Determines whether the instrument is up- or downlooking
        (based on the median value of the ori_ud field).

        Rejects profiles where *ori_ud* does not match the determined
        orientation.
        '''

        mval_ud = np.round(np.ma.median(self.ori_ud))
        if mval_ud == 1:
            self._ori = 'up'
            not_ud = 'down'
        else:
            self._ori = 'down'
            not_ud = 'up'

        wrong_way = self.ori_ud != mval_ud
        
        if wrong_way.any():
            val_tuple = self._remove_cols(np.where(wrong_way)[0], 
                write_to_proc_dict=False, 
                return_val_tuple=True)
            print('Removed %i profiles where the ori_ud field\n'
                  ' indicates that the instrument was facing %sward.'
                  %(val_tuple[0], not_ud))
            self._record_proc_step('colrem_wrongway', (val_tuple[0], 
                not_ud, *val_tuple[1:]))

    #######################################################################


    def _calculate_bin_depths(self, mfile):
        """
        From RDI ADCP matfile loaded as bunch: Calculate 2D depth field
        in meters.

        Ref https://www.bodc.ac.uk/data/documents/series/1014447/
        """

        if self._ori == 'up':
            sign = -1
        elif self_ori == 'down':
            sign = +1

        self.dep = (mfile['AnDepthmm'][np.newaxis, :] * 0.001 
                   + sign * mfile['RDIBin1Mid'] + sign * mfile['RDIBinSize'] 
                   * (mfile['SerBins'][:, np.newaxis] - 1))[0].T
        self.units['dep'] = 'm'
        self.var_desc['dep'] = 'Time-varying depth of bin (m)'

    #######################################################################

    def _extract_config(self, mfile):
        """
        Read configuration info / meta information from a source
        matfile loaded as dictionary. Output a dictionary.
        """
        self._conf_dict = {}
        conf_vars = [
            'RDIFileName', 'RDISystem', 'RDIPingsPerEns',
            'RDISecPerPing', 'RDIEnsDate', 'RDIEnsTime',
            'RDIEnsInterval', 'RDIBin1Mid', 'RDIBinSize']
        for conf_var in conf_vars:
            self._conf_dict[conf_var] = np.squeeze(mfile[conf_var])

    #######################################################################

    def _calculate_time(self, mfile):
        """
        Take time information in file (year, mon, ..) and create array
        of matplotlibtimes. Rounding down to nearest second (not
        concerned with sub-second time precision).
        """
        yr, mo, da = ('SerYear', 'SerMon', 'SerDay')
        hr, mi, se = ('SerHour', 'SerMin', 'SerSec')
        
        t_mpl = np.array([])
        t_datetime = np.array([])
        Nt = len(mfile['SerYear'])

        for nn in np.arange(Nt):
            # Clunky printing due to clunky handling in Jupyter..
            if nn/5 == nn//5:
                print('Calculating time.. (%.0f %%)\r'%((nn+1)/Nt*100), end = '')

            t_ = dt.datetime(int(2000.0) + mfile[yr][nn][0], mfile[mo][nn][0], 
                             mfile[da][nn][0], mfile[hr][nn][0], 
                             mfile[mi][nn][0], mfile[se][nn][0])
            t_datetime = np.append(t_datetime, t_)
            t_mpl = np.append(t_mpl, date2num(t_))
        
        self.t_mpl = t_mpl
        self.t_datetime = t_datetime
        self.units['t_mpl'] = '(Days since 1970-1-1)'
        self.var_desc['t_mpl'] = ('Calculated matplotlib time '
                                    '(days since 1970-1-1)')
        self.var_desc['t_datetime'] = ('Time stamp (matplotlib '
                                        'datetime object)')
        print('Calculating time: -> done.          ')

    #######################################################################
    # Generic functions for masking / chopping / rotating / applying
    # offsets etc
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
        uvc_rot = uvc_orig * np.exp(1j * ang_deg)
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
    # Functions building new data arrays based on existing ones
    #######################################################################

    #######################################################################

    def _calculate_tilt(self):
        """
        Calculate tilt from pitch/roll components.
        See Mantovanelli1 et al 2014 and Woodgate et al 2011.
        """
        P_rad = self.pitch / 180 * np.pi
        R_rad = self.roll / 180 * np.pi
        T_rad = np.arccos(np.sqrt(1 - np.sin(P_rad) ** 2 
                          - np.sin(R_rad) ** 2))
        self.tilt = T_rad * 180 / np.pi
        self.units['tilt'] = 'degrees'
        self.var_desc['tilt'] = 'Tilt (deg)'

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
            raise Exception('Cannot compute shear after we have removed'
                ' individual rows (bins). Either compute shear *before*'
                ' removing bins, or compute the shear yourself taking'
                ' varying depth interval into account.')
        

        self.shu[1:-1] = (0.01 * (self.u[2:] - self.u[:-2]) 
                        / (2 * self.binsize))
        self.shv[1:-1] = (0.01 * (self.v[2:] - self.v[:-2]) 
                        / (2 * self.binsize))
        self.s2 = self.shu ** 2 + self.shv ** 2
        self.units['shu'], self.units['shv'], self.units['s2'] = (
                '1/s', '1/s', '1/s2')

        self._record_proc_step('shear', '')





    #######################################################################
    # PRINTING/HISTORY/DOC FUNCTIONS
    #######################################################################

    def print_proc(self, return_str=False):
        """
        Print a processing history showing the steps applied to the dataset.
        """
        proc_str = ''
        for key in self._proc_dict.keys():
            type_str = self._proc_dict[key][0]
            desc_str = (variables_dictionary.proc_str_dict_rdi[type_str] 
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

    #######################################################################
    # EXPORTING
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
                outvars = variables_dictionary.outvars_rdi_sparse
            else:
                outvars = variables_dictionary.outvars_rdi_full

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
            outvars = variables_dictionary.outvars_rdi_full.copy()
        else:
            outvars = variables_dictionary.outvars_rdi_sparse.copy()

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

    def to_netcdf(self, file_name, sparse = True):
        '''
        Save to a NetCDF file (TBD).
        '''
        print('NetCDF export not implemented yet..')
        pass

    
#######################################################################

def example():
    '''
    Print a quick example of processing steps.
    '''
    
    example_str = '''\
    EXAMPLE CALL FOR PROCESSING:
    ## X-mark: not implemented yet.
    

    d = rdi_adcp.RdiObj(dloc + fn.mat)     # Load from matfile (produced from .000 file using WinADCP).

    d.print_system_info()                  # Print some system parameters (instrument configuration).

    d.remove_ship_time()                   # Chop away deployment and recovery.
    d.apply_depth_offset(3.2)              # Adjust transducer and bin depths 3.2 m *downwards*. 
    
    d.reject_surface_bins()                # Reject bins where the median depth is above the surface.

    d.reject_rows([0])                     # Reject row 0 (bin nearest to transducer). Will prompt y/n.

    d.set_latlon(lat, lon)                 # Set lat and lon

    d.correct_magdec()                     # Correct current direction for magnetic declination.

    --------------------------------------------------------------------------------------------------------
    ## Masking based on criteria (apply the relevant ones and modify the criteria) ##
    ## (Masks will end up as NaNs when exporting to matlab)
    
    d.mask_umax(thr_uvamp = 100)           # Mask entries where u or v exceed 100 cm/s amplitude.
    d.mask_surf_sidelobe() (X)              # Mask entries falling within the estimated range of sidelobe
                                           # interference of the surface.
    d.mask_pg(thr_pg = 75)                 # Mask entries where percent good (PG1 + PG4) is below 75%.
    d.mask_amp_a(thr_amp = 64)             # Mask entries where beam average amplitude is less than 64 db.
    d.mask_errvel(thr_errvel = 50)         # Mask entries where error velocity is greater than 50 cm/s.
    d.mask_ca(thr_cor = 30)                # Mask entries where the mean beam correlation in two or more 
                                           # beams is below 45 counts.
    d.mask_ca_mean(thr_cor = 30)           # Mask entries where the mean beam correlation is below 45 counts.
    d.mask_w(thr_w = 30)                   # Mask entries where the mean vertical is below 30 cm/s.
    d.mask_amp_jump(max_amp_increase = 30) # Masking entries where the beam amplitude of any beam has a jump 
                                           # exceeding 30 dB. 
    d.mask_amp_jump(max_amp_increase = 30, # Same, but also masks all entries *above* such jumps.
                    mask_above = True)
    d.mask_tilt(thr_tilt = 20)             # Mask entries where tile exceeds 20 degrees.
    --------------------------------------------------------------------------------------------------------

    d.print_summary()                      # Print a short summary of the dataset including mean velocities.

    d.calculate_shear()                    # Adds shear (s2, shu, shv)
        
    d.reject_rows(masked_max = 50)         # Reject rows with less than 50% valid (unmasked) entries.

    d.print_summary()                      # Print string showing a summary of the dataset (time/depth ranges,
                                           # mean velocities, etc).
    d.print_proc()                         # Print a processing history listing the individual steps applied  
                                           # to the dataset.

    b = d.to_dict()                        # Export to python Bunch.
    d.to_matfile('out_dir/fn.mat',         # Save as matfile with the typically most important parameters 
                 sparse = True)            # (t, depth, u, v, ..).
    d.to_matfile('out_dir/fn_full.mat')    # Save as pickled python dictionary (all parameters).
    d.to_pickle('out_dir/fn.p')            # Save as pickled python dictionary (all parameters).
    '''

    print(example_str)