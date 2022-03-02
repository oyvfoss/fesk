# Methods for reading Winadcp .mat files the rdi_adcp.RdiObj class

import numpy as np
import datetime as dt
from matplotlib.dates import num2date, date2num

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
            except:
                raise Warning('Issue with variable "%s"/"%s"' % (
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

    self.dep = (self.tdepth[np.newaxis, :]
                + sign * mfile['RDIBin1Mid'] + sign * mfile['RDIBinSize'] 
                * (mfile['SerBins'][0] - 1)[:, np.newaxis])
    self.units['dep'] = 'm'
    self.bin_size = mfile['RDIBinSize'][0][0]
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
    self.d_t = np.median(np.diff(t_mpl))
    
    print('Calculating time: -> done.          ')