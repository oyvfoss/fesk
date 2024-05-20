'''
Various tools to be used on rdi class objects. Some of this is integrated
into rdi_adcp.

Some basic plots:
- UV panels.
- UV distributions with variance ellipses.
- Histograms (1d and by depth) of any parameter.

Some basic utility functions:
- Depth averaging
- Running means in time
- Variance ellipse calculations

May want to expand with basic spectra and EOF modes eventually.
Want to keep external dependencies to a minimum, though..

'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import num2date
from fesk.misc_functions import closest
from fesk import _rdi_defs


### PLOTTING

def uv_panels(d, wlen_days = None, subsamp_days = None, cmap_saturate = 99, 
              save_file = None, return_ax = False):
    '''
    Quick panels of U and V.

    wlen_days: Window length for running mean (default is a window of ~1 day)
    subsamp_days: Subsampling frequency (for reducing the array sizes)
             Default is around 1/4 of wlen.
    cmap_saturate: Saturating the colormap at this percentile. 
    save_file: Enter a valid file path to save the figure.
    '''

    if wlen_days == None:
        wlen = int(np.round(1/d.d_t))
    else:
        wlen = int(np.round(wlen_days/d.d_t))
    if subsamp_days == None:
        subsamp = int(np.round(wlen/4))
    else:
        subsamp = int(np.round(subsamp_days/4))

    # Get LPFed currents
    LP_UV = _runmean_uv(d, wlen, subsamp)
    
    # Strings displayng windows (hours or days)
    if LP_UV['d_t_window'] > 1.5:
        wlen_str = 'Running mean: %.1f days'%LP_UV['d_t_window']
        subsamp_str = 'Subsampled every %.1f days'%LP_UV['d_t_subsamp']
    else:
        wlen_str = 'Running mean: %.1f hours'%(24*LP_UV['d_t_window'])
        subsamp_str = 'Subsampled every %.1f hours'%(24*LP_UV['d_t_subsamp'])

    # Set up figure
    fig = plt.figure(figsize = (12, 7))
    colspan, rowspan = 5, 5
    dims = (rowspan*2 + 2, colspan +1, )
    ax0 = plt.subplot2grid(dims, (0, 0), colspan = colspan,
          rowspan = rowspan)
    ax1 = plt.subplot2grid(dims, (rowspan+1, 0), colspan = colspan,
            rowspan = rowspan, sharex = ax0, sharey = ax0)
    ax2 = plt.subplot2grid(dims, (0, colspan), rowspan = rowspan,
                            sharey = ax0)
    ax3 = plt.subplot2grid(dims, (rowspan+1, colspan), rowspan = rowspan, 
                            sharey = ax0, sharex = ax2)

    cax = plt.subplot2grid(dims, (rowspan*2+1, 1), colspan = colspan-2)


    # Plot parameters
    vmax = np.percentile(np.concatenate(
        (np.abs(LP_UV['u'][~np.isnan(LP_UV['u'])]), 
         np.abs(LP_UV['v'][~np.isnan(LP_UV['v'])]))),
                         cmap_saturate)
    pcm_kws = {'vmax': vmax, 'vmin': -vmax, 'cmap': 'RdBu_r', 
                'shading':'nearest'}

    # Plot panels
    C = ax0.pcolormesh(LP_UV['t'], d.dep_med, LP_UV['u'], **pcm_kws)
    ax1.pcolormesh(num2date(LP_UV['t']), d.dep_med, LP_UV['v'], **pcm_kws)
    yl = ax0.get_ylim()

    # Colorbar
    plt.colorbar(C, cax = cax, label = 'cm/s', extend = 'both', 
        orientation = 'horizontal')

    # Plot mean/STDs
    fillkws = {'color': 'g', 'alpha':0.3, 'label': 'STD'}
    linekws = {'color': 'k', 'alpha':0.9, 'linewidth':2, 
                'label':'mean', 'zorder':2}

    ustd = np.ma.masked_invalid(LP_UV['u']).std(axis=1) 
    umean = np.ma.masked_invalid(LP_UV['u']).mean(axis=1) 
    ax2.plot(umean, d.dep_med, **linekws)
    ax2.fill_betweenx(d.dep_med, umean-ustd, umean+ustd, **fillkws)

    vstd = np.ma.masked_invalid(LP_UV['v']).std(axis=1) 
    vmean = np.ma.masked_invalid(LP_UV['v']).mean(axis=1) 
    ax3.plot(vmean, d.dep_med, **linekws)
    ax3.fill_betweenx(d.dep_med, vmean-vstd, vmean+vstd, **fillkws)


    # Axis and label details
    for axn in [ax0, ax1]:
        axn.set_facecolor('gray')
        axn.set_ylabel('Median depth [m]')

    for axn in [ax2, ax3]:
        axn.grid()
        axn.set_xlabel('cm/s')
        axn.xaxis.tick_top()
        axn.xaxis.set_label_position('top')
        axn.yaxis.tick_right()
        axn.yaxis.set_label_position('right')
        axn.set_ylabel('Median depth [m]')
        axn.legend(handlelength = 1, ncol = 2, bbox_to_anchor = (1.3, 0))

    ax0.set_ylim(yl)
    ax0.invert_yaxis()

    text_pos = yl[0] + 0.03*(yl[0] - yl[1])
    ax0.text(LP_UV['t'][0], text_pos, wlen_str, fontsize = 10, ha = 'left',
            clip_on = False, fontstyle = 'italic')
    ax0.text(LP_UV['t'][-1], text_pos, subsamp_str, fontsize = 10, ha = 'right', 
            clip_on = False, fontstyle = 'italic')

    ax0.set_title('U - Eastward velocity')
    ax1.set_title('V - Northward velocity')

    plt.subplots_adjust(wspace = 0.5, hspace = 5)
    plt.show()

    # Save plot
    if save_file:
        fig.savefig(save_file, bbox_inches = 'tight')

    if return_ax:
        return np.array([ax0, ax1, ax2, ax3])



def plot_ellipse(d, dep0 = None, dep1 = None, lp_days = 5, ax = None, 
                return_ax = True):
    '''
    Plot of u and v components averaged between *dep0* and *dep1* and
    low pass filtered with a running mean of *lp_days*.

    Showing the mean current vector, the low-pass-filtered 
    and subsampled currents, and the semi-major and -minor axes
    of the variance ellipse. 
    '''

    # Depth average u, v
    u_, drange = _depavg_uv(d.u, d.dep_med, dep0=dep0, dep1=dep1, 
                  return_deprange = True)

    v_ = _depavg_uv(d.v, d.dep_med, dep0=dep0, dep1=dep1, 
                  return_deprange = False)

    # Mean vector
    UM, VM = u_.mean(), v_.mean()

    # LPFed 
    wlen = int(np.round(lp_days/d.d_t))
    ULP = np.convolve(u_, np.ones(wlen)/wlen,
                mode = 'valid')[::wlen]
    VLP = np.convolve(v_, np.ones(wlen)/wlen,
                mode = 'valid')[::wlen]

    # Demeaned

    thp, majax, minax = _uv_angle(u_ - u_.mean() ,  v_ - v_.mean())

    if ax == None:
        fig, ax = plt.subplots(figsize = (10, 10))

    ax.set_aspect('equal')

    ax.plot(u_, v_, '.', ms = 1, color = 'Grey', alpha = 0.3, lw = 2,
            zorder = 0)
    ax.plot(u_[-1], v_[-1], '.', ms= 1,color = 'k', alpha =0.5,
            lw = 2, label ='Full')

    ax.plot(ULP, VLP, '.', ms = 3, color = 'b', alpha = 0.5)
    ax.plot(ULP[0], VLP[0], '.', ms = 5, color = 'b', alpha = 0.5, 
        label = '%.1f-day means'%(lp_days), zorder = 0)

    vmaj = np.array([-majax*np.sin(thp), majax*np.sin(thp)])
    umaj = np.array([-majax*np.cos(thp), majax*np.cos(thp)])
    vmin = np.array([-minax*np.sin(thp+np.pi/2),
                        minax*np.sin(thp+np.pi/2)])
    umin = np.array([-minax*np.cos(thp+np.pi/2),
                        minax*np.cos(thp+np.pi/2)])
    ax.plot(UM + umaj , VM + vmaj, '-k', lw = 2, label ='Maj axis')
    ax.plot(UM + umin , VM + vmin, '--k', lw = 2, label ='Min axis')

    ax.quiver(0, 0, UM, VM, 
        color = 'r', scale_units = 'xy', scale = 1, width = 0.03, 
        headlength = 2, headaxislength = 2, alpha = 0.6, 
        label = 'Mean',  edgecolor='k', linewidth = 0.6)

    ax.set_ylabel('v [cm s$^{-1}$]'); ax.set_xlabel('u [cm s$^{-1}$]')
    ax.legend(fontsize= 10, loc =3, handlelength = 1, ncol = 2) 

    ax.set_title('Mean currents, %.0f to %.0f m'%(drange[0], drange[1]))
    ax.grid()
    plt.show()

    if return_ax:
        return ax


def hist(d, varnm, nbins = 50, line = None, return_ax = False):
    '''
    Show histogram(s) of a variable.
    
    For 1D parameters (tdepth, heading, tilt..): 
        Showing a simple histogram
    For 2D parameters (u, amp_a, ca1..)
        Simple histogram + histograms by depth
    
    d: RdiObj
    varnm: Short variable name ('u', 'amp1', 'tilt'..).
    nbins: Number of bins to use for histogram
    line: Plot a vertical line here. Useful for looking at 
          thresholds.
    return_ax: return the axis object
    
    '''
    vl_kws = {'ls':':', 'color': 'k'}

    if getattr(d, varnm).ndim == 1:
        fig, ax = plt.subplots()
        _histogram_param_1d(d, varnm, ax, nbins = nbins)
        if line:
            ax.axvline(line, **vl_kws)
        ax.set_title(varnm, fontsize = 15, fontweight = 'bold')

    elif getattr(d, varnm).ndim == 2:
        fig, ax = plt.subplots(1, 2, sharex = True, figsize = (12, 6))
        _histogram_param_1d(d, varnm, ax[0], nbins = nbins)
        _histogram_param_2d(d, varnm, ax[1], nbins = nbins)
        if line:
            ax[0].axvline(line, **vl_kws)
            ax[1].axvline(line, **vl_kws)
        ax[0].set_title(varnm, fontsize = 15, fontweight = 'bold')
        ax[1].set_title(varnm, fontsize = 15, fontweight = 'bold')

        plt.tight_layout(w_pad = 3)
    
    
    if return_ax:
        return ax



### UTILITY FUNCTIONS

def _runmean_uv(d, wlen, subsamp=1):
    '''
    Create a running mean and subsample to make it easier to do quick
    visualisations etc.

    d: adcpyproc rdi object.
    wlen: Window length for running mean
    subsamp: Subsampling frequency (for reducing the array sizes)
             Default is 1 (no subsampling)

    Returns a dictionary with (t, u, v).
    '''

    LP_UV = {}
    LP_UV['d_t_window'] = wlen * d.d_t
    LP_UV['d_t_subsamp'] = subsamp * d.d_t


    LP_UV['t'] = np.convolve(d.t_mpl, np.ones(wlen)/wlen,
                mode = 'valid')[::subsamp]

    LP_UV['u'] = np.ma.zeros((d.Ndep, len(LP_UV['t'])))
    LP_UV['u'].mask = True
    LP_UV['v'] = LP_UV['u'].copy()

    for nn in np.arange(d.Ndep):
        for key in ['u', 'v']:
            var = getattr(d, key)[nn].copy()
            var[var.mask] = np.nan
            VAR_ = np.convolve(var, 
                np.ones(wlen)/wlen, mode = 'valid')
            LP_UV[key][nn] = VAR_[::subsamp]

    return LP_UV


def _depavg_uv(A, dep, dep0=None, dep1=None, 
                 return_deprange = False):
    '''
    Calculate depth average currents.

    dep0, dep1: Approximate depth limits (will find nearest bin).
    '''
    if dep1:
        dep_ind0 = closest(dep, dep1)
    else:
        dep_ind0 = None
    if dep0:
        dep_ind1 = closest(dep, dep0)+1
    else:
        dep_ind1 = None

    depsl = slice(dep_ind0, dep_ind1)

    Amean = A[depsl].mean(axis = 0)

    if return_deprange:
        return Amean, [dep[depsl][0], dep[depsl][-1]]
    else:
        return Amean




def _uv_angle(u, v):
    '''
    Finds the principal angle angle in [-pi/2, pi/2] where the squares
    of the normal distance  to u,v are maximised. Ref Emery/Thompson
    pp 327.

    Also returns the standard deviation along the semimajor and
    semiminor axes.

    '''
    if u.mean()>1e-7 or v.mean()>1e-7:
        print('Mean of u and/or v is nonzero. Removing mean.')
        u = u - u.mean() ; v = v - v.mean()
    thp = 0.5* np.arctan2(2*np.mean(u*v),(np.mean(u**2)-\
                          np.mean(v**2))) #ET eq 4.3.23b
    uvcr = (u + 1j * v) * np.exp(-1j * thp)
    majax = np.std(uvcr.real)
    minax = np.std(uvcr.imag)

    return thp, majax, minax


### UV PLOT WITH CURRENTS ELLIPSE




def _histogram_param_1d(d, varnm, ax, nbins = 50):
    '''
    Histogram showing the distribution of a parameter - 1D or 2D.
    
    d: A RdiObj
    varnm: Name of the variable ('u', 'amp1', 'tdepth', etc..)
    ax: A figure axis
    nbins: Number of histogram bins
    '''

    VAR_ = getattr(d, varnm) 
    umask = d.u.mask
    
    # Not including invalid entries - give crazy values..
    if varnm in ['u', 'v', 'w', 'errvel']:
        umask = umask[np.abs(VAR_.data) < 1000]
        VAR_ = VAR_[np.abs(VAR_.data) < 1000]
        
    if getattr(d, varnm).ndim == 1:
        VAR_all = VAR_
    elif getattr(d, varnm).ndim == 2:
        VAR_all = VAR_.flatten()
        VAR_val = VAR_[~umask]

    N_all = len(VAR_all)
    col_1 = (1.0, 0.498, 0.055)
    col_2 = (0.122, 0.467, 0.705)

    # Histogram, all entries
    Hargs = {'density': False, }
    H_all, H_bins = np.histogram(VAR_all, bins = nbins, **Hargs)
    Hargs['bins']= H_bins
    H_width = np.ma.median(np.diff(H_bins))

    # Bar plot
    ax.bar(H_bins[:-1], 100*H_all/N_all, width = H_width, align = 'edge', 
            alpha = 0.4, color = col_1, label = 'All')
    
    if getattr(d, varnm).ndim == 2:
        N_val = len(VAR_val)

        # Histogram, valid (unmasked) entries
        H_val, H_val_bins = np.histogram(VAR_val, **Hargs)
        ax.bar(H_bins[:-1], 100*H_val/N_val, width = H_width, align = 'edge', 
            alpha = 0.4, color = col_2, label = 'Valid u')
        ax.legend(loc = 5)
 
    # Cumulative plot
    cumulative = np.concatenate([[0], np.cumsum(100*H_all/N_all)])
    twax = ax.twinx()
    twax.plot(H_bins, cumulative, 'k', clip_on = False)
    twax.set_ylim(0, 105)
    
    # Axis labels
    
    # x label: Long description
    xlab = d.var_desc[varnm]
    if 'From field' in xlab:
        xlab = xlab[:xlab.find('. From field')]
    ax.set_xlabel(xlab)
    
    ax.set_ylabel('Density per bin [%]')
    twax.set_ylabel('Cumulative density [%]')
    

    if getattr(d, varnm).ndim == 2:
        ax.legend(loc = 5)


def _histogram_param_2d(d, varnm, ax, nbins = 50):
    '''
    Histogram showing the distribution of a 2D parameter at each depth.
    
    Showing the distribution only witin "valid" bins - measurements
    where velocities have not been masked.
    
    d: A RdiObj
    varnm: Name of the variable ('u', 'amp1', 'tdepth', etc..)
    ax: A figure axis
    nbins: Number of histogram bins
    '''

    
    # Bins
    varmin = getattr(d, varnm)[~d.u.mask].min()
    varmax = getattr(d, varnm)[~d.u.mask].max()
    bins = varmin + np.arange(nbins+1)*(varmax-varmin)/nbins
    d_bins = np.diff(bins)[0]

    # Histogram at each depth
    H_bydep = np.ma.zeros((d.Ndep, nbins, ))
    H_bydep.mask = True

    for nn in np.arange(d.Ndep):
        VAR_nn = getattr(d, varnm)[nn]
        VAR_nn_valid = VAR_nn[~d.u.mask[nn]]
        H_bydep_, H_bins_ = np.histogram(VAR_nn_valid, bins = bins)
        N_all = len(VAR_nn_valid)
        H_bydep[nn] = H_bydep_*100/N_all
    
    
    # Colorbar
    C = ax.contourf(bins[:-1] + d_bins/2, 
            d.dep_med, H_bydep, 20, cmap = 'cmo.tempo')
    
    plt.colorbar(C, label = 'Density per bin [%]', 
                 orientation = 'horizontal', ax = ax)
    
    # Axis labels
    ax.set_ylabel('Median depth [m]')

    # x label: Long description
    xlab = d.var_desc[varnm]
    if 'From field' in xlab:
        xlab = xlab[:xlab.find('. From field')]
    ax.set_xlabel(xlab)
    
    ax.invert_yaxis()