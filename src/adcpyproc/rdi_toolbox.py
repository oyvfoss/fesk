'''
Various tools to be used on rdi class object ,
objects.

To do here:
- Printing basic stats
- Time series
- Basic profiles
- Direction plots
- Uv panels
- Variance ellipses
- Spectra (uv and shear)
- EOF modes?
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import num2date
from adcpyproc.misc_functions import closest



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
            var = getattr(d, key)[nn]
            var[var.mask] = np.nan
            VAR_ = np.convolve(var, 
                np.ones(wlen)/wlen, mode = 'valid')
            LP_UV[key][nn] = VAR_[::subsamp]

    return LP_UV


def _depavg_uv(A, dep, dep0=None, dep1=None, 
                 return_deprange = False):
    '''
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
    print(depsl)
    if return_deprange:
        return Amean, [dep[depsl][0], dep[depsl][-1]]
    else:
        return Amean



def uv_panels(d, wlen = None, subsamp = None, cmap_saturate = 99, 
              save_file = None):
    '''
    Quick panels of U and V.

    wlen: Window length for running mean (default is a window of ~1 day)
    subsamp: Subsampling frequency (for reducing the array sizes)
             Default is around 1/4 of wlen.
    cmap_saturate: Saturating the colormap at this percentile. 
    save_file: Enter a valid file path to save the figure.
    '''

    if wlen == None:
        wlen = int(np.round(1/d.d_t))
    if subsamp == None:
        subsamp = int(np.round(wlen/4))
    
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
    fig = plt.figure(figsize = (10, 7))
    colspan = 10
    dims = (2, colspan+1)
    ax0 = plt.subplot2grid(dims, (0, 0), colspan = colspan,)
    ax1 = plt.subplot2grid(dims, (1, 0), colspan = colspan,
            sharex = ax0, sharey = ax0)
    cax = plt.subplot2grid((5, colspan+1), (1, colspan), rowspan = 3)

    # Plot parameters
    vmax = np.percentile(np.concatenate(
        (np.abs(LP_UV['u'][~np.isnan(LP_UV['u'])]), 
         np.abs(LP_UV['v'][~np.isnan(LP_UV['v'])]))),
                         cmap_saturate)
    pcm_kws = {'vmax': vmax, 'vmin': -vmax, 'cmap': 'RdBu_r', 
                }

    # Plot
    C = ax0.pcolormesh(LP_UV['t'], d.dep_med, LP_UV['u'], **pcm_kws)
    ax1.pcolormesh(num2date(LP_UV['t']), d.dep_med, LP_UV['v'], **pcm_kws)

    # Colorbar
    plt.colorbar(C, cax = cax, label = 'cm/s', extend = 'both')

    # Axis and label details
    for axn in [ax0, ax1]:
        axn.set_ylabel('Median bin depth [m]')
        axn.set:facecolor('gray')

    ax0.invert_yaxis()
    yl = ax0.get_ylim()
    text_pos = yl[1] - 0.03*(yl[0] - yl[1])
    ax0.text(LP_UV['t'][0], text_pos, wlen_str, fontsize = 10, ha = 'left',
            clip_on = False)
    ax0.text(LP_UV['t'][-1], text_pos, subsamp_str, fontsize = 10, ha = 'right', 
            clip_on = False)

    ax0.set_title('U - Eastward velocity')
    ax1.set_title('V - Northward velocity')

    plt.show()

    # Save plot
    if save_file:
        fig.savefig(save_file, bbox_inches = 'tight')






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



def plot_ellipse(d, dep0 = None, lp_days = 5, dep1 = None, ax = None):
    '''
    '''
    print('depavg..')

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
    u = u_ - u_.mean()  
    v = v_ - v_.mean()


    thp, majax, minax = _uv_angle(u, v)

    if ax == None:
        fig, ax = plt.subplots()


    ax.set_aspect('equal')

    ax.plot(u, v, '.', ms = 1, color = 'Grey', alpha = 0.3, lw = 2,
            zorder = 0)
    ax.plot(u[-1], v[-1], marker= 'o',color = 'k', alpha =0.5,
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
    ax.plot(umaj , vmaj, '-k', lw = 2, label ='Maj axis')
    ax.plot(umin , vmin, '--k', lw = 2, label ='Min axis')

    ax.quiver(0, 0, UM, VM, 
        color = 'r', scale_units = 'xy', scale = 1, width = 0.03, 
        headlength = 2, headaxislength = 2, alpha = 0.6, 
        label = 'Mean',  edgecolor='k', linewidth = 0.6)

    ax.set_ylabel('v [cm s$^{-1}$]'); ax.set_xlabel('u [cm s$^{-1}$]')
    ax.legend(fontsize= 10, loc =3, handlelength = 1, ncol = 2) 

    ax.set_title('Mean currents, %.0f to %.0f m'%(drange[0], drange[1]))
    ax.grid()
    plt.show()