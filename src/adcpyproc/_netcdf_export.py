
import netCDF4
import numpy as np
import os
from adcpyproc.netcdf_formatting._ncattrs_variables import _ncattrs_variables as _ncattrs_vars

def _create_netcdf(d, netcdf_dir, netcdf_file, 
                    attr_file='netcdf_global_attributes.py'):
    '''
    '''
    # If we are not in the specified netcdf output location *netcdf_dir*,
    # we will change directories to move there (moving back after unless we 
    # encounter errors)
    initial_dir = os.getcwd()
    try:
        os.chdir(netcdf_dir)
    except: 
        raise Exception('Unable to change directory to %s.'%netcdf_dir + 
                        'Check that it exists?')

    # Import the variable attribute dictionary
    try:
        from netcdf_global_attributes import ncattrs_global as ncattr_glob
    except:
        dirnow = os.getcwd()
        print(dirnow)
        raise Exception('Cannot find the file in %s in %s. '%(
                        attr_file, netcdf_dir) +
                        'Check that it exists?')


    # Create the ncfile


    # Check if file exists:
    if  os.path.isfile(netcdf_dir + netcdf_file):
        exists_continue = input('File %s%s exists. Overwrite? (y/n): '%(
            netcdf_dir, netcdf_file))
        if exists_continue is not 'y':
            raise Exception('Aborting since we did not want to override.')
        else:
            print('(Deleting previous file and creating new one)')
            os.remove(netcdf_dir + netcdf_file)

    N = netCDF4.Dataset(netcdf_dir + netcdf_file, 'w', format='NETCDF4')
    fill_value = -9999.0

    # Dimensions
    N.createDimension('TIME', d.Nt)
    N.createDimension('nDEPTH', d.Ndep)

    # Time
    N.createVariable('TIME', 'f4', ('TIME',))
    N['TIME'][:] = d.t_mpl#np.array([t_.strftime(t_fmt_nc) for t_ in d.t_datetime])
    N['TIME'].setncatts({
        'axis': 'T', 'long_name':'Time coordinate',
        'standard_name': 'time', 
        'units': 'days since 1970-01-01T00:00:00Z'})

    t_fmt_nc = '%Y-%m-%dT%h:%m:%sZ'
    N.setncattr('time_coverage_start', d.t_datetime[0].strftime(t_fmt_nc))
    N.setncattr('time_coverage_end', d.t_datetime[-1].strftime(t_fmt_nc))



    # Lon / lat
    if d.latlon_single:
        ll_dim = ()
    else:
        ll_dim = ('time',)
    N.createVariable('LATITUDE', 'f4', ll_dim)
    N.createVariable('LONGITUDE', 'f4', ll_dim)
    N['LATITUDE'][:] = d.lat
    N['LONGITUDE'][:] = d.lat

    # Data (u, v, depth)
    for varnm_ in _ncattrs_vars.keys():
        vdict_ = _ncattrs_vars[varnm_]
        N.createVariable(varnm_,
            getattr(d, varnm_).dtype, ('nDEPTH', 'TIME'),
            fill_value=fill_value)
        print(N[varnm_][:].shape, getattr(d, varnm_).shape)
        N[varnm_][:] = getattr(d, varnm_)
        N[varnm_].setncattr('valid_max', N[varnm_][:].max())
        N[varnm_].setncattr('valid_min', N[varnm_][:].min())
        N[varnm_].setncattr('units', getattr(d, 'units')[varnm_])


    # Global attributes (created based on the data)
    date_now = dt.datetime.now().strftime(t_fmt_nc)
    N.setncattr('date_created', date_now)

    if d.latlon_single:
        N.setncattr('geospatial_lon_min', d.lon)
        N.setncattr('geospatial_lon_max', d.lon)
        N.setncattr('geospatial_lat_min', d.lat)
        N.setncattr('geospatial_lat_max', d.lat)
    else:
        N.setncattr('geospatial_lon_min', d.lon.min())
        N.setncattr('geospatial_lon_max', d.lon.max())
        N.setncattr('geospatial_lat_min', d.lat.min())
        N.setncattr('geospatial_lat_max', d.lat.max())
    
    N.setncattr('geospatial_vertical_resolution', d.bin_size)


    # Global attributes (loaded from pre-made file)
    for attrnm in ncattr_glob:
        if ncattr_glob[attrnm] is not None:
            N.setncattr(attrnm, ncattr_glob[attrnm])


    # Global attributes (other)
    N.setncattr('Instrument serial number', str(SN))
    N.setncattr('processing_history', d.print_proc(return_str = True))
    N.setncattr('system_info', d.print_system_info(return_str = True))


    # cd back to the original directory
    os.chdir(initial_dir)
   
    N.close()
    print('Successfully converted to netCDF file: \n%s%s'%(
        netcdf_dir, netcdf_file))