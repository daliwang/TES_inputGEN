
# utility to generate 1D and 2D domain for TESSFA

import os, sys 
import math
import netCDF4 as nc
import numpy as np
from itertools import cycle
from time import process_time
from pyproj import Proj
from pyproj import Transformer
from pyproj import CRS

from datetime import datetime

def calculate_area(xv0, yv0, xv1, yv1, xv2, yv2, xv3, yv3):

    '''#Proj4: +proj=lcc +lon_0=-100 +lat_1=25 +lat_2=60 +k=1 +x_0=0 +y_0=0 +R=6378137 +f=298.257223563 +units=m  +no_defs
    geoxy_proj_str = "+proj=lcc +lon_0=-100 +lat_0=42.5 +lat_1=25 +lat_2=60 +x_0=0 +y_0=0 +R=6378137 +f=298.257223563 +units=m +no_defs"
    geoxyProj = CRS.from_proj4(geoxy_proj_str)
    # EPSG: 4326
    # Proj4: +proj=longlat +datum=WGS84 +no_defs
    lonlatProj = CRS.from_epsg(4326) # in lon/lat coordinates
    Txy2lonlat = Transformer.from_proj(geoxyProj, lonlatProj, always_xy=True)
    Tlonlat2xy = Transformer.from_proj(lonlatProj, geoxyProj, always_xy=True)'''
    
    # Step 1: Calculate the average value (Lat1) of the latitudes of two upper corners
    Lat1 = (yv0 + yv1) / 2

    # Step 2: Calculate the average value (Lat2) of the latitudes of two lower corners
    Lat2 = (yv2 + yv3) / 2

    # Step 3: Calculate the dLat = Lat1-Lat2 in degree 
    dLat = Lat1 - Lat2

    # Step 4: Calculate the average value (Lon1) of the longitude of two left corners
    Lon1 = (xv0 + xv3) / 2

    # Step 5: Calculate the average value (Lon2) of the longitude of two right corners
    Lon2 = (xv1 + xv2) / 2

    # Step 6: Calculate the dLon = Lon1-Lon2 in degree
    dLon = Lon1 - Lon2

    # Step 7: Calculate the area
    area_arcad2 = abs((dLon / 180) * (dLat / 180))

    return area_arcad2


def domain_save_1dTES(output_path, domain_name, total_rows, total_cols, data, lon, lat, XC, YC, fill_value, valid_mask):

    # This function uses (lon, lat) to create a domain for TES project
    # total_rows and total_cols are rows and cols from 2D Forcing dataset
    # lon, lat are 1D arrays
    # XC, YC are grid meshes that created using lon and Lat

    # Get current date
    current_date = datetime.now()
    # Format date to mmddyyyy
    formatted_date = current_date.strftime('%y%m%d')
    '''
    #Proj4: +proj=lcc +lon_0=-100 +lat_1=25 +lat_2=60 +k=1 +x_0=0 +y_0=0 +R=6378137 +f=298.257223563 +units=m  +no_defs
    geoxy_proj_str = "+proj=lcc +lon_0=-100 +lat_0=42.5 +lat_1=25 +lat_2=60 +x_0=0 +y_0=0 +R=6378137 +f=298.257223563 +units=m +no_defs"
    geoxyProj = CRS.from_proj4(geoxy_proj_str)
    # EPSG: 4326
    # Proj4: +proj=longlat +datum=WGS84 +no_defs
    lonlatProj = CRS.from_epsg(4326) # in lon/lat coordinates
    Txy2lonlat = Transformer.from_proj(geoxyProj, lonlatProj, always_xy=True)
    Tlonlat2xy = Transformer.from_proj(lonlatProj, geoxyProj, always_xy=True)
    XC_LCC,YC_LCC = Tlonlat2xy.transform(XC,YC)

    '''
    # add the gridcell IDs. and its x/y indices (will be used to create 2D output)
    total_gridcells = total_rows * total_cols
    grid_ids = np.linspace(0, total_gridcells-1, total_gridcells, dtype=int)
    grid_ids = grid_ids.reshape(total_rows,total_cols)
    grid_xids = np.indices(grid_ids.shape)[1]
    grid_yids = np.indices(grid_ids.shape)[0]

    print("shapes of gridid arrays", grid_ids.shape, grid_xids.shape, grid_yids.shape)
    masked = np.where(valid_mask)
    print("Number of non-NaN and non-fillvalue cells:", masked[0].size)
    print("first elemnt of masked is ", masked[0])
    #create land gridcell mask, area, and landfrac (otherwise 0)

    landmask = np.zeros_like(valid_mask, dtype=int)
    landmask[valid_mask] = 1
    landfrac = landmask.astype(float)*1.0
    print("Number of land cells (sum of landmask):", np.sum(landmask))
    print("Shape of landmask:", landmask.shape)
    
    # area in arcrad2
    '''area_km2 = 16.0 # 4km by 4km resolution
    offset = 2000.0 # distance (m) to the center of the grid
    side_km = math.sqrt(float(area_km2))
    lat[lat==90.0]=lat[lat==90.0]-0.00001
    lat[lat==-90.0]=lat[lat==-90.0]-0.00001
    kmratio_lon2lat = np.cos(np.radians(lat))
    re_km = 6371.22
    yscalar = side_km/(math.pi*re_km/180.0)
    xscalar = side_km/(math.pi*re_km/180.0*kmratio_lon2lat)
    area = xscalar*yscalar'''

    #  need to use the lon-lat bounding box to calculate the area
    #  Assume the resolution is fixed for all the cells
    lon_offset = abs((lon[0]-lon[1])/2)
    lat_offset = abs((lat[0]-lat[1])/2)
    
    # 2d --> 1d, with masked only
    grid_id_arr = grid_ids[masked]
    mask_arr = landmask[masked]
    landfrac_arr = landfrac[masked]
    grid_xids_arr = grid_xids[masked]
    grid_yids_arr = grid_yids[masked]
    XC_arr = XC[masked]
    YC_arr = YC[masked]
    #XC_LCC_arr = XC_LCC[masked]
    #YC_LCC_arr = YC_LCC[masked]
    
    file_name = output_path  + 'domain.lnd.' + domain_name + '.4km.1d.c'+ formatted_date +'.nc'
    print("The domain file is " + file_name)

    # Open a new NetCDF file to write the domain information. For format, you can choose from
    # 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
    w_nc_fid = nc.Dataset(file_name, 'w', format='NETCDF3_64BIT_DATA')
    w_nc_fid.title = '1D domain file for ' + domain_name

    # Create new dimensions of new coordinate system
    x_dim = w_nc_fid.createDimension('x', total_cols)
    y_dim = w_nc_fid.createDimension('y', total_rows)
    dst_var = w_nc_fid.createVariable('x', np.double, ('x'))
    dst_var.units = "degree"
    dst_var.long_name = "x coordinate of projection"
    dst_var.standard_name = "projection_x_coordinate"
    w_nc_fid['x'][...] = np.copy(XC[0,:])
    dst_var = w_nc_fid.createVariable('y', np.double, ('y'))
    dst_var.units = "degree"
    dst_var.long_name = "y coordinate of projection"
    dst_var.standard_name = "projection_y_coordinate"
    w_nc_fid['y'][...] = np.copy(YC[:,0])

    # create the gridIDs, lon, and lat variable
    ni_dim = w_nc_fid.createDimension('ni', grid_id_arr.size)
    nj_dim = w_nc_fid.createDimension('nj', 1)
    nv_dim = w_nc_fid.createDimension('nv', 4)

    w_nc_var = w_nc_fid.createVariable('gridID', np.int32, ('nj','ni'))
    w_nc_var.long_name = 'gridId in the TESSFA domain'
    w_nc_var.decription = "start from #0 at the upper left corner of the domain, covering all land and ocean gridcells" 
    w_nc_fid.variables['gridID'][...] = grid_id_arr

    w_nc_var = w_nc_fid.createVariable('gridXID', np.int32, ('nj','ni'))
    w_nc_var.long_name = 'gridId x in the TESSFA domain'
    w_nc_var.decription = "start from #0 at the upper left corner and from west to east of the domain, with gridID=gridXID+gridYID*x_dim" 
    w_nc_fid.variables['gridXID'][...] = grid_xids_arr

    w_nc_var = w_nc_fid.createVariable('gridYID', np.int32, ('nj','ni'))
    w_nc_var.long_name = 'gridId y in the TESSFA domain'
    w_nc_var.decription = "start from #0 at the upper left corner and from north to south of the domain, with gridID=gridXID+gridYID*y_dim" 
    w_nc_fid.variables['gridYID'][...] = grid_yids_arr
    

    w_nc_var = w_nc_fid.createVariable('xc', np.double, ('nj','ni'))
    w_nc_var.long_name = 'longitude of land gridcell center (GCS_WGS_84), increasing from west to east'
    w_nc_var.units = "degrees_east"
    #w_nc_var.bounds = "xv" ;    
    w_nc_fid.variables['xc'][...] = XC_arr
        
    w_nc_var = w_nc_fid.createVariable('yc', np.double, ('nj','ni'))
    w_nc_var.long_name = 'latitude of land gridcell center (GCS_WGS_84), decreasing from north to south'
    w_nc_var.units = "degrees_north"        
    w_nc_fid.variables['yc'][...] = YC_arr
        
    '''
    # create the XC, YC variable
    w_nc_var = w_nc_fid.createVariable('xc_LCC', np.double, ('nj','ni'))
    w_nc_var.long_name = 'X of land gridcell center (Lambert Conformal Conic), increasing from west to east'
    w_nc_var.units = "m"
    #w_nc_var.bounds = "xv" ;    
    w_nc_fid.variables['xc_LCC'][...] = XC_LCC_arr
        
    w_nc_var = w_nc_fid.createVariable('yc_LCC', np.double, ('nj','ni'))
    w_nc_var.long_name = 'Y of land gridcell center (Lambert Conformal Conic), decreasing from north to south'
    w_nc_var.units = "m"        
    w_nc_fid.variables['yc_LCC'][...] = YC_LCC_arr
    '''
    #
    w_nc_var = w_nc_fid.createVariable('xv', np.double, ('nv','nj','ni'))
    w_nc_var.long_name = 'longitude of land gridcell verticles (GCS_WGS_84), increasing from west to east'
    w_nc_var.units = "degrees_east"
    w_nc_var = w_nc_fid.createVariable('yv', np.double, ('nv','nj','ni'))
    w_nc_var.long_name = 'latitude of land gridcell verticles (GCS_WGS_84), decreasing from north to south'
    w_nc_var.units = "degrees_north"

    xv0,yv0 = (XC_arr-lon_offset,YC_arr + lat_offset)
    w_nc_fid.variables['xv'][0,...] = xv0
    w_nc_fid.variables['yv'][0,...] = yv0
    xv1,yv1 = (XC_arr+ lon_offset,YC_arr+ lat_offset)
    w_nc_fid.variables['xv'][1,...] = xv1
    w_nc_fid.variables['yv'][1,...] = yv1
    xv2,yv2 = (XC_arr+ lon_offset,YC_arr- lat_offset)
    w_nc_fid.variables['xv'][2,...] = xv2
    w_nc_fid.variables['yv'][2,...] = yv2
    xv3,yv3 = (XC_arr- lon_offset,YC_arr- lat_offset)
    w_nc_fid.variables['xv'][3,...] = xv3
    w_nc_fid.variables['yv'][3,...] = yv3

    area_arcad2  = calculate_area(xv0, yv0, xv1, yv1, xv2, yv2, xv3, yv3)

    w_nc_var = w_nc_fid.createVariable('area', np.double, ('nj','ni'))
    w_nc_var.long_name = 'Area of land gridcells'
    w_nc_var.coordinate = 'xc yc' 
    w_nc_var.units = "arcad^2"
    w_nc_fid.variables['area'][...] = area_arcad2

    '''w_nc_var = w_nc_fid.createVariable('area_LCC', np.double, ('nj','ni'))
    w_nc_var.long_name = 'Area of land gridcells (Lambert Conformal Conic)'
    w_nc_var.coordinate = 'xc_LCC yc_LCC' 
    w_nc_var.units = "km^2"
    w_nc_fid.variables['area_LCC'][...] = area_km2'''

    w_nc_var = w_nc_fid.createVariable('mask', np.int32, ('nj','ni'))
    w_nc_var.long_name = 'mask of land gridcells (1 means land)'
    w_nc_var.units = "unitless"        
    w_nc_fid.variables['mask'][...] = mask_arr

    w_nc_var = w_nc_fid.createVariable('frac', np.double, ('nj','ni'))
    w_nc_var.long_name = 'fraction of land gridcell that is active'
    w_nc_var.coordinate = 'xc yc' 
    w_nc_var.units = "unitless" ;        
    w_nc_fid.variables['frac'][...] = landfrac_arr

    '''w_nc_var = w_nc_fid.createVariable('lambert_conformal_conic', np.short)
    w_nc_var.grid_mapping_name = "lambert_conformal_conic"
    w_nc_var.longitude_of_central_meridian = -100. 
    w_nc_var.latitude_of_projection_origin = 42.5         
    w_nc_var.false_easting = 0. 
    w_nc_var.false_northing = 0.
    w_nc_var.standard_parallel = 25., 60.          
    w_nc_var.semi_major_axis = 6378137.
    w_nc_var.inverse_flattening = 298.257223563'''

    w_nc_fid.close()  # close the new file  


def domain_save_2dTES(output_path, domain_name, total_rows, total_cols, data, lon, lat, XC, YC, fill_value, valid_mask):

    # Get current date
    current_date = datetime.now()
    # Format date to mmddyyyy
    formatted_date = current_date.strftime('%y%m%d')

    # add the gridcell IDs.
    total_gridcells = total_rows * total_cols
    grid_ids = np.linspace(0, total_gridcells-1, total_gridcells, dtype=int)
    grid_ids = grid_ids.reshape(total_rows,total_cols)
    grid_xids = np.indices(grid_ids.shape)[1]
    grid_yids = np.indices(grid_ids.shape)[0]

    #create land gridcell mask, area, and landfrac (otherwise 0)
    fill_value = -9999.0  # or get from attrs
    #mask = np.where((~np.isnan(data[0])) & (data[0] != fill_value), 1, 0)
    #mask = np.where(~np.isnan(data[0]), 1, 0)
    mask = np.zeros_like(valid_mask, dtype=int)
    mask[valid_mask] = 1

    #  need to use the lon-lat bounding box to calculate the area
    lon_offset = abs((lon[0]-lon[1])/2)
    lat_offset = abs((lat[0]-lat[1])/2)

    landfrac = mask.astype(float)*1.0
    file_name = output_path  + 'domain.lnd.' + domain_name + '.4km.2d.c'+ formatted_date +'.nc'
    print("The 2D domain file is " + file_name)

    # Open a new NetCDF file to write the domain information. For format, you can choose from
    # 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
    w_nc_fid = nc.Dataset(file_name, 'w', format='NETCDF3_64BIT_DATA')
    w_nc_fid.title = '2D domain file for ' + domain_name

    # create the gridIDs, lon, and lat variable
    x_dim = w_nc_fid.createDimension('ni', total_cols)
    y_dim = w_nc_fid.createDimension('nj', total_rows)
    v_dim = w_nc_fid.createDimension('nv', 4)

    # Create new dimensions of new coordinate system
    dst_var = w_nc_fid.createVariable('ni', np.double, ('ni'))
    dst_var.units = "degree"
    dst_var.long_name = "x coordinate of projection"
    dst_var.standard_name = "projection_x_coordinate"
    w_nc_fid['ni'][...] = np.copy(XC[0,:])
    dst_var = w_nc_fid.createVariable('nj', np.double, ('nj'))
    dst_var.units = "m"
    dst_var.long_name = "y coordinate of projection"
    dst_var.standard_name = "projection_y_coordinate"
    w_nc_fid['nj'][...] = np.copy(YC[:,0])

    w_nc_var = w_nc_fid.createVariable('gridID', np.int32, ('nj','ni'))
    w_nc_var.long_name = 'gridId in the TESSFA domain'
    w_nc_var.decription = "start from #0 at the upper left corner of the domain and row by row, covering all land and ocean gridcells. " \
                          "So gridID=xid+yid*ni" 
    w_nc_fid.variables['gridID'][:] = grid_ids 

    w_nc_var = w_nc_fid.createVariable('xc', np.double, ('nj','ni'))
    w_nc_var.long_name = 'longitude of land gridcell center (GCS_WGS_84), increasing from west to east'
    w_nc_var.units = "degrees_east"
    #w_nc_var.bounds = "xv" ;    
    w_nc_fid.variables['xc'][:] = XC
        
    w_nc_var = w_nc_fid.createVariable('yc', np.double, ('nj','ni'))
    w_nc_var.long_name = 'latitude of land gridcell center (GCS_WGS_84), decreasing from north to south'
    w_nc_var.units = "degrees_north"        
    w_nc_fid.variables['yc'][:] = YC

    w_nc_var = w_nc_fid.createVariable('xv', np.double, ('nv','nj','ni'))
    w_nc_var.long_name = 'longitude of land gridcell verticles (GCS_WGS_84), increasing from west to east'
    w_nc_var.units = "degrees_east"
    w_nc_var = w_nc_fid.createVariable('yv', np.double, ('nv','nj','ni'))
    w_nc_var.long_name = 'latitude of land gridcell verticles (GCS_WGS_84), decreasing from north to south'
    w_nc_var.units = "degrees_north"

    # calculate the bounding box and then the area in arcrad^2 and km2
    xv0,yv0 = (XC-lon_offset,YC + lat_offset)  # upper-left corner
    w_nc_fid.variables['xv'][0,...] = xv0
    w_nc_fid.variables['yv'][0,...] = yv0
    xv1,yv1 = (XC+ lon_offset,YC+ lat_offset)  # upper-right corner
    w_nc_fid.variables['xv'][1,...] = xv1
    w_nc_fid.variables['yv'][1,...] = yv1
    xv2,yv2 = (XC+ lon_offset,YC- lat_offset)  # lower-right corner
    w_nc_fid.variables['xv'][2,...] = xv2
    w_nc_fid.variables['yv'][2,...] = yv2
    xv3,yv3 = (XC- lon_offset,YC- lat_offset)  # lower-left corner
    w_nc_fid.variables['xv'][3,...] = xv3
    w_nc_fid.variables['yv'][3,...] = yv3

    area_arcad2 = calculate_area(xv0, yv0, xv1, yv1, xv2, yv2, xv3, yv3)
    
    w_nc_var = w_nc_fid.createVariable('area', np.double, ('nj','ni'))
    w_nc_var.long_name = 'Area of land gridcells'
    w_nc_var.coordinate = 'xc yc' 
    w_nc_var.units = "arcad^2"
    w_nc_fid.variables['area'][...] = area_arcad2

    w_nc_var = w_nc_fid.createVariable('mask', np.int32, ('nj','ni'))
    w_nc_var.long_name = 'mask of land gridcells (1 means land)'
    w_nc_var.units = "unitless"        
    w_nc_fid.variables['mask'][:] = mask  

    w_nc_var = w_nc_fid.createVariable('frac', np.double, ('nj','ni'))
    w_nc_var.long_name = 'fraction of land gridcell that is active'
    w_nc_var.coordinate = 'xc yc' 
    w_nc_var.units = "unitless" ;        
    w_nc_fid.variables['frac'][:] = landfrac  

    w_nc_fid.close()  # close the new file  


def main():
    args = sys.argv[1:]
    
    if len(sys.argv) != 5  or sys.argv[1] == '--help':  # sys.argv includes the script name as the first argument
        print("Example use: python TES_domainGEN.py <input_path> <input_filename> <output_path> <domain_name>")
        print(" <input_path>: path to the input source data directory")
        print(" <input_filename>: the name of input data")
        print(" <output_path>:  path for the domain output")
        print(" <domain_name>:  output domain name")
        print(" The code uses TESSFA source data to generation 1D and 2D domain files")              
        exit(0)
    
    input_path = args[0]
    if not input_path.endswith("/"): input_path=input_path+'/'
    input_filename = args[1]
    output_path = args[2]
    if not output_path.endswith("/"): output_path=output_path+'/'
    domain_name = args[3]

    #input_path= './'
    #file_name = 'TES_BOT_domain1.nc'
    #output_path = input_path
    timesteps = 1
    var_name = 'FSDS'

    start = process_time()
    full_input = input_path + input_filename

    r_nc_fid = nc.Dataset(full_input, 'r', format='NETCDF3_64BIT_DATA')

    total_cols = r_nc_fid.dimensions['lon'].size
    total_rows = r_nc_fid.dimensions['lat'].size
    total_timesteps = r_nc_fid.dimensions['time'].size
    lon = r_nc_fid['lon']
    lat = r_nc_fid['lat']
    XC, YC = np.meshgrid(lon, lat)  # the array is (y,x) to match the mask
    #YX = np.array([Y, X])
    #print('total timesteps is :' + str(total_timesteps))

    var = r_nc_fid[var_name]  # NetCDF variable object
    data = var[0:timesteps, :, :]  # numpy array
    fill_value = var._FillValue if hasattr(var, '_FillValue') else -9999.0

    print("Unique values in data[0]:", np.unique(data[0]))
    print("Fill value used:", fill_value)

    #masked = np.where((~np.isnan(data[0])) & (~np.isclose(data[0], fill_value)))
    #print("Number of non-NaN and non-fillvalue cells:", masked[0].size)
    #print("first elemnt of masked is ", masked[0])
 
    if np.ma.is_masked(data[0]):
        valid_mask = ~data[0].mask
    else:
        valid_mask = (~np.isnan(data[0])) & (~np.isclose(data[0], fill_value))

    print("Number of non-NaN and non-fillvalue cells:", np.count_nonzero(valid_mask))
 
    print(lon.shape, lat.shape, XC.shape, YC.shape, XC.dtype) 
 
    end = process_time()
    print("Reading " + input_filename + " takes  {}".format(end-start))

    start = process_time()
    domain_save_2dTES(output_path, domain_name, total_rows, total_cols, data, lon, lat, XC, YC, fill_value, valid_mask)
    end = process_time()
    print("Saving 2D domain data takes {}".format(end-start))
    
    start = process_time() 
    domain_save_1dTES(output_path, domain_name, total_rows, total_cols, data, lon, lat, XC, YC, fill_value, valid_mask)
    end = process_time()

    print("Saving 1D domain data takes {}".format(end-start))

if __name__ == '__main__':
    main()
