
### This file uses a shapefile (GIS) to create a list of gridIDs from a large domain
### These gridIDs are used by kiloCraft to create user-defined domain/surfdata/forcing for ELM simulations 

import argparse
import geopandas as gpd
import xarray as xr
import numpy as np
import netCDF4 as nc


from datetime import datetime

def shape2grid(shapefile_path, aoi_name):
    # Get current date
    current_date = datetime.now()
    # Format date to mmddyyyy
    formatted_date = current_date.strftime('%y%m%d')

    # Load the shapefile
    shape = gpd.read_file(shapefile_path)

    print('CRS', shape.crs)
    # Define the target CRS (EPSG:4326, which is the ELM CRS)
    elm_crs = "EPSG:4326"

    # Check and convert the CRS if necessary
    if shape.crs != elm_crs:
        print(f"Shapefile CRS is {shape.crs}. Converting to {elm_crs}...")
        shape = shape.to_crs(elm_crs)
    else:
        print(f"Shapefile is already in {elm_crs}.")

    # save to the gridID  file
    AOI_gridID = f"{aoi_name}_gridID.c{formatted_date}.nc"
    dst = nc.Dataset(AOI_gridID, 'w', format='NETCDF3_64BIT')

    # Load the NetCDF file
    netcdf_file_path = 'domain.lnd.TES_SE.4km.1d.c240827.nc'  # replace with the actual path to your NetCDF file
    ds = xr.open_dataset(netcdf_file_path)

    # Extract the xc and yc coordinates and the gridID
    xc = ds['xc'].values.squeeze()
    yc = ds['yc'].values.squeeze()
    gridID = ds['gridID'].values.squeeze()

    # Create a GeoDataFrame for the grid cells
    grid_cells = gpd.GeoDataFrame({
        'gridID': gridID,
        'geometry': gpd.points_from_xy(xc, yc)
    })

    # Set the same coordinate reference system (CRS) as the shapefile
    grid_cells.set_crs(epsg=4326, inplace=True)

    # Check which grid cells are within the TVA shape
    grid_cells_within_AOI = grid_cells[grid_cells['geometry'].within(shape.unary_union)]

    # Get the list of gridIDs that are inside Tennessee
    grid_ids_within_AOI = grid_cells_within_AOI['gridID'].values

    # Print the resulting grid IDs
    print("Grid IDs within " + aoi_name +":", grid_ids_within_AOI)

    # create the gridIDs, lon, and lat variable
    ni_dim = dst.createDimension('ni', grid_ids_within_AOI.size)
    nj_dim = dst.createDimension('nj', 1)

    gridID_var = dst.createVariable('gridID', np.int32, ('nj','ni'), zlib=True, complevel=5)
    gridID_var.long_name = 'gridId in the TESSFA2 domain'
    gridID_var.decription = "start from #0 at the upper left corner of the domain, covering all land and ocean gridcells" 
    dst.variables['gridID'][...] = grid_ids_within_AOI
    dst.title = aoi_name +' land gridcells in the TESSFA2 domain'

    dst.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process a shapefile to find gridIDs within in the AOI region.')
    parser.add_argument('shapefile_path', type=str, help='shapefile witht the full path')
    parser.add_argument('aoi_name', type=str, help='Area of interest name')

    args = parser.parse_args()
    
    shape2grid(args.shapefile_path, args.aoi_name)





