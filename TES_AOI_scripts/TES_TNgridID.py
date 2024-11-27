import geopandas as gpd
import xarray as xr
import numpy as np
import netCDF4 as nc


from datetime import datetime

# Get current date
current_date = datetime.now()
# Format date to mmddyyyy
formatted_date = current_date.strftime('%y%m%d')

# Load the Tennessee shapefile
tn_shapefile_path = 'NRCSTATE_tn/state_nrcs_a_tn.shp'  # replace with the actual path to your shapefile
tn_shape = gpd.read_file(tn_shapefile_path)

AOI='TN'

# save to the gridID  file
AOI_gridID = str(AOI)+'_gridID.c'+ formatted_date + '.nc'
dst = nc.Dataset(AOI_gridID, 'w', format='NETCDF3_64BIT')

# Load the NetCDF file
netcdf_file_path = '../../entire_domain/domain_surfdata/domain.lnd.TES_SE.4km.1d.c240827.nc'  # replace with the actual path to your NetCDF file
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

# Check which grid cells are within the TN shape
grid_cells_within_tn = grid_cells[grid_cells['geometry'].within(tn_shape.unary_union)]

# Get the list of gridIDs that are inside Tennessee
grid_ids_within_tn = grid_cells_within_tn['gridID'].values

# Print the resulting grid IDs
print("Grid IDs within Tennessee:", grid_ids_within_tn
        )

# create the gridIDs, lon, and lat variable
ni_dim = dst.createDimension('ni', grid_ids_within_tn.size)
nj_dim = dst.createDimension('nj', 1)

gridID_var = dst.createVariable('gridID', np.int32, ('nj','ni'), zlib=True, complevel=5)
gridID_var.long_name = 'gridId in the TESSFA2 domain'
gridID_var.decription = "start from #0 at the upper left corner of the domain, covering all land and ocean gridcells" 
dst.variables['gridID'][...] = grid_ids_within_tn
dst.title = 'TN land gridcells in the TESSFA2 domain'

dst.close()


