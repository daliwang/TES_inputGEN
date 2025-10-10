import sys, os
import netCDF4 as nc
from scipy.interpolate import griddata
import numpy as np
import pandas as pd 
from time import process_time
#from memory_profiler import profile
import datetime

# Get current date
current_date = datetime.datetime.now()

# Format date as yymmdd
date_string = current_date.strftime('%y%m%d')

def main():
    '''
    args = sys.argv[1:]

    if  len(sys.argv) != 3 or sys.argv[1] == '--help':  # sys.argv includes the script name as the first argument
        print("Example use: python NA_surfdataGEN.py <input_path> <domain_name> <mask_file> <output_domain_type>")
        print("<input_path: path to the input surfdata>")
        print("<domain_name: unique name for the project or domain>")
        print("<mask_file: one forcing data>")        
        print(" The code generate TESSFA surfdata from 0.5x0.5 degree globla surfdata")
        print(" The code requires a 0.5x0.5 degree globla surfdata (surfdata.nc)")
        print(" The code requires a domain mask information (surfdata.nc)")        
        exit(0)

    input_path = args[0]
    if not input_path.endswith("/"): input_path=input_path+'/'
    domain_name = args[1]

    full_domain = input_path + domain_name
    mask_file = args[2]

    domain_type = args[3]
    '''

    # Get current date
    current_date = datetime.datetime.now()

    # Format date as yymmdd
    date_string = current_date.strftime('%y%m%d')

    domain_type = '1'

    if domain_type =='1':   domain_name = 'ACCESS/domain.lnd.TES_NORTHACCESS.4km.1d.c250521.nc'     
    if domain_type =='2':   domain_name = 'ACCESS/domain.lnd.TES_NORTHACCESS.4km.2d.c250521.nc' 
        
    mask_file = 'TES_NORTHACCESS_domain.nc'     
    output_file = domain_name.replace('domain.lnd', 'surfdata')

    print("the domain file is " + domain_name)
    print("the output file is " + output_file)
        
    # Only variables listed will be processed

    # nearest neighbor:"double" variables
    Variable_nearest = ['SLOPE', 'TOPO', 'PCT_GLACIER', 'PCT_LAKE', 'STD_ELEV']
    
    # nearest neighbor:"int" variables
    Variable_nearest += ['PFTDATA_MASK','SOIL_COLOR', 'SOIL_ORDER', 'abm']
    # nearest neighbor:"double" variables (added 11/07/22023)
    Variable_nearest += ['EF1_BTR', 'EF1_CRP', 'EF1_FDT', 'EF1_FET', 'EF1_GRS', 'EF1_SHR']
    # nearest neighbor: "double" variables (added 11/10/2023)
    Variable_nearest += ['PCT_SAND', 'PCT_CLAY','ORGANIC' ,'PCT_NAT_PFT', 
            'MONTHLY_LAI', 'MONTHLY_SAI' ,'MONTHLY_HEIGHT_TOP', 'MONTHLY_HEIGHT_BOT']

    # nearest neighbor:"int" variables (gridcell)
    Variable_urban_nearest = ['URBAN_REGION_ID']
    # nearest neighbor:"int" variables (numurbl, gridcell)
    Variable_urban_nearest += ['NLEV_IMPROAD' ]
    # nearest neighbor:"double" variables (numurbl, gridcell)
    Variable_urban_nearest += ['T_BUILDING_MAX', 'T_BUILDING_MIN',
            'WIND_HGT_CANYON','WTLUNIT_ROOF','WTROAD_PERV','THICK_ROOF',
            'THICK_WALL','PCT_URBAN','HT_ROOF','EM_IMPROAD','EM_PERROAD',
            'EM_ROOF','EM_WALL','CANYON_HWR']
    # nearest neighbor:"double" variables (nlevurb, numurbl, gridcell)
    Variable_urban_nearest += ['TK_IMPROAD','TK_ROOF','TK_WALL', 
                                     'CV_IMPROAD', 'CV_ROOF', 'CV_WALL']

    # nearest neighbor:"double" variables (numrad, numurbl, gridcell)
    Variable_urban_nearest += ['ALB_IMPROAD_DIF','ALB_IMPROAD_DIR','ALB_PERROAD_DIF',
            'ALB_PERROAD_DIR','ALB_ROOF_DIF', 'ALB_ROOF_DIR',
            'ALB_WALL_DIF', 'ALB_WALL_DIR']

    Variable_nearest += Variable_urban_nearest

    # linear interpolation of "double" variables
    Variable_linear = ['FMAX', 'Ws', 'ZWT0', 'binfl', 'gdp', 
                    'peatf', 'Ds', 'Dsmax', 'F0', 'LAKEDEPTH',
                   'LANDFRAC_PFT','P3', 'PCT_NATVEG', 'PCT_WETLAND', 
                    'SECONDARY_P', 'OCCLUDED_P', 'LABILE_P']

    # linear: "double" variables (added 11/07/22023)
    Variable_linear += ['APATITE_P', 'PCT_CROP']
    

    # Open the source file
    src = nc.Dataset('surfdata.nc', 'r')
    # points = [y,x] coordinates for src's grid
    src_resx = 0.5
    src_resy = 0.5
    src_lat = src.variables['LATIXY'][...]  
    src_lon = src.variables['LONGXY'][...]
    #src_x,src_y = Tlonlat2xy.transform(src_lon,src_lat)
    src_lon[src_lon<0.0]=360+src_lon[src_lon<0.0]
    
    # Create a new file      
    if os.path.isfile(output_file):
        # If it does, delete it
        os.remove(output_file)
        
    dst = nc.Dataset(output_file, 'w')

    # Copy dimensions
    for name, dimension in src.dimensions.items():
        dst.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))

    # Copy global attributes
    dst.setncatts(src.__dict__)

    # get the fine resolution data and the locations (lat, lon) from the forcing dataset

    # get the fine resolution data and the locations (lat, lon)
    #r_daymet = nc.Dataset('clmforc.Daymet4.1km.TBOT.2014-01.nc', 'r', format='NETCDF4')
    r_mask = nc.Dataset(mask_file, 'r', format='NETCDF4')    
    x_dim = r_mask['lon'][...]  # 1D x-axis
    y_dim = r_mask['lat'][...]  # 1D y-axis
    data = r_mask.variables['FSDS'][0,...]

    grid_ids = np.linspace(0, len(x_dim)*len(y_dim)-1, len(x_dim)*len(y_dim), dtype=int)
    grid_ids = grid_ids.reshape(data.shape)
    grid_xids = np.indices(grid_ids.shape)[1]
    grid_yids = np.indices(grid_ids.shape)[0]

    # setup the bool_mask and mesh of the TESSFA domain (lon, lat)
    if np.ma.is_masked(data):
        bool_mask = ~data.mask
    else:
        bool_mask = (~np.isnan(data)) & (~np.isclose(data, fill_value))
    grid_x, grid_y = np.meshgrid(x_dim,y_dim)
    print(data.shape,bool_mask.shape, grid_x.shape)

    grid_lat = np.copy(grid_y[bool_mask])
    grid_lon = np.copy(grid_x[bool_mask])
    grid_lon[grid_lon<0.0]=360+grid_lon[grid_lon<0.0]
    
    gridcells= len(grid_lon)

    print(data.shape,bool_mask.shape, grid_x.shape, grid_lon.shape, gridcells)
    #del grid_x, grid_y

    # make sure the xc, yc, and area are identical in both domain and surfdata 
    
    r_tessfa = nc.Dataset(domain_name, 'r', format='NETCDF4')
    LONGXY = r_tessfa.variables['xc'][...]  # longitude of gridcell center
    LATIXY = r_tessfa.variables['yc'][...]  # latitude of gridcell center
    AREA = r_tessfa.variables['area'][...]

    '''if domain_type == "2": 
        grid_xc = r_tessfa.variables['xc'][...]  # longitude of gridcell center
        grid_lon[grid_lon < 0.0] = 360 + grid_lon[grid_lon < 0.0]
    
        grid_lat = r_tessfa.variables['yc'][...]  # latitude of gridcell center
        
        gridID = r_tessfa.variables['gridID'][...]  

        print("shape of gridID", gridID.shape, r_tessfa.variables['gridID'].shape)
        
        AREA = r_tessfa.variables['area'][...]  # Area is in (lat, lon)
        y_dim = AREA.shape[0]  # 1D x-axis
        x_dim = AREA.shape[1]  # 1D y-axis

        print(" size and shape of grid_lat",  grid_lat.size, grid_lat.shape)'''
    
    # prepare the data source with points_in_daymet_land  (witn lon, lat)
    idxy = np.nonzero((src_lat<=np.amax(grid_lat)+src_resy) & (src_lat>=np.amin(grid_lat)-src_resy) & \
                      (src_lon<=np.amax(grid_lon)+src_resx) & (src_lon>=np.amin(grid_lon)-src_resx) )

    
    points_in_tessfa_land = {}
    points_in_tessfa_land[0] = src_lat[idxy]
    points_in_tessfa_land[1] = src_lon[idxy]
    points_in_tessfa_land[2] = idxy[0]
    points_in_tessfa_land[3] = idxy[1]
    land_points = len(points_in_tessfa_land[0])
    points=np.zeros((land_points, 2), dtype='double')
    points[:,0] = src_lat[idxy]
    points[:,1] = src_lon[idxy]

    # Create new dimensions for TES domain
    dst.createDimension('lon', x_dim.size)
    dst.createDimension('lat', y_dim.size)

    dst_var = dst.createVariable('lon', np.float64, ('lon'), zlib=True, complevel=5)
    dst_var.units = "degree"
    dst_var.long_name = "longitude coordinate"
    dst_var.standard_name = "x_project_coordinate"
    dst['lon'][...] = np.copy(x_dim)
    dst_var = dst.createVariable('lat', np.float64, ('lat'), zlib=True, complevel=5)
    dst_var.units = "degree"
    dst_var.long_name = "latitude coordinate"
    dst_var.standard_name = "projection_y_coordinate"
    dst['lat'][...] = np.copy(y_dim)
    
#########  the default format should be 2D 
    
    if (domain_type == '1'):
    
        dst_var = dst.createVariable('lon2D', np.float64, ('lat','lon'), zlib=True, complevel=5)
        dst_var.units = "degrees_east"
        dst_var.long_name = "longitude coordinate"
        dst_var.standard_name = "longitude"
        dst['lon2D'][...] = np.copy(grid_x)
        dst_var = dst.createVariable('lat2D', np.float64, ('lat','lon'), zlib=True, complevel=5)
        dst_var.units = "degrees_north"
        dst_var.long_name = "latitude coordinate"
        dst_var.standard_name = "latitude"
        dst['lat2D'][...] = np.copy(grid_y)

        dst.createDimension('gridcell', gridcells)
        
        dst_var = dst.createVariable('gridID', np.int32, ('gridcell'), zlib=True, complevel=5)
        dst_var.long_name = 'gridId in the TES domain'
        dst_var.decription = "start from #0 at the upper left corner of the domain, covering all land and ocean gridcells" 
        dst['gridID'][...] = np.copy(grid_ids[bool_mask]) 
        
        dst_var = dst.createVariable('gridXID', np.int32, ('gridcell'), zlib=True, complevel=5)
        dst_var.long_name = 'gridId x in the TESSFA domain'
        dst_var.decription = "start from #0 at the upper left corner and from west to east of the domain, with gridID=gridXID+gridYID*x_dim" 
        dst.variables['gridXID'][...] = np.copy(grid_xids[bool_mask])
    
        dst_var = dst.createVariable('gridYID', np.int32, ('gridcell'), zlib=True, complevel=5)
        dst_var.long_name = 'gridId y in the TESSFA domain'
        dst_var.decription = "start from #0 at the upper left corner and from north to south of the domain, with gridID=gridXID+gridYID*y_dim" 
        dst.variables['gridYID'][...] = np.copy(grid_yids[bool_mask])

    else:
        dst_var = dst.createVariable('gridID', np.int32, ('lat','lon'), zlib=True, complevel=5)
        dst_var.long_name = 'gridId in the TESSFA domain'
        dst_var.decription = "start from #0 at the upper left corner of the domain, covering all land and ocean gridcells"
        dst.variables['gridID'][...] = np.copy(grid_ids)

        dst_var = dst.createVariable('gridXID', np.int32, ('lat','lon'), zlib=True, complevel=5)
        dst_var.long_name = 'gridId x in the TESSFA domain'
        dst_var.decription = "start from #0 at the upper left corner and from west to east of the domain, with gridID=gridXID+gridYID*x_dim"
        dst.variables['gridXID'][...] = np.copy(grid_xids)

        dst_var = dst.createVariable('gridYID', np.int32, ('lat','lon'), zlib=True, complevel=5)
        dst_var.long_name = 'gridId y in the TESSFA domain'
        dst_var.decription = "start from #0 at the upper left corner and from north to south of the domain, with gridID=gridXID+gridYID*y_dim"
        dst.variables['gridYID'][...] = np.copy(grid_yids)

    count = 0 # record how may 2D layers have been processed 

    # Copy variables
    for name, variable in src.variables.items():
        start = process_time()
        print("Checking on varibale: "+ name + " dimensions: " + str(variable.dimensions))
        
        # Check if the last two dimensions are lsmlat and lsmlon
        if (variable.dimensions[-2:] == ('lsmlat', 'lsmlon')):
            # Determine the interpolation method
            if name in Variable_linear:
                iMethod = 'linear'
            else:
                iMethod = 'nearest'

            print("Working on varibale: "+ name + " dimensions: " + str(variable.dimensions))

            # create variables with the new dimensions

            if variable.datatype == np.int32:
                fill_value = -9999  # or any other value that you want to use to represent missing data
            else:
                fill_value = np.nan

            x = dst.createVariable(name, variable.datatype, variable.dimensions[:-2]+ ('gridcell',), fill_value = fill_value, zlib=True, complevel=5)
                
	    # Copy variable attributes
            dst[name].setncatts(src[name].__dict__)

            # prepare the array for the interpolated result
            f_data1 = np.zeros(gridcells, dtype=variable.datatype)

            # original variable data (source) that need to be interpolated
            o_data=np.zeros(land_points, dtype=variable.datatype)
             
            # Handle variables with two dimensions
            if (len(variable.dimensions) == 2):
                source = src[name][:]
                o_data = source[points_in_tessfa_land[2][:],points_in_tessfa_land[3][:]]
                f_data1 = griddata(points, o_data, (grid_lat, grid_lon), method=iMethod)

                # Special handling for AREA, LONGXY, LATIXY
                if name in ['AREA', 'LONGXY', 'LATIXY']:
                    if domain_type == '1':
                        # For 1D domain, assign directly
                        dst[name][...] = eval(name)
                    else:
                        # For 2D domain, fill only valid land cells
                        f_data = np.full((len(y_dim), len(x_dim)), fill_value, dtype=variable.datatype)
                        f_data[bool_mask] = eval(name)
                        dst[name][...] = f_data
                    count += 1
                    continue  # Skip the rest of the loop for these variables
                
                if (domain_type == '1'):
                    # Assign the interpolated data
                    dst[name][...] = np.copy(f_data1)
                    
                else:
      
                # put the masked data back to the data (with the NA land mask

                    f_data = np.ma.array(np.empty((len(y_dim),len(x_dim)), dtype=variable.datatype), mask=bool_mask, fill_value=fill_value)
                    f_data =  np.where(f_data.mask, f_data, fill_value)
                    f_data[bool_mask]=f_data1 
                    
                # Assign the interpolated data
                    dst[name][...] = np.copy(f_data)
                    
                '''print("o_data, f_data1, f_data, dst: max/min/sum")  
                print(np.nanmax(o_data), np.nanmax(f_data1),np.nanmax(f_data[f_data != -9999]),np.nanmax(dst[name]))
                print(np.nanmin(o_data), np.nanmin(f_data1),np.nanmin(f_data[f_data != -9999]),np.nanmin(dst[name]))   
                print(np.nansum(o_data), np.nansum(f_data1),np.nansum(f_data[f_data != -9999]),np.nansum(dst[name]))'''


                count = count + 1

            # Handle variables with three dimensions
            if (len(variable.dimensions) == 3):
                for index in range(variable.shape[0]):
                    # get all the source data (global)
                    source = src[name][index,:,:]
                    o_data = source[points_in_tessfa_land[2][:],points_in_tessfa_land[3][:]]
                    f_data1 = griddata(points, o_data, (grid_lat, grid_lon), method=iMethod)
                    
                    if (domain_type =='1'):
                        # Assign the interpolated data
                        dst[name][index,...] = np.copy(f_data1)
                    else:

                    # create a mask array to hold the interpolated data

                        f_data = np.ma.array(np.empty((len(y_dim),len(x_dim)), dtype=variable.datatype), mask=bool_mask, fill_value=fill_value)
                        f_data =  np.where(f_data.mask, f_data, fill_value)
                        f_data[bool_mask]=f_data1 
    
                    # Assign the interpolated data to dst.variable
                        dst[name][index,...] = np.copy(f_data)
                        
                    '''print("o_data, f_data1, f_data, dst: max/min/sum")  
                    print(np.nanmax(o_data), np.nanmax(f_data1),np.nanmax(f_data[f_data != -9999]),np.nanmax(dst[name][index,:,:]))
                    print(np.nanmin(o_data), np.nanmin(f_data1),np.nanmin(f_data[f_data != -9999]),np.nanmin(dst[name][index,:,:]))   
                    print(np.nansum(o_data), np.nansum(f_data1),np.nansum(f_data[f_data != -9999]),np.nansum(dst[name][index,:,:]))'''  

                count = count + variable.shape[0]

            # Handle variables with four dimensions
            if (len(variable.dimensions) == 4):
                for index1 in range(variable.shape[0]):
                    for index2 in range(variable.shape[1]):
                        # get all the source data (global)

                        source = src[name][index1, index2,:, :]
                        o_data = source[points_in_tessfa_land[2][:],points_in_tessfa_land[3][:]]
                        
                        f_data1 = griddata(points, o_data, (grid_lat, grid_lon), method=iMethod)                      
                        if (domain_type == '1'):
                            # Assign the interpolated data
                            dst[name][index1,index2,...] = np.copy(f_data1)
                        else:

                        # create a mask array to hold the interpolated data
                            f_data = np.ma.array(np.empty((len(y_dim),len(x_dim)), dtype=variable.datatype), mask=bool_mask, fill_value=fill_value)
                            f_data =  np.where(f_data.mask, f_data, fill_value)
                            f_data[bool_mask]=f_data1 
    
                        # Assign the interpolated data to dst.variable
                            dst[name][index1,index2,...] = np.copy(f_data)

                        #print("o_data, f_data1, f_data, dst: max/min/sum")  
                        #print(np.nanmax(o_data), np.nanmax(f_data1),np.nanmax(f_data[f_data != -9999]),np.nanmax(dst[name][index1,index2,:,:]))
                        #print(np.nanmin(o_data), np.nanmin(f_data1),np.nanmin(f_data[f_data != -9999]),np.nanmin(dst[name][index1,index2,:,:]))   
                        #print(np.nansum(o_data), np.nansum(f_data1),np.nansum(f_data[f_data != -9999]),np.nansum(dst[name][index1,index2,:,:]))  

                    count = count + variable.shape[1]

            end = process_time()
            print("Generating variable: " +name+ " takes  {}".format(end-start))

        else:

            # keep variables with the same dimension
            xerr = dst.createVariable(name, variable.datatype, variable.dimensions, zlib=True, complevel=5)
            # Copy variable attributes
            dst[name].setncatts(src[name].__dict__)
            # Copy the data

            dst[name][...] = src[name][...]
            
            end = process_time()
            print("Copying variable: " +name+ " takes  {}".format(end-start))
            
        if count > 50:
            dst.close()   # output the variable into a file to save memory

            dst = nc.Dataset(output_file, 'a')

            count = 0
        
        print(count)

    # Close the files
    src.close()
    dst.close()
    
if __name__ == '__main__':
    main()