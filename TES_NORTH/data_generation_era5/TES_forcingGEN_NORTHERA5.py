
import os,sys 
import netCDF4 as nc
import numpy as np
from time import process_time
from datetime import datetime, timedelta

def is_leap_year(year):
    """ Check if a given year is a leap year. """
    return year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)

# Get current date
current_date = datetime.now()
# Format date to mmddyyyy
formatted_date = current_date.strftime('%m%d%Y')

def forcing_save_1dTES(input_path, file, var_name, period, time, output_path, dataset_id_arg=None):
    # Open a new NetCDF file to write the data to. For format, you can choose from
    # 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
    source_file = input_path + '/'+ file
    src = nc.Dataset(source_file, 'r', format='NETCDF3_64BIT_DATA')
    
    # Build requested output filename:
    # climforc.<DATASETID>.<resolution:4km>.<dimension:1d>.<variable_name>.<time>.nc
    parts = file.split('.')
    # derive resolution token (e.g., '4km') from any part that endswith('km'); fallback to '4km'
    resolution = None
    for p in parts:
        if p.endswith('km'):
            resolution = p
            break
    if resolution is None:
        resolution = '4km'

    # pick DATASETID: prefer explicit arg; else derive from output_path: .../<DATASETID>/entire_domain/forcing/
    if dataset_id_arg and len(dataset_id_arg.strip()) > 0:
        dataset_id = dataset_id_arg.strip()
    else:
        norm_out = os.path.normpath(output_path)
        dataset_id = os.path.basename(os.path.dirname(os.path.dirname(norm_out)))

    new_filename = f"climforc.{dataset_id}.{resolution}.1d.{var_name}.{period}.nc"
    print('new filename is: ' + new_filename)

    total_rows = src.dimensions['lat'].size
    total_cols = src.dimensions['lon'].size
    total_time = src.dimensions['time'].size
    # remove the leap day 
    if total_time == 232:
            total_time = 224

    #print('total timesteps is :' + str(total_timesteps))
    if time == -1:
        time = total_time

    lat= src['lat'][:]
    lon= src['lon'][:]

    lonxy, latxy = np.meshgrid(lon, lat)

    #print(lonxy.shape)
    '''
    # time in 'days since ' current yyyy-mm-01-01 00:00:00
    data_time = src['time'][0:time] # read (time, y, x) format
    tunit = src.variables['time'].getncattr('units')  # Kao's data is with datetime of leap-year 
    t0=str(tunit.lower()).strip('days since')
    t0=datetime.strptime(t0,'%Y-%m-%d %X')
    iyr = t0.year + np.floor(data_time/365.0)
    data_time0 = datetime.strptime(str(int(iyr[0]))+'-01-01','%Y-%m-%d')
    data_time0 = (data_time0-t0).total_seconds()/86400.0  # days from t0 at the beginning of year
    iday = data_time - data_time0   # now in DOY of leapyear format, will be re-filled below
    imm = np.zeros_like(iday)
    mdoy = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]
    for m in range(1,13):
        tpts = np.where((iday>mdoy[m-1]) & (iday<=mdoy[m]))
        if (len(tpts[0])>0): 
            imm[tpts] = m             # in MM, may be 1 day off for leap-year
            iday[tpts]= iday[tpts]-mdoy[m-1]     # day of current month 
 
    data_time=iday          # in days of current year-month
    tunit = tunit.replace(str(t0.year).zfill(4)+'-', str(int(iyr[0])).zfill(4)+'-')
    tunit = tunit.replace('-'+str(t0.month).zfill(2)+'-', '-'+str(int(imm[0])).zfill(2)+'-')
    if(tunit.endswith(' 00') and not tunit.endswith(' 00:00:00')):
        tunit=tunit+':00:00'
    
    '''

    # Assuming src['time'] is already provided and contains the number of days since a reference date

    # Sample code for time variable transformation:
    time_variable = src['time'][:]  # Read full time variable

    # remove the leap day
    if (len(time_variable)) == 232:
        time_variable = time_variable[0:224]
    tunit = src.variables['time'].getncattr('units')  # Expected format: "days since 1950-01-01 00:00:00"
    print(tunit[10:])
    
    t0 = datetime.strptime(tunit[11:], '%Y-%m-%d %H:%M:%S')  # Parse the reference time
    
    # Initialize an output array for the days of the current month
    days_in_month = np.zeros_like(time_variable, dtype=float)
    months = np.zeros_like(time_variable, dtype=int)  # Store month indices

    # Reference handling for each year
    for i, days_since in enumerate(time_variable):
        # Calculate the full year based on total days since t0
        current_date = t0 + timedelta(days=float(days_since))  # Current date from base date
        '''
        year = t0.year + int(days_since // 365.25)  # Using 365.25 for average leap year calculation
        # Handle the end of year transition
        if (current_date.year > year):
            year += 1
        '''

        # Get the year, month, and the day of the month
        year = current_date.year
        month = current_date.month
        day = current_date.day

        # Store month and day information
        months[i] = month
        #days_in_month[i] = day
        # Instead of storing days, calculate times for 8 intervals per day
        days_in_month[i] = (day-1)  + (1.5 + (i%8) *3 )/24

    # Prepare the output time units in the required format
    new_time_unit = f"days since {int(year)}-{str(month).zfill(2)}-01 00:00:00"

    # Example for how to update the data_time output
    # Match the requested number of timesteps for the 'time' dimension
    data_time = days_in_month[0:time]

    # Get mask and create gridID, latxy, lonxy first
    
    for name, variable in src.variables.items():
        
        # Check if the last two dimensions are lat and lon
        # Get the mask and create gridIDs, latxy, longxy
        
        if (variable.dimensions[-2:] == ('lat', 'lon')):
            
            data = src[name][0:1, :, :]
            #create a land mask
            mask = data[0]    # data is in (time, Y, X) format
            # setup the bool_mask and mesh of the TESSFA domain (lon, lat)
            if np.ma.is_masked(data[0]):
                mask = ~data.mask.squeeze()
            else:
                mask = (~np.isnan(data[0])) & (~np.isclose(data[0], fill_value)).squeeze()
            #mask = np.where(~np.isnan(mask), 1, np.nan)

            print("data.shape and maske.shape are", data.shape,mask.shape)
            print("Number of non-NaN and non-fillvalue cells:", np.count_nonzero(mask))

            #create gridIDs
            total_gridcells = total_rows * total_cols
            grid_ids = np.linspace(0, total_gridcells-1, total_gridcells, dtype=int)
            grid_ids = grid_ids.reshape(total_rows,total_cols)

             # create a flattened list of land gridID
            #grid_ids = np.multiply(mask,grid_ids)
            grid_ids = grid_ids[mask]
            grid_ids = grid_ids[~np.isnan(grid_ids)]            
            print("Number of land cells:", np.count_nonzero(grid_ids))

            grid_ids = grid_ids[~np.isnan(grid_ids)]

            #print(lonxy.shape, latxy.shape, mask.shape)
            latxy_arr = latxy[mask]
            lonxy_arr = lonxy[mask]

            # convert local grid_id_lists into an array
            grid_id_arr = grid_ids
            print("Number of land cells in grid_id_arr:", len(grid_id_arr))

            #lonxy_arr= np.array(lonxy1)
            #latxy_arr= np.array(latxy1)
    
    
    dst_name = output_path + '/'+ new_filename

    # Open a new NetCDF file to write the data to. For format, you can choose from
    # 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
    dst = nc.Dataset(dst_name, 'w', format='NETCDF3_64BIT_DATA')
    dst.title = var_name + '('+period+') creted from '+ input_path +' on ' +formatted_date

    # create the gridIDs, lon, and lat variable
    x = dst.createDimension('time', time)
    x = dst.createDimension('ni', grid_id_arr.size)
    x = dst.createDimension('nj', 1)

    fill_value = np.int32(-9999)
    w_nc_var = dst.createVariable('gridID', np.int32, ('nj','ni'), fill_value=fill_value, zlib=True, complevel=5)
    #  set variable attributes
    w_nc_var.long_name = "gridId in the NA domain" ;
    w_nc_var.decription = "Covers all land and ocean gridcells, with #0 at the upper left corner of the domain" ;
    dst.variables['gridID'][...] = grid_id_arr.reshape(grid_id_arr.size,1)

    # Copy the variables from the source to the target

    for name, variable in src.variables.items():
        #if (name == var_name):
        start = process_time()
        print("Working on varibale: "+ name + " dimensions: " + str(variable.dimensions))
       
        if variable.datatype == np.int32:
            fill_value = -9999  # or any other value that you want to use to represent missing data
        if variable.datatype == np.float32:
            fill_value = np.float32(-9999)  # or any other value that you want to use to represent missing data

        # Check if the last two dimensions are lat and lon and save the data 
        
        if (variable.dimensions[-2:] == ('lat', 'lon') ) :
            
            data = src[name][0:time, :, :]
            
            # extract the data over land gridcells

            landcells = len(grid_id_arr)
            data_arr = data[:, mask]  # shape (time, landcells)
            
            w_nc_var = dst.createVariable(name, np.float32, ('time', 'nj', 'ni'), zlib=True, complevel=5)
            dst.variables[name][:] =data_arr.reshape(time,grid_id_arr.size)
            for attr_name in variable.ncattrs():
                if attr_name != '_FillValue':  # Skip the _FillValue attribute
                    dst[name].setncattr(attr_name, variable.getncattr(attr_name))
                    
        if (name == 'time'):
            dvname = 'time'
            w_nc_var = dst.createVariable(dvname, np.float32, ('time'), fill_value=fill_value, zlib=True, complevel=5)
            dst.variables[dvname][...] = data_time
            for attr_name in variable.ncattrs():
                if 'units' in attr_name:
                    dst[dvname].units = new_time_unit
                elif 'calendar' in attr_name:
                    dst[dvname].calendar = "no_leap"                     
                else:
                    dst[dvname].setncattr(attr_name, variable.getncattr(attr_name))

        if (name == 'lat'):
            dvname = 'LATIXY'
            w_nc_var = dst.createVariable(dvname, np.float64, ('nj','ni'),fill_value=fill_value,  zlib=True, complevel=5)
            dst.variables[dvname][...] = latxy_arr
            for attr_name in variable.ncattrs():
                dst[dvname].setncattr(attr_name, variable.getncattr(attr_name))

        if (name == 'lon'):
            dvname = 'LONGXY'
            w_nc_var = dst.createVariable(dvname, np.float64, ('nj','ni'), fill_value=fill_value, zlib=True, complevel=5)
            dst.variables[dvname][...] = lonxy_arr
            for attr_name in variable.ncattrs():
                dst[dvname].setncattr(attr_name, variable.getncattr(attr_name))

    src.close()  # close the source file 
    dst.close()  # close the new file        
    

'''i
#input_path = "/gpfs/wolf2/cades/cli185/proj-shared/Daymet_GSWP3_4KM_TESSFA/netcdf/2014/TPHWL3Hrly/"
#input_path = "/gpfs/wolf2/cades/cli185/proj-shared/Daymet_GSWP3_4KM_TESSFA/"
#input_path = "/gpfs/wolf2/cades/cli185/proj-shared/Daymet_GSWP3_4KM_TESSFA/netcdf/2014/Solar3Hrly/"
#output_path = "./temp/"
#time = 1

files_nc = get_files(input_path)

for f in files_nc: 
    if (not f.startswith('clmforc')): continue
    parts = f.split('.')
    common_string = parts[0]
    project=parts[1]
    res = parts[2]
    var_name = parts[3]
    period = parts[4]
    print('processing '+ var_name + '(' + period + ') in the file ' + f )
    start = process_time() 
    forcing_save_1dTES(input_path, f, var_name, period, time, output_path)
    end = process_time()
    print("Generating 1D forcing data takes {}".format(end-start))
'''



def main():
    args = sys.argv[1:]

    if (len(args) not in (3,4)) or (args and args[0] == '--help'):
        print("Example use: python TES_forcingGEN_NORTHERA5.py <input_path> <output_path> <time steps> [DATASETID]")
        print(" <input_path>: path to the 1D source data directory")
        print(" <output_path>: path for the 1D AOI forcing data directory")
        print(" <time steps>: timesteps to be processed or -1 (all time series)")
        print(" [DATASETID]: optional dataset id used in output filenames (e.g., Daymet_ERA5_TESSFA_NORTH)")
        exit(0)

    input_path = args[0]
    output_path = args[1]
    time = int(args[2])
    dataset_id_arg = args[3] if len(args) == 4 else None

    # Iterate over all subdirectories in the input directory
    allowed_top_dirs = set(["TPHWL3Hrly", "Solar3Hrly", "Precip3Hrly"])  # ERA5 directories
    for root, dirs, files in os.walk(input_path):
        rel_root = os.path.relpath(root, input_path)
        # At the top-level of input_path: prune traversal to only allowed directories
        if rel_root == '.':
            dirs[:] = [d for d in dirs if d in allowed_top_dirs]
            # still ignore any files at the top-level
            continue
        # For deeper levels: ensure we remain under one of the allowed top-level directories
        top_dir = rel_root.split(os.sep, 1)[0]
        if top_dir not in allowed_top_dirs:
            # prune traversal into disallowed branches (e.g., cpl_bypass_full)
            dirs[:] = []
            continue
        for file in files:
            # Check if the file ends with '.nc'
            if file.endswith('.nc'):
                # Parse the filename into separate substrings
                parts = file.split('.')
                # Source pattern: clmforc.<dataset>.<resolution_km>.<var>.<yyyy-mm>.nc
                var_name = parts[3]
                period = parts[4]         
                print('processing '+ var_name + '(' + period + ') in the file ' + file )
                # Create the corresponding subfolder in the output directory
                new_dir = os.path.join(output_path, os.path.relpath(root, input_path))
                os.makedirs(new_dir, exist_ok=True)
                start = process_time()
                # Copy the file to the new location
                print(root, new_dir)
                forcing_save_1dTES(root, file, var_name, period, time, new_dir, dataset_id_arg)
                end = process_time()
                print("Generating 1D forcing data takes {}".format(end-start))


if __name__ == '__main__':
    main()
    
