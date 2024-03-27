
# AOI_forcingGEN for TESSFA

import os,sys
import netCDF4 as nc
import numpy as np
import pandas as pd
from time import process_time
from datetime import datetime

# Get current date
current_date = datetime.now()
# Format date to mmddyyyy
formatted_date = current_date.strftime('%m-%d-%Y')

def AOI_forcing_save_1d(input_path, file, AOI, AOI_points, output_path):
    # Open a new NetCDF file to write the data to. For format, you can choose from
    # 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
    source_file = input_path + '/'+ file
    print ("Opening source file: ", source_file)
    src = nc.Dataset(source_file, 'r', format='NETCDF4')
    
    #read gridIDs
    grid_ids = src['gridID'][...]    # gridID for all NA
 
    #  
    AOI_idx = np.where(np.in1d(grid_ids, AOI_points))[0]
    
    # create the new_filename
    dst_name = output_path + '/'+ AOI + '_'+file
    print ("Generating AOI file: ", dst_name)
    
    # check if file exists then delete it
    if os.path.exists(dst_name):
        os.remove(dst_name)
    
    # Open a new NetCDF file to write the data to. For format, you can choose from
    # 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
    dst = nc.Dataset(dst_name, 'w', format='NETCDF4')
    dst.title = dst_name +' created from '+ source_file +' on ' +formatted_date

    # Copy the global attributes from the source to the target
    for name in src.ncattrs():
        dst.setncattr(name, src.getncattr(name))

    # Copy the dimensions from the source to the target
    for name, dimension in src.dimensions.items():
        if name != 'ni' and name != 'gridcell':
            dst.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else None))
        else:
            # Update the 'ni' dimension with the length of the list
            #dst.dimensions['ni'].set_length(len(AOI_points))
            if name == 'ni' or name == 'gridcell': 
                ni = dst.createDimension(name, AOI_points.size)

    # Copy the variables from the source to the target
    for name, variable in src.variables.items():
        x = dst.createVariable(name, variable.datatype, variable.dimensions)   
        print(name, variable.dimensions)
        
        if (name != 'lambert_conformal_conic'):
            if variable.dimensions[-1] != 'ni' and variable.dimensions[-1] != 'gridcell':
                dst[name][...] = src[name][...]
            else:
                dst[name][...] = src[name][..., AOI_idx]
           
        # Copy the variable attributes
        for attr_name in variable.ncattrs():
            dst[name].setncattr(attr_name, variable.getncattr(attr_name))
        
    src.close()  # close the source file 
    dst.close()  # close the new file        
    
def get_files(input_path):
    print(input_path)
    files = os.listdir(input_path) 

    files.sort() 

    file_no =0

    files_nc = [f for f in files if (f[-2:] == 'nc')] 
    print("total " + str(len(files_nc)) + " files need to be processed")
    return files_nc

def main():
    args = sys.argv[1:]
    
    if len(sys.argv) != 5  or sys.argv[1] == '--help':  # sys.argv includes the script name as the first argument
        print("Example use: python AOI_forcingGEN.py <input_path> <output_path> <AOI_points_file>")
        print(" <input_path>: path to the 1D source data directory")
        print(" <output_path>:  path for the 1D AOI forcing data directory")
        print(" <AOI_gridID_path>:  path to the AOI_gridID_file")
        print(" <AOI_gridID_file>:  <AOI>_gridID.csv or <AOI>_domain.nc")
        print(" The code uses NA forcing to generation 1D AOI forcing")              
        exit(0)
    
    input_path = args[0]
    if not input_path.endswith("/"): input_path=input_path+'/'
    output_path = args[1]
    if not output_path.endswith("/"): output_path=output_path+'/'
    AOI_gridID_path = args[2]
    if not AOI_gridID_path.endswith("/"): AOI_gridID_path=output_path+'/'
    AOI_gridID_file = args[3]
    AOI=AOI_gridID_file.split("_")[0]

    '''
    #input_path = '/gpfs/wolf2/cades/cli185/proj-shared/Daymet_GSWP3_4KM_TESSFA/ELMinputdata/forcing1D_TES/'
    #input_path = '/gpfs/wolf2/cades/cli185/proj-shared/Daymet_GSWP3_4KM_TESSFA/ELMinputdata/'
    input_path = '/gpfs/wolf2/cades/cli185/proj-shared/Daymet_GSWP3_4KM_TESSFA/ELMinputdata/forcing1D_TES/netcdf/1998/'
    output_path = "./temp/dir2"
    AOI_gridID_path = '/gpfs/wolf2/cades/cli185/proj-shared/Daymet_GSWP3_4KM_TESSFA/ELMinputdata/AOI_data/MOFLUX/AOI_domain1D/'
    AOI_gridID_file = 'MOF21points_domain.lnd.Daymet_GSWP3_TESSFA.4km.1d.c231120.nc'
    AOI=AOI_gridID_file.split("_")[0]
    '''

    AOI_gridID_file = AOI_gridID_path + AOI_gridID_file
    
    if (AOI_gridID_file.endswith('.csv')):
        #AOI_gridcell_file = AOI+'_gridID.csv'  # user provided gridcell IDs
        df = pd.read_csv(AOI_gridID_file, sep=",", skiprows=1, names = ['gridID'])
        #read gridIds
        AOI_points = np.array(df['gridID'])
    elif AOI_gridID_file.endswith('.nc'):
        src = nc.Dataset(AOI_gridID_file, 'r')
        AOI_points = src['gridID'][:]
    else:
        print("Error: Invalid AOI_points_file, see help.")

    print(AOI_gridID_file)
        
    '''files_nc = get_files(input_path)

    for f in files_nc: 
        if (not f.startswith('clmforc')): continue
        print('processing '+ f )
        start = process_time() 
        AOI_forcing_save_1d(input_path, f, AOI, AOI_points, output_path)
        end = process_time()
        print("Generating 1D forcing data for "+AOI+ " domain takes {}".format(end-start))'''

    # Iterate over all subdirectories in the input directory
    for root, dirs, files in os.walk(input_path):
        for file in files:
            # Check if the file ends with '.nc'
            if file.endswith('.nc'):
                # Parse the filename into separate substrings
                parts = file.split('.')
                var_name = parts[3]
                period = parts[4]         
                print('processing '+ var_name + '(' + period + ') in the file ' + file )
                # Create the corresponding subfolder in the output directory
                new_dir = os.path.join(output_path, os.path.relpath(root, input_path))
                os.makedirs(new_dir, exist_ok=True)
                start = process_time()
                # Copy the file to the new location
                print(root, new_dir)
                #forcing_save_1dTES(root, file, var_name, period, time, new_dir)

                start = process_time() 
                AOI_forcing_save_1d(root, file, AOI, AOI_points, new_dir)
                end = process_time()
                print("Generating 1D forcing data for "+AOI+ " domain takes {}".format(end-start))

if __name__ == '__main__':
    main()
