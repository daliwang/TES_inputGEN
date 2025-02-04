import os
import glob
import shutil
import re

path = "../atm_forcing.datm7.uELM_NADaymet.1d.c231120/"

# Check if the directory exists
if os.path.isdir(path):
    # If it exists, remove it
    shutil.rmtree(path)

# Create a new folder
os.makedirs(path)

# Navigate to the specified directory
os.chdir('../domain_surfdata/')

# Compile the regular expression for matching filenames
filename_regex = re.compile(r'(.*)_domain\.lnd\.Daymet_NA\.1km\.1d\.c(.*).nc')

# Initialize variables to hold the AOI_case_name and AOI_case_date
AOI_case_name = None
AOI_case_date = None

# Iterate over all files in the current directory
for filename in os.listdir('.'):
    # If the filename matches the regular expression
    match = filename_regex.match(filename)
    if match:
        # Extract the AOI_case_name and AOI_case_date from the filename
        AOI_case_name, AOI_case_date = match.groups()
        break

os.chdir(path)

# If we found a matching file
if AOI_case_name and AOI_case_date:
    # Execute the shell command
    command = 'ln -s ../domain_surfdata/'+ AOI_case_name + '_domain.lnd.Daymet_NA.1km.1d.c'+ AOI_case_date +'.nc domain.lnd.Daymet4.1km.1d.c231120.nc'
    print(command)
    os.system(command)
else:
    print('No matching file found.')


# Get a list of all files in the input_data directory
files = glob.glob('../forcing/*')

# Loop through the files
for file in files:
    # Check if 'clmforc' is in the file name
    if 'clmforc' in file:
        # Split the file name on '_'
        parts = file.split('_')
        # Use the string after the '_' to create a soft link
        link_name = parts[1]

        command = 'ln -s '+ file + ' ' + link_name
        print(command)
        os.system(command)
        #os.exec(ln -s ../data/file link_name)
        #print("ln -s "+ "../data/"+ file + " " + link_name)


