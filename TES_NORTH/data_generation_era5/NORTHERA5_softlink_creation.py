
## Create forcing soft links for the TESSFA cases

import os
import glob
import shutil

path = "../atm_forcing.datm7.km.1d"

# Check if the directory exists
if os.path.isdir(path):
    # If it exists, remove it
    shutil.rmtree(path)

# Create a new folder
os.makedirs(path)
# Change to the new directory
os.chdir(path)

# Get a list of all files in the forcing directory and its subdirectories
files = glob.glob('../forcing/**/*', recursive=True)

 # Loop through the files
for file in files:
    # Only process actual files
    if not os.path.isfile(file):
        continue
    # Accept both 'clmforc' and 'climforc' in the filename
    base_name = os.path.basename(file)
    if ('clmforc' in base_name) or ('climforc' in base_name):
        # Build the link name to start with "clmforc.Daymet.km" and
        # then append everything from the first occurrence of ".1d" onward
        suffix = ".1d"
        end_index = base_name.find(suffix)
        if end_index == -1:
            continue
        new_link_name = f"clmforc.Daymet.km{base_name[end_index:]}"

        # Create a soft link in the target directory
        link_path = os.path.join(path, new_link_name)

        # Only create the link if it does not already exist
        if not os.path.exists(link_path):
            command = f'ln -s "{file}" "{link_path}"'
            print(command)
            os.system(command)
        else:
            print(f"Link {link_path} already exists, skipping.")

print("Forcing Soft links created successfully.")


"""Create or update domain soft link by auto-discovering domain file."""
domain_glob = "../domain_surfdata/domain.lnd.*.4km.1d*.nc"
domain_candidates = glob.glob(domain_glob)
domain_soft_link = "domain.lnd.Daymet.km.1d.nc"

if not domain_candidates:
    print(f"No domain file found matching: {domain_glob}")
else:
    # pick the most recently modified candidate
    domain_file = max(domain_candidates, key=os.path.getmtime)

    # If a link/file already exists with the target name, remove it first
    if os.path.islink(domain_soft_link) or os.path.exists(domain_soft_link):
        try:
            os.remove(domain_soft_link)
        except OSError:
            pass

    os.symlink(domain_file, domain_soft_link)
    print(f"Soft domain link created: {domain_soft_link} -> {domain_file}")


        
