import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

# Step 1: Read the domain file for the mask
mask_filename = 'mask.nc'  # Replace with your file path if needed
mask_dataset = nc.Dataset(mask_filename)

# Step 2: Get the coordinates and gridID from the mask dataset
x_coords = mask_dataset.variables['x'][:]  # x coordinate values
y_coords = mask_dataset.variables['y'][:]  # y coordinate values
gridID = mask_dataset.variables['gridID'][:].compressed()  # gridID values

# Create an array with the same dimensions as the 2D domain
active_array = np.zeros((len(y_coords), len(x_coords)), dtype=int)

# Flatten the active_array to 1D
flattened_active_array = active_array.flatten()

# Create a set for active grid IDs
active_grid_ids = set(gridID)  # Using a set for faster lookup

# Update the flattened active array based on active grid IDs
for grid_id in gridID:
    if grid_id < flattened_active_array.size:  # Ensure grid_id is within bounds
        flattened_active_array[grid_id] = 1

# Reshape the flattened array back to 2D
reshaped_active_array = flattened_active_array.reshape(active_array.shape)

# Step 3: Read the variable data from the data.nc file
data_filename = 'data.nc'  # Replace with your file path if needed
data_dataset = nc.Dataset(data_filename)

# Loop to allow the user to input multiple variable names until they choose to stop
while True:
    variable_name = input("Enter the variable name to extract (or type 'exit' to stop): ")
    
    if variable_name.lower() == 'exit':
        break  # Exit the loop if the user types 'exit'
    
    try:
        variable_data = data_dataset.variables[variable_name][:]  # Get the data variable
        variable_data = variable_data.compressed()  # Assuming same compression as gridID
        
        # Step 4: Create a masked version of the variable data
        masked_variable_data = np.full(active_array.shape, np.nan)  # Fill with NaN
        for i, grid_id in enumerate(gridID):
            if grid_id < masked_variable_data.size:  # Ensure grid_id is within bounds
                masked_variable_data[np.unravel_index(grid_id, masked_variable_data.shape)] = variable_data[i]

        # Step 5: Plot the data
        # Use a pcolormesh or imshow depending on the type of plot you want

        x_min = 200
        x_max = 500
        y_min = 200
        y_max = 500

        # Create new x and y coordinates for the subdomain
        sub_x_coords = x_coords[x_min-1:x_max]  # Since x_coords is 1-indexed
        sub_y_coords = y_coords[y_min-1:y_max]  # Since y_coords is 1-indexed

        # Select the corresponding data for the subdomain
        sub_data = masked_variable_data[y_min-1:y_max, x_min-1:x_max]

        plt.clf()  # Clear the previous plot
        plt.pcolormesh(sub_x_coords, sub_y_coords, sub_data, shading='auto', cmap='viridis')
        plt.colorbar(label=variable_name)
        plt.title('Plot of {} '.format(variable_name))
        plt.xlabel('X Coordinates')
        plt.ylabel('Y Coordinates')
        plt.show()
        
    except KeyError:
        print(f"Error: Variable '{variable_name}' does not exist in the dataset.")

# Step 6: Close the datasets
mask_dataset.close()
data_dataset.close()

