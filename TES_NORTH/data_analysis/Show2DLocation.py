import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

# Step 1: Read the NetCDF file
filename = 'mask.nc'  # Replace with your file path if needed
dataset = nc.Dataset(filename)

# Step 2: Get the coordinates and gridID from the dataset
x_coords = dataset.variables['x'][:]  # x coordinate values
y_coords = dataset.variables['y'][:]  # y coordinate values
gridID = dataset.variables['gridID'][:].compressed()  # gridID values

print(gridID)

# Create an array with the same dimensions as the 2D domain
# Initialize an array of zeros
active_array = np.zeros((len(y_coords), len(x_coords)), dtype=int)

# Flatten the active_array to 1D
flattened_active_array = active_array.flatten()
print(flattened_active_array.shape)

# Create a set for active grid IDs
active_grid_ids = set(gridID)  # Using a set for faster lookup

# Update the flattened
# active array based on active grid IDs
for grid_id in gridID:
    if grid_id < flattened_active_array.size:  # Ensure grid_id is within bounds
        flattened_active_array[grid_id] = 1

# Reshape the flattened array back to 2D
reshaped_active_array = flattened_active_array.reshape(active_array.shape)

# Plot the 2D map using matplotlib
plt.imshow(reshaped_active_array, cmap='viridis', origin='lower')
plt.colorbar(label='Active Status (1=Active, 0=Inactive)')
plt.title('2D Map of Active Grid IDs')
plt.xlabel('X Coordinates')
plt.ylabel('Y Coordinates')
plt.show()

# Step 4: Plot the array using imshow
plt.figure(figsize=(12, 6))
plt.imshow(reshaped_active_array, origin='lower', extent=(x_coords.min(), x_coords.max(), y_coords.min(), y_coords.max()), cmap='gray')
plt.colorbar(label='Active Cells (1 = Active, 0 = Inactive)')
plt.title('Active Grid Cell Representation')
plt.xlabel('Longitude (degrees)')
plt.ylabel('Latitude (degrees)')
plt.show()

# Close the dataset
dataset.close()
