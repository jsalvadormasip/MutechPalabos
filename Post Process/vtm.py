import pyvista as pv
import numpy as np

# Constants
rho = 1.0  # Assuming constant density for incompressible flow

# Load the VTM file
file_path = 'C:/Users/jordi/Documents/GitHub/MutechPalabos/examples/showCases/jordiPowerFlowCopy/tmp/slice_y_00_00004200.vtm'
dataset = pv.read(file_path)

# Initialize lists to hold the datasets for merging
velocity_datasets_to_merge = []
pressure_datasets_to_merge = []

# Variables to store inlet velocities and modified outlet velocities
inlet_velocities = []
modified_outlet_velocities = []
modified_outlet_areas = []

# Define the number of columns to the left of the outlet
columns_to_left = 10
# print("number of blocks: ", dataset.GetNumberOfBlocks())
# Explore the MultiBlock structure
for i in range(dataset.GetNumberOfBlocks()):
    
    block = dataset.GetBlock(i)

    if block is not None and isinstance(block, pv.MultiBlock):
        print("number of subblocks: ", block.GetNumberOfBlocks())
        for j in range(block.GetNumberOfBlocks()):
            sub_block = block.GetBlock(j)

            if sub_block is not None and isinstance(sub_block, pv.DataSet):
                # Process velocity data
                if 'velocity' in sub_block.array_names:
                    velocity_array = sub_block['velocity']
                    
                    # print(velocity_array)
                    # print(velocity_array.shape)
                    points = sub_block.points
                    # print(points.shape)
                    # Identify inlet based on x-coordinate (leftmost points)
                    # inlet_mask = np.isclose(points[:, 0], points[:, 0].min())

                    # Identify the modified outlet based on the position 10 columns to the left
                    # unique_x = np.unique(points[:, 0])
                    # if len(unique_x) >= columns_to_left + 1:
                    #     modified_outlet_x = unique_x[-columns_to_left - 1]
                    #     modified_outlet_mask = np.isclose(points[:, 0], modified_outlet_x)
                    # else:
                    #     modified_outlet_mask = np.array([False] * len(points))

                    # Extract and store inlet and modified outlet velocities
                    # if inlet_mask.any():
                    #     inlet_velocities.append(sub_block.extract_points(inlet_mask)['velocity'])
                    # if modified_outlet_mask.any():
                    #     extracted_outlet = sub_block.extract_points(modified_outlet_mask)
                    #     modified_outlet_velocities.append(extracted_outlet['velocity'])

                    #     # Store outlet areas for integration
                    #     modified_outlet_areas.append(extracted_outlet.compute_cell_sizes()['Area'])

                    # Process and store components for plotting
                    velocity_magnitude = np.linalg.norm(velocity_array, axis=1)
                    velocity_x_component = velocity_array[:, 0]
                    velocity_y_component = velocity_array[:, 1]
                    velocity_z_component = velocity_array[:, 2]
                    sub_block['velocity_magnitude'] = velocity_magnitude
                    sub_block['velocity_x'] = velocity_x_component
                    sub_block['velocity_y'] = velocity_y_component
                    sub_block['velocity_z'] = velocity_z_component
                    velocity_datasets_to_merge.append(sub_block)

                # Process pressure data
                if 'pressure' in sub_block.array_names:
                    pressure_datasets_to_merge.append(sub_block)

# Merge the velocity datasets
if velocity_datasets_to_merge:
    combined_velocity_dataset = velocity_datasets_to_merge[0]
    for ds in velocity_datasets_to_merge[1:]:
        combined_velocity_dataset = combined_velocity_dataset.merge(ds)
# print(combined_velocity_dataset)
print("cells ", combined_velocity_dataset.cells[0:270])
# print("cells ", combined_velocity_dataset.cells.shape)
# print(combined_velocity_dataset['velocity'].shape)
# print("array names ",combined_velocity_dataset.array_names)
# print(np.min(combined_velocity_dataset.points, axis=0))
# print(np.max(combined_velocity_dataset.points, axis=0))
# print(np.max(combined_velocity_dataset['velocity_x'], axis=0))
# print(np.min(combined_velocity_dataset.cells, axis=0))
# print(np.max(combined_velocity_dataset.cells, axis=0))

# Merge the pressure datasets
if pressure_datasets_to_merge:
    combined_pressure_dataset = pressure_datasets_to_merge[0]
    for ds in pressure_datasets_to_merge[1:]:
        combined_pressure_dataset = combined_pressure_dataset.merge(ds)
# print(combined_pressure_dataset)
# Plot the pressure data
if 'pressure' in combined_pressure_dataset.array_names:
    plotter_pressure = pv.Plotter()
    plotter_pressure.add_mesh(combined_pressure_dataset, scalars='pressure', cmap='viridis')
    plotter_pressure.show()

# Plot the velocity magnitude data
if 'velocity_magnitude' in combined_velocity_dataset.array_names:
    plotter_velocity_magnitude = pv.Plotter()
    plotter_velocity_magnitude.add_mesh(combined_velocity_dataset, scalars='velocity_magnitude', cmap='plasma')
    plotter_velocity_magnitude.show()

for i in range(combined_velocity_dataset['velocity_x'].size):
    if i % 60 ==0 and i >=7200*2:  #inlet
        inlet_velocities.append(combined_velocity_dataset['velocity_x'][i])
print(np.max(inlet_velocities), np.min(inlet_velocities), np.average(inlet_velocities))
# combined_velocity_dataset['velocity_x'][59] = -50000
# Plot the x component of velocity
if 'velocity_x' in combined_velocity_dataset.array_names:
    plotter_velocity_x = pv.Plotter()
    plotter_velocity_x.add_mesh(combined_velocity_dataset, scalars='velocity_x', cmap='coolwarm')
    plotter_velocity_x.show()

# Plot the y component of velocity
if 'velocity_y' in combined_velocity_dataset.array_names:
    plotter_velocity_y = pv.Plotter()
    plotter_velocity_y.add_mesh(combined_velocity_dataset, scalars='velocity_y', cmap='coolwarm')
    plotter_velocity_y.show()

# Plot the z component of velocity
if 'velocity_z' in combined_velocity_dataset.array_names:
    plotter_velocity_z = pv.Plotter()
    plotter_velocity_z.add_mesh(combined_velocity_dataset, scalars='velocity_z', cmap='coolwarm')
    plotter_velocity_z.show()

# # Calculate Drag Force based on inlet and modified outlet velocities
# if inlet_velocities and modified_outlet_velocities:
#     # Calculate the average inlet and outlet velocities
#     inlet_velocity = np.mean(np.vstack(inlet_velocities), axis=0)
#     outlet_velocity = np.vstack(modified_outlet_velocities)

#     # Calculate the drag contribution at each point on the outlet
#     drag_contribution = rho * outlet_velocity * (outlet_velocity - inlet_velocity)

#     # Flatten the outlet areas list
#     modified_outlet_areas = np.hstack(modified_outlet_areas)

#     # Integrate over the outlet surface
#     drag_force = np.sum(drag_contribution * modified_outlet_areas[:, None], axis=0)  # Integrate

#     # Output the drag force (for all directions)
#     print(f"Drag Force (x, y, z): {drag_force}")
# else:
#     print("Inlet or modified Outlet data not found. Unable to calculate drag force.")
