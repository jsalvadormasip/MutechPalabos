import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt

# Parameters
plot_iterations = []
plots = False
area = 19 * 2 * 4
rho = 1.0  # Assuming constant density for incompressible flow
nBlock = 6.0
iterations = range(50, 401, 50)  # Modify the range as needed
file_path_template_input = 'C:/Users/jordi/Documents/GitHub/MutechPalabos/examples/showCases/jordiPowerFlowCopy/tmp/slice_x_01_{:08d}.vtm'
file_path_template_output = 'C:/Users/jordi/Documents/GitHub/MutechPalabos/examples/showCases/jordiPowerFlowCopy/tmp/slice_x_02_{:08d}.vtm'

# Lists to store drag force and drag coefficient for each iteration
drag_forces = []
drag_coefficients = []
iteration_numbers = []

for iteration in iterations:
    file_path_input = file_path_template_input.format(iteration)
    file_path_output = file_path_template_output.format(iteration)

    # Read datasets
    dataset_input = pv.read(file_path_input)
    dataset_output = pv.read(file_path_output)

    # Initialize lists to hold the datasets for merging
    velocity_datasets_to_merge_input = []
    velocity_datasets_to_merge_output = []
    pressure_datasets_to_merge_input = []
    pressure_datasets_to_merge_output = []

    # Process input dataset
    for i in range(dataset_input.GetNumberOfBlocks()):
        block = dataset_input.GetBlock(i)
        if block is not None and isinstance(block, pv.MultiBlock):
            for j in range(block.GetNumberOfBlocks()):
                sub_block = block.GetBlock(j)
                if sub_block is not None and isinstance(sub_block, pv.DataSet):
                    if 'velocity' in sub_block.array_names:
                        velocity_array = sub_block['velocity']
                        velocity_magnitude = np.linalg.norm(velocity_array, axis=1)
                        sub_block['velocity_magnitude'] = velocity_magnitude
                        sub_block['velocity_x'] = velocity_array[:, 0]
                        sub_block['velocity_y'] = velocity_array[:, 1]
                        sub_block['velocity_z'] = velocity_array[:, 2]
                        velocity_datasets_to_merge_input.append(sub_block)
                    if 'pressure' in sub_block.array_names:
                        pressure_datasets_to_merge_input.append(sub_block)

    # Process output dataset
    for i in range(dataset_output.GetNumberOfBlocks()):
        block = dataset_output.GetBlock(i)
        if block is not None and isinstance(block, pv.MultiBlock):
            for j in range(block.GetNumberOfBlocks()):
                sub_block = block.GetBlock(j)
                if sub_block is not None and isinstance(sub_block, pv.DataSet):
                    if 'velocity' in sub_block.array_names:
                        velocity_array = sub_block['velocity']
                        velocity_magnitude = np.linalg.norm(velocity_array, axis=1)
                        sub_block['velocity_magnitude'] = velocity_magnitude
                        sub_block['velocity_x'] = velocity_array[:, 0]
                        sub_block['velocity_y'] = velocity_array[:, 1]
                        sub_block['velocity_z'] = velocity_array[:, 2]
                        velocity_datasets_to_merge_output.append(sub_block)
                    if 'pressure' in sub_block.array_names:
                        pressure_datasets_to_merge_output.append(sub_block)

    # Merge the velocity datasets
    if velocity_datasets_to_merge_input:
        combined_velocity_dataset_input = velocity_datasets_to_merge_input[0]
        for ds in velocity_datasets_to_merge_input[1:]:
            combined_velocity_dataset_input = combined_velocity_dataset_input.merge(ds)
    if velocity_datasets_to_merge_output:
        combined_velocity_dataset_output = velocity_datasets_to_merge_output[0]
        for ds in velocity_datasets_to_merge_output[1:]:
            combined_velocity_dataset_output = combined_velocity_dataset_output.merge(ds)
    # Merge the pressure datasets
    if pressure_datasets_to_merge_input:
        combined_pressure_dataset_input = pressure_datasets_to_merge_input[0]
        for ds in pressure_datasets_to_merge_input[1:]:
            combined_pressure_dataset_input = combined_pressure_dataset_input.merge(ds)
    if pressure_datasets_to_merge_output:
        combined_pressure_dataset_output = pressure_datasets_to_merge_output[0]
        for ds in pressure_datasets_to_merge_output[1:]:
            combined_pressure_dataset_output = combined_pressure_dataset_output.merge(ds)
    # Drag force calculation
    inlet_velocity = np.average(np.array(combined_velocity_dataset_input['velocity_x']))
    outlet_velocity = np.array(combined_velocity_dataset_output['velocity_x'])
    drag_contribution = rho * outlet_velocity * (outlet_velocity - inlet_velocity) * area / outlet_velocity.size
    drag_force = np.sum(drag_contribution)  # Integrate

    # Store results
    domainlength = 40
    chord = 0.4 * domainlength
    span = 4
    drag_coefficient = drag_force / (0.5 * rho * inlet_velocity**2 * chord * span)
    drag_forces.append(drag_force)
    drag_coefficients.append(drag_coefficient)
    iteration_numbers.append(iteration)

    # Optionally plot the data for each iteration
    if iteration in plot_iterations:
        if 'pressure' in combined_pressure_dataset_input.array_names:
            plotter_pressure_input = pv.Plotter()
            plotter_pressure_input.add_mesh(combined_pressure_dataset_input, scalars='pressure', cmap='viridis')
            plotter_pressure_input.show()
        if 'pressure' in combined_pressure_dataset_output.array_names:
            plotter_pressure_output = pv.Plotter()
            plotter_pressure_output.add_mesh(combined_pressure_dataset_output, scalars='pressure', cmap='viridis')
            plotter_pressure_output.show()

        # Plot the velocity magnitude data
        if 'velocity_magnitude' in combined_velocity_dataset_input.array_names:
            plotter_velocity_magnitude_input = pv.Plotter()
            plotter_velocity_magnitude_input.add_mesh(combined_velocity_dataset_input, scalars='velocity_magnitude', cmap='plasma')
            plotter_velocity_magnitude_input.show()
        if 'velocity_magnitude' in combined_velocity_dataset_output.array_names:
            plotter_velocity_magnitude_output = pv.Plotter()
            plotter_velocity_magnitude_output.add_mesh(combined_velocity_dataset_output, scalars='velocity_magnitude', cmap='plasma')
            plotter_velocity_magnitude_output.show()
        # for i in range(combined_velocity_dataset['velocity_x'].size):
        #     if i % 59 ==0 :  #inlet
        #         # inlet_velocities.append(combined_velocity_dataset['velocity_x'][i])
        #         combined_velocity_dataset['velocity_x'][i] = i
        # # print("inlet velocities, max min avg ", np.max(inlet_velocities), np.min(inlet_velocities), np.average(inlet_velocities))
        # for i in range(combined_velocity_dataset['velocity_x'].size):
        #     if (i+10) % 60 ==0 and i <= 7200*2:  #outlet?
        #         outlet_velocities.append(combined_velocity_dataset['velocity_x'][i])
        # print("outlet velocities, max min avg ",np.max(outlet_velocities), np.min(outlet_velocities), np.average(outlet_velocities))
        # combined_velocity_dataset['velocity_x'][59] = -50000
        # Plot the x component of velocity
        if 'velocity_x' in combined_velocity_dataset_input.array_names:
            plotter_velocity_x_input = pv.Plotter()
            plotter_velocity_x_input.add_mesh(combined_velocity_dataset_input, scalars='velocity_x', cmap='coolwarm')
            plotter_velocity_x_input.show()
        if 'velocity_x' in combined_velocity_dataset_output.array_names:
            plotter_velocity_x_output = pv.Plotter()
            plotter_velocity_x_output.add_mesh(combined_velocity_dataset_output, scalars='velocity_x', cmap='coolwarm')
            plotter_velocity_x_output.show()
        # Plot the y component of velocity
        if 'velocity_y' in combined_velocity_dataset_input.array_names:
            plotter_velocity_y_input = pv.Plotter()
            plotter_velocity_y_input.add_mesh(combined_velocity_dataset_input, scalars='velocity_y', cmap='coolwarm')
            plotter_velocity_y_input.show()
        if 'velocity_y' in combined_velocity_dataset_output.array_names:
            plotter_velocity_y_output = pv.Plotter()
            plotter_velocity_y_output.add_mesh(combined_velocity_dataset_output, scalars='velocity_y', cmap='coolwarm')
            plotter_velocity_y_output.show()
        # Plot the z component of velocity
        if 'velocity_z' in combined_velocity_dataset_input.array_names:
            plotter_velocity_z_input = pv.Plotter()
            plotter_velocity_z_input.add_mesh(combined_velocity_dataset_input, scalars='velocity_z', cmap='coolwarm')
            plotter_velocity_z_input.show()
        if 'velocity_z' in combined_velocity_dataset_output.array_names:
            plotter_velocity_z_output = pv.Plotter()
            plotter_velocity_z_output.add_mesh(combined_velocity_dataset_output, scalars='velocity_z', cmap='coolwarm')
            plotter_velocity_z_output.show()

# Plotting drag force and drag coefficient vs iteration
times = np.array(iteration_numbers) * 5/400
plt.figure()
plt.subplot(2, 1, 1)
plt.plot(times, drag_forces, marker='o')
plt.title('Drag Force vs. Iteration')
plt.xlabel('Iteration')
plt.ylabel('Drag Force')

plt.subplot(2, 1, 2)
plt.plot(times, drag_coefficients, marker='o')
plt.title('Drag Coefficient vs. Iteration')
plt.xlabel('Iteration')
plt.ylabel('Drag Coefficient')

plt.tight_layout()
plt.show()
