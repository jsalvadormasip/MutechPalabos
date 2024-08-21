import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt

# Parameters
DragCalculation = False
rho = 1.0  # Assuming constant density for incompressible flow


plot_iterations = [50]

if DragCalculation:

    area = 19 * 2 * 4
    iterations = range(50, 401, 50)  # Modify the range as needed
    file_path_template_input = 'C:/Users/jordi/Documents/GitHub/MutechPalabos/examples/showCases/jordiPowerFlowCopy/tmp/slice_x_01_{:08d}.vtm'
    file_path_template_output = 'C:/Users/jordi/Documents/GitHub/MutechPalabos/examples/showCases/jordiPowerFlowCopy/tmp/slice_x_02_{:08d}.vtm'

    # Lists to store drag force and drag coefficient for each iteration
    drag_forces = []
    drag_coefficients = []
    iteration_numbers = []
if not DragCalculation:
    iterations = [50]
    # file_path_template_input = 'examples/showCases/jordiPowerFlowCopy/tmp/slice_y_00_{:08d}.vtm'
    file_path_template_input = 'C:/Users/jordi/Documents/GitHub/MutechPalabos/examples/showCases/jordiPowerFlowCopy/tmp/slice_y_00_{:08d}.vtm'

    #with cp calc there are only input, no output

for iteration in iterations:
    file_path_input = file_path_template_input.format(iteration)
    if DragCalculation:
        file_path_output = file_path_template_output.format(iteration)

    # Read datasets
    dataset_input = pv.read(file_path_input)
    if DragCalculation:
        dataset_output = pv.read(file_path_output)

    # Initialize lists to hold the datasets for merging
    velocity_datasets_to_merge_input = []
    pressure_datasets_to_merge_input = []
    if DragCalculation:
        velocity_datasets_to_merge_output = []
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
                        sub_block['velocity_x'] = velocity_array[:, ]
                        sub_block['velocity_y'] = velocity_array[:, 1]
                        sub_block['velocity_z'] = velocity_array[:, 2]
                        velocity_datasets_to_merge_input.append(sub_block)
                    if 'pressure' in sub_block.array_names:
                        pressure_datasets_to_merge_input.append(sub_block)
    if DragCalculation:
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

    # Merge the input datasets
    if velocity_datasets_to_merge_input:
        combined_velocity_dataset_input = velocity_datasets_to_merge_input[0]
        for ds in velocity_datasets_to_merge_input[1:]:
            combined_velocity_dataset_input = combined_velocity_dataset_input.merge(ds)
    
    
    if pressure_datasets_to_merge_input:
        combined_pressure_dataset_input = pressure_datasets_to_merge_input[0]
        for ds in pressure_datasets_to_merge_input[1:]:
            combined_pressure_dataset_input = combined_pressure_dataset_input.merge(ds)
    if DragCalculation:
        # Merge the output datasets
        if velocity_datasets_to_merge_output:
            combined_velocity_dataset_output = velocity_datasets_to_merge_output[0]
            for ds in velocity_datasets_to_merge_output[1:]:
                combined_velocity_dataset_output = combined_velocity_dataset_output.merge(ds)
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
if not DragCalculation:
    # print(combined_pressure_dataset_input)
    # print(combined_pressure_dataset_input.points.shape)
    # print("velocity points shape",combined_velocity_dataset_input.points.shape)
    # print(combined_velocity_dataset_input['velocity_x'])
    # print(combined_velocity_dataset_input['velocity_magnitude'])
    # print(combined_pressure_dataset_input['pressure'].shape)
    # Assuming your array is named `arr`
    pressure_points = np.array(combined_pressure_dataset_input.points)  
    velocity_points = np.array(combined_pressure_dataset_input.points)  
    
    # Define the target point
    # Step 1: Read the .dat file
    filename = 'UnderstandingGridRefinement/airfoilcoordinates_clean.dat'  # Replace this with the path to your .dat file

    # Initialize a list to store the coordinates
    coordinates = []

    with open(filename, 'r') as file:
        lines = file.readlines()
        # Step 2: Skip the first line (TITLE OF FILE)
        for line in lines[1:]:
            # Step 3: Split the line into x and z, and convert to float
            x, z = map(float, line.split())
            # Step 4: Assign y = 0 and store the point
            coordinates.append([x, 0, z])

    # Convert the list of coordinates to a NumPy array
    target_points = np.array(coordinates)
    target_points[:,0] -= 0.5
    target_points *= 0.4*2 
    print("target points size", target_points.shape)

    # Now `target_points` contains your x, y, z coordinates
    

    # You can proceed to calculate distances or any other operation with target_points

    
    cellsnumbers = np.array([])
    countzeros = 0
    countbadestimates = 0
    # cellsnumbers2 = np.array([])
    # for i in range(combined_velocity_dataset_input['velocity_magnitude'].size):
    #     if combined_velocity_dataset_input['velocity_magnitude'][i] != 0 and combined_velocity_dataset_input['velocity_magnitude'][i-1] == 0: #or combined_velocity_dataset_input['velocity_magnitude'][i] != 0 and combined_velocity_dataset_input['velocity_magnitude'][i+1] == 0:
    #         cellsnumbers2 = np.append(cellsnumbers2, i) 


    for target_point in target_points:
        # Calculate the Euclidean distance from the target point for each sub-array
        distances = np.linalg.norm(pressure_points - target_point, axis=1)
        
        # Find the index of the minimum distance
        shouldbreak = False
        positions = np.argsort(distances)
        for ixx in positions:
            
            cellsposition = np.where(combined_pressure_dataset_input.cells == ixx)[0] // 9  

            # cellspositionfiltered = cellsposition[0]
            # cellspositionfiltered = combined_pressure_dataset_input['pressure'].size
            for i in range(cellsposition.size):
                if combined_velocity_dataset_input['velocity_magnitude'][int(cellsposition[i])] != 0.:
                    cellspositionfiltered = cellsposition[i]
                    shouldbreak = True
                    break
                if i == cellsposition.size-1:
                    # print("They are all zeros!")
                    countzeros += 1
                    cellspositionfiltered = np.nan
            if shouldbreak:
                break
            # should_break = False
                # for ix in range(15):
                #     if combined_pressure_dataset_input['pressure'][int(cellsposition[i])-ix] != 0:
                #         cellspositionfiltered = cellsposition[i]-ix
                #         should_break = True
                #         break
                #     if combined_pressure_dataset_input['pressure'][int(cellsposition[i])+ix] != 0:
                #         cellspositionfiltered = cellsposition[i]+ix
                #         should_break = True
                #         break
                # if should_break:
                #     break
                # else:
                #     print("too bad estimate")
                #     cellspositionfiltered = combined_pressure_dataset_input['pressure'].size
                #     countbadestimates += 1
                
        # print(cellsposition.size, "cellsposition.size")
        cellsnumbers = np.append(cellsnumbers, cellspositionfiltered) 
        
        # print(cellsposition)
        # print(combined_pressure_dataset_input.cells[0:20])
    print(countzeros, " this many all zeroes. ")
    print(countbadestimates, " this many badestimates. ")
    # print(cellsnumbers)
    print(cellsnumbers.shape)
    # print("cellsnumbers2",cellsnumbers2)
    pressure_distribution = []
    for i in range(cellsnumbers.size):
        # if combined_pressure_dataset_input['pressure'][i] == 0:  #inlet
        #     combined_pressure_dataset_input['pressure'][i] = 10
        # else:
        #     combined_pressure_dataset_input['pressure'][i] = 0
        if not np.isnan(cellsnumbers[i]):
            pressure_distribution.append(combined_pressure_dataset_input['pressure'][int(cellsnumbers[i])])
        else:
            pressure_distribution.append(np.nan)
        # combined_pressure_dataset_input['pressure'][int(cellsnumbers[i])] = i*10
    pressure_distribution = np.array(pressure_distribution)
    print(pressure_distribution.size, "pressure distribution size")
    print("outputing pressure file")
    import numpy as np

    # Select the x and z columns from target_points
    x_values = target_points[:, 0]
    z_values = target_points[:, 2]

    # Combine x, z, and pressure into one array
    output_data = np.column_stack((x_values, z_values, pressure_distribution))

    # Define the filename and title
    filename = 'Post Process/pressureMapAroundAirfoil.dat'
    title = 'x z p'

    # Write the data to a .dat file
    with open(filename, 'w') as file:
        # Write the title
        file.write(f"{title}\n")
        # Write the x, z, pressure data
        np.savetxt(file, output_data, fmt='%f %f %f')

    print(f"Data successfully written to {filename}")

    # Optionally plot the data for each iteration
    if iteration in plot_iterations:
        for ix in range(cellsnumbers.size):
            if not np.isnan(cellsnumbers[ix]):
                combined_pressure_dataset_input['pressure'][int(cellsnumbers[ix])] = 0
        if 'pressure' in combined_pressure_dataset_input.array_names:
            plotter_pressure_input = pv.Plotter()
            plotter_pressure_input.add_mesh(combined_pressure_dataset_input, scalars='pressure', cmap='viridis')
            plotter_pressure_input.show()
        if DragCalculation:
            if 'pressure' in combined_pressure_dataset_output.array_names:
                plotter_pressure_output = pv.Plotter()
                plotter_pressure_output.add_mesh(combined_pressure_dataset_output, scalars='pressure', cmap='viridis')
                plotter_pressure_output.show()

        # Plot the velocity magnitude data
        if 'velocity_magnitude' in combined_velocity_dataset_input.array_names:
            plotter_velocity_magnitude_input = pv.Plotter()
            plotter_velocity_magnitude_input.add_mesh(combined_velocity_dataset_input, scalars='velocity_magnitude', cmap='plasma')
            plotter_velocity_magnitude_input.show()
        if DragCalculation:
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
        if DragCalculation:
            if 'velocity_x' in combined_velocity_dataset_output.array_names:
                plotter_velocity_x_output = pv.Plotter()
                plotter_velocity_x_output.add_mesh(combined_velocity_dataset_output, scalars='velocity_x', cmap='coolwarm')
                plotter_velocity_x_output.show()
        # Plot the y component of velocity
        if 'velocity_y' in combined_velocity_dataset_input.array_names:
            plotter_velocity_y_input = pv.Plotter()
            plotter_velocity_y_input.add_mesh(combined_velocity_dataset_input, scalars='velocity_y', cmap='coolwarm')
            plotter_velocity_y_input.show()
        if DragCalculation:
            if 'velocity_y' in combined_velocity_dataset_output.array_names:
                plotter_velocity_y_output = pv.Plotter()
                plotter_velocity_y_output.add_mesh(combined_velocity_dataset_output, scalars='velocity_y', cmap='coolwarm')
                plotter_velocity_y_output.show()
        # Plot the z component of velocity
        if 'velocity_z' in combined_velocity_dataset_input.array_names:
            plotter_velocity_z_input = pv.Plotter()
            plotter_velocity_z_input.add_mesh(combined_velocity_dataset_input, scalars='velocity_z', cmap='coolwarm')
            plotter_velocity_z_input.show()
        if DragCalculation:
            if 'velocity_z' in combined_velocity_dataset_output.array_names:
                plotter_velocity_z_output = pv.Plotter()
                plotter_velocity_z_output.add_mesh(combined_velocity_dataset_output, scalars='velocity_z', cmap='coolwarm')
                plotter_velocity_z_output.show()
if DragCalculation:
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
