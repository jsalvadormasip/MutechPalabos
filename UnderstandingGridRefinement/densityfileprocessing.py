import numpy as np
import matplotlib.pyplot as plt

def read_airfoil_points(file_path, z_offset,x_offset, chord):
    with open(file_path, 'r') as file:
        lines = file.readlines()[1:]  # Skip the first line (title)
        airfoil_points = np.array([list(map(float, line.split())) for line in lines])
    airfoil_points *= chord  # Multiply all points by 10
    airfoil_points[:, 1] += z_offset  # Add the offset to all z values
    airfoil_points[:, 0] += x_offset  # Add the offset to all x values
    return airfoil_points

def calculate_density(x, z, airfoil_points, density_ranges):
    # Calculate the distance from the point (x, z) to the nearest airfoil point
    distances = np.sqrt((x - airfoil_points[:, 0])**2 + (z - airfoil_points[:, 1])**2)
    min_distance = np.min(distances)
    
    # Determine the density based on the distance ranges
    for (dist_min, dist_max, density) in density_ranges:
        # if x >= int(nx/2):
        #     return 0.0
        # if x>=int(nx/4) and x<=int(nx/2):
        #     return 0.25
        # if x<=int(nx/4) and z>= int(3*nz/4):
        #     return 0.25
        # if x<=int(nx/4) and z<= int(nz/4):
        #     return 0.25
        # if x<=int(nx/8) and z<= int(3*nz/4) and z>=int(0.625*nz):
        #     return 0.5
        # if x<=int(nx/8) and z>= int(nz/4) and z<=int(0.375*nz):
        #     return 0.5
        # if x<=int(nx/4) and x>= int(nx/8) and z>= int(nz/4) and z<= int(nz*3/4):
        #     return 0.5
        # else:
        #     return 0.75
        # if x>=int(nx/4) and x <=int(nx/2) and z<= int(nz*0.625) and z>=int(0.5625*nz):
        #     return 0.75
        # if x>=int(nx/4) and x <=int(nx/2) and z>= int(nz*0.375) and z<=int(0.4375*nz):
        #     return 0.75
        # if x>=int(nx/4) and x <=int(nx/2) and z>= int(nz*0.4375) and z<=int(0.5625*nz):
        #     return 0.75
        # Assuming x and z are within the range of 0 to nx and 0 to nz respectively

        # Initialize density
        # Assuming x and z are within the range of 0 to nx and 0 to nz respectively

        # Initialize density
        # if (0 <= x < 0.25 * nx or 0.75 * nx < x <= nx) or (0 <= z < 0.25 * nz or 0.75 * nz < z <= nz):
        #     return 0

        # # Level 3
        # elif (0.25 * nx <= x < 0.375 * nx or 0.625 * nx < x <= 0.75 * nx) or (0.25 * nz <= z < 0.375 * nz or 0.625 * nz < z <= 0.75 * nz):
        #     return 0.25

        # # Level 4
        # elif (0.375 * nx <= x < 0.4375 * nx or 0.5625 * nx < x <= 0.625 * nx) or (0.375 * nz <= z < 0.4375 * nz or 0.5625 * nz < z <= 0.625 * nz):
        #     return 0.5

        # # Level 5 (center)
        # elif 0.4375 * nx <= x <= 0.5625 * nx and 0.4375 * nz <= z <= 0.5625 * nz:
        #     return 0.75

        # Default return value if no condition is met
        # return 1
            
            

        #if x<=int(boundaries[levelmin]*nx) or x>=int((1-boundaries[levelmin])*nx):
        #     return 0.0
        
        # if x<=int((boundaries[levelmin+1]+boundaries[levelmin])*nx) or x>=int(0.75*nx):
        #     return 0.0
        if offsets:
            # if x <= int(0.25*nx) or x >= int(0.75*nx):
            #     return 0.0
            # if z <= int(0.25*nz) or z >= int(0.75*nz):
            #     return 0.0
            if dist_min <= min_distance < dist_max:
                return density
        if not offsets:
            for ix in range(1,boundariesx.size):
                if x< np.floor(boundariesx[ix]*nx) and x>= np.floor(boundariesx[ix-1]*nx) or z< np.floor(boundariesx[ix]*nz) and z>= np.floor(boundariesx[ix-1]*nz):
                    return (ix - 2) / (levelmax - 2)
                if x>= np.floor((1-boundariesx[ix])*nx) and x< np.floor((1-boundariesx[ix-1])*nx) or z>= np.floor((1-boundariesx[ix])*nz) and z< np.floor((1-boundariesx[ix-1])*nz):
                    return (ix - 2) / (levelmax - 2)

    return 0.0  # Default density if no range matches

def create_density_file(nx, ny, nz, airfoil_points, density_ranges, output_file):
    densities2d = np.zeros((nx, nz))
    for ix in range(nx):
        print(ix/nx*100)
        for iz in range(nz):
                x = ix  # x-coordinate in the grid
                z = iz  # z-coordinate in the grid
                densities2d[ix, iz] = calculate_density(x, z, airfoil_points, density_ranges)
    num_layers = ny
    array_expanded = densities2d[:, np.newaxis, :]
    densities = np.repeat(array_expanded, num_layers, axis=1)
    print(densities.shape)
    
    # Write to the .dat file with the specified format
    with open(output_file, 'w') as file:
        # First line: domain
        file.write(f"{domainx0} {domainx1} {domainy0} {domainy1} {domainz0} {domainz1}\n")
        # Second line: dx
        file.write(f"{dx}\n")
        # Third line: nx, ny, nz
        file.write(f"{nx} {ny} {nz}\n")
        # Fourth line: densities
        for density in densities.flatten(order='C'):  # 'F' means column-major order
            file.write(f"{density} ")
        print("outputed to file")

def sum_fourth_line(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        
        # Check if there are at least four lines in the file
        if len(lines) < 4:
            raise ValueError("The file does not have at least four lines.")
        
        # Get the fourth line
        fourth_line = lines[3]
        
        # Split the line into individual numbers (assuming they are space-separated)
        numbers = list(map(float, fourth_line.split()))
        
        # Sum the numbers
        total_sum = sum(numbers)
        
        # Return the total sum and the number of numbers
        return total_sum, len(numbers)

def plot_numbers_at_intervals(file_path, first_num, interval):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        
        # Check if there are at least four lines in the file
        if len(lines) < 4:
            raise ValueError("The file does not have at least four lines.")
        
        # Get the fourth line
        fourth_line = lines[3]
        
        # Split the line into individual numbers (assuming they are space-separated)
        numbers = list(map(float, fourth_line.split()))
        
        # Select numbers at the specified interval
        selected_numbers = numbers[first_num-1::interval]
        
        # Plot the selected numbers
        plt.figure(figsize=(10, 6))
        plt.plot(selected_numbers, marker='o', linestyle='-', color='b')
        plt.title(f'Numbers from the Fourth Line starting at {first_num} with interval {interval}')
        plt.xlabel('Index')
        plt.ylabel('Value')
        plt.grid(True)
        plt.show()

def plot_y_slice(file_path, y_index, nx, ny, nz):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        
        # Check if there are at least four lines in the file
        if len(lines) < 4:
            raise ValueError("The file does not have at least four lines.")
        
        # Get the fourth line
        fourth_line = lines[3]
        
        # Split the line into individual numbers (assuming they are space-separated)
        densities = np.array(list(map(float, fourth_line.split())))
        
        # Reshape densities to (nx, ny, nz)
        densities = densities.reshape((nx, ny, nz), order='C')
        
        # Extract the y slice
        y_slice = densities[:, y_index, :]
        
        # Plot the y slice as a color map
        plt.figure(figsize=(10, 6))
        plt.imshow(y_slice.T, origin='lower', cmap='viridis', aspect='auto', extent=[0, nx, 0, nz])
        plt.colorbar(label='Density')
        plt.title(f'Y Slice at y={y_index}')
        plt.xlabel('X index')
        plt.ylabel('Z index')
        plt.show()

# Example usage:
# Airfoil points file path
airfoil_file_path = 'UnderstandingGridRefinement/refinedairfoil.dat'  # Replace with your airfoil points file path
#grid dimensions
offsets = True
domainx0 = -7.790592
dx = -1/2000*domainx0
domainx1 = 7.790592
domainy0 = -0.081625/2
domainy1 = 0.081625/2
domainz0 = -7.790592
domainz1 = 7.790592
nx = int(domainx1/dx*2+2)
ny = int(domainy1/dx*2+2)
nz = int(domainz1/dx*2+2)
print("nx ", nx, " ny ", ny, " nz ", nz)
# Define the z offset
z_offset = int(nz/2)  # Replace with the desired offset value
chord = int(0.1/7.790592*nx) #clean chord
x_offset = int(nx/2-1/2*chord)
minborderforchordx = 0.5-chord/2/nx
minborderforchordy = 0.5-chord*0.18/2/nx

levelmax = 14
levelmin = 2
boundariesx = np.zeros(levelmax)
boundariesy = np.zeros(levelmax)
for i in range(0, levelmax+1):
    if i>=2 :
        boundariesx[i:] += (1/(2**i))
    if i>=2 and boundariesy[-1]<minborderforchordy-(1/(2**i)):
        boundariesy[i:] += (1/(2**i))
boundariesx = np.append(boundariesx, 0.5)
print(boundariesx)
print(boundariesy)
# Read airfoil points from file with the z offset
airfoil_points = read_airfoil_points(airfoil_file_path, z_offset,x_offset, chord)

# Define density ranges: (distance_min, distance_max, density)
fine = True
medium = False
coarse = False
if fine:
    density_ranges = [
        (0, 1.58691/200*chord, 1),
        (1.58691/200*chord, 5.5542/200*chord, 0.9),
        (5.5542/200*chord, 15.0757/200*chord, 0.8),
        # (5.5542/200*chord, 10/200*chord, 11/13),
        # (10/200*chord, 15.0757/200*chord, 10/13), #fake 10
        (15.0757/200*chord, 37.2925/200*chord, 0.7),
        # (37.2925/200*chord, 55/200*chord, 8/13), #fake 55
        # (55/200*chord, 88.0737/200*chord, 7/13),
        (37.2925/200*chord, 88.0737/200*chord, 0.6), #fake 55
        # (88.0737/200*chord, 10000*np.sqrt(2)/200*chord, 0.5),
        # (88.0737/200*chord, 200/200*chord, 6/13), #fake 200
        # (200/200*chord, 456.256/200*chord, 5/13),
        (88.0737/200*chord, 456.256/200*chord, 0.5), 
        (456.256/200*chord, 600/200*chord, 0.4/1.0),
        # (600/200*chord, 10000*np.sqrt(2)/200*chord, 0.3),
        (600/200*chord, 1600/200*chord, 0.3/1.0),
        (1600/200*chord, 3750/200*chord, 2/10),
        (3750/200*chord, 9375/200*chord, 1/10),
        (9375/200*chord, 10000*np.sqrt(2)/200*chord, 0.0),

        # Add more ranges as needed
    ]
if medium:
    density_ranges = [
        (0, 1.58691/200*chord, 1),
        (1.58691/200*chord, 5.5542/200*chord, 0.875),
        (5.5542/200*chord, 15.0757/200*chord, 0.75),
        (15.0757/200*chord, 37.2925/200*chord, 0.625),
        (37.2925/200*chord, 88.0737/200*chord, 0.5),
        # (88.0737/200*chord, 10000*np.sqrt(2)/200*chord, 0.5),
        (88.0737/200*chord, 456.256/200*chord, 0.375),
        (456.256/200*chord, 600/200*chord, 0.25),
        # (600/200*chord, 10000*np.sqrt(2)/200*chord, 0.3),
        (600/200*chord, 1600/200*chord, 0.125),
        (1600/200*chord, 10000*np.sqrt(2)/200*chord, 0.0)

        # Add more ranges as needed
    ]
if coarse:
    density_ranges = [
        (0, 1.58691/200*chord, 1),
        (1.58691/200*chord, 5.5542/200*chord, 0.858),
        (5.5542/200*chord, 15.0757/200*chord, 0.715),
        (15.0757/200*chord, 37.2925/200*chord, 0.572),
        (37.2925/200*chord, 88.0737/200*chord, 0.429),
        # (88.0737/200*chord, 10000*np.sqrt(2)/200*chord, 0.5),
        (88.0737/200*chord, 456.256/200*chord, 0.286),
        (456.256/200*chord, 600/200*chord, 0.143),
        # (600/200*chord, 10000*np.sqrt(2)/200*chord, 0.3),
        (600/200*chord, 10000*np.sqrt(2)/200*chord, 0.0)

        # Add more ranges as needed
    ]
# density_ranges = [(i/10, (i+1)/10, 0.007*(i+1)/10) for i in range(1430)]
# density_ranges.append((14.30, 14.31, 1.0))


# Grid dimensions






# Output file
output_file = 'examples/showCases/jordiPowerFlowCopy/density.dat'

create_density_file(nx, ny, nz, airfoil_points, density_ranges, output_file)

# Summing and plotting
#result, num_numbers = sum_fourth_line(output_file)

# print(f"The sum of all numbers in the 4th line is: {result}")
# print(f"The total number of numbers in the 4th line is: {num_numbers}")

# Plot numbers at an interval
first_number = 1
interval = 100
# plot_numbers_at_intervals(output_file, first_number, interval)

# Plot y slice
y_index = int(0.5*ny)  # Choose the y index to plot
plot_y_slice(output_file, y_index, nx, ny, nz)
