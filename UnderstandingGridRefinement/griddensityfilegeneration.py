import numpy as np

def read_airfoil_points(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()[1:]  # Skip the first line (title)
        airfoil_points = np.array([list(map(float, line.split())) for line in lines])
        
        # airfoil_points[:, 1] += 1
        airfoil_points *= 9
        print(airfoil_points)
    return airfoil_points  # Multiply all points by 10

def calculate_density(x, z, airfoil_points, density_ranges):
    # Calculate the distance from the point (x, z) to the nearest airfoil point
    distances = np.sqrt((x - airfoil_points[:, 0])**2 + (z - airfoil_points[:, 1])**2)
    min_distance = np.min(distances)
    
    # Determine the density based on the distance ranges
    for (dist_min, dist_max, density) in density_ranges:
        # print(dist_min, dist_max, density)
        if dist_min <= min_distance < dist_max:
            return density
    return 0.0  # Default density if no range matches

def create_density_file(nx, ny, nz, airfoil_points, density_ranges, output_file):
    densities = np.zeros((nx, ny, nz))
    
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                x = ix  # x-coordinate in the grid
                z = iz  # z-coordinate in the grid
                # print("ey")
                densities[ix, iy, iz] = calculate_density(x, z, airfoil_points, density_ranges)
    
    # Write to the .dat file with the specified format
    with open(output_file, 'w') as file:
        # First line: domain
        file.write("0 10 0 10 0 10\n")
        # Second line: dx
        file.write("1\n")
        # Third line: nx, ny, nz
        file.write(f"{nx} {ny} {nz}\n")
        # Fourth line: densities
        for density in densities.flatten(order='F'):  # 'F' means column-major order
            file.write(f"{density} ")

# Example usage:
# Airfoil points file path
airfoil_file_path = 'naca0018.dat'  # Replace with your airfoil points file path

# Read airfoil points from file
airfoil_points = read_airfoil_points(airfoil_file_path)

# Define density ranges: (distance_min, distance_max, density)
density_ranges = [
    (0, 1, 0.07),
    (1, 2, 0.14),
    (2, 3, 0.21),
    (3, 4, 0.28),
    (4, 5, 0.35),
    (5, 6, 0.42),
    (6, 7, 0.49),
    (7, 8, 0.56),
    (8, 9, 0.63),
    (9, 10, 0.7),
    (10, 11, 0.77),
    (11, 12, 0.84),
    (12, 13, 0.91),
    (13, 14.15, 1),
    # Add more ranges as needed
]

# Grid dimensions
nx, ny, nz = 10, 10, 10

# Output file
output_file = 'density.dat'

create_density_file(nx, ny, nz, airfoil_points, density_ranges, output_file)
