import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpltPath

def read_dat_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        
    title = lines[0].strip()
    points = []
    for line in lines[1:]:
        x, z = map(float, line.split())
        points.append([x, z])
        
    return title, np.array(points)

def write_dat_file(file_path, title, points):
    with open(file_path, 'w') as file:
        file.write(f"{title}\n")
        for point in points:
            file.write(f"{point[0]:.6f} {point[1]:.6f}\n")

def generate_grid_points(contour_points, spacing=0.01):
    min_x, min_z = np.min(contour_points, axis=0)
    max_x, max_z = np.max(contour_points, axis=0)
    
    x_values = np.arange(min_x, max_x, spacing)
    z_values = np.arange(min_z, max_z, spacing)
    
    grid_points = np.array([[x, z] for x in x_values for z in z_values])
    return grid_points

def filter_inside_points(contour_points, grid_points):
    path = mpltPath.Path(contour_points)
    inside = path.contains_points(grid_points)
    return grid_points[inside]

def plot_airfoil(contour_points, interior_points):
    plt.figure(figsize=(10, 6))
    plt.scatter(contour_points[:, 0], contour_points[:, 1], c='blue', s=1, label='Contour')
    # plt.scatter(interior_points[:, 0], interior_points[:, 1], c='red', s=1, label='Interior Points')
    plt.title("Airfoil Contour and Filled Points")
    plt.xlabel("X")
    plt.ylabel("Z")
    plt.axis('equal')
    plt.legend()
    plt.grid(True)
    plt.show()

def main(input_file, output_file, spacing=0.01):
    title, contour_points = read_dat_file(input_file)
    grid_points = generate_grid_points(contour_points, spacing)
    inside_points = filter_inside_points(contour_points, grid_points)
    
    # Combine contour and inside points for output
    all_points = np.vstack([contour_points, inside_points])
    
    write_dat_file(output_file, title, all_points)
    print(f"Filled airfoil saved to {output_file}")
    
    # Plotting the airfoil with the interior points
    plot_airfoil(contour_points, inside_points)

# Usage
input_file = "UnderstandingGridRefinement/airfoilcoordinates_clean.dat"
output_file = "UnderstandingGridRefinement/airfoilcoordinates_clean_filled.dat"
spacing = 0.005  # Adjust spacing as needed
main(input_file, output_file, spacing)
