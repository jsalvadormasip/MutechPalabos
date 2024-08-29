import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

def read_dat_file(file_path):
    # Read the .dat file including the title line
    with open(file_path, 'r') as file:
        title = file.readline().strip()
        data = np.loadtxt(file)
    
    x = data[:, 0]
    z = data[:, 1]
    return title, x, z

def refine_airfoil(x, z, num_points):
    # Create a cubic interpolation function based on the existing points
    interp_func = interp1d(x, z, kind='cubic', fill_value="extrapolate")
    
    # Generate new x values with higher resolution
    x_new = np.linspace(x.min(), x.max(), num_points)
    z_new = interp_func(x_new)
    
    return x_new, z_new

def write_dat_file(file_path, title, x, z):
    # Write the title and the refined x and z coordinates to a new .dat file
    with open(file_path, 'w') as file:
        file.write(title + '\n')
        data = np.column_stack((x, z))
        np.savetxt(file, data, fmt='%.6f')

def plot_airfoil(x, z, x_new, z_new):
    # Plot the original and refined airfoil points
    plt.figure(figsize=(10, 6))
    plt.plot(x, z, 'o', label='Original Points', markersize=4)
    plt.scatter(x_new, z_new, label='Refined Airfoil Points', color='red')
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('z')
    plt.title('Airfoil Refinement')
    plt.grid(True)
    plt.axis('equal')
    plt.show()

def main(input_file, output_file, num_points):
    title, x, z = read_dat_file(input_file)
    x_new, z_new = refine_airfoil(x, z, num_points)
    write_dat_file(output_file, title, x_new, z_new)
    
    # Plot for visualization
    plot_airfoil(x, z, x_new, z_new)

if __name__ == "__main__":
    input_file = "UnderstandingGridRefinement/airfoilupperpart.dat"   # Path to your input .dat file
    output_file = "UnderstandingGridRefinement/airfoilupperpartrefined.dat"  # Path to save the refined .dat file
    num_points = 1000  # Number of points in the refined airfoil line
    
    main(input_file, output_file, num_points)
