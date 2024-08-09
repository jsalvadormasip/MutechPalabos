import numpy as np
import matplotlib.pyplot as plt

# Define the filename
filename = 'Post Process/pressureMapAroundAirfoil.dat'  # Replace this with your actual .dat file name

# Step 1: Read the .dat file
data = np.loadtxt(filename, skiprows=1)  # Skip the title row
ambient_pressure = 0.0
freestream_velocity = 1.0
rho = 1.0
# Step 2: Extract the x, z, and pressure data
x_values = data[:, 0]
z_values = data[:, 1]
pressure = data[:, 2]

minXvalues = np.min(x_values)
x_values -= minXvalues
maxXvalues = np.max(x_values)
x_values /= maxXvalues
z_values /= maxXvalues

pressure = (pressure-ambient_pressure)/(0.5*rho*freestream_velocity**2)
# Step 3: Create the first plot: Pressure vs. X and Z (with color map)
plt.figure(figsize=(10, 6))
plt.scatter(x_values, z_values, c=pressure, cmap='viridis')
plt.colorbar(label='Cp')
plt.xlabel('x/c')
plt.ylabel('thickness/c')
plt.title('Pressure Distribution across X and Z')
plt.grid(True)
plt.show()

# Step 4: Split the data into two halves
mid_index = len(x_values) // 2

x_values_first_half = x_values[:mid_index]
pressure_first_half = pressure[:mid_index]

x_values_second_half = x_values[mid_index:]
pressure_second_half = pressure[mid_index:]

# Step 5: Create the second plot: Pressure vs. X (two lines)
plt.figure(figsize=(10, 6))
plt.plot(x_values_first_half, pressure_first_half, marker='o', linestyle='-', color='b', label='Upper side')
plt.plot(x_values_second_half, pressure_second_half, marker='x', linestyle='-', color='r', label='Lower side')
plt.xlabel('x/c')
plt.ylabel('Cp')
plt.title('Pressure vs. chord (Upper side vs lower side)')
plt.legend()
plt.grid(True)
plt.gca().invert_yaxis()
plt.show()
