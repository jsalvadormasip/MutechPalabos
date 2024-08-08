import pandas as pd
import matplotlib.pyplot as plt

# Step 1: Load the data
# Assuming the file is formatted with whitespace separating the columns
file_path = 'C:/Users/jordi/Documents/GitHub/MutechPalabos/examples/showCases/jordiPowerFlowCopy/tmp/MarlonDesign_total_force.dat'
data = pd.read_csv(file_path, delim_whitespace=True, header=None, names=["time", "lift", "drag", "lateralforce"])

# Step 2: Plot the data
plt.figure(figsize=(10, 6))
plt.plot(data['time'], data['lift'], label='cl')
plt.plot(data['time'], data['drag'], label='cd')
plt.plot(data['time'], data['lateralforce'], label='Lateral Force')

# Adding labels and title
plt.xlabel('Time')
plt.ylabel('Forces')
plt.title('Lift, Drag, and Lateral Force vs Time')
plt.legend()

# Show plot
plt.show()
