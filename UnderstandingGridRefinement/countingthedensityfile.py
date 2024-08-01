import matplotlib.pyplot as plt

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

# Example usage:
file_path = 'density.dat'  # Replace with your file path
result, num_numbers = sum_fourth_line(file_path)

print(f"The sum of all numbers in the 4th line is: {result}")
print(f"The total number of numbers in the 4th line is: {num_numbers}")

# Plot numbers at an interval
first_number = 1
interval = 100
plot_numbers_at_intervals(file_path, first_number, interval)
