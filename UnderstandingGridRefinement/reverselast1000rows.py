def process_dat_file(input_file, output_file):
    with open(input_file, 'r') as file:
        lines = file.readlines()

    # The first line (title) remains unchanged
    title = lines[0]
    
    # Extract the first 1000 rows and keep them intact
    first_1000_rows = lines[1:1001]
    
    # Extract the last 1000 rows and reverse them
    last_1000_rows = lines[1001:2001][::-1]

    # Combine the parts together
    processed_lines = [title] + first_1000_rows + last_1000_rows

    # Write the processed lines to the output file
    with open(output_file, 'w') as file:
        file.writelines(processed_lines)

# Example usage:
input_file = 'UnderstandingGridRefinement/refinedairfoil.dat'
output_file = 'UnderstandingGridRefinement/refinedairfoil2.dat'
process_dat_file(input_file, output_file)
