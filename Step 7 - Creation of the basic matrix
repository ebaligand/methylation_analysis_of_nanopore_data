```python3
import os
import csv

# Directory containing the modified BED files
input_directory = './matrix/bed_modified'

# Directory to save the result matrix
output_directory = './matrix/result_matrix'

# Ensure the output directory exists
os.makedirs(output_directory, exist_ok=True)

# Dictionary to store data from the files
data = {}

# Loop through each file in the directory
for filename in os.listdir(input_directory):
   # Check if the file ends with '.txt'
   if filename.endswith('.txt'):
       filepath = os.path.join(input_directory, filename)
       # Open each file and read its contents line by line
       with open(filepath, 'r') as file:
           # Process each line in the file
           for line in file:
               # Split the line into parts based on whitespace
               parts = line.strip().split()
               col1, col2, col3 = parts[0], parts[1], parts[2]
               # If the first column is not already in the data dictionary, add it
               if col1 not in data:
                   data[col1] = {}
               # Store the value from the third column in the data dictionary
               data[col1][col2] = col3

# Get unique values from columns 1 and 2
unique_col1 = sorted(data.keys())
unique_col2 = sorted({col2 for col1_data in data.values() for col2 in col1_data.keys()})

# Create a matrix dictionary with unique column values as keys
matrix = {col1: {col2: None for col2 in unique_col2} for col1 in unique_col1}

# Populate the matrix dictionary with values from the data dictionary
for col1, col2_data in data.items():
   for col2, col3 in col2_data.items():
       matrix[col1][col2] = col3

# Define the output file path for the CSV matrix
output_filepath = os.path.join(output_directory, 'matrix.csv')

# Write the matrix to a CSV file
with open(output_filepath, 'w', newline='') as csvfile:
   writer = csv.writer(csvfile)
   # Write the header row with column names
   writer.writerow([''] + unique_col2)
   # Write data rows to the CSV file
   for col1 in unique_col1:
       row = [col1] + [matrix[col1][col2] if matrix[col1][col2] is not None else 'NULL' for col2 in unique_col2]
       writer.writerow(row)

# Print a message indicating successful creation of the matrix
print(f'Matrix created successfully and saved to {output_filepath}!')
```
