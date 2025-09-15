import os
import glob


def module_files(filename):
    with open(filename, mode='r') as f:
        return [line.rstrip() for line in f]


def assign_to_dict(directory):
    # Initialize an empty dictionary to store module data
    modules_dict = {}

    # Use glob to find all output_file.* files in the specified directory
    for file_path in glob.glob(os.path.join(directory, "output_file.*")):
        # Extract the module number from the filename
        mod_number = os.path.basename(file_path).split('.')[-1]
        module_key = f'mod{mod_number}'

        # Store the contents of each module in the dictionary
        modules_dict[module_key] = module_files(file_path)

    return modules_dict




