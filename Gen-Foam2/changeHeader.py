import os

folder = "/home/carlo/data/GeN-Foam/GeN-Foam_ffn/GeN-Foam/GeN-Foam/include"

# Use os.walk to iterate over all the files in the folder and its subfolders
for root, dirs, files in os.walk(folder):
    for filename in files:
        # Construct the full path to the file
        filepath = os.path.join(root, filename)

        # Open the file in read mode
        with open(filepath, "r") as file:
            # Read the contents of the file
            file_contents = file.read()

        # Replace the old string with the new string
        file_contents = file_contents.replace("Built on OpenFOAM v2212", "Built on OpenFOAM v2212")

        # Open the file in write mode
        with open(filepath, "w") as file:
            # Write the updated contents to the file
            file.write(file_contents)

