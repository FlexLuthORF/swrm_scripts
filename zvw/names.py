import os
import glob

# Specify the directory path to glob
directory_path = "."

# Glob the directory to get a list of subdirectories
subdirectories = [os.path.basename(d) for d in glob.glob(os.path.join(directory_path, "*")) if os.path.isdir(d)]

# Write the directory names to fofn.txt
with open("samples.txt", "w") as file:
    file.write("\n".join(subdirectories))

print("fofn.txt created successfully.")