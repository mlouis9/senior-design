from paraview.simple import *
import numpy as np
import csv
import os
from paths import AVERAGE_COOLANT_CASE

os.chdir(AVERAGE_COOLANT_CASE)

# Create a new 'Open FOAM Reader'
averageCoolantfoam = OpenFOAMReader(registrationName='averageCoolant.foam', FileName='/home/mlouis9/PythonProjects/senior-design/cfd/3Dpincell/averageCoolant/averageCoolant.foam')
averageCoolantfoam.MeshRegions = ['internalMesh']
averageCoolantfoam.CellArrays = ['T', 'U', 'p', 'p_rgh']

# Create a new 'Clip'
clip1 = Clip(registrationName='Clip1', Input=averageCoolantfoam)
clip1.ClipType = 'Cylinder'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = ['POINTS', 'p']
clip1.Value = 25301432.3984375

# Init the 'Cylinder' selected for 'ClipType'
clip1.ClipType.Center = [0.0, 0.0, 1.4999999999019986]
clip1.ClipType.Axis = [0.0, 0.0, 1.0]
clip1.ClipType.Radius = 0.0056

# Init the 'Plane' selected for 'HyperTreeGridClipper'
clip1.HyperTreeGridClipper.Origin = [0.0, 0.0, 1.4999999999019986]

# Initialize an empty list to store the extracted data file paths
extracted_files = []

# Iterate over z levels from 0 to 3 with a step of 0.1
for z in np.arange(0, 3.1, 0.1):
    # Create a new 'Slice'
    slice1 = Slice(registrationName='Slice1', Input=clip1)
    slice1.SliceType = 'Plane'
    slice1.HyperTreeGridSlicer = 'Plane'
    slice1.SliceOffsetValues = [0.0]
    slice1.PointMergeMethod = 'Uniform Binning'

    # Set the slice origin based on the current z level
    slice1.SliceType.Origin = [0.0, 0.0, z]
    slice1.SliceType.Normal = [0.0, 0.0, 1.0]
    slice1.HyperTreeGridSlicer.Origin = [0.0, 0.0, z]

    # Create a new 'Integrate Variables'
    integrateVariables1 = IntegrateVariables(registrationName='IntegrateVariables1', Input=slice1)
    integrateVariables1.DivideCellDataByVolume = 1

    # Save the extracted data to a CSV file
    extracted_file = f'extracted_data_z{z:.1f}.csv'
    SaveData(extracted_file, proxy=integrateVariables1, FieldAssociation='Cell Data')
    
    # Append the extracted file path to the list
    extracted_files.append(extracted_file)

# Combine the extracted CSV files into a single file with an additional column for z height
if not os.path.exists('extracts'):
    os.mkdir('extracts')

with open('extracts/inner_surface_averaged_data.csv', 'w', newline='') as combined_file:
    writer = csv.writer(combined_file)
    
    for i, file_path in enumerate(extracted_files):
        with open(file_path, 'r') as extracted_file:
            reader = csv.reader(extracted_file)
            
            if i == 0:
                # Write the header row with an additional 'Z' column
                header = next(reader)
                header.insert(0, 'Z')
                writer.writerow(header)
            
            else:
                # Skip the header row for subsequent files
                next(reader)
            
            # Write the data rows with the corresponding z height
            z = np.arange(0, 3.1, 0.1)[i]
            for row in reader:
                row.insert(0, z)
                writer.writerow(row)
# Now remove the individual extract csv's
for z in np.arange(0, 3.1, 0.1):
    os.remove(f'extracted_data_z{z:.1f}.csv')