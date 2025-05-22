# created by egor on 10/18/2024

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

# Read the VTK dump file
reader = vtk.vtkPolyDataReader()
reader.SetFileName('particles_0.vtk')
reader.Update()
data = reader.GetOutput()
point_coordinates = data.GetPoints().GetData()

# points_array is a numpy array that contains coordinates of every primary particle in the aggregate
# coordinates are normalized by primary particle radius
points_array = vtk_to_numpy(point_coordinates)

### DATA POST-PROCESSING EXAMPLE ###

# Compute RMS distance from the center of mass

center_of_mass = points_array.mean(axis=0)
distances = points_array - center_of_mass
rms_distance = np.sqrt(np.tensordot(distances, distances, 2) / distances.shape[0])

print(f'RMS distance from the center of mass is: {rms_distance:.2f} * r_part')