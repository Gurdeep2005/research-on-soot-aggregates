import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib.pyplot as plt
import os
from collections import defaultdict
from shapely.geometry import Point
from shapely.ops import unary_union
import statistics
import math

# Read the VTK dump file
def getPointsArrFromVTKFile(filePath):
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(filePath)
    reader.Update()
    data = reader.GetOutput()
    point_coordinates = data.GetPoints().GetData()

    # points_array is a numpy array that contains coordinates of every primary particle in the aggregate
    # coordinates are normalized by primary particle radius
    points_array = vtk_to_numpy(point_coordinates)
    return points_array


# finds the distance between points p1 and p2
def distanceBetweenPoints(p1, p2):
    p1 = np.array(p1)
    p2 = np.array(p2)
    return np.linalg.norm(p1 - p2)

# Multiple of particle radius; Distance from center of particles to be considered neighbors
neighborThreshold = 2.5

# variables
totalTime = 5e-8
numberOfDumps = 301
particleRadius = 14e-9
numberOfAggregates = 1


# finds convexity over the selected plane with 0, 1, 2 representing the x, y, z, axis
def convexityOverProjectedPlane(particlePoints, axis1, axis2):
    projected_points = [Point(p[axis1], p[axis2]) for p in particlePoints]
    buffered_points = [pt.buffer(1) for pt in projected_points]
    
    union_shape = unary_union(buffered_points)
    
    if union_shape.is_empty:
        return 0
    
    convex_hull = union_shape.convex_hull
    if convex_hull.area != 0:
        convexity = union_shape.area / convex_hull.area
    else: 
        convexity = 0
    return convexity

# getting all particle file paths
root_dir = os.path.expanduser("/home/gurdeep/Desktop/simulationPipeline/aggregate-bank")

particles_by_agg_frac = defaultdict(list)

for dirpath, dirnames, filenames in os.walk(root_dir):
    for filename in filenames:
        full_path = os.path.join(dirpath, filename)
        particles_by_agg_frac["aggregate_1"].append(full_path)

totalConvexity = 0
numAggregates = 100

for files in particles_by_agg_frac[(f"aggregate_1")]:
    file_name = os.path.basename(files)
    pointsArr = getPointsArrFromVTKFile(files)

    totalConvexity += convexityOverProjectedPlane(pointsArr, 0, 1)

average = totalConvexity / numAggregates
print(average)
