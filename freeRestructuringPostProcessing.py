import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib.pyplot as plt
import os
from collections import defaultdict
from shapely.geometry import Point
from shapely.ops import unary_union

# this code looks at restructuring simulation results of 1 aggregate and returns a figure with 3 plots; Radius of Gyration vs Time, Coordination number vs Time, and Area vs Time

# variables
totalTime = 5e-8
numberOfDumps = 301
particleRadius = 14e-9
numberOfAggregates = 1

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

# Compute Radius of Gyration distance
def findRadiusOfGyration(points_array):
    center_of_mass = points_array.mean(axis=0)
    distances = points_array - center_of_mass
    rms_distance = np.sqrt(np.tensordot(distances, distances, 2) / distances.shape[0])
    return rms_distance

# finds the distance between points p1 and p2
def distanceBetweenPoints(p1, p2):
    p1 = np.array(p1)
    p2 = np.array(p2)
    return np.linalg.norm(p1 - p2)

# Multiple of particle radius; Distance from center of particles to be considered neighbors
neighborThreshold = 2.5

# Compute the coordination number

def findCoordinationNumber(points_array):
    totalNeighbors = 0
    for i in range(points_array.shape[0]):
        for j in range(i + 1, points_array.shape[0]):
            if distanceBetweenPoints(points_array[i], points_array[j]) <= neighborThreshold:
                totalNeighbors += 1
    return 2 * totalNeighbors / points_array.shape[0]

# finds the area of the monomer over the selected plane with 0, 1, 2 representing the x, y, z, axis
def findAreaOverProjectedPlane(particlePoints, axis1, axis2):
    circles = []
    for particlePoint in particlePoints:
        circles.append(Point(particlePoint[axis1], particlePoint[axis2]).buffer(1))
    combinedCircles = unary_union(circles)
    totalArea = combinedCircles.area
    return totalArea

# finds the area of the monomer by finding the projected areas over the xy, xz, and yz planes; then averaging them out
def findAreaOfMonomer(particlePoints):
    area = (findAreaOverProjectedPlane(particlePoints, 0, 1) + findAreaOverProjectedPlane(particlePoints, 0, 2) + findAreaOverProjectedPlane(particlePoints, 1, 2))/3
    return area

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

# finds average convexity across XY, XZ, and YZ planes
def averageConvexityOfMonomer(particlePoints):
    convexities = [
        convexityOverProjectedPlane(particlePoints, 0, 1),  # XY plane
        convexityOverProjectedPlane(particlePoints, 0, 2),  # XZ plane
        convexityOverProjectedPlane(particlePoints, 1, 2)   # YZ plane
    ]
    return sum(convexities) / len(convexities)

# get aggregate value list which holds the radius of gyration and coordination number at each timestamp in a tuple
def getAggregateValuesList(aggregateNumber):
    aggregateValuesList = list(range(numberOfDumps))
    for files in particles_by_agg_frac[(f"aggregate_{aggregateNumber}")]:
        file_name = os.path.basename(files)
        timestep = int(file_name.replace("particles_", "").replace(".vtk", ""))
        pointsArr = getPointsArrFromVTKFile(files)
        aggregateValuesList[timestep] = (findRadiusOfGyration(pointsArr), findCoordinationNumber(pointsArr), findAreaOfMonomer(pointsArr)*(particleRadius * particleRadius), averageConvexityOfMonomer(pointsArr))
    return aggregateValuesList

# returns a list of the all values averaged from the 3 aggregates at a specific neck fraction value
def getAverageValuesList():
    aggreagte1Values = getAggregateValuesList(1)

    averageValues = list(range(numberOfDumps))
    for i in range(len(aggreagte1Values)):
        averageValues[i] = []
        for index in range(len(aggreagte1Values[i])):
            averageValues[i].append((aggreagte1Values[i][index]))
    
    return averageValues

# plot values
timeBetweenDump = totalTime / numberOfDumps

xVals = [i * timeBetweenDump for i in range(numberOfDumps)]
def plotAverageValues(label = ""):
    averageValues = getAverageValuesList()

    axis[0][0].plot(xVals, [x[0] for x in averageValues], label=label)

    axis[0][1].plot(xVals, [x[1] for x in averageValues], label=label)

    axis[1][0].plot(xVals, [x[2] for x in averageValues], label=label)

    axis[1][1].plot(xVals, [x[3] for x in averageValues], label=label)


# setting up plot
fig, axis = plt.subplots(2, 2, figsize=(20, 10))
axis[0][0].set_title("Radius of Gyration vs Time")
axis[0][1].set_title("Coordination Number vs Time")
axis[1][0].set_title("Area vs Time")
axis[1][1].set_title("Convexity vs Time")

# ***************************************************************************

inputList = [("~/Downloads/aggregate_1.5e-11", "Force = 1.5e-11"), ("~/Downloads/aggregate_1e-11", "Force = 1e-11"), ("~/Downloads/aggregate_5e-11", "Force = 5e-11"), ("~/Downloads/aggregate_5e-12", "Force = 5e-12"), ("~/Downloads/aggregate_1e-12", "Force = 1e-12"),("~/Downloads/aggregate_1e-10", "Force = 1e-10")]

for value in inputList:
    # getting all particle file paths
    root_dir = os.path.expanduser(value[0])

    particles_by_agg_frac = defaultdict(list)

    for dirpath, dirnames, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename.startswith("particles_") and filename.endswith(".vtk"):
                # parts = dirpath.split(os.sep)
                # if len(parts) >= 3 and parts[-3].startswith("aggregate_"):
                #     agg = parts[-3].split("_")[-1]
                #     key = (f"aggregate_{agg}")
                #     full_path = os.path.join(dirpath, filename)
                #     particles_by_agg_frac[key].append(full_path)
                full_path = os.path.join(dirpath, filename)
                particles_by_agg_frac["aggregate_1"].append(full_path)

    plotAverageValues(value[1])

# ***************************************************************************

axis[0][0].legend()
axis[0][1].legend()
axis[1][0].legend()
axis[1][1].legend()

axis[0][0].set_ylabel("Radius of Gyration")
axis[0][1].set_ylabel("Coordination Number")
axis[1][0].set_ylabel("Area")
axis[1][1].set_ylabel("Convexity")
axis[0][0].set_xlabel("Time (s)")
axis[0][1].set_xlabel("Time (s)")
axis[1][0].set_xlabel("Time (s)")
axis[1][1].set_xlabel("Time (s)")
plt.tight_layout()

plt.savefig("postProcessingPlot.png")