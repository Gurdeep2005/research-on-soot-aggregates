import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib.pyplot as plt
import os
from collections import defaultdict
from shapely.geometry import Point
from shapely.ops import unary_union

# this code looks at restructuring simulation results of 3 aggregates and returns a figure with 3 plots; Radius of Gyration vs Time, Coordination number vs Time, and Area vs Time

# variables
totalTime = 5e-8
numberOfDumps = 500
particleRadius = 14e-9

# getting all particle file paths
root_dir = os.path.expanduser("~/Software/soot-dem/dist/soot-dem/example/linux-tutorial/restructuring_simulations")

particles_by_agg_frac = defaultdict(list)

for dirpath, dirnames, filenames in os.walk(root_dir):
    for filename in filenames:
        if filename.startswith("particles_") and filename.endswith(".vtk"):
            parts = dirpath.split(os.sep)
            if len(parts) >= 3 and parts[-3].startswith("aggregate_") and parts[-2].startswith("frac_"):
                agg = parts[-3].split("_")[-1]
                frac = parts[-2].split("_")[-1]
                key = (f"aggregate_{agg}", f"frac_{frac}")
                full_path = os.path.join(dirpath, filename)
                particles_by_agg_frac[key].append(full_path)

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

# get aggregate value list which holds the radius of gyration and coordination number at each timestamp in a tuple
def getAggregateValuesList(aggregateNumber, fracNumber):
    aggregateValuesList = list(range(numberOfDumps))
    for files in particles_by_agg_frac[(f"aggregate_{aggregateNumber}", f"frac_{fracNumber}")]:
        file_name = os.path.basename(files)
        timestep = int(file_name.replace("particles_", "").replace(".vtk", ""))
        pointsArr = getPointsArrFromVTKFile(files)
        aggregateValuesList[timestep] = (findRadiusOfGyration(pointsArr), findCoordinationNumber(pointsArr), findAreaOfMonomer(pointsArr)*(particleRadius * particleRadius))
    return aggregateValuesList

# returns a list of the all values averaged from the 3 aggregates at a specific neck fraction value
def getAverageValuesList(fracNumber):
    aggreagte0Values = getAggregateValuesList(0, fracNumber)
    aggreagte1Values = getAggregateValuesList(1, fracNumber)
    aggreagte2Values = getAggregateValuesList(2, fracNumber)

    averageValues = list(range(numberOfDumps))
    for i in range(len(aggreagte0Values)):
        averageValues[i] = []
        for index in range(3):
            averageValues[i].append((aggreagte0Values[i][index] + aggreagte1Values[i][index] + aggreagte2Values[i][index])/3)
    
    return averageValues

# plot values
timeBetweenDump = totalTime / numberOfDumps

xVals = [i * timeBetweenDump for i in range(numberOfDumps)]
def plotAverageValues():
    averageValues0 = getAverageValuesList(0)
    averageValues30 = getAverageValuesList(30)
    averageValues70 = getAverageValuesList(70)
    averageValues90 = getAverageValuesList(90)

    axis[0].plot(xVals, [x[0] for x in averageValues0], label="0% Necks")
    axis[0].plot(xVals, [x[0] for x in averageValues30], label="30% Necks")
    axis[0].plot(xVals, [x[0] for x in averageValues70], label="70% Necks")
    axis[0].plot(xVals, [x[0] for x in averageValues90], label="90% Necks")

    axis[1].plot(xVals, [x[1] for x in averageValues0], label="0% Necks")
    axis[1].plot(xVals, [x[1] for x in averageValues30], label="30% Necks")
    axis[1].plot(xVals, [x[1] for x in averageValues70], label="70% Necks")
    axis[1].plot(xVals, [x[1] for x in averageValues90], label="90% Necks")

    axis[2].plot(xVals, [x[2] for x in averageValues0], label="0% Necks")
    axis[2].plot(xVals, [x[2] for x in averageValues30], label="30% Necks")
    axis[2].plot(xVals, [x[2] for x in averageValues70], label="70% Necks")
    axis[2].plot(xVals, [x[2] for x in averageValues90], label="90% Necks")

# setting up plot
fig, axis = plt.subplots(3, figsize=(8, 10))
axis[0].set_title("Radius of Gyration vs Time")
axis[1].set_title("Coordination Number vs Time")
axis[2].set_title("Area vs Time")

plotAverageValues()

axis[0].legend()
axis[1].legend()
axis[2].legend()

axis[0].set_ylabel("Radius of Gyration")
axis[1].set_ylabel("Coordination Number")
axis[2].set_ylabel("Area")
axis[0].set_xlabel("Time (s)")
axis[1].set_xlabel("Time (s)")
axis[2].set_xlabel("Time (s)")
plt.tight_layout()

plt.savefig("postProcessingPlot.png")