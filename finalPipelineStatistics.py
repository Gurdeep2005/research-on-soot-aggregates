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

# this code looks at restructuring simulation results of 1 aggregate and returns a figure with 3 plots; Radius of Gyration vs Time, Coordination number vs Time, and Area vs Time
plt.rcParams.update({'font.size': 20})

# variables
totalTime = 5e-8
numberOfDumps = 101
numberOfDumps += 1
particleRadius = 14e-9
numberOfAggregates = 4

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
# def getAggregateValuesList(aggregateNumber):
#     aggregateValuesList = list(range(numberOfDumps))
#     for files in particles_by_agg_frac[(f"aggregate_{aggregateNumber}")]:
#         file_name = os.path.basename(files)
#         timestep = int(file_name.replace("particles_", "").replace(".vtk", ""))
#         pointsArr = getPointsArrFromVTKFile(files)
#         aggregateValuesList[timestep] = (findRadiusOfGyration(pointsArr), findCoordinationNumber(pointsArr), findAreaOfMonomer(pointsArr)*(particleRadius * particleRadius), convexityOverProjectedPlane(pointsArr, 0, 1))
#     return aggregateValuesList

# returns a list of the all values averaged from the 3 aggregates at a specific neck fraction value
# def getAverageValuesList():
#     aggreagte1Values = getAggregateValuesList(1)

#     averageValues = list(range(numberOfDumps))
#     for i in range(len(aggreagte1Values)):
#         averageValues[i] = []
#         for index in range(len(aggreagte1Values[i])):
#             averageValues[i].append((aggreagte1Values[i][index]))
    
#     return averageValues

# plot values
timeBetweenDump = totalTime / numberOfDumps

# xVals = [i * timeBetweenDump for i in range(numberOfDumps)]
# def plotAverageValues(label = ""):
#     averageValues = getAverageValuesList()

#     axis[0][0].plot(xVals, [x[0] for x in averageValues], label=label)

#     axis[0][1].plot(xVals, [x[1] for x in averageValues], label=label)

#     axis[1][0].plot(xVals, [x[2] for x in averageValues], label=label)

#     axis[1][1].plot(xVals, [x[3] for x in averageValues], label=label)

def populateConvexityLists(dict, sequenceNumber):
    finalConvex, secondToLastConvexity = None, None
    for files in dict[(f"aggregate_1")]:
        file_name = os.path.basename(files)
        timestep = int(file_name.replace("particles_", "").replace(".vtk", ""))
        pointsArr = getPointsArrFromVTKFile(files)

        if(timestep == numberOfDumps-1):
            finalConvex = convexityOverProjectedPlane(pointsArr, 0, 1)
        elif(timestep == numberOfDumps-2):
            secondToLastConvexity = convexityOverProjectedPlane(pointsArr, 0, 1)
    finalConvexity[f"sequence{sequenceNumber}"].append(finalConvex)
    differenceInConvexity[f"sequence{sequenceNumber}"].append(finalConvex-secondToLastConvexity)
    
experimentalConvexity = ((0.73713, 0.099442), (0.880628, 0.040486), (0.56, 0.06), (0.87, 0.04))

labels = ["TEG anchored", "TEG free", "H\u20820 anchored", "H\u20820 free"]
shape = ["o", "s", "v", "D"]
def plotConvexities():
    #axis[1].hist(differenceInConvexity["sequence1"] + differenceInConvexity["sequence2"] + differenceInConvexity["sequence3"] + differenceInConvexity["sequence4"], bins=10)
    
    x = np.linspace(0.5, 0.9, 100)
    y=x
    axis.plot(x, y, color="black", linestyle="--")
    
    x=0
    for sequence in differenceInConvexity.keys():
        mean = statistics.mean(finalConvexity[sequence])
        standardDeviation = statistics.stdev(finalConvexity[sequence])
        CI = 1.96*standardDeviation/math.sqrt(numberOfAggregates)

        axis.errorbar(experimentalConvexity[x][0], mean, yerr=CI, xerr=experimentalConvexity[x][1],label=labels[x], fmt=shape[x])
        x += 1

# setting up plot
fig, axis = plt.subplots(1, figsize=(9, 8))
#axis.set_title("Simulation VS Experimental Final Convexity")
#axis[1].set_title("Change in Convexity")

finalConvexity = {"sequence1":[], "sequence2":[], "sequence3":[], "sequence4":[]}
differenceInConvexity = {"sequence1":[], "sequence2":[], "sequence3":[], "sequence4":[]}

numberOfAggregates = 100
for num in range(1, 101):

    # ***************************************************************************

    inputList = [(f"/home/gurdeep/Desktop/simulationPipeline/simulationSequences/aggregate_{num}/sequence1/spherical_restructuring_anchored/", "A"), (f"//home/gurdeep/Desktop/simulationPipeline/simulationSequences/aggregate_{num}/sequence2/spherical_restructuring_free", "B"), (f"/home/gurdeep/Desktop/simulationPipeline/simulationSequences/aggregate_{num}/sequence3/spherical_restructuring_anchored", "C"), (f"/home/gurdeep/Desktop/simulationPipeline/simulationSequences/aggregate_{num}/sequence4/spherical_restructuring_free", "D")]
    sequenceNum = 1
    for value in inputList:
        # getting all particle file paths
        root_dir = os.path.expanduser(value[0])

        particles_by_agg_frac = defaultdict(list)

        for dirpath, dirnames, filenames in os.walk(root_dir):
            for filename in filenames:
                if filename.startswith("particles_") and (filename.endswith(f"{numberOfDumps-1}.vtk") or filename.endswith(f"{numberOfDumps-2}.vtk")):
                    # parts = dirpath.split(os.sep)
                    # if len(parts) >= 3 and parts[-3].startswith("aggregate_"):
                    #     agg = parts[-3].split("_")[-1]
                    #     key = (f"aggregate_{agg}")
                    #     full_path = os.path.join(dirpath, filename)
                    #     particles_by_agg_frac[key].append(full_path)
                    full_path = os.path.join(dirpath, filename)
                    particles_by_agg_frac["aggregate_1"].append(full_path)

        populateConvexityLists(particles_by_agg_frac, sequenceNum)
        sequenceNum += 1

    # ***************************************************************************
plotConvexities()

axis.legend()
#axis[1].legend()

axis.set_ylabel("Simulation Projected-Convexity")
#axis[1].set_ylabel("Amount of aggregates")
axis.set_xlabel("Experimental Projected-Convexity")
#axis[1].set_xlabel("Convexity Change")
plt.tight_layout()

plt.savefig(f"PipelineResults.pdf")
#plt.show()
