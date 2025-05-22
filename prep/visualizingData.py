import matplotlib.pyplot as plt
import numpy

studentID = numpy.loadtxt("student-data.csv", skiprows=1, delimiter=",", usecols=0)
exam_1_grades = numpy.loadtxt("student-data.csv", skiprows=1, delimiter=",", usecols=1)
exam_2_grades = numpy.loadtxt("student-data.csv", skiprows=1, delimiter=",", usecols=2)

figure, axis = plt.subplots(3, 1, figsize=(10,8))
graph1Axis2 = axis[0].twinx()
graph2Axis2 = axis[1].twinx()

def findMean(sequence):
    total = 0
    amount = 0
    for number in sequence:
        total += number
        amount += 1
    return total/amount

def findStandardDeviation(sequence):
    mean = findMean(sequence)

    added = 0
    amount = 0
    for number in sequence:
        added += (number - mean)**2
        amount += 1
    
    return numpy.sqrt(added/amount)

def findNormalDistribution(sequence, mean, standardDeviation):
    normalDistrbution = list()

    for number in sequence:
        value = (1 / (numpy.sqrt(2 * numpy.pi * standardDeviation**2))) * numpy.exp(-1 * (number - mean)**2 / (2 * standardDeviation**2))
        normalDistrbution.append(value)
    
    return normalDistrbution

exam1Mean = findMean(exam_1_grades)
exam1StandardDeviation = findStandardDeviation(exam_1_grades)
exam2Mean = findMean(exam_2_grades)
exam2StandardDeviation = findStandardDeviation(exam_2_grades)


# finding linear regression of exam2 as a function of exam1
sumOfCrossProducts = 0
sumOfSquaredDeviation = 0
for i in range(0, exam_1_grades.shape[0]):
    sumOfCrossProducts += (exam_1_grades[i] - exam1Mean)*(exam_2_grades[i] - exam2Mean)
    sumOfSquaredDeviation += (exam_1_grades[i] - exam1Mean)**2

slope = sumOfCrossProducts/sumOfSquaredDeviation
yIntercept = exam2Mean - slope * exam1Mean
# hard coded min and max values
min = 35
max = 70
xPoint = list(range(min, max))
yPoint = list(range(min, max))
for i in range(max-min):
    xPoint[i] = i + min
    yPoint[i] = slope*(i + min) + yIntercept

# finding Coefficient of determination

SSTotal = 0
SSRes = 0
for i in range(exam_2_grades.shape[0]):
    SSTotal += (exam_2_grades[i] - exam2Mean)**2
    SSRes += (exam_2_grades[i] - (slope * exam_1_grades[i] + yIntercept))**2

CoefficientOfDetermination = 1 - SSRes/SSTotal

# normal distributions
values = numpy.linspace(20, 80)

# plotting
axis[0].hist(exam_1_grades)
graph1Axis2.plot(values, findNormalDistribution(values, exam1Mean, exam1StandardDeviation), color="red")
axis[0].set_title("Exam 1 Grades")
axis[1].hist(exam_2_grades)
graph2Axis2.plot(values, findNormalDistribution(values, exam2Mean, exam2StandardDeviation), color="red")
axis[1].set_title("Exam 2 Grades")
axis[2].scatter(exam_1_grades, exam_2_grades)
axis[2].plot(xPoint, yPoint, color="red")
axis[2].set_title("Exam 2 Grades as a function of Exam 1 grades")

print(exam1Mean)
print(exam1StandardDeviation)
print(exam2Mean)
print(exam2StandardDeviation)
print(sumOfCrossProducts)
print(sumOfSquaredDeviation)
print(slope)
print(yIntercept)
print(CoefficientOfDetermination)


plt.show()