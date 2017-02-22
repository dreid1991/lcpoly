#Script for calculating radius of gyration from data.txt files

import sys
import math
import numpy as np
from scipy import stats
import copy
import glob

class Point:
    def __init__(self, x, z):
        self.x = x
        self.z = z

def readInput(inFile):
    f = open(inFile, 'r')
    out = []
    contour = []
    for line in f:
        data = map(float, line.strip().split())
        out.append(copy.deepcopy(data))
    return out

def organizeData(inArray):
    out = []
    snap = []
    for data in inArray:
        if (len(data) == 1 and len(snap) != 0):
            out.append(copy.deepcopy(snap))
            snap = []
        snap.append(copy.deepcopy(data))
    return out


def get_xAvg(inArray):
    avg = 0.0
    count = 0.0
    minVal = inArray[0].z
    maxVal = inArray[0].z
    for p in inArray:
        if (p.z > maxVal):
            maxVal = p.z
        if (p.z < minVal):
            minVal = p.z
        avg += p.x
        count += 1.0
    avg /= count
    return [avg, minVal, maxVal]

def getRegression(points):
    xData = np.array([])
    zData = np.array([])
    for p in points:
        xData = np.append(xData, p.x)
        zData = np.append(zData, p.z)
    slope, intercept, r_value, p_value, std_err = stats.linregress(zData,xData)
    return [slope, intercept]


eqTime = int(sys.argv[1])
data = readInput("data.txt")
numChains = int(data[0][0])
numDisks = int(data[1][0])
numTotal = numChains*numDisks
data = data[2 :]
sim = organizeData(data)
R_g = 0.0
count = 0.0
for system in sim:
    time = system[0][0]
    system = system[1 :]
    if time < eqTime:
        continue
    for i in range(numChains):
        x_avg = 0.0
        y_avg = 0.0
        z_avg = 0.0
        for j in range(numDisks):

            m = i*numDisks+j

            x = system[m][0]
            y = system[m][1]
            z = system[m][2]
            x_avg += x
            y_avg += y
            z_avg += z

        x_avg /= float(numDisks)
        y_avg /= float(numDisks)
        z_avg /= float(numDisks)
        avg = 0.0
        for j in range(numDisks):

            m = i*numDisks+j

            x = system[m][0]
            y = system[m][1]
            z = system[m][2]

            avg += (x-x_avg)*(x-x_avg)
            avg += (y-y_avg)*(y-y_avg)
            avg += (z-z_avg)*(z-z_avg)

        avg /= numDisks
        R_g += avg
        count += 1.0


R_g /= count

R_g = math.sqrt(R_g)


output = open('R_g.dat','w')

output.write(str(R_g))

output.close()

