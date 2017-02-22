#Script for calculating persistence length from data.txt file

import sys
import math
import numpy as np
from scipy import stats
import copy
import glob
from scipy.optimize import curve_fit

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

def exp(x,P):
    return np.exp(-x/P)

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
pers = []
count = []
for i in range(numDisks):
    pers.append(0.0)
    count.append(0.0)
for system in sim:
    time = system[0][0]
    system = system[1 :]
    if time < eqTime:
        continue
    for i in range(numChains):
        for j in range(numDisks):
            for k in range(j,numDisks):
                dist = k - j
                m1 = i*numDisks+j
                m2 = i*numDisks+k
                ux1 = system[m1][3]
                uy1 = system[m1][4]
                uz1 = system[m1][5]
                ux2 = system[m2][3]
                uy2 = system[m2][4]
                uz2 = system[m2][5]
                pers[dist] += ux1*ux2+uy1*uy2+uz1*uz2
                count[dist] += 1.0


for i in range(numDisks):
    pers[i] /= count[i]

xData = np.linspace(0,numDisks-1,numDisks)
yData = np.array(pers)

pops, pcov = curve_fit(exp,xData,yData)

output = open('Fits.dat','w')

for i in range(numDisks):
    output.write(str(i) + "\t" + str(pers[i]) + "\t" + str(exp(float(i),pops[0]))  + "\n")

output.close()

output = open("Persistence_Length.dat","w")
output.write(str(pops[0]))
output.close()
