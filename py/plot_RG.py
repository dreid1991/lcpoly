import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import math

bendList = [0.01,0.05,0.07,0.1,0.3,0.5,0.7,0.9,1.0,1.3,1.6,1.8,2.0]
epsList = [0.0,0.1,0.5,1.0,1.5,2.0,5.0,10.0,15.0,20.0,30.0,50.0]

numCurves = len(epsList)
numPoints = len(bendList)

x = np.array(bendList)

curves = np.zeros((numCurves,numPoints))

for i,eps in enumerate(epsList):
    for j,bend in enumerate(bendList):
        os.chdir(str(eps)+"/"+str(bend)+"/1/")
        data = np.loadtxt("R_g.dat")
        curves[i][j] = data
        os.chdir("../../../")

colors = ['blue','green','red','cyan','magenta','yellow','black','brown','pink']

for i in range(len(epsList)):
    plt.errorbar(x,curves[i],marker = 'o', label = "Epsilon = "+ str(epsList[i]))

plt.xlabel("Bending Parameter")
plt.ylabel("Radius of Gyration (LJ $\sigma$)")
x1,x2,y1,y2 = plt.axis()
plt.axis((x1,x2,y1,y2))
legend = plt.legend(loc=2)
for label in legend.get_texts():
    label.set_fontsize('small')


plt.savefig("RG_plot.png")

outfile = open("RG_data.txt","w")

for i in range(numPoints):
    outfile.write(str(bendList[i])+"\t")
    for j in range(numCurves):
        outfile.write(str(curves[j][i])+"\t")
    outfile.write("\n")
outfile.close()
