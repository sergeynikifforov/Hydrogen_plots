#plot the graph of vapor pressure of water on temperature
import matplotlib.pyplot as plt
import math
from scipy.interpolate import CubicSpline
import numpy as np
f = open('test.txt', 'r')
x = list()
y = list()
for line in f:
    a = line.split(' ')
    x.append(float(a[0]))
    y.append(float(a[1]))
plt.plot(x,y)
plt.xlabel('T, C (Temperature)')
plt.ylabel('P, KpA (Saturated vapor pressure)')
plt.show()
