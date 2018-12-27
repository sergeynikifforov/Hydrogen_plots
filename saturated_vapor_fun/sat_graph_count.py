#count the vapor pressure of water on temperature
import matplotlib.pyplot as plt
import math
from scipy.interpolate import CubicSpline
import numpy as np
f = open('test.txt', 'r')
print('Please, input temperature in Celsium')
T = float(input())
print('Please, input the whole pressure in Kpa')
P = float(input())
x = list()
y = list()
for line in f:
    a = line.split(' ')
    x.append(float(a[0]))
    y.append(float(a[1]))
cs = CubicSpline(x, y)
print('The pressure of saturated vapor in KpA:')
print(cs(T))
