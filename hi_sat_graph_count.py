#count the mole fraction of vapor of water on temperature
import matplotlib.pyplot as plt
import math
from scipy.interpolate import CubicSpline
import numpy as np
def hi_sat_pres(P_nas,P):
    if(P_nas >= P):
        return 100.
    else:
        return P_nas/P*100
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
hi = list()
for i in range(len(x)):
    hi.append(hi_sat_pres(x[i],P))
cs = CubicSpline(x, y)
print('The mole fraction of vapor, %:')
print(hi_sat_pres(cs(T),P))
