import matplotlib.pyplot as plt
import math
from scipy import optimize as opt
import cantera as ct

#alpha_new_upd
alpha_new = [1.,0.35,0.43,14.3]
#alpha_new_upd
#alpha_new = [1.,0.35,0.44,6.5]

def from_per_to_alpha(ER,perc):
    return perc*(2*ER+4.76)/(1-perc)
def f_0_ER_H2(ER = 1,alpha = 0):
    return 2*ER/(1+2*ER + alpha + 3.76)
def f_0_ER_O2(ER = 1,alpha = 0):
    return 1/(1+2*ER + alpha + 3.76)
def f_0_ER_H2O(ER = 1,alpha = 0):
    return alpha/(1+2*ER + alpha + 3.76)
def f_0_ER_N2(ER = 1,alpha = 0):
    return 3.76/(1+2*ER + alpha + 3.76)
def f_1_ER_H2(ER = 1,alpha = 0):
    if(ER>1):
        return (2*ER-2)/(2*ER + alpha + 3.76)
    else:
        return 0
def f_1_ER_O2(ER = 1,alpha = 0):
    if(ER>1):
        return 0
    else:
        return (1-ER)/(ER + alpha + 4.76)
def f_1_ER_H2O(ER = 1,alpha = 0):
    if(ER>1):
        return (2+alpha)/(2*ER + alpha + 3.76)
    else:
        return (2*ER+alpha)/(ER + alpha + 4.76)
def f_1_ER_N2(ER = 1,alpha = 0):
    if(ER>1):
        return (3.76)/(2*ER + alpha + 3.76)
    else:
        return (3.76)/(ER + alpha + 4.76)
def obj_func(a,b):
    return (a-b)**2
def f1(value,ER_,P_ = 1,alpha_ = 0):
    k_inf = 4.65e12 * value ** 0.44
    k_0 = 5.75e19 * value ** (-1.4)
    fin = list()
    sum = 0
    fin.append(f_1_ER_H2(ER_,alpha_))
    fin.append(f_1_ER_O2(ER_,alpha_))
    fin.append(f_1_ER_N2(ER_, alpha_))
    fin.append(f_1_ER_H2O(ER_, alpha_))
    for i in range(len(fin)):
        sum += alpha_new[i]*fin[i]
    n = 2.4e19*300*P_/(6.02e23 * value)*sum
    return  (k_inf * k_0 * n)/(k_inf + k_0*n)
def f2(value):
    return 3.52 * math.pow(10,16) * math.pow(value,-0.7) * math.exp(-8590/value)
def f3(ER_,P_= 1,alpha_ = 0):
    min = opt.minimize(lambda x: obj_func(f1(x,ER_,P_,alpha_),f2(x)),900.,method='Nelder-Mead')
    #min = opt.minimize_scalar(lambda x: obj_func(f1(x,ER_,P_,alpha_),f2(x)),900.,method='brent')
    return float(min.x[0])
def f4(ER_,T_0=300,P_=1,alpha_=0):
    gas_new = ct.Solution('gri30.cti')
    gas_new.TPX = T_0,P_*ct.one_atm,'H2:%f, O2:1, N2:3.76, H2O:%f'%(2.*ER_,alpha_)
    gas_new.equilibrate('HP')
    return float(gas_new.T)

print('Please, input the pressure, initial temperature')
P_new, T_new = input().split(' ')
P_new = float(P_new)
T_new = float(T_new)

min_new_1 = list()
min_new_2 = list()
value_H2 = list()
value_H2O = list()

min_new_alpha = opt.minimize(lambda y: obj_func(f3(y,P_new,from_per_to_alpha(1.,y)),f4(1.,T_new,P_new,from_per_to_alpha(1.,y))),0.7,method='Nelder-Mead')
#min_new_alpha = opt.minimize_scalar(lambda y: obj_func(f3(y,P_new,from_per_to_alpha(1.,y)),f4(1.,T_new,P_new,from_per_to_alpha(1.,y))),0.4,method='brent')
alpha_perc_max = int(round(100*round(min_new_alpha.x[0],4)))
perc_new = [i/100. for i in range(alpha_perc_max)]
print(perc_new)
min_new_1_perem = 0.4
min_new_2_perem = 5.0
#temperature of rich/lean limits
#T(ER) only
#out_T = list()
#out_H2 = list()
for i in range(len(perc_new)):
    min_new_1_perem_new = opt.minimize(lambda y: obj_func(f3(y,P_new,from_per_to_alpha(y,perc_new[i])),f4(y,T_new,P_new,from_per_to_alpha(y,perc_new[i]))),min_new_1_perem-min_new_1_perem/3,method='Nelder-Mead')
    #min_new_1_perem_new = opt.minimize_scalar(lambda y: obj_func(f3(y,P_new,from_per_to_alpha(y,perc_new[i])),f4(y,T_new,P_new,from_per_to_alpha(y,perc_new[i]))),min_new_1_perem,method='brent')
    min_new_1.append(min_new_1_perem_new.x[0])
    min_new_1_perem = min_new_1_perem_new.x[0]
    min_new_2_perem_new = opt.minimize(lambda y: obj_func(f3(y,P_new,from_per_to_alpha(y,perc_new[i])),f4(y,T_new,P_new,from_per_to_alpha(y,perc_new[i]))),min_new_2_perem+min_new_1_perem/2,method='Nelder-Mead')
    #min_new_2_perem_new = opt.minimize_scalar(lambda y: obj_func(f3(y,P_new,from_per_to_alpha(y,perc_new[i])),f4(y,T_new,P_new,from_per_to_alpha(y,perc_new[i]))),min_new_2_perem,method='brent')
    min_new_2.append(min_new_2_perem_new.x[0])
    min_new_2_perem = min_new_2_perem_new.x[0]
for i in range(len(min_new_1)):
    value_H2.append(100*f_0_ER_H2(min_new_1[i],from_per_to_alpha(min_new_1[i],perc_new[i])))
    #T(ER) only
    #out_T.append(min_new_1[i])
    #out_H2.append(value_H2[-1])
    value_H2O.append(100*perc_new[i])
for i in range(len(min_new_2)):
    k = len(min_new_2)-i-1
    value_H2.append(100*f_0_ER_H2(min_new_2[k],from_per_to_alpha(min_new_2[k],perc_new[k])))
    #T(ER) only
    #out_T.append(min_new_2[i])
    #out_H2.append(value_H2[-1])
    value_H2O.append(100*perc_new[k])
#373.15
res_x_1 = [0,2.234043,5.265957,7.819149,10.531915,13.404255,15.797872,18.191489,20.744680,23.457447,26.648936,29.680851,32.393617,35.425532,38.617021,41.648936,44.202128,46.755319,49.468085,52.500000,54.893617,57.287234,58.563830,55.053191,49.627660,46.117021,41.808510,37.659574,33.191489,26.329787,20.904255,16.914894,13.404255,8.936170,4.787234,1.9148936,0.159574]
res_y_1 = [75.749208,73.128112,70.063377,67.437755,64.169307,61.321865,58.700770,56.079674,53.241286,50.398370,47.116342,43.626075,40.783160,37.505659,34.010865,30.946129,28.107741,25.269352,22.426437,18.936170,16.315075,13.906745,11.530104,10.140335,9.443187,9.117248,9.026709,8.718877,8.845631,8.401992,8.555908,8.456315,8.768674,8.682662,8.800362,8.881847,8.931643]
#473.15
#res_x_1 = [0,19.148936,38.617021,59.3617021,61.436170,59.042553,39.414894,19.468085,0]
#res_y_1 = [80.221819,58.392938,36.776822,13.635129,11.023087,8.750566,7.605251,6.681756,6.378451]
#375
#res_x_1 = [0,1.812081,3.691275,5.167785,7.046980,8.791946,10.402685,12.013423,13.892617,15.503356,16.979866,18.590604,20.335570,22.214765,23.691275,25.302013,26.510067,28.120805,29.731544,31.073825,32.684564,34.295302,35.771812,36.979866,38.590604,39.798658,41.006711,42.348993,43.557047,44.765101,45.973154,47.181208,48.255033,49.463087]
#res_x_2 = [48.389262,46.510067,44.765101,42.885906,41.275168,39.395973,37.516778,35.771812,34.026846,32.013423,30.402685,28.255033,26.510067,24.765101,22.885906,21.006711,18.993288,17.114094,15.369127,13.489933,11.610738,9.731543,7.852349,6.107382,4.093960,2.483221,1.275168,0]
#res_y_1 = [70,69.230769,68.461538,67.846154,66.769231,66,65.076923,64.307692,63.384615,62.461538,61.538461,60.615385,59.538461,58.153846,57.076923,56,55.076923,53.692308,52.461538,51.384615,50,48.461538,47.076923,45.692308,44.153846,42.769231,41.230769,39.538461,37.846154,36.461538,34.923077,33.230769,31.692307,30]
#res_y_2 = [20.923077,20.153846,19.384615,18.615384,18.153846,17.538461,16.923077,16.461538,15.846154,15.230769,14.769231,14.307692,13.846154,13.538461,13.076923,12.769231,12.461538,12.153846,11.692308,11.538461,11.230769,10.923077,10.615384,10.461538,10.153846,10,9.846154,9.846154]
#un_sq_x = [48.657718,56.442953,56.577181,48.791946,48.657718]
#un_sq_y = [29.538461,29.538461,21.384615,21.384615,29.538461]
#T(ER)
#plt.plot(out_H2, out_T, label='initial temperature {} K'.format(T_new))
plt.plot(value_H2O, value_H2, label='our model with temperature {} K'.format(T_new))
plt.plot(res_x_1, res_y_1, label='theoretical with temperature {} K'.format(T_new),dashes = [2,2],color = 'red')
#375 only
#plt.plot(res_x_2, res_y_2, label='theoretical with temperature {} K'.format(T_new),dashes = [2,2],color = 'red')
#375 only
#plt.plot(un_sq_x, un_sq_y, label='Uncertainty defined by stability boundary')
#T(ER) only
plt.xlabel('H$_{2}$,%')
plt.ylabel('T,K')
#all but T(ER)
#plt.xlabel('H$_{2}$O,%')
#plt.ylabel('H$_{2}$,%')
plt.legend()
plt.show()
