import matplotlib.pyplot as plt
import math
from scipy import optimize as opt
import cantera as ct

alpha_new = [1.,0.35,0.43,14.3]

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
    min = opt.minimize(lambda x: obj_func(f1(x,ER_,P_,alpha_),f2(x)),800.,method='Nelder-Mead')
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
for i in range(len(perc_new)):
    min_new_1_perem_new = opt.minimize(lambda y: obj_func(f3(y,P_new,from_per_to_alpha(y,perc_new[i])),f4(y,T_new,P_new,from_per_to_alpha(y,perc_new[i]))),min_new_1_perem-min_new_1_perem/3,method='Nelder-Mead')
    #min_new_1_perem_new = opt.minimize_scalar(lambda y: obj_func(f3(y,P_new,from_per_to_alpha(y,perc_new[i])),f4(y,T_new,P_new,from_per_to_alpha(y,perc_new[i]))),min_new_1_perem,method='brent')
    min_new_1.append(min_new_1_perem_new.x[0])
    min_new_1_perem = min_new_1_perem_new.x[0]
    min_new_2_perem_new = opt.minimize(lambda y: obj_func(f3(y,P_new,from_per_to_alpha(y,perc_new[i])),f4(y,T_new,P_new,from_per_to_alpha(y,perc_new[i]))),min_new_2_perem+min_new_1_perem,method='Nelder-Mead')
    #min_new_2_perem_new = opt.minimize_scalar(lambda y: obj_func(f3(y,P_new,from_per_to_alpha(y,perc_new[i])),f4(y,T_new,P_new,from_per_to_alpha(y,perc_new[i]))),min_new_2_perem,method='brent')
    min_new_2.append(min_new_2_perem_new.x[0])
    min_new_2_perem = min_new_2_perem_new.x[0]
for i in range(len(min_new_1)):
    value_H2.append(f_0_ER_H2(min_new_1[i],from_per_to_alpha(min_new_1[i],perc_new[i])))
    value_H2O.append(perc_new[i])
for i in range(len(min_new_2)):
    k = len(min_new_2)-i-1
    value_H2.append(f_0_ER_H2(min_new_2[k],from_per_to_alpha(min_new_2[k],perc_new[k])))
    value_H2O.append(perc_new[k])
plt.plot(value_H2O, value_H2)
plt.show()
