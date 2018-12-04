import matplotlib.pyplot as plt
import math
import numpy as np
from scipy import optimize as opt
import cantera as ct


alpha_new = [2.5,1.,1.,16.,1.2]
def sign_func(array):
    if((np.max(array)>0 and np.min(array)>0) or (np.max(array)<0 and np.min(array)<0)):
        return True
    else:
        return False
def min_val(a,b):
    if(a>b):
        return b
    else:
        return a
def from_per_to_alpha(ER, perc_alpha = 0, perc_beta = 0):
    return perc_alpha*(4.76+2*ER)/(1-perc_alpha-perc_beta)
def from_per_to_beta(ER, perc_alpha = 0, perc_beta = 0):
    return perc_beta*(4.76+2*ER)/(1-perc_alpha-perc_beta)
def f_0_ER_H2(ER = 1,alpha = 0, beta = 0):
    return 2*ER/(4.76+2*ER + alpha + beta)
def f_0_ER_O2(ER = 1,alpha = 0, beta = 0):
    return 1/(4.76+2*ER + alpha + beta)
def f_0_ER_H2O(ER = 1,alpha = 0, beta = 0):
    return alpha/(4.76+2*ER + alpha + beta)
def f_0_ER_N2(ER = 1,alpha = 0, beta = 0):
    return 3.76/(4.76+2*ER + alpha + beta)
def f_0_ER_CO(ER = 1,alpha = 0, beta = 0):
    return beta/(4.76+2*ER + alpha + beta)
def f_0_ER_CO2(ER = 1,alpha = 0, beta = 0):
    return 0.
def f_1_ER_H2(ER = 1,alpha = 0, beta = 0):
    if(ER>1):
        return (2*ER-2)/(2*ER + alpha + 3.76 + beta)
    else:
        return 0.
def f_1_ER_O2(ER = 1,alpha = 0, beta = 0):
    if(ER>1):
        return 0.
    else:
        if(beta <= (2-2*ER)):
            return (1-ER-beta/2)/(ER + beta/2 + alpha + 4.76)
        else:
            return 0.
def f_1_ER_H2O(ER = 1,alpha = 0, beta = 0):
    if(ER>1):
        return (2+alpha)/(2*ER + beta + alpha + 3.76)
    else:
        if(beta <= (2-2*ER)):
            return (2*ER + alpha)/(ER + beta/2 + alpha + 4.76)
        else:
            return (alpha + 2*ER)/(2*ER + alpha + beta + 3.76)
def f_1_ER_N2(ER = 1,alpha = 0, beta = 0):
    if(ER>1):
        return (3.76)/(2*ER + alpha + 3.76 + beta)
    else:
        if(beta <= (2-2*ER)):
            return ((3.76)/(ER + beta/2 + alpha + 4.76))
        else:
            return (3.76/(2*ER + alpha + beta + 3.76))
def f_1_ER_CO(ER = 1,alpha = 0, beta = 0):
    if(ER>1):
        return beta/(2*ER + alpha + 3.76 + beta)
    else:
        if(beta <= (2-2*ER)):
            return 0.
        else:
            return ((beta-2+2*ER)/(2*ER + alpha + beta + 3.76))
def f_1_ER_CO2(ER = 1,alpha = 0, beta = 0):
    if(ER>1):
        return 0.
    else:
        if(beta <= (2-2*ER)):
            return (beta/(ER + beta/2 + alpha + 4.76))
        else:
            return ((2-2*ER)/(2*ER + alpha + beta + 3.76))
def obj_func(a,b):
    return (a-b)**2
def f1(value,ER_,P_ = 1,alpha_ = 0, beta_ = 0):
    k_inf = 4.65e12 * value ** 0.44
    k_0 = 5.75e19 * value ** (-1.4)
    fin = list()
    sum = 0
    fin.append(f_1_ER_H2(ER_, alpha_, beta_))
    fin.append(f_1_ER_O2(ER_, alpha_, beta_))
    fin.append(f_1_ER_N2(ER_, alpha_, beta_))
    fin.append(f_1_ER_H2O(ER_, alpha_, beta_))
    fin.append(f_1_ER_CO(ER_, alpha_, beta_))
    for i in range(len(fin)):
        sum += alpha_new[i]*fin[i]
    n = 2.4e19*300*P_/(6.02e23 * value)*sum
    return  (k_inf * k_0 * n)/(k_inf + k_0*n)
def f2(value):
    return 3.52 * math.pow(10,16) * math.pow(value,-0.7) * math.exp(-8590/value)
def f3(ER_,P_= 1,alpha_ = 0, beta_ = 0):
    min = opt.minimize(lambda x: obj_func(f1(x,ER_,P_,alpha_,beta_),f2(x)),800.,method='Nelder-Mead')
    return float(min.x[0])
def f4(ER_,T_0=300,P_=1,alpha_= 0,beta_ = 0):
    gas_new = ct.Solution('gri30.cti')
    gas_new.TPX = T_0,P_*ct.one_atm,'H2:%f, O2:1, N2:3.76, H2O:%f, CO:%f'%(2.*ER_,alpha_,beta_)
    gas_new.equilibrate('HP')
    return float(gas_new.T)


print('Please, input the pressure, percentage of water from zero to one,percentage of carbon monoxide from zero to one,initial temperature')
P_new, perc_new_alpha, perc_new_beta ,T_new = input().split(' ')
perc_new_alpha = float(perc_new_alpha)
perc_new_beta = float(perc_new_beta)
P_new = float(P_new)
T_new = float(T_new)

assert (perc_new_alpha <= 1. and perc_new_alpha >= 0.), "incorrect percentage of water"
assert (perc_new_beta <= 1. and perc_new_beta >= 0. ), "incorrect percentage of carbon monoxide"
assert ((perc_new_alpha+perc_new_beta) <= 1.) , "incorrect summary percentage"



ER = [i/10. for i in range(0,100)]
#HP_solution
gas = ct.Solution('gri30.cti')
H2conc = [value for value in ER]
T_HP = list()
for i in range(len(H2conc)):
    gas.TPX = T_new, P_new*ct.one_atm, 'H2:%f, O2:1, N2:3.76, H2O:%f, CO:%f'%(2*H2conc[i],from_per_to_alpha(H2conc[i],perc_new_alpha,perc_new_beta),from_per_to_beta(H2conc[i],perc_new_alpha,perc_new_beta))
    gas.equilibrate('HP')
    T_HP.append(float(gas.T))


row1_1 = [i for i in range(900,1200)]
ans_row1_1 = [f1(val,ER[10],P_new,alpha_=from_per_to_alpha(ER[10],perc_new_alpha,perc_new_beta),beta_ = from_per_to_beta(ER[10],perc_new_alpha,perc_new_beta)) for val in row1_1]

row2 = [i for i in range(900,1200)]
ans_row2 = [f2(val) for val in row2]

min_val = []
for i in range(len(ER)):
    min = opt.minimize(lambda x: obj_func(f1(x,ER[i],P_new,alpha_=from_per_to_alpha(ER[i],perc_new_alpha,perc_new_beta),beta_ = from_per_to_beta(ER[10],perc_new_alpha,perc_new_beta)),f2(x)),900.,method='Nelder-Mead')
    min_val.append(min.x[0])


control_arr = list()
for i in range(len(T_HP)):
    control_arr.append(T_HP[i]-min_val[i])
min_new_1 = list()
min_new_2 = list()
if(sign_func(control_arr)):
    min_new_1.append(-1)
    min_new_2.append(-1)
    print('There is no fuel-lean limit or fuel-rich limit')
    print('The efficiency coefficient for H2 is {}, for O2 is {}, for N2 is {}, for H2O is {}, for CO is {}'.format(alpha_new[0],alpha_new[1],alpha_new[2],alpha_new[3],alpha_new[4]))
else:
    #for i in range(len(alpha)):
    bouds_new_1 = [(0,res)]
    bouds_new_2 = [(res,np.max(ER))]
    fun = lambda y: obj_func(f3(y,P_new,from_per_to_alpha(y,perc_new_alpha,perc_new_beta),from_per_to_beta(y,perc_new_alpha,perc_new_beta)),f4(y,T_new,P_new,from_per_to_alpha(y,perc_new_alpha,perc_new_beta),from_per_to_beta(y,perc_new_alpha,perc_new_beta))) if y>=0 else np.Inf
    min_new_2.append(opt.minimize(fun,res-res/3,method='Nelder-Mead'))
    min_new_1.append(opt.minimize(fun,res+res/3,method='Nelder-Mead'))
    print('fuel-lean limit is', min_new_1[0].x[0],'with temperature:',f3(min_new_1[0].x[0],alpha_=from_per_to_alpha(min_new_1[0].x[0],perc_new_alpha,perc_new_beta),beta_=from_per_to_beta(min_new_1[0].x[0],perc_new_alpha,perc_new_beta)))
    print('fuel-rich limit is', min_new_2[0].x[0],'with temperature:',f3(min_new_2[0].x[0],alpha_=from_per_to_alpha(min_new_2[0].x[0],perc_new_alpha,perc_new_beta),beta_=from_per_to_beta(min_new_2[0].x[0],perc_new_alpha,perc_new_beta)))
    res_initial_H2_0 = f_0_ER_H2( min_new_1[0].x[0],from_per_to_alpha(min_new_1[0].x[0],perc_new_alpha,perc_new_beta),from_per_to_beta(min_new_1[0].x[0],perc_new_alpha,perc_new_beta))
    res_initial_H2_1 = f_0_ER_H2( min_new_2[0].x[0],from_per_to_alpha(min_new_2[0].x[0],perc_new_alpha,perc_new_beta),from_per_to_beta(min_new_2[0].x[0],perc_new_alpha,perc_new_beta))
    res_final_H2_0 = f_1_ER_H2( min_new_1[0].x[0],from_per_to_alpha(min_new_1[0].x[0],perc_new_alpha,perc_new_beta),from_per_to_beta(min_new_1[0].x[0],perc_new_alpha,perc_new_beta))
    res_final_H2_1 = f_1_ER_H2( min_new_2[0].x[0],from_per_to_alpha(min_new_2[0].x[0],perc_new_alpha,perc_new_beta),from_per_to_beta(min_new_2[0].x[0],perc_new_alpha,perc_new_beta))
    print('The initial percentage of H2 for fuel-lean limit is {:.4f}'.format((round(res_initial_H2_0,4)*100)))
    print('The initial percentage of H2 for fuel-rich limit is {:.4f}'.format((round(res_initial_H2_1,4)*100)))
    print('The final percentage of H2 for fuel-lean limit is {:.4f}'.format((round(res_final_H2_0,4)*100)))
    print('The final percentage of H2 for fuel-rich limit is {:.4f}'.format((round(res_final_H2_1,4)*100)))
    print('The efficiency coefficient for H2 is {}, for O2 is {}, for N2 is {}, for H2O is {}, for CO is {}'.format(alpha_new[0],alpha_new[1],alpha_new[2],alpha_new[3],alpha_new[4]))
     #plot
plt.subplot(211)
plt.plot(row1_1, ans_row1_1, label='k1')
plt.plot(row2, ans_row2, label='k2')
plt.xlabel('T')
plt.ylabel('k1 and k2, ER=1.0')
plt.legend()
plt.subplot(212)
plt.plot(ER, min_val,label='Crossover')
plt.plot(ER, T_HP,label='T_HP')
plt.xlabel('ER')
plt.ylabel('T_cross,T_HP')
plt.legend()
plt.show()
