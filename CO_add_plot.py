import matplotlib.pyplot as plt
import math
from scipy import optimize as opt
import cantera as ct
import numpy as np

#alpha_new = [2.5,1.,1.,16.,1.2]

def sign_func(array):
    counter = 0
    for i in range(len(array)-1):
        if((array[i]>=0 and array[i+1]<0) or (array[i]>0 and array[i+1]<=0) or (array[i]<=0 and array[i+1]>0) or (array[i]<0 and array[i+1]>=0)):
            counter+=1
    if(counter == 0):
        return 0
    elif(counter == 1):
        return 1
    else:
        return 2
def func_calc_T_HP(ER_,T_new_,P_new_,perc_new_alpha_,perc_new_beta_):
    #HP_solution
    gas = ct.Solution('gri30.cti')
    T_HP = list()
    for i in range(len(ER_)):
        gas.TPX = T_new_, P_new_*ct.one_atm, 'H2:%f, O2:1, N2:3.76, H2O:%f, CO:%f'%(2*ER_[i],from_per_to_alpha(ER_[i],perc_new_alpha_,perc_new_beta_),from_per_to_beta(ER_[i],perc_new_alpha_,perc_new_beta_))
        gas.equilibrate('HP')
        T_HP.append(float(gas.T))
    return T_HP
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
def f1(value,ER_,P_ = 1,alpha_ = 0, beta_ = 0, alpha_new_out = 16.):
    k_inf = 4.65e12 * value ** 0.44
    k_0 = 5.75e19 * value ** (-1.4)
    #old_version
    #fin = list()

    #new_version
    fin = 0
    sum = 0
    #old_version

    #fin.append(f_1_ER_H2(ER_, alpha_, beta_))
    #fin.append(f_1_ER_O2(ER_, alpha_, beta_))
    #fin.append(f_1_ER_N2(ER_, alpha_, beta_))
    #fin.append(f_1_ER_H2O(ER_, alpha_, beta_))
    #fin.append(f_1_ER_CO(ER_, alpha_, beta_))

    #new version
    fin = f_1_ER_H2O(ER_, alpha_, beta_)

    #old_version
    #for i in range(len(fin)):
    #    sum += alpha_new_*fin[i]
    #new_version
    sum = 1. + fin*(alpha_new_out - 1.)
    n = 2.4e19*300*P_/(6.02e23 * value)*sum
    return  (k_inf * k_0 * n)/(k_inf + k_0*n)
def f2(value):
    return 3.52 * math.pow(10,16) * math.pow(value,-0.7) * math.exp(-8590/value)
def f3(ER_,P_= 1,alpha_ = 0, beta_ = 0, alpha_new_out = 16.):
    min = opt.minimize(lambda x: obj_func(f1(x,ER_,P_,alpha_,beta_,alpha_new_out),f2(x)),900.,method='Nelder-Mead')
    return float(min.x[0])
def f4(ER_,T_0=300,P_=1,alpha_= 0,beta_ = 0):
    gas_new = ct.Solution('gri30.cti')
    gas_new.TPX = T_0,P_*ct.one_atm,'H2:%f, O2:1, N2:3.76, H2O:%f, CO:%f'%(2.*ER_,alpha_,beta_)
    gas_new.equilibrate('HP')
    return float(gas_new.T)

print('Please, input the pressure, initial temperature')
alpha_new_out_value = 16.
P_new, T_new = input().split(' ')
P_new = float(P_new)
T_new = [float(T_new)]
final_x = list()
final_y = list()

perc_beta_input = [0.0, 0.1, 0.3]
T_HP_final = list()
for i in range(len(perc_beta_input)):
    perc_new_beta = perc_beta_input[i]
    ER = [i/10. for i in range(0,200)]
    perc_new_beta = [float(perc_new_beta)]
    lim_up = int(100 - 100*perc_new_beta[0]) - 1
    perc_new_alpha = [0]
    assert (perc_new_beta[0] <= 1. and perc_new_beta[0] >= 0. ), "incorrect percentage of carbon monoxide"
    for i in range(len(perc_new_alpha)):
        assert (perc_new_alpha[i] <= 1. and perc_new_alpha[i] >= 0.), "incorrect percentage of water"
        assert ((perc_new_alpha[i]+perc_new_beta[0]) < 1.) , "incorrect summary percentage"
    T_HP_final.append(func_calc_T_HP(ER,T_new[0],P_new,perc_new_alpha[0],perc_new_beta[0]))
    res_final =list()
for i in range(len(perc_beta_input)):
    H2_out_plot = list()
    for j in range(len(ER)):
        H2_out_plot.append(100*f_0_ER_H2(ER[j],from_per_to_alpha(ER[j],perc_new_alpha[0],perc_beta_input[i]),from_per_to_beta(ER[j],perc_new_alpha[0],perc_beta_input[i])))
    #plt.plot(H2_out_plot,T_HP_final[i],label='CO:{}'.format(perc_beta_input[i]))
    plt.plot(ER,T_HP_final[i],label='CO:{}'.format(perc_beta_input[i]))
    #plt.scatter(final_y[i], final_x[i],label='CO:{}'.format(perc_beta_input[i]), s = 2)
    plt.xlabel('H$_{2}$')
    plt.ylabel('T$_{add}$,K')
    #plt.xlabel('H$_{2}O$,%')
    #plt.ylabel('H$_{2}$,%')
plt.legend()
plt.show()
