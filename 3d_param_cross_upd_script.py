#defines critical point where 2 roots pass into 1 root (dry mixture)
import matplotlib.pyplot as plt
import math
from scipy import optimize as opt
import cantera as ct
import numpy as np

alpha_new = [2.5,1.,1.,16.,1.2]


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
    min = opt.minimize(lambda x: obj_func(f1(x,ER_,P_,alpha_,beta_),f2(x)),900.,method='Nelder-Mead')
    return float(min.x[0])
def f4(ER_,T_0=300,P_=1,alpha_= 0,beta_ = 0):
    gas_new = ct.Solution('gri30.cti')
    gas_new.TPX = T_0,P_*ct.one_atm,'H2:%f, O2:1, N2:3.76, H2O:%f, CO:%f'%(2.*ER_,alpha_,beta_)
    gas_new.equilibrate('HP')
    return float(gas_new.T)


ER = [i/10. for i in range(0,200)]
print('Please, input the initial temperature, initial pressure')
T_new, P_new = input().split(' ')
P_new = float(P_new)
T_new = [float(T_new)]
perc_new_alpha = [0.]
perc_new_beta = [i/100. for i in range(100)]
for i in range(len(perc_new_beta)):
    perc_new_beta[i] = float(perc_new_beta[i])
assert (perc_new_alpha[0] <= 1. and perc_new_alpha[0] >= 0.), "incorrect percentage of water"
for i in range(len(perc_new_beta)):
    assert (perc_new_beta[i] <= 1. and perc_new_beta[i] >= 0. ), "incorrect percentage of carbon monoxide"
    assert ((perc_new_alpha[0]+perc_new_beta[i]) < 1.) , "incorrect summary percentage"

T_HP_final = list()
for i in range(len(perc_new_beta)):
    T_HP_final.append(func_calc_T_HP(ER,T_new[0],P_new,perc_new_alpha[0],perc_new_beta[i]))
res_final =list()
for i in range(len(perc_new_beta)):
    T_HP_zip = zip(T_HP_final[i],ER)
    T_HP_zip = sorted(T_HP_zip)
    res = T_HP_zip[-1][1]
    res_final.append(res)
min_val_final = list()
for j in range(len(perc_new_beta)):
    min_val = []
    for i in range(len(ER)):
        min = opt.minimize(lambda x: obj_func(
                                                f1( x, ER[i], P_new, alpha_=from_per_to_alpha( ER[i], perc_new_alpha[0], perc_new_beta[j] ), beta_ = from_per_to_beta( ER[i], perc_new_alpha[0] ,perc_new_beta[j] ) ),
                                                f2( x ) ),900.,method='Nelder-Mead')
        min_val.append(min.x[0])
    min_val_final.append(min_val)
control_arr_final = list()
for j in range(len(perc_new_beta)):
    control_arr = list()
    for i in range(len(T_HP_final[0])):
        control_arr.append(T_HP_final[j][i]-min_val_final[j][i])
    control_arr_final.append(control_arr)
for j in range(len(perc_new_beta)-1):
    if(sign_func(control_arr_final[j]) == 2 and sign_func(control_arr_final[j+1]) == 1):
        beta_per = perc_new_beta[j]
        beta_per_index = j
        break
#####
min_new_1 = list()
min_new_2 = list()
fun = lambda y: obj_func(
                        f3( y ,P_new , from_per_to_alpha( y, perc_new_alpha[0], beta_per ), from_per_to_beta( y, perc_new_alpha[0], beta_per ) ),
                        f4( y, T_new[0], P_new, from_per_to_alpha( y, perc_new_alpha[0],beta_per ), from_per_to_beta( y, perc_new_alpha[0] , beta_per ) )
                        ) if y>=0 else np.Inf
min_new_2.append(opt.minimize(fun,res_final[j]-res_final[j]/3,method='Nelder-Mead'))
min_new_1.append(opt.minimize(fun,res_final[j]+res_final[j]/3,method='Nelder-Mead'))

print('fuel-rich limit is', min_new_1[0].x[0],'with temperature:',f4(min_new_1[0].x[0],T_new[0],alpha_=from_per_to_alpha(min_new_1[0].x[0],perc_new_alpha[0],beta_per),beta_=from_per_to_beta(min_new_1[0].x[0],perc_new_alpha[0],beta_per)))
print('fuel-lean limit is', min_new_2[0].x[0],'with temperature:',f4(min_new_2[0].x[0],T_new[0],alpha_=from_per_to_alpha(min_new_2[0].x[0],perc_new_alpha[0],beta_per),beta_=from_per_to_beta(min_new_2[0].x[0],perc_new_alpha[0],beta_per)))
res_initial_H2_0 = f_0_ER_H2( min_new_1[0].x[0],from_per_to_alpha(min_new_1[0].x[0],perc_new_alpha[0],beta_per),from_per_to_beta(min_new_1[0].x[0],perc_new_alpha[0],beta_per))
res_initial_H2_1 = f_0_ER_H2( min_new_2[0].x[0],from_per_to_alpha(min_new_2[0].x[0],perc_new_alpha[0],beta_per),from_per_to_beta(min_new_2[0].x[0],perc_new_alpha[0],beta_per))
res_final_H2_0 = f_1_ER_H2( min_new_1[0].x[0],from_per_to_alpha(min_new_1[0].x[0],perc_new_alpha[0],beta_per),from_per_to_beta(min_new_1[0].x[0],perc_new_alpha[0],beta_per))
res_final_H2_1 = f_1_ER_H2( min_new_2[0].x[0],from_per_to_alpha(min_new_2[0].x[0],perc_new_alpha[0],beta_per),from_per_to_beta(min_new_2[0].x[0],perc_new_alpha[0],beta_per))
print('The initial percentage of H2 for fuel-lean limit is {:.4f}'.format((round(res_initial_H2_0,4)*100)))
print('The initial percentage of H2 for fuel-rich limit is {:.4f}'.format((round(res_initial_H2_1,4)*100)))
print('The final percentage of H2 for fuel-lean limit is {:.4f}'.format((round(res_final_H2_0,4)*100)))
print('The final percentage of H2 for fuel-rich limit is {:.4f}'.format((round(res_final_H2_1,4)*100)))
plt.plot(ER, min_val_final[beta_per_index],label='Crossover with percentage of carbon monoxide {}'.format(beta_per))
plt.plot(ER, T_HP_final[beta_per_index],label='Ad. temp. with percentage of carbon monoxide {}'.format(beta_per),dashes = [2,2])
plt.xlabel('ER')
plt.ylabel('T$_{cross}$, T$_{adiab}$ , K')
plt.legend()
plt.show()
