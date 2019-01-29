import matplotlib.pyplot as plt
import math
from scipy import optimize as opt
import cantera as ct
import numpy as np

#alpha_new = [2.5,1.,1.,16.,1.2]

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
def f1(value,ER_,P_ = 1,alpha_ = 0, alpha_new_out = 16.):
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
    fin = f_1_ER_H2O(ER_, alpha_)

    #old_version
    #for i in range(len(fin)):
    #    sum += alpha_new_*fin[i]
    #new_version
    sum = 1. + fin*(alpha_new_out - 1.)
    n = 2.4e19*300*P_/(6.02e23 * value)*sum
    return  (k_inf * k_0 * n)/(k_inf + k_0*n)
def f2(value):
    return 3.52 * math.pow(10,16) * math.pow(value,-0.7) * math.exp(-8590/value)
def f3(ER_,P_= 1,alpha_ = 0, alpha_new_out = 16.):
    min = opt.minimize(lambda x: obj_func(f1(x,ER_,P_,alpha_,alpha_new_out),f2(x)),900.,method='Nelder-Mead')
    #min = opt.minimize_scalar(lambda x: obj_func(f1(x,ER_,P_,alpha_),f2(x)),900.,method='brent')
    return float(min.x[0])
def f4(ER_,T_0=300,P_=1,alpha_=0):
    gas_new = ct.Solution('gri30.cti')
    gas_new.TPX = T_0,P_*ct.one_atm,'H2:%f, O2:1, N2:3.76, H2O:%f'%(2.*ER_,alpha_)
    gas_new.equilibrate('HP')
    return float(gas_new.T)
def func_calc_T_HP(ER_,T_new_,P_new_,perc_new_alpha_):
    #HP_solution
    gas = ct.Solution('gri30.cti')
    T_HP = list()
    for i in range(len(ER_)):
        gas.TPX = T_new_, P_new_*ct.one_atm, 'H2:%f, O2:1, N2:3.76, H2O:%f'%(2*ER_[i],from_per_to_alpha(ER_[i],perc_new_alpha_))
        gas.equilibrate('HP')
        T_HP.append(float(gas.T))
    return T_HP

print('Please, input the pressure ')
alpha_new_out_value = 16.
P_new = input()
P_new = float(P_new)
#T_new = [373., 473., 573., 673., 773., 813.,863., 873., 973., 1073.]
T_new =[973.]
final_x = list()
final_y = list()
T_HP_final = list()
H2_out_plot_final = list()
for i in range(len(T_new)):
    ER = [i/10. for i in range(0,200)]
    perc_new_alpha = [0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7]
    for j in range(len(perc_new_alpha)):
        T_HP_final.append(func_calc_T_HP(ER,T_new[i],P_new,perc_new_alpha[j]))
    res_final =list()
    for m in range(len(perc_new_alpha)):
        H2_out_plot = list()
        for j in range(len(ER)):
            H2_out_plot.append(100*f_0_ER_H2(ER[j],from_per_to_alpha(ER[j],perc_new_alpha[m])))
        H2_out_plot_final.append(H2_out_plot)
        #plt.plot(H2_out_plot,T_HP_final[i],label='H2O:{}'.format(perc_new_alpha[i]))
        plt.plot(H2_out_plot,T_HP_final[i],label='initial_temperature:{}, $\epsilon$: {}, H2O: {}'.format(T_new[0],alpha_new_out_value,perc_new_alpha[m]))
    #plt.scatter(final_y[i], final_x[i],label='CO:{}'.format(perc_beta_input[i]), s = 2)
    plt.xlabel('H$_{2}$,%')
    plt.ylabel('T$_{add}$,K')
    #plt.xlabel('H$_{2}O$,%')
    #plt.ylabel('H$_{2}$,%')
f = open('H2O_add_alpha_16.txt', 'w')
f.write(str(T_new[0]) + '\n')
for i in range(len(perc_new_alpha)):
    f.write(str(perc_new_alpha[i]) + '\n')
    for j in range(len(H2_out_plot_final[i])):
        f.write(str(H2_out_plot_final[i][j]) + ' ' + str(T_HP_final[i][j]) + '\n')
f.close()
plt.legend()
plt.show()
