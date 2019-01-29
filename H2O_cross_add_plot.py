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

print('Please, input the pressure')
P_new = input().split(' ')
print('Please, input the initial temperature')
P_new = float(P_new)
for i in range(len(T_new)):
    T_new[i] = float(T_new[i])

alpha_new_out_value_new = 10.
min_new_1 = list()
min_new_2 = list()
#value_H2 = list()
#value_H2O = list()
for i in range(len(T_new)):
    #temperature of rich/lean limits
    #T(ER) only
    out_T_r = list()
    out_T_l = list()
    out_H2_r = list()
    out_H2_l = list()
    T_HP_final = list()
    for i in range(len(T_new)):
        T_HP_final.append(func_calc_T_HP(ER,T_new[i],P_new,perc_new_alpha[0],perc_new_beta[0]))
    res_final =list()
    for i in range(len(T_new)):
        T_HP_zip = zip(T_HP_final[i],ER)
        T_HP_zip = sorted(T_HP_zip)
        res = T_HP_zip[-1][1]
        res_final.append(res)
    min_val_final = list()
    for j in range(len(T_new)):
        min_val = []
        for i in range(len(ER)):
            min = opt.minimize(lambda x: obj_func(
                                                f1( x, ER[i], P_new, alpha_=from_per_to_alpha( ER[i], perc_new_alpha[0], alpha_new_out_value_new)),
                                                f2( x ) ),900.,method='Nelder-Mead')
            min_val.append(min.x[0])
        min_val_final.append(min_val)
    control_arr_final = list()
    for j in range(len(T_new)):
        control_arr = list()
        for i in range(len(T_HP_final[0])):
            control_arr.append(T_HP_final[j][i]-min_val_final[j][i])
        control_arr_final.append(control_arr)
    min_new_1_final = list()
    min_new_2_final = list()
    #print('The efficiency coefficient for H2 is {}, for O2 is {}, for N2 is {}, for H2O is {}, for CO is {}'.format(alpha_new[0],alpha_new[1],alpha_new[2],alpha_new[3],alpha_new[4]))
    for j in range(len(T_new)):
        min_new_1 = list()
        min_new_2 = list()
        if(sign_func(control_arr_final[j]) == 0):
            print('There is no fuel-lean limit or fuel-rich limit')
        elif(sign_func(control_arr_final[j]) == 1):
            min_new_1 = list()

            fun = lambda y: obj_func(
                                      f3(  y, P_new, from_per_to_alpha( y, perc_new_alpha[0]), alpha_new_out_value_new ),
                                      f4(  y, T_new[j], P_new, from_per_to_alpha( y, perc_new_alpha[0]), alpha_new_out_value_new )
                                    ) if y>=0 else np.Inf

            min_new_1.append(  opt.minimize( fun, res_final[j]+res_final[j]/3, method='Nelder-Mead' )  )

            print('There is no fuel-lean limit')
            print('fuel-rich limit is', min_new_1[0].x[0],'with temperature:',f4(min_new_1[0].x[0],T_new[j],alpha_=from_per_to_alpha(min_new_1[0].x[0],perc_new_alpha[0])))
            res_initial_H2_0 = f_0_ER_H2( min_new_1[0].x[0],from_per_to_alpha(min_new_1[0].x[0],perc_new_alpha[0]))
            res_final_H2_0 = f_1_ER_H2( min_new_1[0].x[0],from_per_to_alpha(min_new_1[0].x[0],perc_new_alpha[0]))
            print('The initial percentage of H2 for fuel-rich limit is {:.4f}'.format((round(res_initial_H2_0,4)*100)))
            print('The final percentage of H2 for fuel-rich limit is {:.4f}'.format((round(res_final_H2_0,4)*100)))
        else:

            fun = lambda y: obj_func(
                                        f3( y, P_new, from_per_to_alpha( y, perc_new_alpha[0]),from_per_to_beta( y, perc_new_alpha[0]] ) ),
                                        f4( y, T_new[j], P_new, from_per_to_alpha( y, perc_new_alpha[0] ),from_per_to_beta( y, perc_new_alpha[0] ) )
                                        ) if y>=0 else np.Inf
            min_new_2.append(opt.minimize(fun,res_final[j]-res_final[j]/3,method='Nelder-Mead'))
            min_new_1.append(opt.minimize(fun,res_final[j]+res_final[j]/3,method='Nelder-Mead'))

            print('fuel-rich limit is', min_new_1[0].x[0],'with temperature:',f4(min_new_1[0].x[0],T_new[j],alpha_=from_per_to_alpha(min_new_1[0].x[0],perc_new_alpha[0],perc_new_beta[0]))))
            print('fuel-lean limit is', min_new_2[0].x[0],'with temperature:',f4(min_new_2[0].x[0],T_new[j],alpha_=from_per_to_alpha(min_new_2[0].x[0],perc_new_alpha[0],perc_new_beta[0]))))
            res_initial_H2_0 = f_0_ER_H2( min_new_1[0].x[0],from_per_to_alpha(min_new_1[0].x[0],perc_new_alpha[0],perc_new_beta[0]))
            res_initial_H2_1 = f_0_ER_H2( min_new_2[0].x[0],from_per_to_alpha(min_new_2[0].x[0],perc_new_alpha[0],perc_new_beta[0]))
            res_final_H2_0 = f_1_ER_H2( min_new_1[0].x[0],from_per_to_alpha(min_new_1[0].x[0],perc_new_alpha[0],perc_new_beta[0]))
            res_final_H2_1 = f_1_ER_H2( min_new_2[0].x[0],from_per_to_alpha(min_new_2[0].x[0],perc_new_alpha[0],perc_new_beta[0]))
            print('The initial percentage of H2 for fuel-lean limit is {:.4f}'.format((round(res_initial_H2_1,4)*100)))
            print('The initial percentage of H2 for fuel-rich limit is {:.4f}'.format((round(res_initial_H2_0,4)*100)))
            print('The final percentage of H2 for fuel-lean limit is {:.4f}'.format((round(res_final_H2_1,4)*100)))
            print('The final percentage of H2 for fuel-rich limit is {:.4f}'.format((round(res_final_H2_0,4)*100)))
        min_new_1_final.append(min_new_1)
        min_new_2_final.append(min_new_2)
for i in range(len(T_new)):
    plt.plot(ER, min_val_final, label='our model with temperature {} K, $\epsilon$: {}'.format(T_new[i],alpha_new_out_value_new), s = 4)
plt.xlabel('H$_{2}$O,%')
plt.ylabel('H$_{2}$,%')
plt.legend()
plt.show()
