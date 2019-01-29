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
#alpha_new_out_value_new = [10.,16]
alpha_new_out_value_new = [16.]
P_new, T_new = input().split(' ')
P_new = float(P_new)
T_new = [float(T_new)]
final_x = list()
final_y = list()
perc_beta_input = [0.0, 0.02, 0.04, 0.06, 0.08, 0.1]
final_x_2 = list()
final_y_2 = list()
graph_max_x = list()
graph_max_y = list()

for m in range(len(alpha_new_out_value_new)):
    alpha_new_out_value = alpha_new_out_value_new[m]
    #perc_beta_input = [0.2, 0.3]
    for i in range(len(perc_beta_input)):
        perc_new_beta = perc_beta_input[i]
        ER = [i/10. for i in range(0,200)]
        perc_new_beta = [float(perc_new_beta)]
        lim_up = int(100 - 100*perc_new_beta[0]) - 1
        perc_new_alpha = [float(i/100.) for i in range(0,lim_up)]
        assert (perc_new_beta[0] <= 1. and perc_new_beta[0] >= 0. ), "incorrect percentage of carbon monoxide"
        for i in range(len(perc_new_alpha)):
            assert (perc_new_alpha[i] <= 1. and perc_new_alpha[i] >= 0.), "incorrect percentage of water"
            assert ((perc_new_alpha[i]+perc_new_beta[0]) < 1.) , "incorrect summary percentage"
        T_HP_final = list()
        for i in range(len(perc_new_alpha)):
            T_HP_final.append(func_calc_T_HP(ER,T_new[0],P_new,perc_new_alpha[i],perc_new_beta[0]))
        res_final =list()
        for i in range(len(perc_new_alpha)):
            T_HP_zip = zip(T_HP_final[i],ER)
            T_HP_zip = sorted(T_HP_zip)
            res = T_HP_zip[-1][1]
            res_final.append(res)
        min_val_final = list()
        for j in range(len(perc_new_alpha)):
            min_val = []
            for k in range(len(ER)):
                min = opt.minimize(lambda x: obj_func(
                                                        f1( x, ER[k], P_new, alpha_=from_per_to_alpha( ER[k], perc_new_alpha[j] , perc_new_beta[0] ), beta_ = from_per_to_beta( ER[k] , perc_new_alpha[j] , perc_new_beta[0] ) , alpha_new_out = alpha_new_out_value) ,
                                                        f2( x ) ),900.,method='Nelder-Mead')
                min_val.append(min.x[0])
            min_val_final.append(min_val)
        control_arr_final = list()
        for j in range(len(perc_new_alpha)):
            control_arr = list()
            for i in range(len(T_HP_final[0])):
                control_arr.append(T_HP_final[j][i]-min_val_final[j][i])
            control_arr_final.append(control_arr)
        min_new_1 = list()
        min_new_2 = list()
        alpha_out_r = list()
        alpha_out_l = list()
        out_T_r = list()
        out_T_l = list()
        out_H2_r = list()
        out_H2_l = list()
        value_H2 = list()
        value_H2O_l = list()
        value_H2O_r = list()
        #print('The efficiency coefficient for H2 is {}, for O2 is {}, for N2 is {}, for H2O is {}, for CO is {}'.format(alpha_new[0],alpha_new[1],alpha_new[2],alpha_new[3],alpha_new[4]))
        for j in range(len(perc_new_alpha)):
            if(sign_func(control_arr_final[j]) == 1):
                fun = lambda y: obj_func(
                                          f3(  y, P_new, from_per_to_alpha( y, perc_new_alpha[j], perc_new_beta[0] ), from_per_to_beta( y, perc_new_alpha[j], perc_new_beta[0]  ) , alpha_new_out_value ),
                                          f4(  y, T_new[0], P_new, from_per_to_alpha( y, perc_new_alpha[j], perc_new_beta[0] ), from_per_to_beta( y, perc_new_alpha[j], perc_new_beta[0] )  )
                                        ) if y>=0 else np.Inf

                min_new_1_perem_new =  opt.minimize( fun, res_final[j]+res_final[j]/3, method='Nelder-Mead' )
                min_new_1.append(float(min_new_1_perem_new.x[0]))
                alpha_out_l.append(perc_new_alpha[j])
            elif(sign_func(control_arr_final[j]) == 2):

                fun = lambda y: obj_func(f3(y,P_new,from_per_to_alpha(y,perc_new_alpha[j],perc_new_beta[0]),from_per_to_beta(y,perc_new_alpha[j],perc_new_beta[0]),alpha_new_out_value),f4(y,T_new[0],P_new,from_per_to_alpha(y,perc_new_alpha[j],perc_new_beta[0]),from_per_to_beta(y,perc_new_alpha[j],perc_new_beta[0]))) if y>=0 else np.Inf
                min_new_2_perem_new =  opt.minimize( fun, res_final[j]-res_final[j]/2, method='Nelder-Mead' )
                min_new_2.append(float(min_new_2_perem_new.x[0]))
                alpha_out_l.append(perc_new_alpha[j])
                min_new_1_perem_new =  opt.minimize( fun, res_final[j]+res_final[j]/2, method='Nelder-Mead' )
                min_new_1.append(float(min_new_1_perem_new.x[0]))
                alpha_out_r.append(perc_new_alpha[j])


        for i in range(len(min_new_1)):
            val_new = from_per_to_beta(min_new_1[i],alpha_out_l[i],perc_new_beta[0])
            value_H2.append(100*f_0_ER_H2(min_new_1[i],from_per_to_alpha(min_new_1[i],alpha_out_l[i],perc_new_beta[0]),val_new))
            value_H2O_l.append(100*alpha_out_l[i])
            #T(ER) only
            out_T_l.append(f4(min_new_1[i],T_new[0],P_new,alpha_out_l[i],val_new))
            out_H2_l.append(value_H2[-1])
        for i in range(len(min_new_2)):
            k = len(min_new_2)-i-1
            val = from_per_to_beta(min_new_2[k],alpha_out_r[k],perc_new_beta[0])
            value_H2.append(100*f_0_ER_H2(min_new_2[k],from_per_to_alpha(min_new_2[k],alpha_out_r[k],perc_new_beta[0]),val))
            value_H2O_r.append(100*alpha_out_r[k])
            #T(ER) only
            out_T_r.append(f4(min_new_2[k],T_new[0],P_new,alpha_out_r[i],val))
            out_H2_r.append(value_H2[-1])
        perem_x = out_H2_l
        for i in range(len(out_H2_r)):
            perem_x.append(out_H2_r[i])
        perem_y = value_H2O_l
        for i in range(len(value_H2O_r)):
            perem_y.append(value_H2O_r[i])
        final_x.append(perem_x)
        final_y.append(perem_y)
        final_res_val = sorted(zip(final_y[-1],final_x[-1]))
        print(final_res_val[-1][0],final_res_val[-1][1])
        print(np.max(final_y[-1]),alpha_new_out_value)
        graph_max_x.append(final_res_val[-1][0])
        graph_max_y.append(final_res_val[-1][1])
for i in range(len(perc_beta_input)):
    plt.scatter(final_y[i], final_x[i],label='CO:{}, $\epsilon$: {}'.format(perc_beta_input[i], alpha_new_out_value_new[0]), s = 2)
plt.plot(graph_max_x, graph_max_y, label = 'max_points')
plt.xlabel('H$_{2}O$,%')
plt.ylabel('H$_{2}$,%')
plt.legend()
plt.show()
