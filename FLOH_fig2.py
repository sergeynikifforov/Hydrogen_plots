import matplotlib.pyplot as plt
import math
from scipy import optimize as opt
import cantera as ct
import numpy as np

alpha_new = [2.5,1.,1.,16.,1.2]


def f_x_calc(ER = 1,alpha = 0, beta = 0):
    return 2*ER/(2*ER + beta)

def f_y_calc(ER = 1,alpha = 0, beta = 0):
    return (2*ER+beta)/(2*ER + beta + 4.76 + alpha)

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
    '''
    if((np.max(array)>0 and np.min(array)>0) or (np.max(array)<0 and np.min(array)<0)):
        return 0
    elif(array[0]>0):
        return 1
    else:
        return 2
    '''

def func_calc_T_HP(ER_,T_new_,P_new_,perc_new_alpha_,perc_new_beta_):
    #HP_solution
    gas = ct.Solution('gri30.cti')
    T_HP = list()
    for i in range(len(ER_)):
        gas.TPX = T_new_, P_new_*ct.one_atm, 'H2:%f, O2:1, N2:3.76, H2O:%f, CO:%f'%( 2*ER_[i],
                                                                                     from_per_to_alpha(ER_[i],perc_new_alpha_,perc_new_beta_),
                                                                                     from_per_to_beta(ER_[i],perc_new_alpha_,perc_new_beta_)
                                                                                   )
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

#termination rate constant
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

#branching rate constant
def f2(value):
    return 3.52 * math.pow(10,16) * math.pow(value,-0.7) * math.exp(-8590/value)

#Tcross
def f3(ER_,P_= 1,alpha_ = 0, beta_ = 0):
    min = opt.minimize(lambda x: obj_func(  f1( x, ER_, P_, alpha_, beta_ ), f2(x) )  ,900. ,method='Nelder-Mead' )
    return float(min.x[0])

#Tad
def f4(ER_,T_0=300,P_=1,alpha_= 0,beta_ = 0):
    gas_new = ct.Solution('gri30.cti')
    gas_new.TPX = T_0,P_*ct.one_atm,'H2:%f, O2:1, N2:3.76, H2O:%f, CO:%f'%(2.*ER_,alpha_,beta_)
    gas_new.equilibrate('HP')
    return float(gas_new.T)

#def run():

ER = [i/10. for i in range(0,200)]
perc_new_beta = [i/100. for i in range(0,100)]
perc_new_alpha = [0.0] #[0.]
P_new = 1.

print('initial temperature')
T_new = input()
T_new = [float(T_new)]

assert (perc_new_alpha[0] <= 1. and perc_new_alpha[0] >= 0.), "incorrect percentage of water"

for i in range(len(perc_new_beta)):
    assert (perc_new_beta[i] <= 1. and perc_new_beta[i] >= 0. ), "incorrect percentage of carbon monoxide"
    assert ((perc_new_alpha[0]+perc_new_beta[i]) < 1.) , "incorrect summary percentage"

#list of adiabatic curves
T_HP_final = list()

for i in range(len(perc_new_beta)):
    T_HP_final.append(  func_calc_T_HP( ER, T_new[0], P_new, perc_new_alpha[0], perc_new_beta[i] )  )

res_final =list()

#find maxT of adatabic curve
for i in range(len(perc_new_beta)):
    T_HP_zip = zip(T_HP_final[i],ER)
    T_HP_zip = sorted(T_HP_zip)
    res = T_HP_zip[-1][1]
    res_final.append(res)

min_val_final = list()

#construct Tcross curves
for j in range(len(perc_new_beta)):
    min_val = []
    for i in range(len(ER)):
        min = opt.minimize( lambda x: obj_func(
                                               f1( x,
                                                   ER[i],
                                                   P_new,
                                                   alpha_=from_per_to_alpha(ER[i],perc_new_alpha[0],perc_new_beta[j]),
                                                   beta_ = from_per_to_beta(ER[i],perc_new_alpha[0],perc_new_beta[j])
                                                 ),
                                               f2(x)
                                              ),
                             900.,
                             method='Nelder-Mead'
                          )
        min_val.append(min.x[0])
    min_val_final.append(min_val)

control_arr_final = list()

#difference Tad-Tcross
for j in range(len(perc_new_beta)):
    control_arr = list()
    for i in range(len(T_HP_final[0])):
        control_arr.append(T_HP_final[j][i]-min_val_final[j][i])
    control_arr_final.append(control_arr)

min_new_1_final = list()       #rich limit
min_new_2_final = list()       #lean limit

perc_new_beta_min1 = list()
perc_new_beta_min2 = list()

for j in range(len(perc_new_beta)):
    if(sign_func(control_arr_final[j]) == 1):

        min_new_1 = list()

        fun = lambda y: obj_func(
                                  f3(  y, P_new, from_per_to_alpha( y, perc_new_alpha[0], perc_new_beta[j] ), from_per_to_beta( y, perc_new_alpha[0], perc_new_beta[j] )  ),
                                  f4(  y, T_new[0], P_new, from_per_to_alpha( y, perc_new_alpha[0], perc_new_beta[j] ), from_per_to_beta( y, perc_new_alpha[0], perc_new_beta[j] )  )
                                ) if y>=0 else np.Inf

        min_new_1.append(  opt.minimize( fun, res_final[j]+res_final[j]/3, method='Nelder-Mead' )  )

        #arrays of intersection points
        min_new_1_final.append(min_new_1[0].x[0])
        perc_new_beta_min1.append(perc_new_beta[j])


    elif(sign_func(control_arr_final[j]) == 2):

        min_new_1 = list()
        min_new_2 = list()

        fun = lambda y: obj_func(
                                  f3(  y, P_new, from_per_to_alpha( y, perc_new_alpha[0], perc_new_beta[j] ), from_per_to_beta( y, perc_new_alpha[0], perc_new_beta[j] )  ),
                                  f4(  y, T_new[0], P_new, from_per_to_alpha( y, perc_new_alpha[0], perc_new_beta[j] ), from_per_to_beta( y, perc_new_alpha[0], perc_new_beta[j] )  )
                                ) if y>=0 else np.Inf

        min_new_2.append(  opt.minimize( fun, res_final[j]-res_final[j]/3 , method='Nelder-Mead' )  )
        min_new_1.append(  opt.minimize( fun, res_final[j]+res_final[j]/3, method='Nelder-Mead' )  )

        #arrays of intersection points
        min_new_1_final.append(min_new_1[0].x[0])
        min_new_2_final.append(min_new_2[0].x[0])
        perc_new_beta_min1.append(perc_new_beta[j])
        perc_new_beta_min2.append(perc_new_beta[j])
        #print(min_new_1_final[-1]-min_new_2_final[-1],perc_new_beta[j] )

x_new = list()
y_new = list()
z_new = list()

for j in range(len(min_new_1_final)):
    alpha_fin = from_per_to_alpha( min_new_1_final[j], perc_new_alpha[0], perc_new_beta_min1[j] )
    beta_fin  = from_per_to_beta( min_new_1_final[j], perc_new_alpha[0], perc_new_beta_min1[j] )

    #ER of intersection point
    val_X = f_x_calc( min_new_1_final[j], alpha_fin, beta_fin )
    #T of intersection point
    val_Y = f_y_calc( min_new_1_final[j], alpha_fin, beta_fin )

    x_new.append(val_X)
    y_new.append(val_Y)
    z_new.append(perc_new_beta_min1[j])


x2_new = list()
y2_new = list()
z2_new = list()

for j in range(len(min_new_2_final)):
    alpha_fin = from_per_to_alpha( min_new_2_final[j], perc_new_alpha[0], perc_new_beta_min2[j] )
    beta_fin  = from_per_to_beta( min_new_2_final[j], perc_new_alpha[0], perc_new_beta_min2[j] )

    #ER of intersection point
    val_X = f_x_calc( min_new_2_final[j], alpha_fin, beta_fin )
    #T of intersection point
    val_Y = f_y_calc( min_new_2_final[j], alpha_fin, beta_fin )

    x2_new.append(val_X)
    y2_new.append(val_Y)
    z2_new.append(perc_new_beta_min2[j])



#plotting
value_  = zip(x_new,  y_new,  z_new)
value2_ = zip(x2_new, y2_new, z2_new)
value_  = sorted(value_)  #sort by x values
value2_ = sorted(value2_) #sort by x values

for i in range(len(min_new_1_final)):
    x_new[i] = value_[i][0]*100
    y_new[i] = value_[i][1]*100

for i in range(len(min_new_2_final)):
    x2_new[i] = value2_[i][0]*100
    y2_new[i] = value2_[i][1]*100


#data    visualistaion
#150
#x_exp = [0, 2.63584,6.19139,25.73025,50.39820,75.28581,100]
#y_exp = [76.35838,78.90173,78.90173,78.30058,78.34682,78.30058,79.36416]
#18
x_exp = [0,2.75658,6.08478,25.61721,50.26718,75.14965,100]
y_exp = [72.33526,75.24856,75.06358,74.231213,73.63006,73.39884,74.41618]
#300
#x_exp = [0,2.46114,6.23250,25.78805,50.23250,75.34361,100]
#y_exp = [76.26590,80.61272,80.38150,80.38150,80.38150,80.38150,80.56647]
plt.plot(x_new,y_new,label='experimental (18 degrees Celsius)')
plt.plot(x_exp,y_exp,label='theoretical data (18 degrees Celsius)')
plt.ylim(66,100)
plt.xlabel('percentage of hydrogen in fuel, %')
plt.ylabel('percentage of fuel in mixture, %')
plt.legend()
plt.show()

#150
#x2_exp = [0,5.83475,25.51280,50.17858,75.079545,100]
#y2_exp = [12.14035,9.824561,6.84211,5.017544,3.75439,3.05263]
#18
x2_exp = [0,1.62800,6.10891,25.33421,50.20121,75.31472,100]
y2_exp = [13.57895,12.456140,11.47368,7.85965,5.64912,4.31578,3.75438]
#300
#x2_exp = [0,25.65994,50.12198,75.25813,100]
#y2_exp = [9.96491,4.94737,3.43860,2.73684,2.21053]
plt.plot(x2_new,y2_new,label='experimental data (18 degrees Celsius)')
plt.plot(x2_exp,y2_exp,label='theoretical data (18 degrees Celsius)')
plt.xlabel('percentage of hydrogen in fuel, %')
plt.ylabel('percentage of fuel in mixture, %')
plt.legend()
plt.show()

#data output
f1 = open("min_new_1_final.txt","w+")
for j in range(len(x_new)):
    f1.write(str(value_[j][0]) + ' ' + str(value_[j][1])+' ' + str(value_[j][2])+'\n')
f1.close()

f2 = open("min_new_2_final.txt","w+")
for j in range(len(x2_new)):
    f2.write(str(value2_[j][0]) + ' ' + str(value2_[j][1])+' ' + str(value_[j][2])+'\n')
f2.close()
