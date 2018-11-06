import matplotlib.pyplot as plt
import math
from scipy import optimize as opt
import cantera as ct

alpha_new = [1.,0.35,0.43,14.3]


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
    return float(min.x[0])

def f4(ER_,P_=1,alpha_=0):
    gas_new = ct.Solution('gri30.cti')
    gas_new.TPX = 300,P_*ct.one_atm,'H2:%f, O2:1, N2:3.76, H2O:%f'%(2.*ER_,alpha_)
    gas_new.equilibrate('HP')
    return float(gas_new.T)

print('Please, input the pressure')
P_new = input().split(' ')
P_new = float(P_new[0])

alpha = [100]

ER = [i/10. for i in range(0,100)]
#HP_solution
gas = ct.Solution('gri30.cti')
H2conc = [2*value for value in ER]
T_HP = list()
for i in range(len(H2conc)):
    gas.TPX = 300,P_new*ct.one_atm,'H2:%f, O2:1, N2:3.76'%(H2conc[i])
    gas.equilibrate('HP')
    T_HP.append(float(gas.T))


row1_1 = [i for i in range(900,1500)]
ans_row1_1 = [f1(val,ER[10],P_new,alpha_=alpha[0]) for val in row1_1]

row2 = [i for i in range(900,1500)]
ans_row2 = [f2(val) for val in row2]

min_val = []
for i in range(len(ER)):
    min = opt.minimize(lambda x: obj_func(f1(x,ER[i],P_new,alpha_=alpha[0]),f2(x)),900.,method='Nelder-Mead')
    min_val.append(min.x[0])

min_new_1 = list()
min_new_2 = list()
for i in range(len(alpha)):
    min_new_1.append(opt.minimize(lambda y: obj_func(f3(y,P_new,alpha[i]),f4(y,P_new,alpha[i])),0.2100,method='Nelder-Mead'))
    min_new_2.append(opt.minimize(lambda y: obj_func(f3(y,P_new,alpha[i]),f4(y,P_new,alpha[i])),8.0,method='Nelder-Mead'))
print('fuel-lean limit is' , min_new_1[0].x[0],'with temperature:',f3(min_new_1[0].x[0],alpha_=alpha[0]))
print('fuel-rich limit is' , min_new_2[0].x[0],'with temperature:',f3(min_new_2[0].x[0],alpha_=alpha[0]))
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
'''
value_H2O = list()
value_H2 = list()
new_arr = list()
for i in range(len(alpha)):
    new_arr.append([val/1000 for val in range(int(round(min_new_1[i].x[0]*1000)),int(round(1000*min_new_2[i].x[0])))])
for i in range(len(alpha)):
    value_H2O.append([f_0_ER_H2O(value,alpha[i]) for value in new_arr[i]])
    value_H2.append([f_0_ER_H2(value,alpha[i]) for value in new_arr[i]])
for i in range(len(alpha)):
    plt.plot(value_H2O[i], value_H2[i])
plt.legend()
plt.show()
'''