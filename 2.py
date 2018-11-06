import matplotlib.pyplot as plt
import math
from scipy import optimize as opt
import cantera as ct

alpha = [1.,0.35,0.43,14.3]


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
def f1(value,ER):
    k_inf = 4.65e12 * value ** 0.44
    k_0 = 5.75e19 * value ** (-1.4)
    fin = list()
    sum = 0
    fin.append(f_1_ER_H2(ER))
    fin.append(f_1_ER_O2(ER))
    fin.append(f_1_ER_N2(ER))
    fin.append(f_1_ER_H2O(ER))
    for i in range(len(fin)):
        sum += alpha[i]*fin[i]
    n = 2.4e19*300/(6.02e23 * value)*sum
    return  (k_inf * k_0 * n)/(k_inf + k_0*n)
def f2(value):
    return 3.52 * math.pow(10,16) * math.pow(value,-0.7) * math.exp(-8590/value)

ER = [i/100. for i in range(20,30)]
#HP_solution
gas = ct.Solution('gri30.cti')
H2conc = [2*value for value in ER]
T_HP = list()
for i in range(len(H2conc)):
    gas.TPX = 300,1*ct.one_atm,'H2:%f, O2:1, N2:3.76'%(H2conc[i])
    gas.equilibrate('HP')
    T_HP.append(float(gas.T))

row1_1 = [i for i in range(900,1200)]
ans_row1_1 = [f1(val,ER[0]) for val in row1_1]
'''
row1_2 = [i for i in range(900,1200)]
ans_row1_2 = [f1(val,ER[1]) for val in row1_2]

row1_3 = [i for i in range(900,1200)]
ans_row1_3 = [f1(val,ER[2]) for val in row1_3]
'''
#new_row1 = [1/val for val in row1]
#new_ans_row1 = [math.log(val) for val in ans_row1]

row2 = [i for i in range(900,1200)]
ans_row2 = [f2(val) for val in row2]

#new_row2 = [1/val for val in row2]
#new_ans_row2 = [math.log(val) for val in ans_row2]

#for i in range(len(ans_row1)):
    #print(ans_row1[i],ans_row2[i],row1[i])
#fig = plt.figure()
#ax = fig.add_subplot(1, 1, 1)
#ax.plot(row1, ans_row1, 'k1', row1, ans_row2, 'k2')
min_val = []
for i in range(len(ER)):
    min = opt.minimize(lambda x: obj_func(f1(x,ER[i]),f2(x)),900.,method='Nelder-Mead')
    min_val.append(min.x[0])
#plt.subplot(311)
print(min_val[0])
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
min_1 = opt.minimize(lambda x: obj_func(f1(x,ER[0]),f2(x)),900.,method='Nelder-Mead')
min_2 = opt.minimize(lambda x: obj_func(f1(x,ER[1]),f2(x)),900.,method='Nelder-Mead')
min_3 = opt.minimize(lambda x: obj_func(f1(x,ER[2]),f2(x)),900.,method='Nelder-Mead')
print(min_1.x[0])
print(min_2.x[0])
print(min_3.x[0])
'''
'''
plt.subplot(311)
plt.plot(row1_1, ans_row1_1, label='k1')
plt.plot(row2, ans_row2, label='k2')
plt.xlabel('T')
plt.ylabel('k1 and k2, ER=0.5')
plt.subplot(312)
plt.plot(row1_2, ans_row1_2, label='k1')
plt.plot(row2, ans_row2, label='k2')
plt.xlabel('T')
plt.ylabel('k1 and k2, ER=1.0')
plt.subplot(313)
plt.plot(row1_3, ans_row1_3, label='k1')
plt.plot(row2, ans_row2, label='k2')
plt.xlabel('T')
plt.ylabel('k1 and k2, ER=1.5')
#plt.yscale('log')
plt.title("Crossover Temperature")
#fig.tight_layout()
plt.legend()
plt.show()
'''
