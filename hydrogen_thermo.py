import cantera as ct
import matplotlib.pyplot as plt

ER = [i/10. for i in range(0,40)]
gas = ct.Solution('gri30.cti')
H2conc = [2*value for value in ER]
T_HP = list()
for i in range(len(H2conc)):
    gas.TPX = 300,1*ct.one_atm,'H2:%f, O2:1, N2:3.76'%(H2conc[i])
    gas.equilibrate('HP')
    T_HP.append(float(gas.T))
plt.plot(ER, T_HP, label='k1')
plt.xlabel('ER')
plt.ylabel('T_HP')
plt.legend()
plt.show()
