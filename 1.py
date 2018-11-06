import matplotlib.pyplot as plt
import math

def f1o(ER):
    return 2*ER/(2*ER + 4.76)
def f2o(ER):
    return 1/(2*ER + 4.76)
def f3o(ER):
    return 3.76/(2*ER + 4.76)
def f4o(ER):
    return 0;
def f1k(ER):
    return (2*ER-2)/(2*ER+3.76) if (ER>=1) else 0
def f2k(ER):
    return 0 if (ER>=1) else (1-ER)/(ER+4.76)
def f3k(ER):
    return 3.76/(2*ER+3.76) if (ER>=1) else 3.76/(ER+4.76)
def f4k(ER):
    return 2/(2*ER+3.76) if (ER>=1) else (2*ER)/(ER+4.76)

ER_value = [0.01*i for i in range(0,1000)]

value_h2_old = [f1o(val) for val in ER_value]
value_o2_old = [f2o(val) for val in ER_value]
value_n2_old = [f3o(val) for val in ER_value]
value_h2o_old = [f4o(val) for val in ER_value]

value_h2_new = [f1k(val) for val in ER_value]
value_o2_new = [f2k(val) for val in ER_value]
value_n2_new = [f3k(val) for val in ER_value]
value_h2o_new = [f4k(val) for val in ER_value]

plt.subplot(211)
plt.plot(ER_value, value_h2_old, label='H$_{2}$')
plt.plot(ER_value, value_o2_old, label='O$_{2}$')
plt.plot(ER_value, value_n2_old, label='N$_{2}$')
plt.plot(ER_value, value_h2o_old, label='H$_{2}$O')
plt.xlabel('ER')
plt.ylabel('[C]')
plt.title("Before")
plt.subplot(212)
plt.plot(ER_value, value_h2_new, label='H$_{2}$')
plt.plot(ER_value, value_o2_new, label='O$_{2}$')
plt.plot(ER_value, value_n2_new, label='N$_{2}$')
plt.plot(ER_value, value_h2o_new, label='H$_{2}$O')
plt.xlabel('ER')
plt.ylabel('[C]')
plt.title("After")
plt.legend()
plt.show()
