import numpy as np
#Сelsius
T_cr = 374+273.15
r = 9685.7142

# P: 0.1 MpA .. 22 MpA
def T_boil(P):
    return (179.47*P**(0.2391))+273.15

def K_val(r,T):
     return r/(4.567*T_cr*(1-T/T_cr)**0.38)
     #return (1-T/T_cr)**0.38
#A const in Miller's correlation
def A_value(k,T_b):
    res1 = 0.607*k*(4*T_cr/T_b - (T_b/T_cr)**2)
    res2 = 1.448*k*(T_cr/T_b-T_b/T_cr) - 2.88081
    return res1-res2

#B const in Miller's correlation
def B_value(k):
    return 0.98 * k * T_cr

#C1 const in Miller's correlation
def C1_value(k):
    return -1.448 * k / T_cr

#C2 const in Miller's correlation
def C2_value(k):
    return 0.607 * k / (T_cr**2)
#P_saturated calculation
def P_sat(T,A,B,C1,C2):
    res = A - B/T + C1*T + C2*(T**2)
    print('A,B/T,C1*T,C2*(T**2)')
    print(A,B/T,C1*T,C2*(T**2))
    return 10**(res)


print('please enter temperature in Kelvins')
T = float(input())
T = T
print('please enter pressure in MpA')
P_full = float(input())
k_ = K_val(r,T)
print('k_')
print(k_)
T_boil_ = T_boil(P_full)
print('T_boil_')
print(T_boil_)
A_ = A_value(k_,T_boil_)
B_ =  B_value(k_)
C1_ = C1_value(k_)
C2_ = C2_value(k_)
print('A_,B_,C1_,C2_')
print(A_,B_,C1_,C2_)
P_sat_ = P_sat(T,A_,B_,C1_,C2_)
print('P_sat_')
print(P_sat_)
