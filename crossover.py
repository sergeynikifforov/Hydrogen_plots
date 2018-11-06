import math
def P_saturated_count(T,T_ref,T_boil,r_T):
    k = r_T / (4.567 * T_ref * math.pow((1 - T/T_ref),0.38))
    A = 0.607*k*abs(4*T_ref/T_boil - (T_boil/T_ref)**2) - 1.448*k*abs(T_ref/T_boil - T_boil/T_ref) + 2.88081
    B = 0.980 * k * T_ref
    C_1 =  - 1.448 * k / T_ref
    C_2 = 0.607 * k  / T_ref**2
    P_saturated = math.exp(A - B/T + C_1*T + C_2*T**2)
    return P_saturated
def f_0_ER_H2(ER,alpha = 0):
    return 2*ER/(1+2*ER + alpha + 3.76)
def f_0_ER_O2(ER,alpha = 0):
    return 1/(1+2*ER + alpha + 3.76)
def f_0_ER_H2O(ER,alpha = 0):
    return alpha/(1+2*ER + alpha + 3.76)
def f_0_ER_N2(ER,alpha = 0):
    return 3.76/(1+2*ER + alpha + 3.76)
def f_1_ER_H2(ER,alpha = 0):
    if(ER>1):
        return (2*ER-2)/(2*ER + alpha + 3.76)
    else:
        return 0
def f_1_ER_O2(ER,alpha = 0):
    if(ER>1):
        return 0
    else:
        return (1-ER)/(ER + alpha + 4.76)
def f_1_ER_H2O(ER,alpha = 0):
    if(ER>1):
        return (2+alpha)/(2*ER + alpha + 3.76)
    else:
        return (2*ER+alpha)/(ER + alpha + 4.76)
def f_1_ER_N2(ER,alpha = 0):
    if(ER>1):
        return (3.76)/(2*ER + alpha + 3.76)
    else:
        return (3.76)/(ER + alpha + 4.76)
def alpha_pred(P,ER,T,T_ref_upd,T_boil_upd,r_T_upd):
    P_sat = P_saturated_count(T,T_ref_upd,T_boil_upd,r_T_upd)
    return (2*ER+4.76)*P_sat/(P-P_sat)
