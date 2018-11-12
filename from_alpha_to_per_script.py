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


print('Please, enter the ER and alpha')
alpha_new_H2O, ER_new_H2O = input().split(' ')
alpha_new_H2O = float(alpha_new_H2O)
ER_new_H2O = float(ER_new_H2O)

res_initial_H2O = f_0_ER_H2O(ER_new_H2O,alpha_new_H2O)
res_final_H2O =  f_1_ER_H2O(ER_new_H2O,alpha_new_H2O)
res_initial_H2 = f_0_ER_H2(ER_new_H2O,alpha_new_H2O)
res_final_H2 =  f_1_ER_H2(ER_new_H2O,alpha_new_H2O)
res_initial_O2 = f_0_ER_O2(ER_new_H2O,alpha_new_H2O)
res_final_O2 =  f_1_ER_O2(ER_new_H2O,alpha_new_H2O)
res_initial_N2 = f_0_ER_N2(ER_new_H2O,alpha_new_H2O)
res_final_N2 =  f_1_ER_N2(ER_new_H2O,alpha_new_H2O)


print(f"The initial percentage of H2O is {(round(res_initial_H2O,4)*100):.4f}")
print(f"The initial percentage of H2 is {(round(res_initial_H2,4)*100):.4f}")
print(f"The initial percentage of O2 is {(round(res_initial_O2,4)*100):.4f}")
print(f"The initial percentage of N2 is {(round(res_initial_N2,4)*100):.4f}")
print(f"The final percentage of H2O is {(round(res_final_H2O,4)*100):.4f}")
print(f"The final percentage of H2 is {(round(res_final_H2,4)*100):.4f}")
print(f"The final percentage of O2 is {(round(res_final_O2,4)*100):.4f}")
print(f"The final percentage of N2 is {(round(res_final_N2,4)*100):.4f}")
