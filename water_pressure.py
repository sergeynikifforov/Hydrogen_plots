'''
import cantera as ct

mix = ct.Solution('liquidvapor.cti')
#initial_temp,initial_press =[float(val) for val in input().split()]
mix.TPX = 300,ct.one_atm,'H2O:1'

mix.equilibrate('HP')

print(mix.P_sat/ct.one_atm)
'''
import cantera as ct

liquid = ct.Solution('water.cti', 'liquid_water')
solid = ct.Solution('water.cti', 'ice')
gas = ct.Solution('h2o2.cti')

for T in [280, 400]:
    gas.TPX = 300, 101325, 'AR:1.0, H2O:0.0'
    mix = ct.Mixture([(gas, 0.01), (liquid, 0.99), (solid, 0.0)])
    mix.T = T
    mix.P = 101325
    mix.equilibrate('TV', solver='gibbs')
    print(gas['H2O'].X*mix.P)
    print('T = {}'.format(T))
    #mix()
