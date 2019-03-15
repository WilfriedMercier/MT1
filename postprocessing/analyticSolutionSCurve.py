# coding: utf-8
# !/usr/bin/python3

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import f90nml

# Physical constant
G = 6.67408e-11
c = 299792458.
alpha = 0.1
solarMass = 1.98847e30

# Load simulation name and path
if (len(sys.argv) > 2):
    print('Too much argument')
    sys.exit()
elif (len(sys.argv) == 1):
    print('Use default setting')
    paramFile = 'physicalsetting'
else:
    paramFile = sys.argv[1]


# Read the parameter
param = f90nml.read(paramFile)

path = param['physicalsetting']['path']
simulationName = param['physicalsetting']['simulationName']

# Move to simulation file
os.chdir(path + '/' + simulationName)
if not os.path.exists('plot'):
    os.makedirs('plot')
os.chdir('plot')
if not os.path.exists('sCourbe'):
    os.makedirs('sCourbe')
os.chdir('sCourbe')

# Read simulation value
blackHoleMass = param['physicalsetting']['blackHoleMass']
accretionRate0 = param['physicalsetting']['accretionRate0']
rmax = param['physicalsetting']['rmax']
nSPaceStep = param['physicalsetting']['nSPaceStep']

# Define space step
rs = 2.*G*blackHoleMass/(c*c)
rmax = rs*rmax
rmin = rs*3.
xmax = np.sqrt(rmax/rs)
xmin = np.sqrt(rmin/rs)
x = np.zeros((nSPaceStep))
for i in range(np.size(x)):
    x[i] = xmin + ((xmax-xmin)/nSPaceStep)*(float(i)+0.5)
r = (x*rs)**2.
accretionRate = np.logspace(np.log10(accretionRate0*1e-12), np.log10(accretionRate0*1e12), num=100)


# generate mesh
r, accretionRate = np.meshgrid(r, accretionRate)

# Calulate T and Sigma
f = 1.-(1./np.sqrt(r/(3.*rs)))
T = 1.4e4*(alpha**-0.2)*((accretionRate/1e13)**0.3)*((blackHoleMass/solarMass)**0.25)*((r/1e8)**-0.75)*(f**0.3)
Sigma = 52.*(alpha**-0.8)*((accretionRate/1e13)**0.7)*((blackHoleMass/solarMass)**0.25)*((r/1e8)**-0.75)*(f**0.7)

for i in range(nSPaceStep):
    plt.figure()
    plt.plot(np.log10(Sigma[i, :]), np.log10(T[i, :]))
    plt.title('Courbe en S au rayon ' + str(i))
    plt.xlabel('log($\Sigma$)')
    plt.ylabel('log($T$)')
    plt.savefig('Analytic_S-'+str(i)+'.pdf')
    plt.close('all')
