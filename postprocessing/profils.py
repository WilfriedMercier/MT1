#coding: utf-8
#!/usr/bin/python3

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

#Load simulation name and path
if (len(sys.argv) > 3) :
    print('Too many arguments')
    sys.exit()
elif (len(sys.argv) == 1):
    print('Use default path and default simulation name')
    path = '.'
    simulationName = 'test'
elif (len(sys.argv) == 2):
    print('Use default path')
    path = '.'
    simulationName = sys.argv[1]
else:
    path = sys.argv[2]
    simulationName = sys.argv[1]


#Move to simulation file
os.chdir(path + '/' + simulationName)
print(os.getcwd())

#Choose snapshot
snapnumber = input('Which snapshot ? ')
snapshot = int(snapnumber)

zeros = '0'*(6 - len(str(snapshot)))
chaine = zeros + str(snapshot)   #Nom du fichier, c'est le i+1ème snap
rawData = np.loadtxt('snapshot_' + chaine + '.dat', skiprows = 2)

with open('snapshot_' + chaine + '.dat') as f:
    out = f.readline()
time = (float(out.split()[-1])) #Extraction du temps et append() dans la liste time

print(snapshot, time)

n = rawData.shape[0]
radius = np.zeros((n))
radius = rawData[:,0]
height = np.zeros(n)
temperature = np.zeros(n)
P = np.zeros(n)
Pgaz = np.zeros(n)
Prad = np.zeros(n)
beta = np.zeros(n)
sigma = np.zeros(n)
cs = np.zeros(n)
nu = np.zeros(n)
v = np.zeros(n)
accretionRate = np.zeros(n)
Qp = np.zeros(n)
Qm = np.zeros(n)
Qadv = np.zeros(n)
Cv = np.zeros(n)
Fz = np.zeros(n)
kff = np.zeros(n)
ke = np.zeros(n)
epsilonff = np.zeros(n)
tauff = np.zeros(n)

height = rawData[:,1]
temperature = rawData[:,2]
P = rawData[:,3]
Pgaz = rawData[:,4]
Prad = rawData[:,5]
beta = rawData[:,6]
sigma = rawData[:,7]
cs = rawData[:,8]
nu = rawData[:,9]
v = rawData[:,10]
accretionRate = rawData[:,11]
Qp = rawData[:,12]
Qm = rawData[:,13]
Qadv = rawData[:,14]
Cv = rawData[:,15]
Fz = rawData[:,16]
kff = rawData[:,17]
ke = rawData[:,18]
epsilonff = rawData[:,19]
tauff = rawData[:,20]

plt.figure()
plt.plot(radius,Qp-Qm)
plt.title('deltaQ at a given time')
plt.xlabel('radius (m)')
plt.ylabel('deltaQ (J/s)')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.savefig('deltaQ'+ str(snapnumber) + '.pdf')

# 
# plt.figure()
# plt.plot(radius,height)
# plt.title('Height as a function of radius and time')
# plt.xlabel('radius')
# plt.ylabel('height')
# 
# plt.figure()
# plt.plot(radius,temperature)
# plt.title('temperature as a function of radius and time')
# plt.xlabel('radius')
# plt.ylabel('temperature')
# 
# plt.figure()
# plt.plot(radius,P)
# plt.title('P as a function of radius and time')
# plt.xlabel('radius')
# plt.ylabel('P')
# 
# plt.figure()
# plt.plot(radius,Pgaz)
# plt.title('Pgaz as a function of radius and time')
# plt.xlabel('radius')
# plt.ylabel('Pgaz')
# 
# plt.figure()
# plt.plot(radius,Prad)
# plt.title('Prad as a function of radius and time')
# plt.xlabel('radius')
# plt.ylabel('Prad')
# 
# plt.figure()
# plt.plot(radius,beta)
# plt.title('beta as a function of radius and time')
# plt.xlabel('radius')
# plt.ylabel('beta')
# 
# plt.figure()
# plt.plot(radius,sigma)
# plt.title('sigma as a function of radius and time')
# plt.xlabel('radius')
# plt.ylabel('sigma')
# 
# plt.figure()
# plt.plot(radius,cs)
# plt.title('cs as a function of radius and time')
# plt.xlabel('radius')
# plt.ylabel('cs')
# 
# plt.figure()
# plt.plot(radius,nu)
# plt.title('nu as a function of radius and time')
# plt.xlabel('radius')
# plt.ylabel('nu')
# 
# plt.figure()
# plt.plot(radius,v)
# plt.title('v as a function of radius and time')
# plt.xlabel('radius')
# plt.ylabel('v')
# 
plt.figure()
plt.plot(radius,accretionRate)
plt.title('accretionRate at a given time')
plt.xlabel('radius (m)')
plt.ylabel('accretionRate (kg/s)')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.savefig('accretionRate'+ str(snapnumber) + '.pdf')
# 
# plt.figure()
# plt.plot(radius,Qp)
# plt.title('Qp as a function of radius and time')
# plt.xlabel('radius')
# plt.ylabel('Qp')
# 
# plt.figure()
# plt.plot(radius,Qm)
# plt.title('Qm as a function of radius and time')
# plt.xlabel('radius')
# plt.ylabel('Qm')
# 
# plt.figure()
# plt.plot(radius,Qp-Qm)
# plt.title('Qp-Qm as a function of radius and time')
# plt.xlabel('radius')
# plt.ylabel('Qp-Qm')
# 
# plt.figure()
# plt.plot(radius,Qadv)
# plt.title('Qadv as a function of radius and time')
# plt.xlabel('radius')
# plt.ylabel('Qadv')
# 
# plt.figure()
# plt.plot(radius,Cv)
# plt.title('Cv as a function of radius and time')
# plt.xlabel('radius')
# plt.ylabel('Cv')
# 
# plt.figure()
# plt.plot(radius,Fz)
# plt.title('Fz as a function of radius and time')
# plt.xlabel('radius')
# plt.ylabel('Fz')
# 
# plt.figure()
# plt.plot(radius,kff)
# plt.title('kff as a function of radius and time')
# plt.xlabel('radius')
# plt.ylabel('kff')
# 
# plt.figure()
# plt.plot(radius,ke)
# plt.title('ke as a function of radius and time')
# plt.xlabel('radius')
# plt.ylabel('ke')
# 
# plt.figure()
# plt.plot(radius,epsilonff)
# plt.title('epsilonff as a function of radius and time')
# plt.xlabel('radius')
# plt.ylabel('epsilonff')
# 
# plt.figure()
# plt.plot(radius,tauff)
# plt.title('tauff as a function of radius and time')
# plt.xlabel('radius')
# plt.ylabel('tauff')
