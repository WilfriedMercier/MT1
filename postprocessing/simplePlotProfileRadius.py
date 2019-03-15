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

#Load time
snapnumber = input('Which snapshot ? ')
# while (type(snapnumber) != <class 'int'>) :
#     radius = input('Integer please. Type again : ')
snapshot = int(snapnumber)

#Load main data
zeros = '0'*(6 - len(str(snapshot)))
chaine = zeros + str(snapshot)   #Nom du fichier, c'est le i+1ème snap
rawData = np.loadtxt('snapshot_' + chaine + '.dat', skiprows = 2)

with open('snapshot_' + chaine + '.dat') as f:
    out = f.readline()
time = float(out.split()[-1]) #Extraction du temps et append() dans la liste time


n = rawData.shape[0]
radius = np.zeros((n))
radius = rawData[:,0]
height = np.zeros((n))  #(ntimeStep,n)
temperature = np.zeros((n))
P = np.zeros((n))
Pgaz = np.zeros((n))
Prad = np.zeros((n))
beta = np.zeros((n))
sigma = np.zeros((n))
cs = np.zeros((n))
nu = np.zeros((n))
v = np.zeros((n))
accretionRate = np.zeros((n))
Qp = np.zeros((n))
Qm = np.zeros((n))
Qadv = np.zeros((n))
Cv = np.zeros((n))
Fz = np.zeros((n))
kff = np.zeros((n))
ke = np.zeros((n))
epsilonff = np.zeros((n))
tauff = np.zeros((n))

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


print (len(radius))
#Create directory for plot
if not os.path.exists('plot'):
    os.makedirs('plot')
os.chdir('plot')


def build_window():
   f = plt.figure()
   ax = f.add_subplot(111)
   ax.yaxis.set_ticks_position('both')
   ax.xaxis.set_ticks_position('both')
   ax.yaxis.label.set_size(16)
   ax.xaxis.label.set_size(16)
   ax.title.set_size(16)
   ax.tick_params(which='both', direction='in', labelsize=12)
   ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))

build_window()
plt.plot(radius,height)
plt.title('Height as a function of radius and time')
plt.xlabel('radius [m]')
plt.ylabel('height [m]')
plt.tight_layout()
plt.savefig('ProfileRheight.pdf')

build_window()
plt.plot(radius,temperature)
plt.title('temperature as a function of radius and time')
plt.xlabel('radius [m]')
plt.ylabel('temperature [K]')
plt.tight_layout()
plt.savefig('ProfileRtemperature.pdf')

build_window()
plt.figure()
plt.plot(radius,P)
plt.title('P as a function of radius and time')
plt.xlabel('radius')
plt.ylabel('P')
plt.savefig('ProfileRP.pdf')

build_window()
plt.figure()
plt.plot(radius,Pgaz)
plt.title('Pgaz as a function of radius and time')
plt.xlabel('radius')
plt.ylabel('Pgaz')
plt.savefig('ProfileRPgaz.pdf')

build_window()
plt.figure()
plt.plot(radius,Prad)
plt.title('Prad as a function of radius and time')
plt.xlabel('radius')
plt.ylabel('Prad')
plt.savefig('ProfileRPrad.pdf')

build_window()
plt.figure()
plt.plot(radius,beta)
plt.title('beta as a function of radius and time')
plt.xlabel('radius')
plt.ylabel('beta')
plt.savefig('ProfileRbeta.pdf')

build_window()
plt.figure()
plt.plot(radius,sigma)
plt.title('sigma as a function of radius and time')
plt.xlabel('radius [m]')
plt.ylabel('sigma $[kg.m^{-2}]$')
plt.tight_layout()
plt.savefig('ProfileRsigma.pdf')

build_window()
plt.figure()
plt.plot(radius,cs)
plt.title('cs as a function of radius and time')
plt.xlabel('radius')
plt.ylabel('cs')
plt.savefig('ProfileRcs.pdf')

build_window()
plt.figure()
plt.plot(radius,nu)
plt.title('nu as a function of radius and time')
plt.xlabel('radius')
plt.ylabel('nu')
plt.savefig('ProfileRnu.pdf')

build_window()
plt.figure()
plt.plot(radius,v)
plt.title('v as a function of radius and time')
plt.xlabel('radius')
plt.ylabel('v')
plt.savefig('ProfileRv.pdf')

build_window()
plt.figure()
plt.plot(radius,accretionRate)
plt.title('accretionRate as a function of radius and time')
plt.xlabel('radius')
plt.ylabel('accretionRate')
plt.savefig('ProfileRaccretionRate.pdf')

build_window()
plt.figure()
plt.plot(radius,Qp)
plt.title('Qp as a function of radius and time')
plt.xlabel('radius')
plt.ylabel('Qp')
plt.savefig('ProfilReQp.pdf')

build_window()
plt.figure()
plt.plot(radius,Qm)
plt.title('Qm as a function of radius and time')
plt.xlabel('radius')
plt.ylabel('Qm')
plt.savefig('ProfileRQm.pdf')

build_window()
plt.figure()
plt.plot(radius,Qp-Qm)
plt.title('Qp-Qm as a function of radius and time')
plt.xlabel('radius')
plt.ylabel('Qp-Qm')
plt.savefig('ProfileRQp-Qm.pdf')

build_window()
plt.figure()
plt.plot(radius,np.log10(np.abs(Qp-Qm)))
plt.title('Log(Qp-Qm) as a function of radius and time')
plt.xlabel('radius')
plt.ylabel('log(Qp-Qm)')
plt.savefig('ProfileRLog_Qp-Qm.pdf')

build_window()
plt.figure()
plt.plot(radius,Qadv)
plt.title('Qadv as a function of radius and time')
plt.xlabel('radius')
plt.ylabel('Qadv')
plt.savefig('ProfileRQadv.pdf')

build_window()
plt.figure()
plt.plot(radius,Cv)
plt.title('Cv as a function of radius and time')
plt.xlabel('radius')
plt.ylabel('Cv')
plt.savefig('ProfileRCv.pdf')

build_window()
plt.figure()
plt.plot(radius,Fz)
plt.title('Fz as a function of radius and time')
plt.xlabel('radius')
plt.ylabel('Fz')
plt.savefig('ProfileRFz.pdf')

build_window()
plt.figure()
plt.plot(radius,kff)
plt.title('kff as a function of radius and time')
plt.xlabel('radius')
plt.ylabel('kff')
plt.savefig('ProfileRkff.pdf')

build_window()
plt.figure()
plt.plot(radius,ke)
plt.title('ke as a function of radius and time')
plt.xlabel('radius')
plt.ylabel('ke')
plt.savefig('ProfileRke.pdf')

build_window()
plt.figure()
plt.plot(radius,epsilonff)
plt.title('epsilonff as a function of radius and time')
plt.xlabel('radius')
plt.ylabel('epsilonff')
plt.savefig('ProfileRepsilonff.pdf')

build_window()
plt.figure()
plt.plot(radius,tauff)
plt.title('tauff as a function of radius and time')
plt.xlabel('radius')
plt.ylabel('tauff')
plt.savefig('ProfileRtauff.pdf')



# plt.figure()
# plt.plot(radius,height)
