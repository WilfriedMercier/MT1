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
snapnumber = input('How many snapshots do you have ? ')
# while (type(snapnumber) != <class 'int'>) :
#     radius = input('Integer please. Type again : ')
snapshots = int(snapnumber)

time = []

#Load main data
for i in range(0,snapshots) :
    zeros = '0'*(6 - len(str(i+1)))
    chaine = zeros + str(i+1)   #Nom du fichier, c'est le i+1ème snap
    rawData = np.loadtxt('snapshot_' + chaine + '.dat', skiprows = 2)

    with open('snapshot_' + chaine + '.dat') as f:
        out = f.readline()
    time.append(float(out.split()[-1])) #Extraction du temps et append() dans la liste time

    if i == 0:
        n = rawData.shape[0]
        radius = np.zeros((n))
        radius = rawData[:,0]
        height = np.zeros((snapshots,n))  #(ntimeStep,n)
        temperature = np.zeros((snapshots,n))
        P = np.zeros((snapshots,n))
        Pgaz = np.zeros((snapshots,n))
        Prad = np.zeros((snapshots,n))
        beta = np.zeros((snapshots,n))
        sigma = np.zeros((snapshots,n))
        cs = np.zeros((snapshots,n))
        nu = np.zeros((snapshots,n))
        v = np.zeros((snapshots,n))
        accretionRate = np.zeros((snapshots,n))
        Qp = np.zeros((snapshots,n))
        Qm = np.zeros((snapshots,n))
        Qadv = np.zeros((snapshots,n))
        Cv = np.zeros((snapshots,n))
        Fz = np.zeros((snapshots,n))
        kff = np.zeros((snapshots,n))
        ke = np.zeros((snapshots,n))
        epsilonff = np.zeros((snapshots,n))
        tauff = np.zeros((snapshots,n))
        omega = np.zeros((snapshots,n))
    #k = snapshots-i-1 #reverse time to plot from bottom to up while increasing time
    k=i
    height[k,:] = rawData[:,1]
    temperature[k,:] = rawData[:,2]
    P[k,:] = rawData[:,3]
    Pgaz[k,:] = rawData[:,4]
    Prad[k,:] = rawData[:,5]
    beta[k,:] = rawData[:,6]
    sigma[k,:] = rawData[:,7]
    cs[k,:] = rawData[:,8]
    nu[k,:] = rawData[:,9]
    v[k,:] = rawData[:,10]
    accretionRate[k,:] = rawData[:,11]
    Qp[k,:] = rawData[:,12]
    Qm[k,:] = rawData[:,13]
    Qadv[k,:] = rawData[:,14]
    Cv[k,:] = rawData[:,15]
    Fz[k,:] = rawData[:,16]
    kff[k,:] = rawData[:,17]
    ke[k,:] = rawData[:,18]
    epsilonff[k,:] = rawData[:,19]
    tauff[k,:] = rawData[:,20]
    omega[k,:] = rawData[:,21]

print (len(time))
print (len(radius))

#Test the conservation of the mass
massConservation = np.zeros((sigma.shape[0]-2, sigma.shape[1]-2))
for j in range(1,snapshots-1) :
    for k in range(1,n-1) :
        massConservation[j-1,k-1] = (1./radius[k])*(radius[k+1]*sigma[j,k+1]*v[j,k+1]-radius[k-1]*sigma[j,k-1]*v[j,k-1])/(np.abs(radius[k+1]-radius[k-1]))
        massConservation[j-1,k-1] = massConservation[j-1,k-1] + (sigma[j+1,k]-sigma[j-1,k])/(np.abs(time[j+1]-time[j-1]))

cineticMConservation = np.zeros((sigma.shape[0]-2, sigma.shape[1]-4))
for j in range(1,snapshots-1) :
    for k in range(2,n-2) :
        tTerm = radius[k]*radius[k]*omega[j,k]*(sigma[j+1,k]-sigma[j-1,k])/(np.abs(time[j+1]-time[j-1]))
        fOrderTerm = (1./radius[k])*((radius[k+1]**3.)*omega[j,k+1]*sigma[j,k+1]*v[j,k+1]-(radius[k-1]**3.)*omega[j,k-1]*sigma[j,k-1]*v[j,k-1])/(np.abs(radius[k+1]-radius[k-1]))
        sOrderTerm = (1./radius[k])*((radius[k+1]**3.)*((omega[j,k+2]-omega[j,k])/(np.abs(radius[k+2]-radius[k])))*sigma[j,k+1]*v[j,k+1]-(radius[k-1]**3.)*((omega[j,k]-omega[j,k-2])/(np.abs(radius[k]-radius[k-2])))*sigma[j,k-1]*v[j,k-1])/(np.abs(radius[k+1]-radius[k-1]))
        cineticMConservation[j-1,k-2] = tTerm+fOrderTerm-sOrderTerm

print(np.max(abs(massConservation)))
print(np.max(abs(cineticMConservation)))

#Create directory for plot
if not os.path.exists('plot'):
    os.makedirs('plot')
os.chdir('plot')

bins = 256
plt.figure()
plt.contourf(radius,time,height,bins)
plt.title('Height as a function of radius and time')
plt.xlabel('radius (m)')
plt.ylabel('time (s)')
cbar = plt.colorbar()
cbar.set_label('height (m)')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.tight_layout()
plt.savefig('height.pdf')

plt.figure()
plt.contourf(radius,time,temperature,bins)
plt.title('temperature as a function of radius and time')
plt.xlabel('radius (m)')
plt.ylabel('time (s)')
cbar = plt.colorbar()
cbar.set_label('temperature')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.tight_layout()
plt.savefig('temperature.pdf')

plt.figure()
plt.contourf(radius,time,P,bins)
plt.title('P as a function of radius and time')
plt.xlabel('radius (m)')
plt.ylabel('time (s)')
cbar = plt.colorbar()
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
cbar.set_label('P')
plt.tight_layout()
plt.savefig('P.pdf')

plt.figure()
plt.contourf(radius,time,Pgaz,bins)
plt.title('Pgaz as a function of radius and time')
plt.xlabel('radius (m)')
plt.ylabel('time (s)')
cbar = plt.colorbar()
cbar.set_label('Pgaz')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.tight_layout()
plt.savefig('Pgaz.pdf')

plt.figure()
plt.contourf(radius,time,Prad,bins)
plt.title('Prad as a function of radius and time')
plt.xlabel('radius (m)')
plt.ylabel('time (s)')
cbar = plt.colorbar()
cbar.set_label('Prad')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.tight_layout()
plt.savefig('Prad.pdf')

plt.figure()
plt.contourf(radius,time,beta,bins)
plt.title('beta as a function of radius and time')
plt.xlabel('radius (m)')
plt.ylabel('time (s)')
cbar = plt.colorbar()
cbar.set_label('beta')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.tight_layout()
plt.savefig('beta.pdf')

plt.figure()
plt.contourf(radius,time,sigma,bins)
plt.title('sigma as a function of radius and time')
plt.xlabel('radius (m)')
plt.ylabel('time (s)')
cbar = plt.colorbar()
cbar.set_label('sigma')
plt.tight_layout()
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.savefig('sigma.pdf')

plt.figure()
plt.contourf(radius,time,cs,bins)
plt.title('cs as a function of radius and time')
plt.xlabel('radius (m)')
plt.ylabel('time (s)')
cbar = plt.colorbar()
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
cbar.set_label('cs')
plt.tight_layout()
plt.savefig('cs.pdf')

plt.figure()
plt.contourf(radius,time,nu,bins)
plt.title('nu as a function of radius and time')
plt.xlabel('radius (m)')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ylabel('time (s)')
cbar = plt.colorbar()
cbar.set_label('nu')
plt.tight_layout()
plt.savefig('nu.pdf')

plt.figure()
plt.contourf(radius,time,v,bins)
plt.title('v as a function of radius and time')
plt.xlabel('radius (m)')
plt.ylabel('time (s)')
cbar = plt.colorbar()
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
cbar.set_label('v')
plt.tight_layout()
plt.savefig('v.pdf')

plt.figure()
plt.contourf(radius,time,accretionRate,bins)
plt.title('accretionRate as a function of radius and time')
plt.xlabel('radius (m)')
plt.ylabel('time (s)')
cbar = plt.colorbar()
cbar.set_label('Taux d\'accrétion $(kg.s^{-1})$')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.tight_layout()
plt.savefig('accretionRate.pdf')

plt.figure()
plt.contourf(radius,time,Qp,bins)
plt.title('Qp as a function of radius and time')
plt.xlabel('radius (m)')
plt.ylabel('time (s)')
cbar = plt.colorbar()
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
cbar.set_label('Qp')
plt.tight_layout()
plt.savefig('Qp.pdf')

plt.figure()
plt.contourf(radius,time,Qm,bins)
plt.title('Qm as a function of radius and time')
plt.xlabel('radius (m)')
plt.ylabel('time (s)')
cbar = plt.colorbar()
cbar.set_label('Qm')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.tight_layout()
plt.savefig('Qm.pdf')

plt.figure()
plt.contourf(radius,time,Qp-Qm,bins)
plt.title('Qp-Qm as a function of radius and time')
plt.xlabel('radius (m)')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ylabel('time (s)')
cbar = plt.colorbar()
cbar.set_label('Qp-Qm')
plt.tight_layout()
plt.savefig('Qp-Qm.pdf')

plt.figure()
plt.contourf(radius,time,np.log10(np.abs(Qp-Qm)),bins)
plt.title('Log(Qp-Qm) as a function of radius and time')
plt.xlabel('radius (m)')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ylabel('time (s)')
cbar = plt.colorbar()
cbar.set_label('log(Qp-Qm)')
plt.tight_layout()
plt.savefig('log(Qp-Qm).pdf')

plt.figure()
plt.contourf(radius,time,Qadv,bins)
plt.title('Qadv as a function of radius and time')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.xlabel('radius (m)')
plt.ylabel('time (s)')
cbar = plt.colorbar()
cbar.set_label('Qadv')
plt.tight_layout()
plt.savefig('Qadv.pdf')

plt.figure()
plt.contourf(radius,time,Cv,bins)
plt.title('Cv as a function of radius and time')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.xlabel('radius (m)')
plt.ylabel('time (s)')
cbar = plt.colorbar()
cbar.set_label('Cv')
plt.tight_layout()
plt.savefig('Cv.pdf')

plt.figure()
plt.contourf(radius,time,Fz,bins)
plt.title('Fz as a function of radius and time')
plt.xlabel('radius (m)')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ylabel('time (s)')
cbar = plt.colorbar()
cbar.set_label('Fz')
plt.tight_layout()
plt.savefig('Fz.pdf')

plt.figure()
plt.contourf(radius,time,kff,bins)
plt.title('kff as a function of radius and time')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.xlabel('radius (m)')
plt.ylabel('time (s)')
cbar = plt.colorbar()
cbar.set_label('kff')
plt.tight_layout()
plt.savefig('kff.pdf')

plt.figure()
plt.contourf(radius,time,ke,bins)
plt.title('ke as a function of radius and time')
plt.xlabel('radius (m)')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ylabel('time (s)')
cbar = plt.colorbar()
cbar.set_label('ke')
plt.tight_layout()
plt.savefig('ke.pdf')

plt.figure()
plt.contourf(radius,time,epsilonff,bins)
plt.title('epsilonff as a function of radius and time')
plt.xlabel('radius (m)')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ylabel('time (s)')
cbar = plt.colorbar()
cbar.set_label('epsilonff')
plt.tight_layout()
plt.savefig('epsilonff.pdf')

plt.figure()
plt.contourf(radius,time,tauff,bins)
plt.title('tauff as a function of radius and time')
plt.xlabel('radius (m)')
plt.ylabel('time (s)')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
cbar = plt.colorbar()
cbar.set_label('tauff')
plt.tight_layout()
plt.savefig('tauff.pdf')



# plt.figure()
# plt.contourf(radius,time,height)
