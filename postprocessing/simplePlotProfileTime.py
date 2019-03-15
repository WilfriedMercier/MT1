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
radChoice = int(input('Wich radius ?'))-1
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

print (len(time))
print (len(radius))
#Create directory for plot
if not os.path.exists('plot'):
    os.makedirs('plot')
os.chdir('plot')

plt.figure()
plt.plot(time,height[:,radChoice])
plt.title('Height as a function of radius and time')
plt.xlabel('time')
plt.ylabel('height')
plt.savefig('Profileheight.pdf')

plt.figure()
plt.plot(time,temperature[:,radChoice])
plt.title('temperature as a function of radius and time')
plt.xlabel('time')
plt.ylabel('temperature')
plt.savefig('Profiletemperature.pdf')

plt.figure()
plt.plot(time,P[:,radChoice])
plt.title('P as a function of radius and time')
plt.xlabel('time')
plt.ylabel('P')
plt.savefig('ProfileP.pdf')

plt.figure()
plt.plot(time,Pgaz[:,radChoice])
plt.title('Pgaz as a function of radius and time')
plt.xlabel('time')
plt.ylabel('Pgaz')
plt.savefig('ProfilePgaz.pdf')

plt.figure()
plt.plot(time,Prad[:,radChoice])
plt.title('Prad as a function of radius and time')
plt.xlabel('time')
plt.ylabel('Prad')
plt.savefig('ProfilePrad.pdf')

plt.figure()
plt.plot(time,beta[:,radChoice])
plt.title('beta as a function of radius and time')
plt.xlabel('time')
plt.ylabel('beta')
plt.savefig('Profilebeta.pdf')

plt.figure()
plt.plot(time,sigma[:,radChoice])
plt.title('sigma as a function of radius and time')
plt.xlabel('time')
plt.ylabel('sigma')
plt.savefig('Profilesigma.pdf')

plt.figure()
plt.plot(time,cs[:,radChoice])
plt.title('cs as a function of radius and time')
plt.xlabel('time')
plt.ylabel('cs')
plt.savefig('Profilecs.pdf')

plt.figure()
plt.plot(time,nu[:,radChoice])
plt.title('nu as a function of radius and time')
plt.xlabel('time')
plt.ylabel('nu')
plt.savefig('Profilenu.pdf')

plt.figure()
plt.plot(time,v[:,radChoice])
plt.title('v as a function of radius and time')
plt.xlabel('time')
plt.ylabel('v')
plt.savefig('Profilev.pdf')

plt.figure()
plt.plot(time,accretionRate[:,radChoice])
plt.title('accretionRate as a function of radius and time')
plt.xlabel('time')
plt.ylabel('accretionRate')
plt.savefig('ProfileaccretionRate.pdf')

plt.figure()
plt.plot(time,Qp[:,radChoice])
plt.title('Qp as a function of radius and time')
plt.xlabel('time')
plt.ylabel('Qp')
plt.savefig('ProfileQp.pdf')

plt.figure()
plt.plot(time,Qm[:,radChoice])
plt.title('Qm as a function of radius and time')
plt.xlabel('time')
plt.ylabel('Qm')
plt.savefig('ProfileQm.pdf')

plt.figure()
plt.plot(time,Qp[:,radChoice]-Qm[:,radChoice])
plt.title('Qp-Qm as a function of radius and time')
plt.xlabel('time')
plt.ylabel('Qp-Qm')
plt.savefig('ProfileQp-Qm.pdf')

plt.figure()
plt.plot(time,np.log10(np.abs(Qp[:,radChoice]-Qm[:,radChoice])))
plt.title('Log(Qp-Qm) as a function of radius and time')
plt.xlabel('time')
plt.ylabel('log(Qp-Qm)')
plt.savefig('ProfileLog_Qp-Qm.pdf')

plt.figure()
plt.plot(time,Qadv[:,radChoice])
plt.title('Qadv as a function of radius and time')
plt.xlabel('time')
plt.ylabel('Qadv')
plt.savefig('ProfileQadv.pdf')

plt.figure()
plt.plot(time,Cv[:,radChoice])
plt.title('Cv as a function of radius and time')
plt.xlabel('time')
plt.ylabel('Cv')
plt.savefig('ProfileCv.pdf')

plt.figure()
plt.plot(time,Fz[:,radChoice])
plt.title('Fz as a function of radius and time')
plt.xlabel('time')
plt.ylabel('Fz')
plt.savefig('ProfileFz.pdf')

plt.figure()
plt.plot(time,kff[:,radChoice])
plt.title('kff as a function of radius and time')
plt.xlabel('time')
plt.ylabel('kff')
plt.savefig('Profilekff.pdf')

plt.figure()
plt.plot(time,ke[:,radChoice])
plt.title('ke as a function of radius and time')
plt.xlabel('time')
plt.ylabel('ke')
plt.savefig('Profileke.pdf')

plt.figure()
plt.plot(time,epsilonff[:,radChoice])
plt.title('epsilonff as a function of radius and time')
plt.xlabel('time')
plt.ylabel('epsilonff')
plt.savefig('Profileepsilonff.pdf')

plt.figure()
plt.plot(time,tauff[:,radChoice])
plt.title('tauff as a function of radius and time')
plt.xlabel('time')
plt.ylabel('tauff')
plt.savefig('Profiletauff.pdf')



# plt.figure()
# plt.plot(time,height)
