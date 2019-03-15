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
radius = input('Which radius number ? ')
radiusnumber = int(radius)
snapnumber = int(input('How many snapshots ? '))
time = []
Qp = np.zeros(snapnumber)
Qm = np.zeros(snapnumber)
deltaQ = np.zeros(snapnumber)
accretionrate = np.zeros(snapnumber)
temperature = np.zeros(snapnumber)
sigma = np.zeros(snapnumber)

for i in range(0,snapnumber):
    zeros = '0'*(6 - len(str(i+1)))
    chaine = zeros + str(i+1)   #Nom du fichier, c'est le i+1ème snap
    rawData = np.loadtxt('snapshot_' + chaine + '.dat', skiprows = 2)

    with open('snapshot_' + chaine + '.dat') as f:
        out = f.readline()
    time.append(float(out.split()[-1])) #Extraction du temps et append() dans la liste time

    Qp[i] = rawData[radiusnumber,12]
    Qm[i] = rawData[radiusnumber,13]
    deltaQ[i] = Qp[i] - Qm[i]
    temperature[i] = rawData[radiusnumber,2]
    accretionrate[i] = rawData[radiusnumber,11]
    sigma[i] = rawData[radiusnumber,7]

plt.figure()
plt.plot(time,deltaQ,'r*-')
plt.title('bilan thermique pour le rayon numéro ' + str(radiusnumber))
plt.xlabel('temps (s)')
plt.ylabel('$Q_+ - Q_-$')
plt.savefig('deltaQ.pdf')

plt.figure()
plt.plot(time,sigma,'r*-')
plt.title('densité au rayon numéro' + str(radiusnumber))
plt.xlabel('temps (s)')
plt.ylabel('$\Sigma$')
plt.savefig('density.pdf')

plt.show()

plt.figure()
plt.plot(time,accretionrate)
plt.title('accretionrate as a function of radius and time')
plt.xlabel('temps (s)')
plt.ylabel('accretionrate')
plt.savefig('accretionr-time.pdf')


plt.figure()
plt.plot(time,temperature,'ro-')
plt.title('temperature as a function of radius and time')
plt.xlabel('temps (s)')
plt.ylabel('temperature')
plt.savefig('temp-time.pdf')


