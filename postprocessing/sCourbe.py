#coding: utf-8
#!/usr/bin/python3

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

#Load simulation name and path
if (len(sys.argv) > 3) :
    print('Too much argument')
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
print(os.getcwd())
os.chdir(path + '/' + simulationName)

radius = int(input('Which sCourbe would you want to show ? Please type an integer : '))
#while (type(radius) != <class 'int'>) :
#    radius = input('Integer please : ')

i,j,k,x,T,sigma,Qp,Qm = np.loadtxt('ScurveMesh.dat', unpack=True)

i = i.astype(int)
j = j.astype(int)
k = k.astype(int)

S = []
for p in range(len(i)):  #How to know which value corresponds to the radius
    if (i[p] == radius):
        S.append(p)
S.sort()
imin = S[0]
imax = S[-1]

mshape = (int(np.max(j[S])),int(np.max(k[S])))


#Focus on the range of values between imin and imax that corresponds to this input radius
#x = x[imin:(imax+1)]
T = np.reshape(T[imin:(imax+1)], mshape)
sigma = np.reshape(sigma[imin:(imax+1)], mshape)
Qp = np.reshape(Qp[imin:(imax+1)], mshape)
Qm = np.reshape(Qm[imin:(imax+1)], mshape)


deltaQ = Qp-Qm

if not os.path.exists('plot'):
    os.makedirs('plot')
os.chdir('plot')
if not os.path.exists('sCourbe'):
    os.makedirs('sCourbe')
os.chdir('sCourbe')


plt.figure()
plt.contourf(np.log10(sigma),np.log10(T),np.log10(np.abs(deltaQ)))
plt.title('Courbe en S au rayon ' + str(radius))
plt.xlabel('log($\Sigma$)')
plt.ylabel('log($T$)')
cbar = plt.colorbar()
cbar.set_label('$Q_+ - Q_-$')
plt.savefig('S-'+str(radius)+'.pdf')





#A venir : tracer toutes les courbes en S à la fois sans demander quel radius
