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
os.chdir(path + '/' + simulationName)

print('Read file')

i,j,k,x,T,sigma,Qp,Qm,TauEff = np.loadtxt('ScurveMesh.dat', unpack=True)

i = i.astype(int)
j = j.astype(int)
k = k.astype(int)

time_step = int(input('How many snapshots do you have ? '))
space_step= int(input('How many space steps were used ? '))

forme=(space_step,time_step)

sigm=np.zeros(forme)
temperature=np.zeros(forme)

for q in range(0,time_step-1):
    zeros = '0'*(6 - len(str(q+1)))
    chaine = zeros + str(q+1)   #Nom du fichier, c'est le i+1ème snap
    rawData = np.loadtxt('snapshot_' + chaine + '.dat', skiprows = 2)

    temperature[:,q] = rawData[:,2]
    sigm[:,q] = rawData[:,7]

print(temperature)
print(sigm)


if not os.path.exists('plot'):
    os.makedirs('plot')
os.chdir('plot')
if not os.path.exists('sCourbe'):
    os.makedirs('sCourbe')
os.chdir('sCourbe')

plt.rcParams.update({'font.size': 14})

for l in np.unique(i):
    if ((l== 1) or (l%10 == 0)):
        print('Plot',l)
        S = []
        for p in range(len(i)):  #How to know which value corresponds to the radius
            if (i[p] == l):
                S.append(p)

        mshape = (int(np.max(j[S])),int(np.max(k[S])))


        #Focus on the range of values between imin and imax that corresponds to this input radius
        #x = x[imin:(imax+1)]
        Tl = np.zeros(mshape)
        Taueffl = np.zeros(mshape)
        sigmal = np.zeros(mshape)
        Qpl = np.zeros(mshape)
        Qml = np.zeros(mshape)

        Tl[j[S]-1,k[S]-1] = T[S]
        Taueffl[j[S]-1,k[S]-1] = TauEff[S]
        sigmal[j[S]-1,k[S]-1] = sigma[S]
        Qpl[j[S]-1,k[S]-1] = Qp[S]
        Qml[j[S]-1,k[S]-1] = Qm[S]

        deltaQ = Qpl-Qml


        plt.figure()
        plt.contourf(sigmal,Tl,np.log10(np.abs(deltaQ)),256)
        plt.plot(sigm[l-1,0],temperature[l-1,0],'bo')
        for g in range(1, time_step-2):
            plt.plot(sigm[l-1,g],temperature[l-1,g],'r+')
        plt.plot(sigm[l-1,-1],temperature[l-1,0],'ko')
        plt.title('Courbe en S au rayon ' + str(l))
        plt.xlabel('$\Sigma\ [kg.m^{-2}]$')
        plt.ylabel('T $[K]$')
        plt.yscale('log')
        plt.xscale('log')
        cbar = plt.colorbar()
        cbar.set_label('log($|Q_+ - Q_-|$) $[W.kg^{-1}]$')
        plt.tight_layout()
        plt.savefig('SLog-'+str(l)+'.pdf')

        plt.figure()
        plt.contourf(sigmal,Tl,np.log10(Taueffl),256)
        plt.plot(sigm[l-1,0],temperature[l-1,0],'bo')
        for g in range(1, time_step-2):
        	plt.plot(sigm[l-1,g],temperature[l-1,g],'r+')
        plt.plot(sigm[l-1,-1],temperature[l-1,0],'ko')
        plt.title('Courbe en S au rayon ' + str(l))
        plt.xlabel('$\Sigma\ [kg.m^{-2}]$')
        plt.ylabel('T $[K]$')
        plt.yscale('log')
        plt.xscale('log')
        cbar = plt.colorbar()
        cbar.set_label('log($|Q_+ - Q_-|$) $[W.kg^{-1}]$')
        plt.tight_layout()
        plt.savefig('SLogFrame-'+str(l)+'.pdf')

        plt.figure()
        plt.contourf(sigmal,Tl,deltaQ,256)
        plt.plot(sigm[l-1,0],temperature[l-1,0],'bo')
        for g in range(1, time_step-2):
        	plt.plot(sigm[l-1,g],temperature[l-1,g],'r+')
        plt.plot(sigm[l-1,-1],temperature[l-1,0],'ko')
        plt.title('Courbe en S au rayon ' + str(l))
        plt.xlabel('$\Sigma\ [kg.m^{-2}]$')
        plt.ylabel('T $[K]$')
        plt.yscale('log')
        plt.xscale('log')
        cbar = plt.colorbar()
        cbar.set_label('$Q_+ - Q_-$ $[W.kg^{-1}]$')
        plt.tight_layout()
        plt.savefig('S-'+str(l)+'.pdf')
        plt.close('all')
