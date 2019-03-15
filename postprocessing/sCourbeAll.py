#coding: utf-8
#!/usr/bin/python3

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mlt

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

iFrame,jFrame,lowSigma,highSigma,lowTemp,highTemp = np.loadtxt('ScurveFrame.dat', unpack=True)

i = i.astype(int)
j = j.astype(int)
k = k.astype(int)
iFrame = iFrame.astype(int)
jFrame = jFrame.astype(int)

if not os.path.exists('plot'):
    os.makedirs('plot')
os.chdir('plot')
if not os.path.exists('sCourbe'):
    os.makedirs('sCourbe')
os.chdir('sCourbe')

for l in np.unique(i):
    print('Plot',l)
    S = []
    SFrame = []
    for p in range(len(i)):  #How to know which value corresponds to the radius
        if (i[p] == l):
            S.append(p)
    for p in range(len(iFrame)):  #How to know which value corresponds to the radius
        if (iFrame[p] == l):
            SFrame.append(p)

    mshape = (int(np.max(j[S])),int(np.max(k[S])))
    mshapeFrame = int(np.max(jFrame[SFrame]))

    #Focus on the range of values between imin and imax that corresponds to this input radius
    #x = x[imin:(imax+1)]
    Tl = np.zeros(mshape)
    Taueffl = np.zeros(mshape)
    sigmal = np.zeros(mshape)
    Qpl = np.zeros(mshape)
    Qml = np.zeros(mshape)
    lowSigmal = np.zeros(mshapeFrame)
    highSigmal = np.zeros(mshapeFrame)
    lowTempl = np.zeros(mshapeFrame)
    highTempl = np.zeros(mshapeFrame)

    Tl[j[S]-1,k[S]-1] = T[S]
    Taueffl[j[S]-1,k[S]-1] = TauEff[S]
    sigmal[j[S]-1,k[S]-1] = sigma[S]
    Qpl[j[S]-1,k[S]-1] = Qp[S]
    Qml[j[S]-1,k[S]-1] = Qm[S]
    lowSigmal[jFrame[SFrame]-1] = lowSigma[SFrame]
    highSigmal[jFrame[SFrame]-1] = highSigma[SFrame]
    lowTempl[jFrame[SFrame]-1] = lowTemp[SFrame]
    highTempl[jFrame[SFrame]-1] = highTemp[SFrame]

    deltaQ = Qpl-Qml

    plt.rcParams.update({'font.size': 14})
    plt.figure()
    plt.contourf(sigmal,Tl, np.log10(np.abs(deltaQ)), 256)
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
    plt.contourf(sigmal, Tl, np.log10(np.abs(deltaQ)), 256)
    plt.plot(lowSigmal, lowTempl, color='black', marker='x', linestyle='none')
    plt.plot(highSigmal, highTempl, color='black', marker='+', linestyle='none')
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
    plt.contourf(sigmal,Tl,np.log10(Taueffl),256)
    plt.title('Courbe en S au rayon ' + str(l))
    plt.xlabel('$\Sigma\ [kg.m^{-2}]$')
    plt.ylabel('T $[K]$')
    plt.yscale('log')
    plt.xscale('log')
    cbar = plt.colorbar()
    cbar.set_label('log($\tau_{eff}$)')
    plt.tight_layout()
    plt.savefig('TauEff-'+str(l)+'.pdf')

    plt.figure()
    plt.contourf(sigmal,Tl,deltaQ,256)
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
