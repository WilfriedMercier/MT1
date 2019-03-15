# coding: utf-8
#!/usr/bin/python3

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# Load simulation name and path
if (len(sys.argv) > 3):
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

# Move to simulation file
os.chdir(path + '/' + simulationName)

print('Read file')

i, j, x, T, sigma, Qp, Qm, TauEff = np.loadtxt('SMeshProfileTemperature.dat', unpack=True)

i = i.astype(int)
j = j.astype(int)


if not os.path.exists('plot'):
    os.makedirs('plot')
os.chdir('plot')
if not os.path.exists('sCourbe'):
    os.makedirs('sCourbe')
os.chdir('sCourbe')

l = int(input('Wich radius ?'))

print('Plot', l)
S = []
for p in range(len(i)):  # How to know which value corresponds to the radius
    if (i[p] == l):
        S.append(p)


mshape = int(np.max(j[S]))

# Focus on the range of values between imin and imax that corresponds to this input radius
# x = x[imin:(imax+1)]
Tl = np.zeros(mshape)
Taueffl = np.zeros(mshape)
sigmal = np.zeros(mshape)
Qpl = np.zeros(mshape)
Qml = np.zeros(mshape)

Tl[j[S] - 1,] = T[S]
Taueffl[j[S] - 1,] = TauEff[S]
sigmal[j[S] - 1] = sigma[S]
Qpl[j[S] - 1] = Qp[S]
Qml[j[S] - 1] = Qm[S]


deltaQ = Qpl - Qml

plt.figure()
plt.loglog(Tl, np.abs(deltaQ))
plt.title('Profil température pour une densité de ' + str(sigmal[0]) + ' au rayon ' + str(l))
plt.xlabel('T')
plt.ylabel('Q_+ - Q_-')
plt.savefig('ProfilT-' + str(l) + '.pdf')

plt.figure()
plt.semilogx(Tl, deltaQ)
plt.title('Profil température pour une densité de ' + str(sigmal[0]) + ' au rayon ' + str(l))
plt.xlabel('T')
plt.ylabel('Q_+ - Q_-')
plt.ylim(-1e12,1e12)
plt.savefig('UnlogProfilT-' + str(l) + '.pdf')

plt.figure()
plt.loglog(Tl, np.abs(deltaQ))
plt.title('Profil profondeur optique pour une densité de ' + str(sigmal[0]) + ' au rayon ' + str(l))
plt.xlabel('T')
plt.ylabel('Taueff')
plt.savefig('ProfilTTauEff-' + str(l) + '.pdf')
