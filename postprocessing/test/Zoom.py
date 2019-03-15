import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as colors

f = plt.figure()
ax = f.add_subplot(111)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.tick_params(which='both', direction='in')


#Data from dichotomy
X, Y, FUNC, TAU  = np.genfromtxt("../../test/dichospaceZoom", unpack=True, skip_header=1)

sX = np.size(np.unique(X))
sY = np.size(np.unique(Y))

gridX = np.log10(X.reshape(sX, sY))
gridY = np.log10(Y.reshape(sX, sY))

grid = FUNC.reshape(sX, sY)
gridabs = np.abs(grid)

gridTau = TAU.reshape(sX, sY)
gridTau = np.abs(gridTau)

X0 = X[TAU==1]
Y0 = Y[TAU==1]

#Plot Q+-Q-
p = ax.pcolormesh(gridX, gridY, gridabs, cmap=cm.coolwarm, norm=colors.LogNorm(10**14, gridabs.max()))
#p = ax.pcolormesh(gridX, gridY, gridabs, cmap=cm.coolwarm, vmin=10**14, vmax=gridabs.max())
cbar = plt.colorbar(p, label=r'$Q^+ - Q^-$', ax=ax)
cbar.set_label(r'$\mid \delta Q \mid (\rm{W/kg})$', size=16)

#plt.plot(np.log10(X0), np.log10(Y0), label=r"$\tau_{\rm{eff}} = 1$")

ax.set_xlabel(r"$log_{10}(\rm{S/m^2})$", fontsize=16)
ax.set_ylabel(r"$log_{10}(\rm{T/K})$", fontsize=16)

#ax.legend(loc=2, prop={'size': 12})

#plt.show()

fig = plt.figure()
ax = fig.gca(projection='3d')

ax.zaxis.set_rotate_label(False)
ax.set_zlabel('label text', rotation=90)

p = ax.plot_surface(gridX, gridY, grid, cmap=cm.coolwarm, antialiased=False)
ax.set_xlabel(r"$log_{10}(\rm{S/m^2})$", fontsize=16)
ax.set_ylabel(r"$log_{10}(\rm{T/K})$", fontsize=16)
ax.set_zlabel(r"$\delta Q (\rm{W/kg})$", fontsize=16)
#cbar = plt.colorbar(p, ax=ax)
#cbar.set_label(r'$\mid \delta Q \mid (\rm{W/kg})$', size=16)

plt.show()

