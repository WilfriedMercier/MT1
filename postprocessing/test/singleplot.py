import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as colors

f = plt.figure()
ax = f.add_subplot(111)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.tick_params(which='both', direction='in')

#x = np.arange(-101, 0, 0.2)
#x1 = np.arange(0, 101, 0.2)
#y = np.arange(-101, 101, 0.2)
#x2 = np.arange(-200, 200, 0.5)
#y2 = x2[:]

#
#array = np.meshgrid(x, y)
#array1 = np.meshgrid(x1, y)
#array2 = np.meshgrid(x2, y2)
#
#func = array[:][0]**2 - array[:][1]**2
#funcfun = np.array(array1[:][0])*0
#
#func2 = array[:][0] - array[:][1]
#func3 = np.array((array2[:][0])**2 + (array2[:][1])**2) - 30000
#
#print(func3[func3==0])
#
#
#plt.xlabel("X")
#plt.ylabel("Y")

#Full analytical function
#c = plt.contour(array[:][0], array[:][1], func)
#p = plt.pcolormesh(array[:][0], array[:][1], func, cmap=cm.coolwarm)
#p = plt.pcolormesh(array2[:][0], array2[:][1], func3, cmap=cm.coolwarm, )
#plt.colorbar(p, label="Analytical values of F")
#plt.clim(-20000,20000)

#p1 = plt.pcolormesh(array1[:][0], array1[:][1], funcfun, cmap=cm.coolwarm)

#Data from dichotomy
X, Y, FUNC, TAU  = np.genfromtxt("../../test/dichospace", unpack=True, skip_header=1)

sX = np.size(np.unique(X))
sY = np.size(np.unique(Y))

gridX = np.log10(X.reshape(sX, sY))
gridY = np.log10(Y.reshape(sX, sY))

grid = FUNC.reshape(sX, sY)
grid = np.abs(grid)

gridTau = TAU.reshape(sX, sY)
gridTau = np.abs(gridTau)

X0 = X[TAU==1]
Y0 = Y[TAU==1]

#plot end of lower branch
xho = [np.log10(X.min()), np.log10(X.max())]
yho = [7.053, 7.053]
xvert = [3.83, 3.83]
yvert = [np.log10(Y.min()), np.log10(Y.max())]

ax.plot(xho, yho, "--", color="black")
ax.plot(xvert, yvert, "--", color="black")

plt.text(4.5, 7.25, r"$T_{\rm{end}}$", size=14)
plt.text(3.5, 7.6, r"$S_{\rm{end}}$", size=14)


#Plot Q+-Q-
p = ax.pcolormesh(gridX, gridY, grid, cmap=cm.coolwarm, norm=colors.LogNorm(10**13, 10**26))
cbar = plt.colorbar(p, label=r'$Q^+ - Q^-$', ax=ax)
cbar.set_label(r'$\mid \delta Q \mid (\rm{W/kg})$', size=16)

#p = ax.pcolormesh(gridX, gridY, gridTau, cmap=cm.coolwarm, norm=colors.LogNorm(vmin=gridTau.min(), vmax=gridTau.max()))
#cbar = plt.colorbar(p, ax=ax)
#cbar.set_label(r'$\tau_{\rm{eff}}$', size=16)

#plt.plot(np.log10(X0), np.log10(Y0), label=r"$\tau_{\rm{eff}} = 1$")

ax.set_xlabel(r"$log_{10}(\rm{S/m^2})$", fontsize=16)
ax.set_ylabel(r"$log_{10}(\rm{T/K})$", fontsize=16)

Xdic, Ydic = np.genfromtxt("../../test/dichodataGoodBranch", unpack=True, usecols=(0,1))
print(Xdic)
ax.plot(np.log10(Xdic), np.log10(Ydic), "x", color='black', label="Zeros found via dichotomy")
ax.legend(loc=2, prop={'size': 12})

plt.show()

