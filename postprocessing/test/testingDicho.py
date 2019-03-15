import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as colors

f = plt.figure()
ax = f.add_subplot(221)
bx = f.add_subplot(222)
cx = f.add_subplot(223)
dx = f.add_subplot(224)
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
#sX, sY           = np.genfromtxt("../../test/dichospace", unpack=True, )

sX = np.size(np.unique(X))
sY = np.size(np.unique(Y))

grid = FUNC.reshape(sX, sY)
grid = np.abs(grid)

gridTAU = TAU.reshape(sX, sY)
gridTAU = np.abs(gridTAU)

gridX = np.log10(X.reshape(sX, sY))
gridY = np.log10(Y.reshape(sX, sY))

print(np.size(gridX), np.size(gridY), np.size(grid))

#Plot Q+-Q-
c = ax.contour(gridX, gridY, grid)
p = ax.pcolormesh(gridX, gridY, grid, cmap=cm.coolwarm, norm=colors.LogNorm(vmin=grid.min(), vmax=grid.max()))
plt.colorbar(p, label=r'$Q^+ - Q^-$', ax=ax)
#plt.clim(-1000,1000)

ax.set_xlabel(r"$log_{10}(S/m^2)$")
ax.set_ylabel(r"$log_{10}(T/K)$")

#Plot TauEff
c2 = bx.contour(gridX, gridY, gridTAU)
p2 = bx.pcolormesh(gridX, gridY, gridTAU, cmap=cm.coolwarm,  norm=colors.LogNorm(vmin=gridTAU.min(), vmax=gridTAU.max()))
plt.colorbar(p2, label=r'$\tau_{eff}$', ax=bx)

bx.set_xlabel(r"$log_{10}(S/m^2)$")
bx.set_ylabel(r"$log_{10}(T/K)$")

#Xdic, Ydic, taueff, i = np.genfromtxt("../../test/dichodata", unpack=True, usecols=(0,1,4,5))
Xdic, Ydic, taueff, i = np.genfromtxt("../../test/dichodataGoodBranch", unpack=True, usecols=(0,1,4,5))
print(Xdic)
ax.plot(np.log10(Xdic), np.log10(Ydic), "x", color='black', label="Zeros found via dichotomy")

cx.plot(np.log10(Xdic), taueff, ".", color='black', label=r"$\tau_{eff} (\delta Q = 0)$")
cx.set_xlabel(r"$log_{10}(S/m^2)$")
cx.set_ylabel(r"$\tau_{eff}$")

dx.plot(np.log10(Xdic), i, ".", color='black', label=r"Dimension of the dichotomy (0=X, 1=Y)")
dx.set_xlabel(r"$log_{10}(S/m^2)$")

ax.legend()
bx.legend()
cx.legend()
dx.legend()
plt.show()

#plt.savefig('dichoOnSimpleFunction.pdf')

