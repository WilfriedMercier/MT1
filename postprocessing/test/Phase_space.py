import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as colors

#fig = plt.figure()
#ax = fig.gca(projection='3d')
#ax.yaxis.set_ticks_position('both')
#ax.xaxis.set_ticks_position('both')
#ax.tick_params(which='both', direction='in')

f = plt.figure()
ax = f.add_subplot(111)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.tick_params(which='both', direction='in', labelsize=20)

#Data from dichotomy
name = "SCurveFromFullDicho"
for i in range(35, 36, 1):
   if (i<10):
      fullname = name + "00000" + str(i)
   elif ( (i>9) and (i<100) ):
      fullname = name + "0000" + str(i)
   else:
      fullname = name + "000" + str(i)

   X, Y, FUNC = np.genfromtxt("../../test/dichospace", unpack=True, skip_header=1, usecols=(0,1,2))

   sX = np.size(np.unique(X))
   sY = np.size(np.unique(Y))

   FUNC1 = np.where(FUNC<=0, np.nan, FUNC)
   gridX = np.log10(X.reshape(sX, sY))
   gridY = np.log10(Y.reshape(sX, sY))

   grid = FUNC1.reshape(sX, sY)
   grid = np.abs(grid)

   FUNC2 = np.where(FUNC>0, np.nan, FUNC)

   grid2 = FUNC2.reshape(sX, sY)
#   grid2 = np.abs(grid2)

   print(grid.min(), grid.max(), grid2.min(), grid2.max())

   #Plot Q+-Q-
   p = ax.pcolormesh(gridX, gridY, grid, cmap=cm.coolwarm, norm=colors.LogNorm(np.nanmin(grid), np.nanmax(grid)))
   cbaxes = f.add_axes([0.91, 0.52, 0.01, 0.35])
   cbar = plt.colorbar(p, label=r'$Q^+ - Q^-$', cax=cbaxes, shrink=0.5)
   cbar.set_label(r'$ \rm{Positive}\ \ \delta Q (\rm{W/kg})$', size=24)
   plt.tick_params(axis='both', which='major', labelsize=20)

   p2 = ax.pcolormesh(gridX, gridY, grid2, cmap=cm.BrBG, norm=colors.SymLogNorm(linthresh=0.03, vmin=np.nanmin(grid2), vmax=np.nanmax(grid2)))
   cbaxes2 = f.add_axes([0.91, 0.13, 0.01, 0.35])
   cbar2 = plt.colorbar(p2, label=r'$Q^+ - Q^-$', cax=cbaxes2, shrink=0.5)
   cbar2.set_label(r'$ \rm{Negative}\ \ \delta Q (\rm{W/kg})$', size=24)
#   cbaxes2.invert_yaxis()
   plt.tick_params(axis='both', which='major', labelsize=20)

#   X, Y = np.genfromtxt("../../test/" + fullname + ".dat", unpack=True, skip_header=1, usecols=(0,1))
#   p = ax.plot(np.log10(X), np.log10(Y), ".", color="black")

ax.set_xlabel(r"$log_{10}(\rm{S/m^2})$", fontsize=24)
ax.set_ylabel(r"$log_{10}(\rm{T/K})$", fontsize=24)
#ax.set_zlabel(r"$\delta Q (\rm{W/kg})$", fontsize=16)
#ax.zaxis.set_rotate_label(False)
#ax.set_zlabel('Radius index', rotation=90)

plt.show()
