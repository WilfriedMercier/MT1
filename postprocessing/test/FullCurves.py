import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as colors

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.tick_params(which='both', direction='in')

f1 = plt.figure()
ax1 = f1.add_subplot(111)
ax1.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('both')
ax1.tick_params(which='both', direction='in')

#Data from dichotomy
name = "SCurveFromFullDicho"
for i in range(1, 30, 1):
   if (i<10):
      fullname = name + "00000" + str(i)
   elif ( (i>9) and (i<100) ):
      fullname = name + "0000" + str(i)
   else:
      fullname = name + "000" + str(i)

   X, Y = np.genfromtxt("../../test/" + fullname + ".dat", unpack=True, skip_header=1, usecols=(0,1))
   p = ax.plot(np.log10(X), np.log10(Y), X*0+i, linestyle='-')

   ax.set_xlabel(r"$log_{10}(\rm{S/m^2})$", fontsize=16)
   ax.set_ylabel(r"$log_{10}(\rm{T/K})$", fontsize=16)
   ax.zaxis.set_rotate_label(False)
   ax.set_zlabel('Radius index', rotation=90, fontsize=16)

   ax1.set_xlabel(r"$log_{10}(\rm{S/m^2})$", fontsize=16)
   ax1.set_ylabel(r"$log_{10}(\rm{T/K})$", fontsize=16)
   ax1.plot(np.log10(X), np.log10(Y), linestyle=':')

ax1.quiver(3.5, 6.75, 0.4, -0.6, angles='xy', scale_units='xy', scale=1)
ax1.text(3.82, 6.33, r'$x\nearrow$', size=16)

plt.draw()
plt.show()
