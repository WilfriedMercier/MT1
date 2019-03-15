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

fig = plt.figure()
ax = fig.gca(projection='3d')

ax.zaxis.set_rotate_label(False)
ax.set_zlabel('label text', rotation=90)

ax.set_xlabel(r"$log_{10}(\rm{S_{\rm{max}}/m^2})$", fontsize=16)
ax.set_ylabel(r"$log_{10}(\rm{T_{\rm{max}}/K})$", fontsize=16)
ax.set_zlabel(r"$\rm{Radius index}$", fontsize=16)

N, Xdic, Ydic = np.genfromtxt("../../test/MaxSlopePoint", unpack=True, usecols=(0,2,4))
print(Xdic)
ax.plot(np.log10(Xdic), np.log10(Ydic), N, "x", color='black', label="Maximum slope point on the lower branch")
ax.legend(loc=2, prop={'size': 12})

plt.show()

