"""
.. versionadded:: 1.1.0
   This demo depends on new features added to contourf3d.
"""

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
a = np.linspace(0.,10.,100)
b = np.linspace(0.,30.,100)
a1 = 10
b1 = 20


fig = plt.figure()
ax = fig.gca(projection='3d')
u = np.linspace(0.,10.,100)
X, Y = np.meshgrid(u,u)
Z = np.sin(3*X) + 3*np.sin(Y) + 50
print 'x: ',X
print 'y: ',Y
print 'z: ',Z
ax.plot_surface(X, Y, Z, rstride=8, cstride=8, alpha=0.3)
cset = ax.contourf(X, Y, Z, zdir='z', offset=-100, cmap=cm.coolwarm)
cset = ax.contourf(X, Y, Z, zdir='x', offset=-40, cmap=cm.coolwarm)
cset = ax.contourf(X, Y, Z, zdir='y', offset=40, cmap=cm.coolwarm)

ax.set_xlabel('X')
ax.set_xlim(-15, 15)
ax.set_ylabel('Y')
ax.set_ylim(-15, 15)
ax.set_zlabel('Z')
ax.set_zlim(0, 100)

plt.show()