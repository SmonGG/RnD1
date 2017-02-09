import numpy as np
import scipy.stats as spy
from astropy import units
from astropy import constants
import matplotlib.pyplot as plt

mass = 10. * units.GeV/constants.c**2
maxw = spy.maxwell.rvs(size=1000)* units.km/units.s
vel = np.linspace(0.,800.,1000)*units.km/units.s
temp = 1 * units.Kelvin#((0.5 * mass * (vel**2))/constants.k_B)


def veldist(m,v,T):
	k = constants.k_B
	F_v_ = (np.sqrt((m/(2.*np.pi*k*T))**3.)) * (4.*np.pi*(v**2.)) * (np.exp(-((m*(v**2.))/(2.*k*T))))
	return F_v_

print 'Mass: ', mass
print 'Velocity: ', min(vel), 'to ', max(vel)
#print 'Temperature: ', min(temp.cgs), 'to ', max(temp.cgs)
Fv = veldist(mass, vel, temp)
print 'F(v): ', min(Fv), 'to ', max(Fv)

plt.plot(maxw,Fv)
#plt.hist(Fv, bins=100., alpha=0.3)
plt.show()