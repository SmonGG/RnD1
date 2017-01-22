from scipy.stats import rv_continuous
import scipy.stats as spy
import numpy as np
from astropy import units
from astropy import constants
import matplotlib.pyplot as plt

class dmvdf_gen(rv_continuous):
	"Dark Matter VPDF"
	def pdf(self, m, v, T):
	#	m = m * units.GeV / constants.c**2
	#	v = v * units.km / units.s
	#	T = T * units.Kelvin
		k = constants.k_B 	#Boltzmann Constant
		return (np.sqrt((m/(2.*np.pi*k*T))**3)) * (4.*np.pi*(v**2)) * (np.exp(-(m*(v**2))/(2.*k*T))) / 600000

	def rvs(self, size, m, T):
	#	m = m * units.GeV / constants.c**2
		v = spy.maxwell.rvs(size=size) * units.km/units.s
	#	T = T * units.Kelvin
		k = constants.k_B 	#Boltzmann Constant
		return v#(np.sqrt((m/(2.*np.pi*k*T))**3)) * (4.*np.pi*(v**2)) * (np.exp(-(m*(v**2))/(2.*k*T))) / 600000

dmvdf = dmvdf_gen(name= 'dmvdf')

v_min = 0 * units.km/units.s
v_max = 5 * units.km/units.s
vrange = np.linspace(v_min,v_max,num=1000)
m=2 * units.GeV/constants.c**2
T=300 * units.Kelvin

dist = dmvdf.rvs(10000,m,T)
theo = dmvdf.pdf(m,vrange,T)
print dist
plt.hist(dist,normed=1,alpha=0.4,bins=500)
plt.plot(vrange, theo, lw=2)
plt.xlabel('Velocity (km/s)')
plt.ylabel('F(v)')
plt.grid(True)
plt.title('Velocity Distribution Function')
plt.show()
print "Program Terminated."