from scipy.stats import rv_continuous
import scipy.stats as spy
import numpy as np
from astropy import units
from astropy import constants
import matplotlib.pyplot as plt

v_min = 0 * units.km/units.s
v_max = 50 * units.km/units.s
##correction = 6e5
vrange = np.linspace(v_min,v_max,num=1000)
m=1 * units.GeV/constants.c**2
T=230 * units.Kelvin
num = 1000##int(input('Enter the number of random samples: '))

class dmvdf_gen(rv_continuous):
	"Dark Matter VPDF"
	def pdf(self, m, v, T):
	#	m = m * units.GeV / constants.c**2
	#	v = v * units.km / units.s
	#	T = T * units.Kelvin
		k = constants.k_B 	#Boltzmann Constant
		return (np.sqrt((m/(2.*np.pi*k*T))**3)) * (4.*np.pi*(v**2)) * (np.exp(-(m*(v**2))/(2.*k*T)))## / correction

	def rvs(self, size, m, T):
	#	m = m * units.GeV / constants.c**2
		v = spy.maxwell.rvs(size=size) * units.km/units.s
	#	T = T * units.Kelvin
		k = constants.k_B 	#Boltzmann Constant
		return v#(np.sqrt((m/(2.*np.pi*k*T))**3)) * (4.*np.pi*(v**2)) * (np.exp(-(m*(v**2))/(2.*k*T))) / 600000

dmvdf = dmvdf_gen(name= 'dmvdf')


dist = dmvdf.rvs(num,m,T)
theo = dmvdf.pdf(m,vrange,T)
##plt.hist(dist,normed=1,alpha=0.4,bins=num/20)
plt.plot(vrange, theo, lw=2)
plt.plot(vrange,dmvdf.pdf(m,vrange,T*10),lw=2)
plt.plot(vrange,dmvdf.pdf(m,vrange,T*100),lw=2)

plt.xlabel('Velocity (km/s)')
plt.ylabel('F(v)')
plt.grid(True)
plt.title('Velocity Distribution Function')
plt.show()
print "Program Terminated."