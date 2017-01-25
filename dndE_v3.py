from scipy.stats import rv_continuous
import scipy.stats as spy
import numpy as np
from astropy import units
from astropy import constants
import matplotlib.pyplot as plt

# #DM Parameters
class DMclass(object):
	def __init__(self,mass,v_min,v_max,temp):
		self.mass=mass * units.GeV / constants.c**2
		self.v_min=v_min * units.km / units.s
		self.v_max=v_max * units.km / units.s
		self.temp=temp * units.Kelvin

mass = 10 ##(float(input("Enter the mass of dark matter (GeV): ")))
v_min = 0.
v_max = 5.
temp = 230.

DM = DMclass(mass,v_min,v_max,temp)
# #Velocity limitations
def lim(vdf, vmin=DM.v_min, vmax=DM.v_max):
	vdf = vdf[vdf > vmin]
	vdf = vdf[vdf < vmax]
	return vdf
# #
# #DM VDF generation
class dmvdf_gen(rv_continuous):
	"Dark Matter VPDF"
	def pdf(self, m, v, T):
		k = constants.k_B 	#Boltzmann Constant
		return (np.sqrt((m/(2.*np.pi*k*T))**3)) * (4.*np.pi*(v**2)) * (np.exp(-(m*(v**2))/(2.*k*T))) / 600000

	def rvs(self, size, m, T, mu=0, sigma=0.5):
		v = spy.maxwell.rvs(size=size, loc=mu, scale=sigma) * units.km/units.s
		lim(v)
		k = constants.k_B 	#Boltzmann Constant
		return v, (np.sqrt((m/(2.*np.pi*k*T))**3)) * (4.*np.pi*(v**2)) * (np.exp(-(m*(v**2))/(2.*k*T))) / 600000

dmvdf = dmvdf_gen(name= 'dmvdf')
# #note: rvs function returns a duple. when assigning a variable that uses the rvs command, make sure that you split it: eg. test_v, test_DF = rvs(...)
# #Noise Emulator
def noisy(vdf, mu=0, sigma=0.5):
	noise = np.random.normal(mu, sigma, len(vdf))
	soln = vdf * (1 + 0.2*noise)
	return soln
# #	
num = 100 ##int(input('Enter the number of random samples: '))
vrange = np.linspace(DM.v_min,DM.v_max,num=500)

print 'mass: ', DM.mass
print 'v_min: ', DM.v_min
print 'v_max: ', DM.v_max
print 'temp: ', DM.temp
print 'vrange: ', vrange[0], 'to', vrange[-1]

theoretical = dmvdf.pdf(DM.mass, vrange, DM.temp)
print 'theoretical max: ', max(theoretical)
test_v, test_DF = dmvdf.rvs(num, DM.mass, DM.temp)
print 'velocity min: ', min(test_v)
print 'velocity max: ', max(test_v)
print 'max value: ', max(test_DF)
nlvl= ((max(test_DF) / max(theoretical)) - 1) * 100
print 'noise level: ', nlvl,'%'
test_DF = noisy(test_DF)

plt.plot(vrange, theoretical)
plt.scatter(test_v, test_DF)
plt.show()