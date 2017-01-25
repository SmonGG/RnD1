from scipy.stats import rv_continuous
import scipy.stats as spy
import numpy as np
from astropy import units
from astropy import constants
import matplotlib.pyplot as plt



#	def pdf(self, m, v, T):
#		k = constants.k_B 	#Boltzmann Constant
#		return (np.sqrt((m/(2.*np.pi*k*T))**3)) * (4.*np.pi*(v**2)) * (np.exp(-(m*(v**2))/(2.*k*T))) / 600000

#	def rvs(self, size, m, T, mu=0, sigma=0.5):
#		v = spy.maxwell.rvs(size=size, loc=mu, scale=sigma) * units.km/units.s
#		lim(v)
#		k = constants.k_B 	#Boltzmann Constant
#		return v, (np.sqrt((m/(2.*np.pi*k*T))**3)) * (4.*np.pi*(v**2)) * (np.exp(-(m*(v**2))/(2.*k*T))) / 600000

# #DM Parameters
class DMclass(object):
	def __init__(self,mass,v,temp):
		self.mass=mass * units.GeV / constants.c**2
		self.v=v * units.km / units.s
		if temp.cgs.unit == units.Kelvin:
			self.temp = temp
		else:
			self.temp=temp * units.Kelvin

class HALOclass(object):
	def __init__(self,v_min,v_max):
		self.v_min=v_min * units.km / units.s
		self.v_max=v_max * units.km / units.s

mass = 10.
v = 1.

v_min = 1.
v_max = 800.
temp = ((0.5 * mass * (v * units.km/units.s)**2)/constants.k_B).cgs *units.g

DM = DMclass(mass,v,temp)
MW = HALOclass(v_min,v_max)

print DM.mass ##(float(input("Enter the mass of dark matter (GeV): ")))
print DM.temp
print DM.v
print MW.v_min
print MW.v_max

# DM halo properties
print "Effective temp for DM particles moving at ",DM.v, ":", temp.to(units.Kelvin)
a = np.sqrt(constants.k_B * temp/mass)
print "Normalisation velocity ",a
scale = a

m = DM.mass
k_B = constants.k_B
T = DM.temp
v = DM.v
norm = np.sqrt((m/(2.*np.pi*k_B*T))**3) * (4.*np.pi*(v**2)) * (np.exp(-(m*(v**2))/(2.*k_B*T)))
print "test ", norm.to(units.s/units.km)

# #Velocity limitations
def lim(vdf, vmin=MW.v_min, vmax=MW.v_max):
	vdf = vdf[vdf > vmin]
	vdf = vdf[vdf < vmax]
	return vdf
# 

class dmdv_newgen(rv_continuous):
	" Units out are inverse velocity, 2.02799e-11 s/cm for v=1 km/s "
	def _pdf(self, v, m, T):

		m=m * units.GeV / constants.c**2
		v=v * units.km / units.s
		T=T * units.Kelvin		
		k_B = constants.k_B
		out = np.sqrt((m/(2.*np.pi*k_B*T))**3) * (4.*np.pi*(v**2)) * (np.exp(-(m*(v**2))/(2.*k_B*T)))
		out = out.cgs
		out = out.to(units.s/units.km)
		return out.value

# Instantiation
my_dmdv = dmdv_newgen(name='dmdv', a=0.)

print "Cumulative Distribution Function "
print "Out to 100 km/s not many particles ", my_dmdv.cdf(100., m=DM.mass, T=DM.temp)
print "Out to 300 km/s about half of all particles ", my_dmdv.cdf(300., m=DM.mass, T=DM.temp)
print "Out to 600 km/s should include most of all particles ", my_dmdv.cdf(600., m=DM.mass, T=DM.temp)

vrange = 10.*(np.linspace(MW.v_min,MW.v_max,num=10000.))

rv = my_dmdv.pdf(vrange, m=DM.mass.value, T=DM.temp.value, loc=0, scale=scale)
plt.plot(vrange,rv,'-',color='b')

v_values = my_dmdv.rvs(m=DM.mass.value, T=DM.temp.value, size=2, loc=0, scale=scale)
plt.hist(v_values,2,color='b',normed=1)
print v_values
plt.show()

#note: rvs function returns a duple. when assigning a variable that uses the rvs command, make sure that you split it: eg. test_v, test_DF = rvs(...)
# #Noise Emulator
def noisy(vdf, mu=0, sigma=0.5):
	noise = np.random.normal(mu, sigma, len(vdf))
	soln = vdf * (1 + 0.2*noise)
	return soln
# #	
#>num = 10000 ##int(input('Enter the number of random samples: '))

print 'mass: ', DM.mass
print 'v_min: ', MW.v_min
print 'v_max: ', MW.v_max
print 'temp: ', DM.temp.cgs
print 'vrange: ', vrange[0], 'to', vrange[-1]
plt.show()

#>theoretical = dmvdf.pdf(DM.mass, vrange, DM.temp)
#>print 'theoretical max: ', max(theoretical.to(units.s/units.km))
#>test_v, test_KE = dmvdf.rvs(num, DM.mass, DM.temp)
#>print 'velocity min: ', min(test_v)
#>print 'velocity max: ', max(test_v)
#>print 'max value: ', max(test_KE)
#>nlvl= ((max(test_KE) / max(theoretical)) - 1) * 100
#>print 'noise level: ', nlvl,'%'
#>test_KE = noisy(test_v)

#>plt.plot(vrange, theoretical)
#plt.scatter(test_v, test_DF)
#>plt.hist(test_v, bins=500, normed=1,alpha=0.5,color='g')
#>plt.hist(test_KE, bins=500, normed=1,alpha=0.5,color='b')
#>plt.show()