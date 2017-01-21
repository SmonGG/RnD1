import numpy as np
import matplotlib.pyplot as plt 
import scipy.stats as spy
from astropy import constants
from astropy import units

#DM Parameters
class DMclass(object):
	def __init__(self,mass,v_min,v_max,temp):
		self.mass=mass * units.GeV / constants.c**2
		self.v_min=v_min * units.km / units.s
		self.v_max=v_max * units.km / units.s
		self.temp=temp * units.Kelvin

mass = 0.1 #(float(input("Enter the mass of dark matter (GeV): ")))
v_min = 0
v_max = 5
temp = 1.

DM = DMclass(mass,v_min,v_max,temp)

num = int(input("Enter the number of random samples: "))

#Velocity distribution function
def vpdf(m,v,T):
	k = constants.k_B #Boltzmann Constant
	##kT = (m*(v**2))/2 	#kT ~ KE 
	ans = (np.sqrt((m/(2.*np.pi*k*T))**3)) * (4.*np.pi*(v**2)) * (np.exp(-(m*(v**2))/(2.*k*T)))
	return ans

#Kinetic Energy distribution function
def epdf(m,v,T):
	KE = (m*(v**2))/2 *units.m**2/units.km**2
	print(KE.unit)
	k = constants.k_B #Boltzmann Constant
	kT = (m*(v**2))/2 	#kT ~ KE 
	ans = (np.sqrt((m/(2.*np.pi*k*T))**3)) * (4.*np.pi*(KE**2)) * (np.exp(-(m*(KE**2))/(2.*k*T)))
	return ans

nbins = 500
#Theoretical plot between v_min and v_max, resolution of 500 points
vrange = np.linspace(DM.v_min,DM.v_max,num=nbins)
theo = vpdf(DM.mass, vrange, DM.temp)
#>KEtheo = epdf(DM.mass, vrange, DM.temp)
#plt.plot(vrange, theo)
#>plt.plot(vrange, KEtheo)

#Noise emulator. DF = distribution function, amt = noise amount based on a percentage of DF.
def noisy(DF, amt):
	NDF = (spy.expon.rvs(size=len(DF))-spy.expon.rvs(size = len(DF))) * 0.05*max(DF)
	DF = DF + NDF
	return DF

mu, sigma = 0., 0.5
mxw = spy.maxwell.rvs(loc=mu, scale=sigma, size=num)*units.km/units.s

#Maxwell sampling, limited between v_min and v_max
mxw = spy.maxwell.rvs(size=num)*units.km/units.s
mxw = mxw[mxw > DM.v_min]
mxw = mxw[mxw < DM.v_max]

s = np.random.normal(mu, sigma, len(mxw))

mxw_noise = mxw * (1. + s)
soln = noisy(vpdf(DM.mass,mxw,DM.temp).value, 0.5)
#plt.scatter(mxw, soln, alpha=0.5)
plt.plot(vrange,spy.maxwell.pdf(vrange), 'k-', lw=2)

#hist, bin_edges = np.hist(mxw,range=(DM.v_min,DM.v_max),bins=nbins,density=True)
plt.hist(mxw, nbins, normed=1, alpha=0.2)
plt.hist(mxw_noise, nbins, normed=1, alpha=0.2)
#Main
plt.xlabel('Velocity (km/s)')
plt.ylabel('F(v)')
plt.grid(True)
plt.title('Velocity Distribution Function')
plt.show()
print "Program Terminated."