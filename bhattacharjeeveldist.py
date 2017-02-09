import numpy as np
import scipy.stats as spy
from scipy.stats import rv_continuous
from astropy import units
from astropy import constants
import matplotlib.pyplot as plt

def fn_epsilon(x):
    '''For use in the fn_VDF_sol funciton'''
    k = -1.47
    eps1 = (1+x)**k
    eps2 = np.exp(-1.*(x**(1-k)))
    epsT = eps1 * eps2
    return epsT

def normalise(VDF):
    '''Normalise the Distribution'''
    #Normalising large sample sizes can take time.
    norm = [float(i)/sum(VDF) for i in VDF]
    #print 'Normalising factor: ', np.mean(VDF) / np.mean(norm)
    return norm

class maxwbolt_gen(rv_continuous):
	'''Maxwell-Boltzmann Velocity Probability Distribution Function (from Kuhlen et al. 2009)'''
	def _pdf(self, x):
		x = 230. * units.km/units.s
		m = 100. * units.GeV/constants.c**2
		k = constants.k_B   #Boltzmann Constant
		T = ((0.5 * m * (x**2))/k).cgs
		return ((np.sqrt((m/(2.*np.pi*k*T))**3)) * (4.*np.pi*(x**2)) * (np.exp(-1*((m*(x**2))/(2.*k*T))))).cgs
maxboltgen = maxwbolt_gen(name='maxboltgen', a=0)

class bhattachatjee_gen(rv_continuous):
	'''Bhattachatjee Velocity Probability Distribution Function'''
	def _pdf(self, x):
		vel_0 = 339. * units.km / units.s
		vel_max = 516. * units.km / units.s
		beta = (x**2) / (vel_0**2)
		beta_max = (vel_max**2) / (vel_0**2)

		vdf1 = 4. * np.pi * (x**2)
		vdf2 = fn_epsilon(beta.value) - fn_epsilon(beta_max.value)
		vdfT = vdf1 * vdf2
		return vdfT


vdfgen = bhattachatjee_gen(name='vdfgen', a=0)


def fn_VDF_sol(velocity):
    '''Create the Bhattachatjee Velocity Probability Distribution'''
    vel_0 = 339. * units.km / units.s
    vel_max = 516. * units.km / units.s
    beta = (vel**2) / (vel_0**2)
    beta_max = (vel_max**2) / (vel_0**2)

    vdf1 = 4. * np.pi * (velocity**2)
    vdf2 = fn_epsilon(beta) - fn_epsilon(beta_max)
    vdfT = vdf1 * vdf2

def runplots():
	vel = np.linspace(25.,800.,1000.) * units.km / units.s
	maxw = spy.maxwell.pdf(vel,loc=0, scale=162.634559673)
	bcje = vdfgen._pdf(vel)
	maxbolt = maxboltgen.pdf(vel)
	maxw_max = vel[maxw.argmax()]
	bcje_max = vel[bcje.argmax()]
	maxbolt_max = vel[maxbolt.argmax()]

	#mxfit_parameters = spy.maxwell.fit(bcje)
	#print mxfit_parameters
	#mxfit = spy.maxwell.pdf(vel, *mxfit_parameters)

	print '---'
	print 'Running Plots...'
	print 'Processing... 25%'
	plt.plot(vel,normalise(maxw), label='Maxwell-Boltzmann')
	plt.plot([maxw_max.value, maxw_max.value], [0, max(normalise(maxw))], color='black', alpha=0.5)

	print 'Processing... 50%'
	plt.plot(vel,normalise(bcje.value), label='Bhattacharjee 2013')
	plt.plot([bcje_max.value, bcje_max.value], [0, max(normalise(bcje.value))], color='black', alpha=0.5)

	print 'Processing... 75%'
	plt.plot(vel,normalise(maxbolt), label='Kuhlen 2009')
	plt.plot([maxbolt_max.value, maxbolt_max.value], [0, max(normalise(maxbolt))], color='black', alpha=0.5)

	print 'Processing... Done!'
	print '---'
	print 'Maxwell-Boltzmann Maximum: ', maxw_max
	print 'Bhattacharjee 2013 Maximum: ', bcje_max
	print 'Kuhlen et al. 2009 Maximum: ', maxbolt_max

	plt.xlabel('Dark Matter Velocity (km/s)')
	plt.ylabel('F(v) (normalized)')
	plt.grid(True)
	plt.legend(loc='best', frameon=True)
	plt.title('DM Velocity Distribution Function Comparison')
	plt.show()
	return True

def rvsplots():
	vel = np.linspace(25.,800.,1000.) * units.km / units.s
	randomvel = (np.random.random(size=1000.) * 800)*units.km / units.s

	print '---'
	print 'Running random variables...'
	print 'Processing... 33%'
	maxw_rvs = spy.maxwell.rvs(size=1000, scale=162.634559673)
	#plt.hist(maxw_rvs, bins=1000, normed=1)

	print 'Processing... 66%'
	bcje_rvs = normalise(vdfgen.pdf(vel)) * vel
	#plt.scatter(randomvel, bcje_rvs)
	testE = maxw_rvs
	finalE = testE * bcje_rvs
	plt.scatter(vel, finalE)
	print bcje_rvs

	print 'Processing... Done!'
	print 'Maxwell rvs: ', min(maxw_rvs), 'to', max(maxw_rvs)
	print 'Bhattacharjee rvs: ', min(bcje_rvs), 'to', max(bcje_rvs)

	plt.xlabel('Dark Matter Velocity (km/s)')
	plt.ylabel('F(v) (normalized)')
	plt.grid(True)
	plt.legend(loc='best', frameon=True)
	plt.title('DM Velocity Distribution Function Comparison')
	plt.show()
	return True

runplots()
rvsplots()