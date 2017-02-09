from scipy.stats import maxwell
import numpy as np
import astropy as ap
import scipy.stats as spy
import astropy.units as units
import astropy.constants as constants
import matplotlib.pyplot as plt
from scipy.stats import rv_continuous

# #Dark Matter Parameters
class E_recoil_DMClass(object):
	def __init__(self, mass, vel, temp):
		self.mass=mass * units.GeV/constants.c**2
		self.vel=vel * units.km/units.s

mass = 10.
vel = 230.
temp = (0.5 * mass * (vel**2)/constants.k_B.value)
DM = E_recoil_DMClass(mass, vel, temp)

DM.temp = ((0.5 * DM.mass * (DM.vel**2)/constants.k_B)).cgs
DM.temp.to(units.Kelvin)
# #
# #DM Halo Parameters
class E_recoil_HALOClass(object):
	def __init__(self, vmin, vmax, vrange, mrange):
		self.vmin=vmin * units.km/units.s
		self.vmax=vmax * units.km/units.s
		self.vrange=vrange
		self.mrange=mrange

vmin = 0.
vmax = 800.
vrange = np.linspace(vmin, vmax, num=100.) * units.km/units.s 
mrange = np.linspace(1., 10000., num=100.) * units.GeV/constants.c**2
MW = E_recoil_HALOClass(vmin, vmax, vrange, mrange)
# #
# #DM VDF generation
class E_recoil_dmvdf_gen(rv_continuous):
	"Dark Matter VPDF"
	def pdf(self, m, v, T):
		k = constants.k_B 	#Boltzmann Constant
		return (np.sqrt((m/(2.*np.pi*k*T))**3)) * (4.*np.pi*(v**2)) * (np.exp(-(m*(v**2))/(2.*k*T))) /2500

	def rvs(self, size, m, T, mu=0, sigma=0.5):
		v = spy.maxwell.rvs(size=size, loc=mu, scale=sigma) * units.km/units.s
		k = constants.k_B 	#Boltzmann Constant
		return v, (np.sqrt((m/(2.*np.pi*k*T))**3)) * (4.*np.pi*(v**2)) * (np.exp(-(m*(v**2))/(2.*k*T))) / 600000

dmvdf = E_recoil_dmvdf_gen(name= 'dmvdf')
# #note: rvs function returns a duple. when assigning a variable that uses the rvs command, make sure that you split it: eg. test_v, test_DF = rvs(...)
# #Recoil Energy
def recoilE(AtomMass, DMmass=DM.mass, DMvel=MW.vrange, angle=90):
	KE = (0.5 * DMmass * constants.c**2 * ((DMvel/constants.c)**2)).to(units.GeV)
	r = (4 * DMmass * AtomMass) / ((DMmass + AtomMass)**2)
	soln = (KE*r*(1-np.cos(angle))) /2
	return soln.to(units.keV)
# #
# #Array Multiplier
#>def arraymul(array1, array2):


# #Atomic Masses
sodium = 22.989769 * (0.932 * units.GeV/constants.c**2)
iodine = 126.90447 * (0.932 * units.GeV/constants.c**2)
# #

print "Vel: ", MW.vmin, 'to ', MW.vmax
print "Resonant DM mass for sodium atom (GeV): ", sodium*constants.c**2
print "Resonant DM mass for iodine atom (GeV): ", iodine*constants.c**2

theo_VPDF = dmvdf.pdf(DM.mass, MW.vrange, DM.temp)

sodium_energies = recoilE(sodium)
iodine_energies = recoilE(iodine)
VPDFxER_sodium = theo_VPDF*sodium_energies
VPDFxER_iodine = theo_VPDF*iodine_energies
ERPDF_NaI = VPDFxER_sodium+VPDFxER_iodine
sodium_M_energies = recoilE(sodium, DMmass = np.linspace(1.,10000.,100000) * units.GeV/constants.c**2, DMvel = DM.vel)
iodine_M_energies = recoilE(iodine, DMmass = np.linspace(1.,10000.,100000) * units.GeV/constants.c**2, DMvel = DM.vel)
#plt.plot(MW.vrange, sodium_energies, label='Sodium')
#plt.plot(MW.vrange, iodine_energies, label='Iodine')
#plt.show()
#plt.plot(np.linspace(1.,10000.,100000) * units.GeV/constants.c**2, sodium_M_energies, label='Sodium M')
#plt.plot(np.linspace(1.,10000.,100000) * units.GeV/constants.c**2, iodine_M_energies, label='Iodine M')

plt.plot(MW.vrange, theo_VPDF, label='Theoretical VPDF')
plt.show()
plt.plot(MW.vrange, VPDFxER_sodium, label='Sodium')
plt.plot(MW.vrange, VPDFxER_iodine, label='Iodine')
plt.plot(MW.vrange, ERPDF_NaI, label='Sodium Iodide')
plt.xlabel('Dark Matter Velocity (km/s)')
plt.ylabel('Nuclear Recoil Energy (keV)')
plt.grid(True)
plt.title('Nuclear Energy Recoil Probability Distribution')
plt.legend(loc='best', frameon=True)
plt.show()
# #3D contour plot
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm

fig = plt.figure()
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(MW.vrange,MW.mrange*(constants.c**2))
Z = (recoilE(sodium, DMvel=X)+recoilE(iodine, DMvel=X)) + (recoilE(sodium, DMmass=Y/constants.c**2, DMvel=230.*(units.km/units.s))+recoilE(iodine, DMmass=Y/constants.c**2, DMvel=230.*(units.km/units.s)))
ax.plot_surface(X, Y, Z, rstride=8, cstride=8, alpha=0.3)
cset = ax.contourf(X, Y, Z, zdir='z', offset=-10, cmap=cm.coolwarm)
cset = ax.contourf(X, Y, Z, zdir='x', offset=-200, cmap=cm.coolwarm)
cset = ax.contourf(X, Y, Z, zdir='y', offset=-3000, cmap=cm.coolwarm)

ax.set_xlabel('DM velocity (km/s)')
ax.set_xlim(-200, 1000)
ax.set_ylabel('DM mass (GeV)')
ax.set_ylim(-3000, 13000)
ax.set_zlabel('Recoil Energy (keV)')
ax.set_zlim(-10, 200)

# #
#plt.xlabel('Dark Matter Velocity (km/s)')
#plt.ylabel('Recoil Energy (keV)')
plt.grid(True)
plt.legend(loc='best', frameon=True)
plt.title('Recoil Energy responses for NaI with varying DM mass and velocity')
plt.show()