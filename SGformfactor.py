import astropy.units as units
import scipy.special
import astropy.constants as constants
import numpy as np
import matplotlib.pyplot as plt

def mtmtransfer(M, E):
	""" The momentum transfer q """

	q = np.sqrt(2.*M*E)
	#q = np.sqrt(2.*M.to('GeV') * E.to('keV') * 1e-6) / (0.197*(units.GeV/constants.c**2) * units.fm)

	return (q/constants.hbar).to('fm**-1')

def tpf_val(elementsymbol='I'):
	""" Take 2 Parameter Fermi coefficient from Duda, Kemper and Gondolo 2007 """

	if elementsymbol.lower() == 'i':
		coeff = {'t': 2.3*units.fm, 'c':5.5931*units.fm, 'r': 4.749*units.fm, 'a': 0.523*units.fm }	
	elif elementsymbol.lower() == 'na':
		## Sodium
		coeff = {'t': 2.3*units.fm, 'c':2.9393*units.fm, 'r': 2.994*units.fm, 'a': 0.523*units.fm }
	else:
		print "What element have you asked for?", elementsymbol

	return coeff

def formfactor(M,E):
	""" For now just use Helm Form factor """

	#if detvar.upper() == 'SABRE':
	#	coeff = tpf_val('Na')
	#	coeff = tpf_val('I')
	#else:
	#	print "Create other detector technologies "

	# Eqn 31 from Duda, Kemper and Gondolo 2007 form factor for Two-parameter Fermi 
	#radial_rho = rho_c / (np.exp( (coeff['r'] - coeff['c'])/coeff['a']) +1.)

	## Helm form factor form DarkSUSY 4.1 code 

	# Eqn 19 from Duda, Kemper and Gondolo 2007 
	R = (0.91 * (((M*constants.c**2)/units.GeV))**(1./3.) + 0.3)*units.fm
	s = 1. * units.fm
	R_one = np.sqrt(R**2. - 5.*s**2.)

	q = mtmtransfer(M,E)
	# Eqn 9 from Duda, Kemper and Gondolo 2007 
	QR = (q*R_one).value
	QR_norm = (1.e-6*R_one).value
	bessel = scipy.special.j1(QR)
	f1 = (3. * bessel/QR )*2.
	f2 = np.exp(-1. * q**2. * s**2.)
	'''	print 'average f1: ', np.mean(f1)
	print 'average f2: ', np.mean(f2)
	print 'Average Fq_sq: ', np.mean(f1 * f2)'''
	Fq_sq = (3. * bessel/QR )**2. * np.exp(-1. * q**2. * s**2.)
	Fq_norm = Fq_sq[0]#(3. * bessel/QR )**2. * np.exp(-1. * q**2. * s**2.)
	#plt.plot(q, Fq_sq)
	print Fq_sq / Fq_norm
	plt.semilogy(q, Fq_sq/Fq_norm)
	plt.xlim(0,1.5)
	plt.ylim(1.e-8,1.)

	return Fq_sq

def recoilE(AtomMass, DMmass, DMvel, costheta = -1.):
    """Calculates and returns recoil energy for DM (mass, velocity) colliding with atom (mass in GeV/c**2) at scattering angle (degrees)."""
    ## Gives DM particles random angles between 0 and 360 degrees.
    #costheta = (np.random.random(size=len(DMvel))*2.)-1.
    ## check that atomic mass units are correct for calc (GeV s2 / m2)
    if AtomMass.unit == 'g / mol':
        AtomMass = AtomMass * units.mol / units.g
        AtomMass = AtomMass * (0.932 * units.GeV/constants.c**2.)
    else:
        print '***units error - recoilE()***'

    KE = (0.5 * DMmass * DMvel**2).to(units.GeV)
    r = (4. * DMmass * AtomMass) / ((DMmass + AtomMass)**2.)
    E_recoil = (KE*r*(1.-costheta))/2.
    return E_recoil.to(units.keV)

TwoParameterFermi = tpf_val(elementsymbol='I')
# in g/mol because the recoilE() function converts it to GeV s2/m2
sodium = 22.989769 * units.g / units.mol
iodine = 126.90447 * units.g / units.mol
Na_RE = recoilE(sodium, 100.*units.GeV/constants.c**2, np.linspace(25.,800.,1000.)*units.km/units.s)
print 'Average Sodium RE: ', np.mean(Na_RE)
I_RE = recoilE(iodine, 100.*units.GeV/constants.c**2, np.linspace(25.,800.,10000.)*units.km/units.s)
print 'Average Iodine RE: ', np.mean(I_RE)
# in GeV s2/m2 because those are the appropriate units for the formfactor() function.
sodium = 22.989769 * (0.932 * units.GeV/constants.c**2.)
iodine = 126.90447 * (0.932 * units.GeV/constants.c**2.)
Na_Form = formfactor(sodium, Na_RE)
#print 'Sodium Form Factor: ',np.min(Na_Form)
I_Form = formfactor(iodine, I_RE)
#print 'Iodine Form Factor: ',np.min(I_Form)

plt.xlabel('q[fm-1]')
plt.ylabel('|F(q)|2')
plt.grid(True)
plt.legend(loc='best', frameon=True)
plt.title('Form Factors of Sodium (blue) and Iodine (green)')
plt.show()