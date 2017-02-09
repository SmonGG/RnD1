import astropy.units as units
import scipy.special
import astropy.constants as constants
import numpy as np

def mtmtransfer(M, E):
	""" The momentum transfer q """

	q = np.sqrt(2.*M.to('GeV s**2 / m**2') * E.to('keV') * 1e-6) / (0.197*(units.GeV/constants.c**2) * units.fm)

	return (q/constants.c).to('fm**-1')

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
	R = (0.91 * (M.to('GeV s**2 / m**2').value)**(1./3.) + 0.3)*units.fm
	s = 1. * units.fm
	R_one = np.sqrt(R**2 - 5.*s**2)

	q = mtmtransfer(M,E)
	# Eqn 9 from Duda, Kemper and Gondolo 2007 
	QR = q*R_one
	bessel = scipy.special.j1(QR)
	Fq_sq = (3. * bessel/QR )*2. * np.exp(-1. * q**2. * s**2.)

	return Fq_sq