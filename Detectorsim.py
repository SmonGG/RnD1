# -*- coding: utf-8 -*-

"""Routine for calculating DM collision rates with a generic direct detection experiment."""

__author__ = 'Alan Duffy & Simon Goode'
__email__ = 'mail@alanrduffy.com, simongoode95@gmail.com'
__version__ = '0.1.0'

import astropy as ap
from astropy import units
from astropy import constants
import scipy.stats as spy
import scipy.special
from scipy.stats import rv_continuous
import matplotlib.pyplot as plt
import numpy as np

class Detector:
        def __init__(self, name, atomic_mass, mass_of_detector):
            """Class for variables of Detector properties"""
            self.name=name
            self.atomic_num=atomic_mass 
            self.atomic_mass=atomic_mass * units.g / units.mol
            self.mass_of_detector=mass_of_detector * units.kg

class DM:
        def __init__(self, density, velocity, mass, cross_section):
            """Class for variables of DM particle"""
            self.density=density * units.GeV / constants.c**2 / units.cm**3 #units.msolMass / units.mpc**3
            self.velocity=velocity * units.km / units.s
            self.mass=mass * units.GeV / constants.c**2
            self.cross_section=cross_section * units.pbarn

def fn_epsilon(x):
	'''For use in the Bhattacharjee Distribution Funciton'''
	k = -1.47
	eps1 = (1+x)**k
	eps2 = np.exp(-1.*(x**(1-k)))
	epsT = eps1 * eps2
	return epsT

class gen_Bhattacharjee(rv_continuous):
	'''Generate the Velocity Probability Distribution Function used in Bhattacharjee 2013'''
	def _pdf(self, x):
		vel_0 = 339. * units.km / units.s
		vel_max = 516. * units.km / units.s
		vpdf = (4.*np.pi*(x**2)) * (fn_epsilon(((x**2)/(vel_0**2)).value)-fn_epsilon(((vel_max**2)/(vel_0**2)).value))
		return vpdf - min(vpdf)

Bhattacharjee = gen_Bhattacharjee(name='Bhattacharjee', a=0)

def normalise(VDF):
    '''Normalise the Distribution'''
    nsum = max(VDF)
    norm = [float(i)/nsum for i in VDF]
    return norm

def noisy(vdf, mu=0, sigma=0.1):
    noise = np.random.normal(mu, sigma, len(vdf))
    soln = vdf * (1. + noise)
    return soln

def velocityPDF(name, VDF):
	if name.lower() == 'bhattacharjee':
		V_PDF = normalise(Bhattacharjee.pdf(VDF))
	elif name.lower() == 'maxwell':
		V_PDF = normalise(spy.maxwell.pdf(VDF, scale=162.6))
	else:
		print 'What velocity distribution do you want? You said: ', name
	return V_PDF

def calcrate(det,dm_p,vdf):
    """ Calculate flythough rates per second """
    rate = (dm_p.velocity*(dm_p.density/dm_p.mass))*vdf
    return rate#*len(vdf)

def calcsigma(det,dm_p):
    """ Modify WIMP-nucleon cross section for collision cross-section with detector material """

    Na_mass_reduced = det.atomic_num[0] * (dm_p.mass/constants.m_p) / (det.atomic_num[0] + dm_p.mass/constants.m_p)
    I_mass_reduced = det.atomic_num[1] * (dm_p.mass/constants.m_p) / (det.atomic_num[1] + dm_p.mass/constants.m_p)
    
    Na_sigma = det.atomic_num[0]**2. * Na_mass_reduced**2. * dm_p.cross_section
    I_sigma = det.atomic_num[1]**2. * I_mass_reduced**2. * dm_p.cross_section

    return Na_sigma, I_sigma

def calccol(sigma,rate,det):
    
    Na_colrate = sigma[0] * rate * constants.N_A*(det.mass_of_detector/det.atomic_mass[0])
    I_colrate = sigma[1] * rate * constants.N_A*(det.mass_of_detector/det.atomic_mass[1])
    return Na_colrate, I_colrate

def recoilE(AtomMass, DMmass, DMvel, costheta = -1.):
    """Calculates and returns recoil energy for DM (mass, velocity) colliding with atom (mass in GeV/c**2) at scattering angle (degrees)."""
    ## Gives DM particles random angles between 0 and 360 degrees.
    costheta = (np.random.random(size=len(DMvel))*2.)-1.
    ## check that atomic mass units are correct for calc (GeV s2 / m2)
    if AtomMass.unit == 'g / mol':
        AtomMass = AtomMass * units.mol / units.g
        AtomMass = AtomMass * (0.932 * units.GeV/constants.c**2.)
    else:
        print '***units error - recoilE()***'

    KE = (0.5 * DMmass * constants.c**2. * ((DMvel/constants.c)**2.)).to(units.GeV)
    r = (4. * DMmass * AtomMass) / ((DMmass + AtomMass)**2.)
    E_recoil = (KE*r*(1.-costheta))/2.
    return E_recoil.to(units.keV)

def mtmtransfer(M, E):
    """ The momentum transfer q """
    q = np.sqrt(2.*M*E)
    return (q/constants.hbar).to('fm**-1')

def formfactor(M, E):
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

    Fq_sq = (3. * bessel/QR )**2. * np.exp(-1. * q**2. * s**2.)
    Fq_norm = max(Fq_sq)
    FF = Fq_sq/Fq_norm
    if max(FF) == 1.0:
        print 'Form Factor Success!', min(FF), 'to', 1.0
    else:
        print 'Form Factor Error! ', max(FF)
    return FF

def quenchingfactor(keVnr):
    """Calculate Quenching Factor for nuclear recoil energies and return electron equivalent energies (keVee)"""
    quench = 0.
    if keVnr >= 0. and keVnr <= 10.:
        quench = ((0.4 * keVnr) + 12.)/100.
    elif keVnr >= 10. and keVnr <=1000.:
        quench = ((keVnr / 8.) + 14.75)/100.
    else:
        print 'keVnr value out of QF range - returning 0...'
    return quench

def QFcuts(keVee):
	##Set the detector's lowest detection energy here (keVee)
    lowerlim = 1.
    arraylen = len(keVee)
    keVee = keVee[keVee.value > lowerlim]
    cuts = 100.*(arraylen-len(keVee))/arraylen
    print 'Lower limit set to', lowerlim, 'keVee. Datapoints cut: ', cuts, '%'
    return keVee

def extraplots(atm_mass, veldist):
	'''Showing the differences between 10GeV and 100GeV WIMP recoil energies'''
	Na10=recoilE(atm_mass[0], 10.*units.GeV/constants.c**2, veldist)
	Na100=recoilE(atm_mass[0], 100.*units.GeV/constants.c**2, veldist)
	I10=recoilE(atm_mass[1], 10.*units.GeV/constants.c**2, veldist)
	I100=recoilE(atm_mass[1], 100.*units.GeV/constants.c**2, veldist)
	plt.scatter(veldist, Na10, color='red', alpha=0.3, label='Sodium-10GeV WIMP')
	plt.scatter(veldist, Na100, color='orange', alpha=0.3, label='Sodium-100GeV WIMP')
	plt.scatter(veldist, I10, color='blue', alpha=0.3, label='Iodine-10GeV WIMP')
	plt.scatter(veldist, I100, color='purple', alpha=0.3, label='Iodine-100GeV WIMP')
	'''plot parameters'''
	plt.xlabel('Dark Matter Velocity (km/s)')
	plt.ylabel('Unaltered Nuclear Recoil Energy (keVnr)')
	plt.grid(True)
	plt.legend(loc='best', frameon=True)
	plt.title('Comparison of 10GeV and 100GeV WIMP Nuclear Recoil Energies')
	plt.xlim(min(veldist.value)-max(veldist.value)/40,max(veldist.value)+max(veldist.value)/40)
	plt.ylim(-max(I100.value)/10,max(I100.value)+(max(I100.value)/10))
	plt.show()	

	return True

def run(detvar, dmvar='CDM'):
    """ Calculate rates and plot """

    if isinstance(detvar, dict):
        ## user is providing own parameters, i.e. det = ['name' : 'SABRE', 'atomic_mass' : 150., 'mass_of_detector' : 50.]
        det = Detector(detvar['name'],detvar['atomic_mass'],detvar['mass_of_detector'])
    elif detvar == 'SABRE':
    	AMsodium = 22.989769
    	AMiodine = 126.90447
        det = Detector('SABRE', [AMsodium,AMiodine], 50.)
        ## Note: detvar 'atomic_mass' has 2 values, [0] for Sodium and [1] for Iodine.
    else:
        print "What detector do you want? You said ", detvar
        ## Add a return crash statement

    if isinstance(dmvar, dict):
        ## user is providing own parameters, i.e. particle = ['density' : 0.3 'GeV/cm^3', 'velocity' : 230., 'mass' : 100., 'cross_section' : 1e-5]
        dm_p = DM(dmvar['density'], dmvar['velocity'], dmvar['mass'], dmvar['cross_section'])
    elif dmvar == 'CDM':
    	##Generate [samplenum] random velocities between DMvmin and DMvmax (in km/s)
    	DMvmin = 25.
    	DMvmax = 800.
        samplenum = 100000
        randomvels = ((DMvmax-DMvmin)*np.random.random(size=samplenum)) + DMvmin
        
        dm_p = DM(0.3, randomvels, 100., 1e-5)
    else:
        print "What particle is this? You passed ", dmvar

    print ' '
    print '*** WIMP MASS: ', dm_p.mass*constants.c**2, '/c**2 ***'
    print ' '

    print '---VELOCITY DISTRIBUTION---'
    V_PDF = 'Bhattacharjee'
    print 'Velocity Distribution used: ', V_PDF
    Th_PDF = velocityPDF(V_PDF, np.linspace(DMvmin,DMvmax,samplenum))
    V_PDF = velocityPDF(V_PDF, dm_p.velocity)
    '''plot theoretical VPDF'''
    plt.plot(np.linspace(DMvmin,DMvmax,samplenum), Th_PDF, color='red', lw=4., label='Theoretical')
    '''add noise'''
    V_PDF = noisy(V_PDF)
    print 'Most probable velocity: ', dm_p.velocity[np.asarray(V_PDF).argmax()]
    '''plot velocity sample'''
    plt.scatter(dm_p.velocity, V_PDF, color='blue', alpha=0.1, label='Velocity Sample')
    '''plot parameters'''
    plt.xlabel('Dark Matter Velocity (km/s)')
    plt.ylabel('F(v) (normalized)')
    plt.grid(True)
    plt.legend(loc='best', frameon=True)
    plt.title('Dark Matter Velocity Sample and Theoretical Distribution')
    plt.xlim(DMvmin-DMvmax/40,DMvmax+DMvmax/40)
    plt.ylim(-max(V_PDF)/10,max(V_PDF)+(max(V_PDF)/10))
    plt.show()
    print ' '

    print '---NUMBER DENSITY---'
    print "Number density of Sodium atoms per unit mass: ", (constants.N_A / det.atomic_mass[0])
    print "Number density of Iodine atoms per unit mass: ", (constants.N_A / det.atomic_mass[1])
    print "Number density of DM particles per unit volume: ", (dm_p.density / dm_p.mass).cgs
    print ' '

    print '---FLYTHROUGH RATE---'
    rate = calcrate(det, dm_p, V_PDF)
    print "Average Flythough rate per sec: ", np.ndarray.mean(rate.cgs)
    print "Average Flythough rate per year: ", np.ndarray.mean((rate * units.yr).cgs)
    '''plot flythrough rates (1/s) for the velocity sample'''
    plt.scatter(dm_p.velocity, rate.cgs, color='green', alpha=0.1, label='Flythrough rates of velocity sample')
    '''plot parameters'''
    plt.xlabel('Dark Matter Velocity (km/s)')
    plt.ylabel('Flythrough rate (1/s)')
    plt.grid(True)
    plt.legend(loc='best', frameon=True)
    plt.title('Dark Matter Flythrough rates (1/s) convolved with Velocity Probability Distribution')
    plt.xlim(DMvmin-DMvmax/40,DMvmax+DMvmax/40)
    plt.ylim(-max((rate.cgs).value)/10,max((rate.cgs).value)+(max((rate.cgs).value)/10))
    plt.show()
    print ' '

    print '---COLLISION CROSS-SECTIONS---'
    sigma = calcsigma(det, dm_p)
    print "Collision WIMP-nucleon cross-section per nucleon: ",dm_p.cross_section
    print "Collision cross-section per WIMP-Sodium atom pair: ",sigma[0].cgs
    print "Collision cross-section per WIMP-Iodine atom pair: ",sigma[1].cgs
    print ' '

    print '---THEORETICAL COLLISION RATES---'
    col_rate = calccol(sigma,rate,det)
    print ' -Sodium-'
    print "Average rate per sec: ",np.ndarray.mean(col_rate[0].cgs)
    print "Average rate per hour: ",np.ndarray.mean((col_rate[0] * units.hr).cgs)
    print "Average rate per day: ",np.ndarray.mean((col_rate[0] * units.day).cgs)
    print "Average rate per year: ",np.ndarray.mean((col_rate[0] * units.yr).cgs)
    print ' -Iodine-'
    print "Average rate per sec: ",np.ndarray.mean(col_rate[1].cgs)
    print "Average rate per hour: ",np.ndarray.mean((col_rate[1] * units.hr).cgs)
    print "Average rate per day: ",np.ndarray.mean((col_rate[1] * units.day).cgs)
    print "Average rate per year: ",np.ndarray.mean((col_rate[1] * units.yr).cgs)
    '''Plot Theoretical count rates (per year) from velocity sample'''
    plt.scatter(dm_p.velocity, (col_rate[0] * units.yr).cgs, color='orange', alpha=0.1, label='Theoretical Sodium collision rates')
    plt.scatter(dm_p.velocity, (col_rate[1] * units.yr).cgs, color='purple', alpha=0.1, label='Theoretical Iodine collision rates')
    '''plot parameters'''
    plt.xlabel('Dark Matter Velocity (km/s)')
    plt.ylabel('Theoretical Collision rate (1/yr)')
    plt.grid(True)
    plt.legend(loc='best', frameon=True)
    plt.title('Theoretical Dark Matter Collision rates for Sodium and Iodine nuclei (1/yr)')
    plt.xlim(DMvmin-DMvmax/40,DMvmax+DMvmax/40)
    plt.ylim(-max(((col_rate[1] * units.yr).cgs).value)/10,max(((col_rate[1] * units.yr).cgs).value)+(max(((col_rate[1] * units.yr).cgs).value)/10))
    plt.show()
    print ' '

    print '---NUCLEAR RECOIL ENERGIES---'
    '''Calculate Sodium and Iodine Nuclear Recoil Energies'''
    Na_ER = recoilE(det.atomic_mass[0], dm_p.mass, dm_p.velocity)
    I_ER = recoilE(det.atomic_mass[1], dm_p.mass, dm_p.velocity)
    '''Plot Sodium and Iodine Nuclear Recoil Energies'''
    plt.hist(Na_ER, color='orange', alpha = 0.2, label='Sodium Unaltered Recoil Energies', bins=500)
    plt.hist(I_ER, color='purple', alpha = 0.2, label='Iodine Unaltered Recoil Energies', bins=500)
    '''plot parameters'''
    plt.xlabel('Nuclear Recoil Energy (keVnr)')
    plt.ylabel('Entries')
    plt.grid(True)
    plt.legend(loc='best', frameon=True)
    plt.title('Distribution of Unaltered Nuclear Recoil Energies for Sodium and Iodine')
    plt.xlim(-0.02*max(I_ER.value),max(I_ER.value)+0.02*max(I_ER.value))
    plt.show()
    
    '''Calculate Sodium and Iodine Helm Form Factors'''
    print ' -Evaluating Helm Form Factors...'
    Na_FF = formfactor(dm_p.mass, Na_ER)
    I_FF = formfactor(dm_p.mass, I_ER)
    '''Plot Sodium and Iodine Helm Form Factors'''
    plt.scatter(Na_ER, Na_FF, color='red', alpha=0.1, label='Sodium Helm Form Factors')
    plt.scatter(I_ER, I_FF, color='blue', alpha=0.1, label='Iodine Helm Form Factors')
    '''plot parameters'''
    plt.xlabel('Unaltered Nuclear Recoil Energy (keVnr)')
    plt.ylabel('Helm Form Factor (|F(v)|^2)')
    plt.grid(True)
    plt.legend(loc='best', frameon=True)
    plt.title('Helm Form Factors for Sodium and Iodine from velocity sample')
    plt.xlim(min(I_ER.value)-max(I_ER.value)/40,max(I_ER.value)+max(I_ER.value)/40)
    plt.ylim(0.000001,1,)
    plt.yscale('log')
    plt.show()

    print ' -'
    '''Multiply Nuclear Recoil Energies with their respective Helm Form Factors'''
    Na_ER *= Na_FF
    I_ER *= I_FF
    print "Sodium Average (keVnr): ", np.ndarray.mean(Na_ER)
    print "Iodine Average (keVnr): ", np.ndarray.mean(I_ER)
    '''Plot True Nuclear Recoil Energies (histogram)'''
    plt.hist(Na_ER, color='orange', alpha = 0.2, label='Sodium True Recoil Energies', bins=500)
    plt.hist(I_ER, color='purple', alpha = 0.2, label='Iodine True Recoil Energies', bins=500)
    '''plot parameters'''
    plt.xlabel('Nuclear Recoil Energy (keVnr)')
    plt.ylabel('Entries')
    plt.grid(True)
    plt.legend(loc='best', frameon=True)
    plt.title('Distribution of True Nuclear Recoil Energies for Sodium and Iodine')
    plt.xlim(-0.02*max(I_ER.value),max(I_ER.value)+0.02*max(I_ER.value))
    plt.show()

    '''Plot True Nuclear Recoil Energies for velocity sample'''
    plt.scatter(dm_p.velocity, Na_ER, alpha=0.1, color='red', label='Sodium True Recoil Energies')
    plt.scatter(dm_p.velocity, I_ER, alpha=0.1, color='blue', label='Iodine True Recoil Energies')
    '''plot parameters'''
    plt.xlabel('Dark Matter Velocity (km/s)')
    plt.ylabel('True Nuclear Recoil Energy (keVnr)')
    plt.grid(True)
    plt.legend(loc='best', frameon=True)
    plt.title('True Nuclear Recoil Energies for Sodium and Iodine from velocity sample')
    plt.xlim(DMvmin-DMvmax/40,DMvmax+DMvmax/40)
    plt.ylim(-max(I_ER.value)/10,max(I_ER.value)+(max(I_ER.value)/10))
    plt.show()
    print ' '

    print '---QUENCHING FACTORS---'
    print '*** TODO: Find / create a more accurate quenching factor formula! ***'
    '''Calculate Quenching factors for Sodium and Iodine Nuclear Recoil Energies'''
    Na_QFactor = list(map(lambda x: quenchingfactor(Na_ER[x].value), np.arange(0,len(Na_ER),1)))
    I_QFactor = list(map(lambda x: quenchingfactor(I_ER[x].value), np.arange(0,len(I_ER),1)))
    print 'Average Quenching Factor for Sodium Nuclei: ', np.mean(Na_QFactor) * 100.,'%'
    print 'Average Quenching Factor for Iodine Nuclei: ', np.mean(I_QFactor) * 100.,'%'
    '''Plot Quenching Factors for Sodium and Iodine from velocity sample'''
    plt.scatter(Na_ER, Na_QFactor, color='blue', alpha=0.1, label='Sodium Quenching Factors')
    plt.scatter(I_ER, I_QFactor, color='red', alpha=0.1, label='Iodine Quenching Factors')
    '''plot parameters'''
    plt.xlabel('Nuclear Recoil Energy (keVnr)')
    plt.ylabel('Quenching Factor')
    plt.grid(True)
    plt.legend(loc='best', frameon=True)
    plt.title('Quenching Factors for Sodium and Iodine from Nuclear Recoil Energy sample')
    plt.xlim(-max(I_ER.value)/40,max(I_ER.value)+max(I_ER.value)/40)
    plt.ylim(min(I_QFactor)-max(I_QFactor)/20,max(I_QFactor)+(max(I_QFactor)/20))
    plt.show()
    print ' '

    print '---ELECTRON EQUIVALENT ENERGIES---'
    Na_EEenergies = Na_ER * Na_QFactor
    I_EEenergies = I_ER * I_QFactor
    print 'Sodium Average (keVee): ', np.mean(Na_EEenergies)
    print 'Iodine Average (keVee): ', np.mean(I_EEenergies)
    '''Plot the full set of Electron Equivalent energies - to be compared to the cuts below'''
    plt.hist(Na_EEenergies, color='black', alpha = 0.05, label='Cut Values', bins=500)
    plt.hist(I_EEenergies, color='black', alpha = 0.05, bins=500)
    print ' '

    print '---DATA POINT CUTS---'
    '''Cut data points that are below the lower detection limit of the detector'''
    Na_Cutenergies = QFcuts(Na_EEenergies)
    I_Cutenergies = QFcuts(I_EEenergies)
    Na_cutfactor = float(len(Na_Cutenergies))/len(dm_p.velocity)
    I_cutfactor = float(len(I_Cutenergies))/len(dm_p.velocity)
    '''Plot the electron equivalent energies of sodium and iodine from Nuclear Recoil Energy sample'''
    plt.hist(Na_Cutenergies, color='orange', alpha = 0.2, label='Sodium Electron Equivalent Energies', bins=500)
    plt.hist(I_Cutenergies, color='purple', alpha = 0.2, label='Iodine Electron Equivalent Energies', bins=500)
    '''plot parameters'''
    plt.xlabel('Electron Equivalent Energy (keVee)')
    plt.ylabel('Entries')
    plt.grid(True)
    plt.legend(loc='best', frameon=True)
    plt.title('Distribution of Electron Equivalent Energies of Sodium and Iodine from the nuclear recoil energy sample')
    plt.xlim(-0.02*max(I_Cutenergies.value),max(I_Cutenergies.value)+0.02*max(I_Cutenergies.value))
    plt.show()
    print ' '

    print '---FINALISED COLLISION RATES---'
    print ' -Sodium'
    print 'Average per sec: ', np.ndarray.mean(((col_rate[0] * units.s).cgs) * Na_cutfactor)
    print 'Average per min: ', np.ndarray.mean(((col_rate[0] * units.minute).cgs) * Na_cutfactor)
    print 'Average per hour: ', np.ndarray.mean(((col_rate[0] * units.hour).cgs) * Na_cutfactor)
    print 'Average per year: ', np.ndarray.mean(((col_rate[0] * units.year).cgs) * Na_cutfactor)
    print ' -Iodine'
    print 'Average per sec: ', np.ndarray.mean(((col_rate[1] * units.s).cgs) * I_cutfactor)
    print 'Average per min: ', np.ndarray.mean(((col_rate[1] * units.minute).cgs) * I_cutfactor)
    print 'Average per hour: ', np.ndarray.mean(((col_rate[1] * units.hour).cgs) * I_cutfactor)
    print 'Average per year: ', np.ndarray.mean(((col_rate[1] * units.year).cgs) * I_cutfactor)
    print ' '

    print '---'
    print 'Showing Extra Plots...'
    extraplots(det.atomic_mass, dm_p.velocity)
    return True

run('SABRE')