# -*- coding: utf-8 -*-

"""Routine for calculating DM collision rates with a generic direct detection experiment."""

__author__ = 'Alan Duffy & Simon Goode'
__email__ = 'mail@alanrduffy.com, simongoode95@gmail.com'
__version__ = '0.1.2'

import astropy as ap
from astropy import units
from astropy import constants
import scipy.stats as spy
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import rv_continuous
import scipy.special

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

class recoildetection_dmvdf_gen(rv_continuous):
    "Dark Matter VPDF"
    def pdf(self, m, v, T):
        k = constants.k_B   #Boltzmann Constant
        return (np.sqrt((m/(2.*np.pi*k*T))**3)) * (4.*np.pi*(v**2)) * (np.exp(-(m*(v**2))/(2.*k*T)))

    def rvs(self, size, m, T, mu=0, sigma=0.5):
        v = spy.maxwell.rvs(size=size, loc=mu, scale=sigma) * units.km/units.s
        k = constants.k_B   #Boltzmann Constant
        return v, (np.sqrt((m/(2.*np.pi*k*T))**3)) * (4.*np.pi*(v**2)) * (np.exp(-(m*(v**2))/(2.*k*T))) / 600000

dmvdf = recoildetection_dmvdf_gen(name= 'dmvdf')
# #note: rvs function returns a duple. when assigning a variable that uses the rvs command, make sure that you split it: eg. test_v, test_DF = rvs(...)

def calcsigma(det,dm_p):
    """ Modify WIMP-nucleon cross section for collision cross-section with detector material """

    Na_mass_reduced = det.atomic_num[0] * (dm_p.mass/constants.m_p) / (det.atomic_num[0] + dm_p.mass/constants.m_p)
    I_mass_reduced = det.atomic_num[1] * (dm_p.mass/constants.m_p) / (det.atomic_num[1] + dm_p.mass/constants.m_p)
    NaI_mass_reduced = det.atomic_num[2] * (dm_p.mass/constants.m_p) / (det.atomic_num[2] + dm_p.mass/constants.m_p)
    
    Na_sigma = det.atomic_num[0]**2. * Na_mass_reduced**2. * dm_p.cross_section
    I_sigma = det.atomic_num[1]**2. * I_mass_reduced**2. * dm_p.cross_section
    NaI_sigma = det.atomic_num[2]**2. * NaI_mass_reduced**2. * dm_p.cross_section

    return Na_sigma, I_sigma, NaI_sigma

def calcrate(det,dm_p):
    """ Calculate flythough rates per second """
    avg_vel = np.ndarray.mean(dm_p.velocity)
    print 'Average velocity of sample: ', avg_vel
    rate = dm_p.velocity*(dm_p.density/dm_p.mass)

    return rate

def calccol(sigma,rate,det):
    
    Na_colrate = sigma[0] * rate * constants.N_A*(det.mass_of_detector/det.atomic_mass[0])
    I_colrate = sigma[1] * rate * constants.N_A*(det.mass_of_detector/det.atomic_mass[1])
    NaI_colrate = sigma[2] * rate * constants.N_A*(det.mass_of_detector/det.atomic_mass[2])
    return Na_colrate, I_colrate, NaI_colrate

def recoilE(AtomMass, DMmass, DMvel, costheta = -1.):
    """Calculates and returns recoil energy for DM (mass, velocity) colliding with atom (mass in GeV/c**2) at scattering angle (degrees)."""
    ## Gives DM particles random angles between 0 and 360 degrees.
    costheta = (np.random.random(size=len(DMvel))*2.)-1.
    ## check that atomic mass units are correct for calc (GeV s2 / m2)
    if AtomMass.unit != 'GeV s2 / m2':
        AtomMass = AtomMass * units.mol / units.g
        AtomMass = AtomMass * (0.932 * units.GeV/constants.c**2.)
    else:
        print '***units error - recoilE()***'

    KE = (0.5 * DMmass * constants.c**2. * ((DMvel/constants.c)**2.)).to(units.GeV)
    r = (4. * DMmass * AtomMass) / ((DMmass + AtomMass)**2.)
    E_recoil = (KE*r*(1.-costheta))/2.
    return E_recoil.to(units.keV)

def noisy(vdf, mu=0, sigma=0.5):
    noise = np.random.normal(mu, sigma, len(vdf))
    soln = vdf * (1. + 0.2*noise)
    return soln

def quenchingfactor(keVnr):
    """Calculate Quenching Factor for nuclear recoil energies and return electron equivalent energies (keVee)"""
    quench = 0.
    if keVnr >= 0. and keVnr <= 10.:
        quench = ((0.4 * keVnr) + 12.)/100.
    elif keVnr >= 10. and keVnr <=1000.:
        quench = ((keVnr / 8.) + 14.75)/100.
    else:
        print 'keVnr value out of QF range - returning 0...'
    ##print 'Quenching Factor: ', quench * 100.,'%'
    return quench
    #keVee = keVnr * quench
    #return keVee

def QFcuts(keVee):
    lowerlim = 2.
    arraylen = len(keVee)
    keVee = keVee[keVee > lowerlim]
    cuts = 100.*(arraylen-len(keVee))/arraylen
    print 'Lower limit set to', lowerlim, 'keVee. Datapoints cut: ', cuts, '%'
    return keVee

def energyrate(dm_p, det):
    """Accepts time period, calculates time-dependent velocity, graphs expected modulation"""
    time=int(input("Please Enter the Observing Period (in days): "))
    t=np.arange(0, time)
    v_esc=(232+15*np.cos((2*np.pi)*(((t)-152.5)/365.25)))*(units.km/units.s)
    E_esc=recoilE(det.atomic_mass[2], dm_p.mass, v_esc)
    rate_1 = dm_p.cross_section*v_esc*(dm_p.density/dm_p.mass)*(ap.constants.N_A/det.atomic_mass[2])*det.mass_of_detector
    rate_2 = dm_p.cross_section*E_esc*(dm_p.density/dm_p.mass)*(ap.constants.N_A/det.atomic_mass[2])*det.mass_of_detector
    plt.plot(t, rate_1*det.mass_of_detector)
    plt.ylabel('Counts (per day)')
    plt.xlabel('Time (days)')
    plt.title('Predicted Modulation')
    plt.grid(True)
    plt.show()

    plt.plot(t, rate_2*det.mass_of_detector)
    plt.ylabel('Energy of Counts (per day)(keVnr)')
    plt.xlabel('Time (days)')
    plt.title('Predicted Modulation')
    plt.grid(True)
    plt.show()

def mtmtransfer(M, E):
    """ The momentum transfer q """

    q = np.sqrt(2.*M*E)
    #q = np.sqrt(2.*M.to('GeV') * E.to('keV') * 1e-6) / (0.197*(units.GeV/constants.c**2) * units.fm)

    return (q/constants.hbar).to('fm**-1')

def formfactor(M, E):
    """ For now just use Helm Form factor """

    #if detvar.upper() == 'SABRE':
    #   coeff = tpf_val('Na')
    #   coeff = tpf_val('I')
    #else:
    #   print "Create other detector technologies "

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
    ''' print 'average f1: ', np.mean(f1)
    print 'average f2: ', np.mean(f2)
    print 'Average Fq_sq: ', np.mean(f1 * f2)'''
    Fq_sq = (3. * bessel/QR )**2. * np.exp(-1. * q**2. * s**2.)
    Fq_norm = max(Fq_sq)#(3. * bessel/QR )**2. * np.exp(-1. * q**2. * s**2.)
    FF = Fq_sq/Fq_norm
    if max(FF) == 1.0:
        print 'Form Factor Success!', min(FF), 'to', max(FF)
    else:
        print 'Form Factor Error! ', max(FF)
    return FF

def run(detvar, dmvar='CDM'):
    """ Calculate rates and plot """

    if isinstance(detvar, dict):
        ## user is providing own parameters, i.e. det = ['name' : 'SABRE', 'atomic_mass' : 150., 'mass_of_detector' : 50.]
        det = Detector(detvar['name'],detvar['atomic_mass'],detvar['mass_of_detector'])
    elif detvar == 'SABRE':
        det = Detector('SABRE', [22.989769,126.90447,149.89], 50.)
        ## Note: detvar 'atomic_mass' has 3 values, [0] for Sodium, [1] for Iodine and [2] for NaI. All are in AMU.
    else:
        print "What detector do you want? You said ", detvar
        ## Add a return crash statement

    if isinstance(dmvar, dict):
        ## user is providing own parameters, i.e. particle = ['density' : 0.3 'GeV/cm^3', 'velocity' : 230., 'mass' : 100., 'cross_section' : 1e-5]
        dm_p = DM(dmvar['density'], dmvar['velocity'], dmvar['mass'], dmvar['cross_section'])
    elif dmvar == 'CDM':
        samplenum = 10000.
        dm_p = DM(0.3, spy.maxwell.rvs(size=samplenum, loc=0, scale=143.75), 100., 1e-5)
    else:
        print "What particle is this? You passed ", dmvar

    ## Check numbers
    print "Number density of Na atoms per unit mass ", (constants.N_A / det.atomic_mass[0])
    print "Number density of I atoms per unit mass ", (constants.N_A / det.atomic_mass[1])
    print "Number density of NaI atoms per unit mass ", (constants.N_A / det.atomic_mass[2])
    print "Number density of DM particles per unit volume ", (dm_p.density / dm_p.mass).cgs

    rate = calcrate(det, dm_p)
    print "Average Flythough rate per s ", np.ndarray.mean(rate.cgs)
    print "Average Flythough rate per year ", np.ndarray.mean((rate * units.yr).cgs)

    sigma = calcsigma(det, dm_p) 
    print "Collision WIMP-nucleon cross section per nucelon ",dm_p.cross_section
    print "Colllision-cross section per WIMP-Na atom pair ",sigma[0].cgs
    print "Colllision-cross section per WIMP-I atom pair ",sigma[1].cgs
    print "Colllision-cross section per WIMP-NaI atom pair ",sigma[2].cgs

    col_rate = calccol(sigma,rate,det) * V_PDF
    print "Average Na Collision rate per s ",np.ndarray.mean(col_rate[0].cgs)
    print "Average Na Collision rate per hour ",np.ndarray.mean((col_rate[0] * units.hr).cgs)
    print "Average Na Collision rate per day ",np.ndarray.mean((col_rate[0] * units.day).cgs)
    print "Average Na Collision rate per year ",np.ndarray.mean((col_rate[0] * units.yr).cgs)

    print "Average I Collision rate per s ",np.ndarray.mean(col_rate[1].cgs)
    print "Average I Collision rate per hour ",np.ndarray.mean((col_rate[1] * units.hr).cgs)
    print "Average I Collision rate per day ",np.ndarray.mean((col_rate[1] * units.day).cgs)
    print "Average I Collision rate per year ",np.ndarray.mean((col_rate[1] * units.yr).cgs)

    print "Average NaI Collision rate per s ",np.ndarray.mean(col_rate[2].cgs)
    print "Average NaI Collision rate per hour ",np.ndarray.mean((col_rate[2] * units.hr).cgs)
    print "Average NaI Collision rate per day ",np.ndarray.mean((col_rate[2] * units.day).cgs)
    print "Average NaI Collision rate per year ",np.ndarray.mean((col_rate[2] * units.yr).cgs)

    Na_ER = recoilE(det.atomic_mass[0], dm_p.mass, dm_p.velocity)
    I_ER = recoilE(det.atomic_mass[1], dm_p.mass, dm_p.velocity)
    NaI_ER = recoilE(det.atomic_mass[2], dm_p.mass, dm_p.velocity)
    print "Average recoil energy of Sodium nuclei (keVnr): ", np.ndarray.mean(Na_ER)
    print "Average recoil energy of Iodine nuclei (keVnr): ", np.ndarray.mean(I_ER)
    print "Average recoil energy of NaI nuclei (keVnr): ", np.ndarray.mean(NaI_ER)

    Na_QFactor = list(map(lambda x: quenchingfactor(Na_ER[x].value), np.arange(0,len(Na_ER),1)))
    I_QFactor = list(map(lambda x: quenchingfactor(I_ER[x].value), np.arange(0,len(I_ER),1)))
    NaI_QFactor = list(map(lambda x: quenchingfactor(NaI_ER[x].value), np.arange(0,len(NaI_ER),1)))
    print 'Average Quenching Factor for Sodium Nuclei: ', np.mean(Na_QFactor) * 100.,'%'
    print 'Average Quenching Factor for Iodine Nuclei: ', np.mean(I_QFactor) * 100.,'%'
    print 'Average Quenching Factor for NaI Nuclei: ', np.mean(NaI_QFactor) * 100.,'%'

    Na_FFenergies = Na_ER * formfactor(dm_p.mass, Na_ER)
    I_FFenergies = Na_ER * formfactor(dm_p.mass, I_ER)
    NaI_FFenergies = NaI_ER * formfactor(dm_p.mass, NaI_ER)

    Na_QFenergies = Na_FFenergies * Na_QFactor
    I_QFenergies = I_FFenergies * I_QFactor
    NaI_QFenergies = NaI_FFenergies * NaI_QFactor

    print "Average electron equivalent energy of Sodium quenching factor (keVee): ", np.mean(Na_QFenergies)
    print "Average electron equivalent energy of Iodine quenching factor (keVee): ", np.mean(I_QFenergies)
    print "Average electron equivalent energy of NaI quenching factor (keVee): ", np.mean(NaI_QFenergies)
#energyrate(dm_p, det)

#plt.hist(dm_p.velocity, normed = 1, bins = samplenum/10, alpha = 0.4, color='blue', label='Velocity')
    plt.xlabel('Dark Matter Velocity (km/s)')
    plt.ylabel('Normalised Entries')
    plt.grid(True)
    plt.legend(loc='best', frameon=True)
    plt.title('Distribution of Velocities')
#plt.show()

    plt.hist(Na_ER, normed = 1, bins = samplenum/10, alpha = 0.2, color='blue', label='Na keVnr')
    plt.hist(I_ER, normed = 1, bins = samplenum/10, alpha = 0.2, color='red', label='I keVnr')
    plt.xlabel('Nuclear Recoil Energy (keVnr)')
    plt.ylabel('Normalised Entries')
    plt.grid(True)
    plt.legend(loc='best', frameon=True)
    plt.title('Distribution of Nuclear recoil energies')
    plt.show()

    Na_QFenergies=np.array(Na_QFenergies*units.keV)
    I_QFenergies=np.array(I_QFenergies*units.keV)
    NaI_QFenergies=np.array(NaI_QFenergies*units.keV)
#plt.hist(QFenergies, bins=samplenum/10, alpha=0.2, color='black', label='Removed')
    Na_QFenergies = QFcuts(Na_QFenergies)
    I_QFenergies = QFcuts(I_QFenergies)
    NaI_QFenergies = QFcuts(NaI_QFenergies)

    print '---'
    print 'Final Average Sodium Collision rates per year: ',np.ndarray.mean(((col_rate[0] * units.yr).cgs)*float(len(Na_QFenergies))/len(dm_p.velocity))
    print 'Final Average Iodine Collision rates per year: ',np.ndarray.mean(((col_rate[1] * units.yr).cgs)*float(len(I_QFenergies))/len(dm_p.velocity))
    print 'Final Average NaI Collision rates per year: ',np.ndarray.mean(((col_rate[2] * units.yr).cgs)*float(len(NaI_QFenergies))/len(dm_p.velocity))

    plt.scatter(Na_ER, Na_QFactor, label='Na Quenching Factor', alpha=0.4, color='blue')
    plt.scatter(I_ER, I_QFactor, label='I Quenching Factor', alpha=0.4, color='red')
    plt.xlabel('Nuclear Recoil Energy (keVnr)')
    plt.ylabel('Quenching Factor')
    plt.grid(True)
    plt.legend(loc='best', frameon=True)
    plt.title('Distribution of Quenching Factors used')
    plt.show()

    plt.hist(Na_QFenergies, normed=1, bins = samplenum/10, alpha = 0.2, color='blue', label='Na keVee')
    plt.hist(I_QFenergies, normed=1, bins = samplenum/10, alpha = 0.2, color='red', label='I keVee')
    plt.xlabel('Electron Equivalent Energy (keVee)')
    plt.ylabel('Normalised Entries')
    plt.grid(True)
    plt.legend(loc='best', frameon=True)
    plt.title('Distribution of Electron Equivalent energies')
    plt.show()

    return rate

run('SABRE')