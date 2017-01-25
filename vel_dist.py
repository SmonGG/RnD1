from scipy.stats import maxwell
import numpy as np
import astropy as ap
import astropy.units as units
import astropy.constants as constants
import matplotlib.pyplot as plt

# velocity range
vrange = 10.**(np.linspace(1.,3.,num=100000.)) * units.km/units.s

# scipy states that Maxwell.pdf(x) = sqrt(2/pi) x^2 exp(-x^2 / 2) 
# but full equation in the help file is with loc and scale arrays as
# Maxwell.pdf(x, loc, scale) = Maxwell.pdf(y) / scale 
# where y = (x - loc) / scale redfines the input array x to y

# we want the full Maxwell velocity distribution given by
# a = sqrt(kT/m) with sqrt(2/pi) v^2 exp(-v^2 / (2a^2)) / a^3

# so if we say scale = a and loc = 0
# then y = v / scale = v / a and hence 
# sqrt(2/pi) (v/a)^2 exp(-v^2 / (2a^2)) / a

# DM particle mass
mass = 10. * units.GeV / constants.c**2.

# kT ~ 0.5 mv^2 = 0.5 * mass * (230km/s)**2.
temp = (0.5 * mass * (230. * units.km/units.s)**2)/constants.k_B
print "Effective temp for DM particles moving at 230km/s ", temp.to(units.Kelvin)
a = np.sqrt(constants.k_B * temp/mass)
print "Normalisation velocity ",a
scale = a

rv = maxwell.pdf(vrange, loc=0, scale=scale)
plt.plot(vrange,rv,'-',color='b')
plt.show()

v_values = maxwell.rvs(size=10000, loc=0, scale=scale)
plt.hist(v_values,1000,color='b',normed=1)
plt.show()


# Apply the v_values to the rate calculation from dmdetect
# plot rate vs the BINNED v_values from plot.hist
# plot rate * maxwell.pdf(BINNED v_values, loc=0, scale=scale) vs BINNED v_values