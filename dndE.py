import numpy as np
import matplotlib.pyplot as plt 
import scipy.stats as spy
from scipy.integrate import quad

#DM Parameters
m = float(input("Enter the mass of dark matter (GeV): "))*1e9
num = int(input("Enter the number of random samples: "))
v_min = 0 #Minimum DM velocity
v_esc = 20 #Galactic escape velocity (max DM velocity)

#Automatically limit VDF?
autolimit = "ask"
	#yes: automatically limit velocity distribution
	#ask: the program will ask the user
	#no: disable the limit feature

#Energy distribution function
def epdf(m,v):
	kT = (m*np.power(v,2))/2 	#kT ~ KE 
	g1 = 2 * np.pi * kT
	g2 = 4 * np.pi * np.power(v, 2)
	g3 = (m*v)/(2*kT)
	g4 = np.power(np.e, -(g3))
	ans = np.sqrt(np.power((m/g1), 3)) * g2 * g4
	return ans

#Theoretical plot
theo = epdf(m, np.arange(0.01,20,0.1))
plt.plot(np.arange(0.01,20,0.1), theo)

#Maxwellian sampling
mxw = spy.maxwell.rvs(size=num)
soln = epdf(m, mxw)
plt.scatter(mxw, soln, alpha=0.5)

#Lognorm sampling
lgn = spy.lognorm.rvs(2,size=num)
soln2 = epdf(m, lgn)
plt.scatter(lgn, soln2, alpha=0.5, color='r')

#Plot unity
##def integrand(v, m):
##	kT = (m*np.power(v,2))/2
##	return np.sqrt(np.power((m/2*np.pi*kT), 3))*(4*np.pi*np.power(v,2))*(np.power(np.e,(-(m*v)/(2*kT))))
##I=quad(integrand, 0.1, np.inf, args=(m))

from astropy import units
from astropy import constants

mxw = spy.maxwell.rvs(size=5*num) * units.km / units.s**2.

# Select velocities with a mask to be fancy
mxw = ma.masked_where(mxw > v_min AND mxw < v_max, mxw)

# OR by hand
mxw = mxw[mxw > v_min]
mxw = mxw[mxw < v_max]

from numpy import random

tmp = random.randint(0,num,dtype='int64')

if len(mxw) < len(tmp):
	do limit
else:
	mxw = mxw[tmp]
	
	


#Fit samples within velocity limits
def limit(v_array, v_min, v_max, type='t_mxw'):
	i = 0
	reps = 0
	replimit = 10000
	global m
	global v_esc
	if type == t_mxw:
		replacement = spy.maxwell.rvs(size=1)
	elif type == t_lgn:
		replacement = spy.lognorm.rvs(2,size=1)
	else:
		print ("ERROR - TYPE NOT RECOGNISED!")
		pass
	
	while i < len(v_array) and reps < replimit:
		if v_array[i] < v_min or v_array[i] > v_max:
			v_array[i] = epdf(m, replacement)
			reps = reps + 1
		else:
			i = i+1
	
	if reps == replimit or reps > replimit:
		print("Replacement limit (",replimit,") reached. Consider changing PDF!")
	else:
		print("Velocity distribution was successfully limited!")
		print("Limits: ", v_min, "to", v_max, "with", reps, "replacements.")

#Global variable determines if VDF limits are applied
def autolimit(autolim):
	global autolimit
	if autolim == "no":
		pass
	elif autolim == "yes":
		limit(mxw, v_min, v_esc)
		limit(lgn, v_min, v_esc)
	elif autolim == "ask":
		user_response = input("Would you like to limit your VDFs? (y/n): ")
		while user_response != "y" and user_response != "n":
			user_response = input("Please enter yes (y) or no (n): ")
		if user_response == "y":
			limit(mxw, v_min, v_esc, t_mxw)
			limit(lgn, v_min, v_esc, t_lgn)
		else:
			pass
	else:
		print("WARNING! ",autolim, "is not valid!")
		user_response = input("Would you like to limit your VDFs? (y/n): ")
		while user_response != "y" and user_response != "n":
			user_response = input("Please enter yes (y) or no (n): ")
		if user_response == "y":
			limit(mxw, v_min, v_esc, mxw)
			limit(lgn, v_min, v_esc, lgn)
		else:
			pass		
#Main
autolimit(autolimit)
print("Showing", num, "samples from maxwellian and lognorm distributions")
##print("Area under curve: ", I)
plt.show()
print("Program terminated.")