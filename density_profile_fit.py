import pylab as plb
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
import numpy as np
from scipy.optimize import minimize
from scipy.optimize import curve_fit
import emcee
import os
import sys

# read 
data = np.loadtxt( 'densidad_bins.dat', dtype='str') 

rho = data[:,1]
d   = data[:,0]
err = data[:,2]

rho    = np.array( map(float, rho) )
d    = np.array( map(float, d) )
err    = np.array( map(float, err) )



min_chi2 = 10**9
# --------------------------------------------- simple power law ---------------------------------------------
# using "minimize"
Y = np.log10(rho)
def chi_2_single(x):
	A = x[0]
	n = x[1]	

	c = 8

	return np.sum( ( Y - ( A + n*np.log10(d/c) ) )**2 )  

x0 = np.array([1.3,-2])
# 		method: Nelder-Mead,Powell,CG,BFGS,Newton-CG,L-BFGS-B,TNC,COBYLA,SLSQP,dogleg,trust-ncg
res = minimize( chi_2_single, x0, method='SLSQP')

A_min = res.x[0]
n_min = res.x[1]
min_chi2 = chi_2_single(res.x)
print "\n\n Parametros sin quiebre (A, n, chi2): \t", A_min, n_min, min_chi2, '\n'

# using "curve_fit"

def log_rho(R, A, n): # A = log(rho_sol)
	R_sol = 8
	return A + n*np.log10( R/R_sol )

popt,pcov = curve_fit(log_rho, d, Y, p0=[ 0, -1]) 
print popt
print pcov

# plotting
fig, ax = plt.subplots()
ax = plt.subplot(1,1,1);
plt.errorbar( d,rho, yerr=err)
plt.plot(d, (10**A_min)*((d/8)**n_min))
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim(10, 250)
ax.set_ylim(10**-6, 1)
ax.set_xticks([10,20,30,50,100,200,250])
ax.set_xticklabels(['10', '20', '30','50','100','200','250'])



# --------------------------------------------- broken power law ---------------------------------------------
# using "minimize"
min_chi2 = 10**9
mini = np.min(d)
maxi = np.max(d)

rango = 0.1
r_b = mini+rango
i=1

# from min(d) to max(d) with steps of 0.1 in possible r_break values
while r_b <= (maxi-rango):
#while i < len(d)-1:
	#r_b = d[i]
	#print r_b
	d_1   = d[ d <= r_b ]
	density_1 = rho[ d <= r_b ]

	d_2   = d[ d >= r_b ]
	density_2 = rho[ d >= r_b ]

	c = 8	

	Y1 = np.log10(density_1)
	Y2 = np.log10(density_2)
	def chi_2(x):
		A2 = x[0]
		n1 = x[1]
		n2 = x[2]
		
		return np.sum( ( Y1 - ( A2 + (n2-n1)*np.log10(r_b/c) + n1*np.log10(d_1/c)) )**2 ) + np.sum( ( Y2 - ( A2 + n2*np.log10(d_2/c) ) )**2 )

	x0 = np.array([2,-3,-5.4])
	res = minimize( chi_2, x0, method='SLSQP')
	
	A2 = res.x[0]
	n1 = res.x[1]
	n2 = res.x[2]		
	A1 = A2 + (n2-n1)*np.log10(r_b/c)


	valor_chi2 = chi_2(res.x)
	if valor_chi2 < min_chi2:
		min_chi2  = valor_chi2
		min_res   = res.x
		r_break   = r_b

		d_1_min = np.linspace(np.min(d), r_break, 1000)
		d_2_min = np.linspace(r_break, np.max(d), 1000)
	
	r_b = r_b+0.1
	i=i+1

# the A1, A2, n1, n2 that minimize chi2
A2 = min_res[0]
n1 = min_res[1]
n2 = min_res[2]


A1 = A2 + (n2-n1)*np.log10(r_break/c)


print '\n Parametros con quiebre (r_break, A1, A2, n1, n2, chi2): ', r_break, A1, A2, n1, n2, min_chi2, '\n'

rho_0_1 = 10**A1
rho_0_2 = 10**A2

plt.xlabel(r'$d$ [kpc]')
plt.ylabel(r'$\rho$ [kpc$^{-3}$]')

#plt.xlim(10, 110)
#plt.ylim(0.0001, 2)

#ax.set_yscale('log')
#ax.set_xscale('log')

#plt.plot(d, density, 'ro')
#plt.plot(np.log10(d), np.log10(rho), 'ro')
#ax.errorbar( d, density, yerr = (err, err), ecolor='b', linestyle="None")
#ax.set_xticks([10,20,30,50,100,200,300])
#ax.set_xticklabels(['10', '20', '30','50','100','200','300'])


eje_Y_1_min = rho_0_1*(d_1_min/8)**n1
eje_Y_2_min = rho_0_2*(d_2_min/8)**n2

plt.plot(d_1_min, eje_Y_1_min, 'r' )
plt.plot(d_2_min, eje_Y_2_min, 'r' )

#print d_1_min, d_2_min
#plt.plot(np.log10(d_1), np.log10(eje_Y_1), 'r' )
#plt.plot(np.log10(d_2), np.log10(eje_Y_2), 'r' )

#ax_0 = plt.subplot(2,1,1)
#plt.errorbar( d,rho, yerr=err)
#plt.plot(d, (10**A_min)*((d/8)**n_min))
#plt.plot(d_1_min, eje_Y_1_min, 'r' )
#plt.plot(d_2_min, eje_Y_2_min, 'r' )
#ax_0.set_xlim(10, 250)
#ax_0.set_ylim(10**-6, 1)

plt.show()
plt.close()


# using "curve_fit"

# assuming the r_break found before
def log_rho_broken(R, A2, n1, n2):
	mask = (R <= r_break)
	invmask = np.invert(mask)
    	c = 8
	rho = np.zeros(np.shape(R))

	rho[invmask] = A2 + n2*(R[invmask] / c)
	rho[mask] = A2 + (n2-n1)*np.log10(r_break/c) + n1*np.log10(R[mask]/c)

	return rho

popt,pcov = curve_fit(log_rho_broken, d, Y, p0=[2,-3,-5.4]) 
print r_break, popt[0] + (popt[2]-popt[1])*np.log10(r_break/c), popt
print pcov

# with no assumptions
def log_rho_broken2(R, R_BREAK, A2, n1, n2):
	mask = (R <= R_BREAK)
	invmask = np.invert(mask)
    	c = 8
	rho = np.zeros(np.shape(R))

	rho[invmask] = A2 + n2*(R[invmask] / c)
	rho[mask] = A2 + (n2-n1)*np.log10(R_BREAK/c) + n1*np.log10(R[mask]/c)

	return rho

popt,pcov = curve_fit(log_rho_broken2, d, Y, p0=[25, 2,-3,-5.4]) 
print popt
print pcov


