import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# simple power law profile
def profile1(r, rho0, n0):
    
    r0 = 10
    return rho0 * (r / r0)**n0

# broken power law profile, forced continuity
def profile2(r, rho0, n0, r1, n1):
    
    mask = (r >= r1)
    invmask = np.invert(mask)
    
    r0 = 10
    rho = np.zeros(np.shape(r))
    rho[invmask] = rho0 * (r[invmask] / r0)**n0
    rho[mask] = rho0 / r0**n0 * r[mask]**n1 * r1**(n0 - n1)

    return rho

# read input data
(d, rho, e_rho) = np.loadtxt("densidad_bins.dat").transpose()

# start plotting
fig, ax = plt.subplots()
ax.errorbar(d, rho, yerr = e_rho, label = "HiTS RR Lyrae")

# try to find values manually to define initial iteration and then comment
rs = np.linspace(12, 250, 1000)
p0s = [(0.07, -2, 25, -4), (0.8, -3, 50, -6), (0.5, -1, 17, -4), (0.7, -1, 30, -7)]
for i in p0s:
    ax.plot(rs, profile2(rs, i[0], i[1], i[2], i[3]), 'r', label = "Initial solution")
    (popt, pcov) = curve_fit(profile2, d, rho, sigma = e_rho, p0 = i)
    ax.plot(rs, profile2(rs, popt[0], popt[1], popt[2], popt[3]), 'k', label = "Broken power law fit")
   
    print "Broken power law best-fitting parameters (rho0(r0=10), n0, r1, n1):", popt

# fit simple power law
(popt, pcov) = curve_fit(profile1, d, rho, sigma = e_rho, p0 = (0.3, -4))
ax.plot(rs, profile1(rs, popt[0], popt[1]), 'gray', label = "Simple power law fit")

print "Simple power law best-fitting parameters (rho0(r0=10), n0):", popt

# save figure
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim(10, 300)
ax.set_xlabel("Distance [kpc]")
ax.set_ylabel("Density [kpc^-3]")
plt.legend(loc = 1, fontsize = 9, framealpha = 0.5)
plt.savefig("bestfit.png")

