from colossus.cosmology import cosmology
from colossus.lss import mass_function
import numpy as np
import matplotlib.pyplot as plt

# set cosmology, z vals, and mass vals for plots
cosmo = cosmology.setCosmology('planck18',persistence='r')
z = [5.0,10.0]
lowLim = 5.0
highLim = 13
M = 10**np.arange(lowLim,highLim,0.1)

# set up plot
plt.figure()
plt.title("Halo Mass Function")
plt.xlabel(r"$M_{200m}$")
plt.ylabel(r"$\frac{dn}{d\ln (M)}$")
plt.loglog()
plt.xlim(10**(lowLim),10**(highLim))
plt.ylim(1E-7,1E3)

# generate halo mass function for lambda CDM
for i in range(len(z)):
    mfunc = mass_function.massFunction(M, z[i], mdef = '200m', model = 'despali16', q_out = 'dndlnM')
    plt.plot(M, mfunc, '-', label = 'z (CDM) = %.1f' % (z[i]))
plt.legend()
plt.savefig("./mass_func_plots/z5_z10_haloMassFunc.png")

'''
# now use truncated power spectrum to generate halo mass function
# to give massFunction() the P(k) we want to use, we use ps_args to pass the path of
# the power spectrum to the matterPowerSpectrum() function
for i in range(len(z)):
    mfunc = mass_function.massFunction(M, z[i], mdef = '200m', model = 'despali16', q_out = 'dndlnM', ps_args = {'model' : 'eisenstein98', 'path' : '/Users/alexstewart/Documents/GitHub/halo_formation_project/axino2keV06_pk.dat'})
    plt.plot(M, mfunc, '-', label = 'z (WDM) = %.1f' % (z[i]))
plt.title("Halo Mass Function")
plt.xlabel(r"$M_{200m}$")
plt.ylabel(r"$\frac{dn}{d\ln (M)}$")
plt.loglog()
plt.xlim(10**(lowLim),10**(highLim))
plt.ylim(1E-7,1E3)
plt.legend()
plt.show()

# testing stuff to see if I can get the power spectrum read from the files bc above code isn't working
from colossus.utils import utilities

file_path = '/Users/alexstewart/Documents/GitHub/halo_formation_project/axino2keV06_pk.dat'
ps_args_table = dict(model = 'mytable', path = file_path)

k_table = 10**np.linspace(0.1, 2.2, 200)
print(k_table)
P_table = cosmo.matterPowerSpectrum(k_table, **ps_args_table)
print(P_table)

plt.figure()
plt.loglog()
plt.xlabel('k(h / Mpc)')
plt.ylabel('P(k)')
plt.plot(k_table, P_table, '-')
plt.show()

# just open the power spectrum and see if I can read it
import csv
with open('/Users/alexstewart/Documents/GitHub/halo_formation_project/axino2keV06_pk.dat') as f:
    kData = [line.split()[0] for line in f]

with open('/Users/alexstewart/Documents/GitHub/halo_formation_project/axino2keV06_pk.dat') as g:
    PData = [line.split()[1] for line in g]
    
print(kData[4:])
print(PData[4:])
'''
