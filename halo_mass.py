from colossus.cosmology import cosmology
from colossus.lss import mass_function
import numpy as np
import matplotlib.pyplot as plt

# set cosmology, z vals, and mass vals for plots
cosmo = cosmology.setCosmology('planck18',persistence='r')
z = [5.0,10.0]
lowLim = 3.3
highLim = 13
M = 10**np.arange(lowLim,highLim,0.01)

# set up plot
plt.figure()
'''
plt.title("Halo Mass Function")
plt.xlabel(r"$M_{200m}$")
plt.ylabel(r"$\frac{dn}{d\ln (M)}$")
'''

# generate halo mass function for lambda CDM
mfunc_CDM_z5 = mass_function.massFunction(M, 5, mdef = '200m', model = 'despali16', q_out = 'dndlnM')
mfunc_CDM_z10 = mass_function.massFunction(M, 10, mdef = '200m', model = 'despali16', q_out = 'dndlnM')

# now use truncated power spectrum to generate halo mass function
# to give massFunction() the P(k) we want to use, we use ps_args to pass the path of
# the power spectrum to the matterPowerSpectrum() function
mfunc_WDM_z5 = mass_function.massFunction(M, 5, mdef = '200m', model = 'despali16', q_out = 'dndlnM', ps_args = {'model' : 'mytable', 'path' : '/Users/alexstewart/Documents/GitHub/halo_formation_project/axino20keV_2k03_pk_log.txt'})
mfunc_WDM_z10 = mass_function.massFunction(M, 10, mdef = '200m', model = 'despali16', q_out = 'dndlnM', ps_args = {'model' : 'mytable', 'path' : '/Users/alexstewart/Documents/GitHub/halo_formation_project/axino20keV_2k03_pk_log.txt'})

mfunc_HDM_z5 = mass_function.massFunction(M, 5, mdef = '200m', model = 'despali16', q_out = 'dndlnM', ps_args = {'model' : 'mytableHDM', 'path' : '/Users/alexstewart/Documents/GitHub/halo_formation_project/hot_dark_matter_pk_log.txt'})

mfunc_HDM_z10 = mass_function.massFunction(M, 10, mdef = '200m', model = 'despali16', q_out = 'dndlnM', ps_args = {'model' : 'mytableHDM', 'path' : '/Users/alexstewart/Documents/GitHub/halo_formation_project/hot_dark_matter_pk_log.txt'})

#plt.plot(M, mfunc_WDM_z5, '-', label = 'WDM')
#plt.plot(M, mfunc_CDM_z5, '-', label = 'CDM')
plt.xlim(10**(lowLim),10**(highLim))
plt.ylim(1E-7,1E7)
plt.xlabel(r"$M_{200m}$")
plt.ylabel(r"$\frac{dn}{d\ln (M)}$")
plt.title("Halo Mass Functions (z = 10)")
plt.plot(M, mfunc_CDM_z10, '-', label = 'CDM')
plt.plot(M, mfunc_HDM_z10, '-', label = 'HDM')
plt.legend()
plt.loglog()
plt.show()
plt.clf()
plt.subplot(121)
plt.plot(M, mfunc_CDM_z10, '-', label = 'CDM')
plt.plot(M, mfunc_WDM_z10, '-', label = 'WDM')
plt.loglog()
plt.xlim(10**(lowLim),10**(highLim))
plt.ylim(1E-7,1E7)
plt.xlabel(r"$M_{200m}$")
plt.ylabel(r"$\frac{dn}{d\ln (M)}$")
plt.title("Halo Mass Functions (z = 10)")
plt.legend()
plt.subplot(122)
plt.plot(M, mfunc_CDM_z5, '-', label = 'CDM')
plt.plot(M, mfunc_WDM_z5, '-', label = 'WDM')
plt.loglog()
plt.xlim(10**(lowLim),10**(highLim))
plt.ylim(1E-7,1E7)
plt.xlabel(r"$M_{200m}$")
plt.ylabel(r"$\frac{dn}{d\ln (M)}$")
plt.title("Halo Mass Functions (z = 5)")
plt.legend()
plt.show()
plt.clf()

plt.xlim(10**(lowLim),10**(highLim))
plt.ylim(1E-7,1E7)
plt.xlabel(r"$M_{200m}$")
plt.ylabel(r"$\frac{dn}{d\ln (M)}$")
plt.title("Halo Mass Functions (z = 10)")
plt.plot(M, mfunc_CDM_z10, '-', label = 'CDM')
plt.plot(M, mfunc_HDM_z10, '-', label = 'HDM')
plt.legend()
plt.show()


# testing stuff to see if I can get the power spectrum read from the files bc above code isn't working

from colossus.utils import utilities

file_path = '/Users/alexstewart/Documents/GitHub/halo_formation_project/axino20keV_2k03_pk_log.txt'
ps_args_table = dict(model = 'mytable', path = file_path)

k_table = 10**np.linspace(-4.9, 3, 500)
P_table = cosmo.matterPowerSpectrum(k_table, **ps_args_table)

Pk = cosmo.matterPowerSpectrum(k_table)

file_path2 = '/Users/alexstewart/Documents/GitHub/halo_formation_project/hot_dark_matter_pk_log.txt'
ps_args_table_HDM = dict(model = 'mytableHDM', path = file_path2)
PkHDM = cosmo.matterPowerSpectrum(k_table, **ps_args_table_HDM)

plt.figure()
plt.loglog()
plt.title("Power Spectra for Plank 18 Data and WDM")
plt.xlabel('k(h / Mpc)')
plt.ylabel('P(k)')
plt.plot(k_table, P_table, '-', label = 'Axino 20keV WDM')
plt.plot(k_table, Pk, '-', label='Planck 18')
plt.plot(k_table, PkHDM, '-', label='HDM')
plt.legend()
plt.show()


'''
# Used this code to compute the log(k) log(P(k)) stuff
with open('/Users/alexstewart/Documents/GitHub/halo_formation_project/axino20keV_2k03_pk.dat') as f:
    kData = [line.split()[0] for line in f]

with open('/Users/alexstewart/Documents/GitHub/halo_formation_project/axino20keV_2k03_pk.dat') as g:
    PData = [line.split()[1] for line in g]
    
kData = kData[4:]
PData = PData[4:]

kDataArrTemp = np.array(kData).astype(float)
PDataArrTemp = np.array(PData).astype(float)

kDataArr = np.zeros(len(kDataArrTemp))
PDataArr = np.zeros(len(PDataArrTemp))

for i in range(0,len(kDataArrTemp)):
    kDataArr[i] = np.log10(kDataArrTemp[i])
    if kDataArrTemp[i] > 1:
        PDataArr[i] = np.log10(10**-10)
    else:
        PDataArr[i] = np.log10(PDataArrTemp[i])

dataArr = np.vstack((kDataArr, PDataArr)).T

np.savetxt("hot_dark_matter_pk_log.txt",dataArr,delimiter='\t',newline='\n')
'''
