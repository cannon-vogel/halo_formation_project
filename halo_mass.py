from colossus.cosmology import cosmology
from colossus.lss import mass_function
import numpy as np
import matplotlib.pyplot as plt

# set cosmology, z vals, and mass vals for plots
cosmo = cosmology.setCosmology('planck18',persistence='r')
z = [5.0,10.0]
lowLim = 1
highLim = 13
M = 10**np.arange(lowLim,highLim,0.01)

# set up plot
plt.figure()
plt.title("Halo Mass Function")
plt.xlabel(r"$M_{200m}$")
plt.ylabel(r"$\frac{dn}{d\ln (M)}$")

# generate halo mass function for lambda CDM
for i in range(len(z)):
    mfunc1 = mass_function.massFunction(M, z[i], mdef = '200m', model = 'bocquet16', q_out = 'dndlnM')
    plt.plot(M, mfunc1, '-', label = 'z (CDM) = %.1f' % (z[i]))
plt.legend()
plt.savefig("./plots/z5_z10_haloMassFunc.png")


# now use truncated power spectrum to generate halo mass function
# to give massFunction() the P(k) we want to use, we use ps_args to pass the path of
# the power spectrum to the matterPowerSpectrum() function
for i in range(len(z)):
    mfunc2 = mass_function.massFunction(M, z[i], mdef = '200m', model = 'bocquet16', q_out = 'dndlnM', ps_args = {'model' : 'mytable', 'path' : '/Users/alexstewart/Documents/GitHub/halo_formation_project/axino2keV06_pk_log.txt'})
    plt.plot(M, mfunc2, '-', label = 'z (WDM) = %.1f' % (z[i]))
plt.loglog()
#plt.xlim(10**(lowLim),10**(highLim))
#plt.ylim(1E-7,1E5)
plt.legend()
plt.show()

for i in range(0,len(mfunc2)):
    print(mfunc2[i])


'''
# testing stuff to see if I can get the power spectrum read from the files bc above code isn't working
from colossus.utils import utilities

file_path = '/Users/alexstewart/Documents/GitHub/halo_formation_project/axino2keV06_pk_log.txt'
ps_args_table = dict(model = 'mytable', path = file_path)

k_table = 10**np.linspace(-4.0, 2.2, 200)
P_table = cosmo.matterPowerSpectrum(k_table, **ps_args_table)

Pk = cosmo.matterPowerSpectrum(k_table)

plt.figure()
plt.loglog()
plt.xlabel('k(h / Mpc)')
plt.ylabel('P(k)')
plt.plot(k_table, P_table, '-', label = 'Axino 2keV WDM')
plt.plot(k_table, Pk, '-', label='Planck 18')
plt.legend()
plt.show()
'''

'''
# Used this code to compute the log(k) log(P(k)) stuff
with open('/Users/alexstewart/Documents/GitHub/halo_formation_project/axino2keV06_pk.dat') as f:
    kData = [line.split()[0] for line in f]

with open('/Users/alexstewart/Documents/GitHub/halo_formation_project/axino2keV06_pk.dat') as g:
    PData = [line.split()[1] for line in g]
    
kData = kData[4:]
PData = PData[4:]

kDataArr = np.log10(np.array(kData).astype(float))
PDataArr = np.log10(np.array(PData).astype(float))

dataArr = np.vstack((kDataArr, PDataArr)).T

np.savetxt("axino2keV06_pk_log.txt",dataArr,delimiter='\t',newline='\n')
'''
