#%%
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


itnum = str(sys.argv[1])

#%% Read in input spectrum
df = pd.read_table('DATA/InputSpectrum_1D.dat', header=0)
inSpectrum = df.values


#%%
plt.clf()
plt.figure(figsize=(8, 6))
plt.plot(inSpectrum[:,0], inSpectrum[:,3])
plt.title(r'Input Intensity Spectrum')
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.xlabel(r'$\omega$')
plt.grid(which='major')
#plt.show()
plt.tight_layout()
plt.savefig('figs/inSpectrum.png')
######################################


#%% Read in plasma density
df = pd.read_table('DATA/Ne_iter_5_Zpos_10630nm.dat', header=0)
nedat = df.values.astype(np.float)
if np.isnan(nedat).any() != True:
    nedat[:,1] = np.where(nedat[:,1] < 1e-8, 1e-8, nedat[:,1])


#%%
plt.clf()
plt.figure(figsize=(8, 6))
plt.semilogy(nedat[:,0], nedat[:,1])
plt.title(r'Electron density')
plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
plt.xlabel(r'$t$ [s]')
plt.grid(which='major')
#plt.show()
plt.margins(x=0)
plt.tight_layout()
plt.savefig('figs/Ne.png')

#%%
idx1 = int(nedat.shape[0] / 2 - nedat.shape[0] / 20)
idx2 = int(nedat.shape[0] / 2 + nedat.shape[0] / 20)
plt.clf()
plt.figure(figsize=(8, 6))
plt.semilogy(nedat[idx1:idx2,0], nedat[idx1:idx2,1])
plt.title(r'Electron density')
plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
plt.xlabel(r'$t$ [s]')
plt.grid(which='major')
#plt.show()
plt.margins(x=0)
plt.tight_layout()
plt.savefig('figs/Ne_zoomed.png')
######################################


#%% Read in forward-propagating pulse
df = pd.read_table('DATA/Efield_iteration_' + itnum + '_Transmitted.dat', header=0)
eFieldT = df.values


#%%
plt.clf()
plt.figure(figsize=(8, 6))
plt.plot(eFieldT[:,0], eFieldT[:,1])
plt.title(r'Transmitted pulse')
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.xlabel(r'$t$ [s]')
plt.grid(which='major')
#plt.show()
plt.tight_layout()
plt.savefig('figs/Et_transmitted.png')

#%%
idx1 = int(eFieldT.shape[0] / 2 - eFieldT.shape[0] / 20)
idx2 = int(eFieldT.shape[0] / 2 + eFieldT.shape[0] / 20)
plt.clf()
plt.figure(figsize=(8, 6))
plt.plot(eFieldT[idx1:idx2,0], eFieldT[idx1:idx2,1])
plt.title(r'Forward-propagating pulse')
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.xlabel(r'$t$ [s]')
plt.grid(which='major')
#plt.show()
plt.tight_layout()
plt.savefig('figs/Et_transmitted_zoomed.png')
######################################


#%% Read in backward-propagating pulse
df = pd.read_table('DATA/Efield_iteration_' + itnum + '_Reflected.dat', header=0)
eFieldR = df.values


#%%
plt.clf()
plt.figure(figsize=(8, 6))
plt.plot(eFieldR[:,0], eFieldR[:,1])
plt.title(r'Backward-propagating pulse')
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.xlabel(r'$t$ [s]')
plt.grid(which='major')
plt.tight_layout()
#plt.show()
plt.savefig('figs/Et_reflected.png')

#%%
idx1 = int(eFieldR.shape[0] / 2 - eFieldR.shape[0] / 20)
idx2 = int(eFieldR.shape[0] / 2 + eFieldR.shape[0] / 20)
plt.clf()
plt.figure(figsize=(8, 6))
plt.plot(eFieldR[idx1:idx2,0], eFieldR[idx1:idx2,1])
plt.title(r'Backward-propagating pulse')
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.xlabel(r'$t$ [s]')
plt.grid(which='major')
#plt.show()
plt.tight_layout()
plt.savefig('figs/Et_reflected_zoomed.png')
######################################


#%% Read in spectrum of transmitted pulse
df = pd.read_table('DATA/Spectrum_iteration_' + itnum + '_Transmitted.dat', header=0)
eSpectrumT = df.values


#%%
plt.clf()
plt.figure(figsize=(8, 6))
plt.plot(eSpectrumT[:,0], eSpectrumT[:,3])
plt.title(r'Transmitted pulse spectrum')
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.xlabel(r'$\omega$')
plt.grid(which='major')
#plt.show()
plt.tight_layout()
plt.savefig('figs/Ew_transmitted.png')



#%% Read in spectrum of reflected pulse
df = pd.read_table('DATA/Spectrum_iteration_' + itnum + '_Reflected.dat', header=0)
eSpectrumR = df.values


#%%
plt.clf()
plt.figure(figsize=(8, 6))
plt.plot(eSpectrumR[:,0], eSpectrumR[:,3])
plt.title(r'Reflected pulse spectrum')
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.xlabel(r'$\omega$')
plt.grid(which='major')
#plt.show()
plt.tight_layout()
plt.savefig('figs/Ew_reflected.png')

######################################


#%% Overlay transmitted and reflected spectra
halfR = int(eSpectrumR.shape[0]/2)
halfT = int(eSpectrumT.shape[0]/2)

plt.clf()
plt.figure(figsize=(8, 6))
plt.semilogy(eSpectrumT[:halfT,0], eSpectrumT[:halfT,3], '-b', 
    eSpectrumR[:halfR,0], eSpectrumR[:halfR,3], '-r',
    #inSpectrum[:halfR,0], inSpectrum[:halfR,3], '--k',
    )
plt.title(r'Transmitted and reflected pulse spectra')
plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
plt.legend(['Transmitted', 'Reflected'])
plt.xlabel(r'$\omega$')
plt.ylabel(r'$|E(\omega)|$')
plt.grid(which='major')
#plt.show()
plt.margins(x=0)
plt.tight_layout()
plt.savefig('figs/Ew_both.png')

# %%
