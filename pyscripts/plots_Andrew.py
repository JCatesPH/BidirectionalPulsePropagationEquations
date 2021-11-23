#%%
import sys
import glob
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#%matplotlib widget

#%%
if len(sys.argv) > 1:
    pathhead = str(sys.argv[1])
else:
    pathhead = 'DATA'
if len(sys.argv) > 2:
    itnum = str(sys.argv[2])
else:
    itnum = '5'

# Set values if run in interactive mode (VSCode)
if hasattr(sys, 'ps1'):
    print("Interactive mode detected..")
    #pathhead = '../DATA_Drude_10umOmeg0_FixedPointComp/5umL_NoFixedPoint_1e-6Noise'
    pathhead = '../DATA1'
    itnum = '12'


CLIGHT = 299792458
#plt.ion()
######################################
#%% Compute and print the preformed plasma parameters 
lamb0 = 10.0e-6
omeg0 = 2 * np.pi * CLIGHT / lamb0 

epsilon_0 = 8.85418782e-12
charge_e = 1.602176634e-19
mass_e = 9.10938356e-31

numE = 9.0e24
omegPlasma = np.sqrt(charge_e**2*numE/(epsilon_0*mass_e))

print("For {:} electrons, the plasma frequency is: {:.8e}".format(numE, omegPlasma))
print("Wavelength of {:.2e} [m] gives frequecy {:.2e}".format(lamb0, omeg0))
print("Approximate plasma skin depth: {:.2e}".format(CLIGHT/omegPlasma))
print("{}".format(omegPlasma/omeg0))

#%% Load in the simulation parameters
df = pd.read_table(pathhead + '/SimParameters.dat', 
    header=None, names=['Variable', 'Value', 'Description']
    )

# Remove trailing whitespace on variable names
df['Variable']=df['Variable'].str.strip()

#%% Parse out desired values (Intensity, Central freq., and pulse duration)
iIntensity = df.loc[df['Variable'] == 'I_0'].values[0,1]
omeg0 = df.loc[df['Variable'] == 'omega_0'].values[0,1]
taup = df.loc[df['Variable'] == 'tau'].values[0,1]
nAtoms = df.loc[df['Variable'] == 'num_atoms'].values[0,1]
sourceThickness = df.loc[df['Variable'] == 'LHSsourceLayerThickness'].values[0,1]
sampleThickness = df.loc[df['Variable'] == 'sampleLayerThickness'].values[0,1]
freqUpperCutoff = int(df.loc[df['Variable'] == 'freqUpperCutoff'].values[0,1])
freqLowerCutoff = int(df.loc[df['Variable'] == 'freqLowerCutoff'].values[0,1])

labstr = r'$I_0={:4.2f}$ [TW/cm$^2$], $\omega_0={:4.2f}$ [THz], $\tau_p={:4.2f}$ [fs]'.format(iIntensity*1e-16, omeg0*1e-12, taup*1e15)

######################################

#%% Load in the structure file
df = pd.read_table(pathhead + '/Structure.dat',
        index_col=False
    )
structure = df.values

#%% Plot structure indices
plt.figure(1, figsize=(8, 6))
plt.clf()
plt.plot(structure[:,0], structure[:,2], label=r'n_0')
plt.plot(structure[:,0], structure[:,3], label=r'n_2')
plt.title(r'Background Refractive Index')
plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
plt.xlabel(r'$z$ [m]')
plt.grid(which='major')
plt.legend()
plt.tight_layout()
plt.savefig(pathhead + '/figs/Structure.png')
#plt.show()

######################################


#%% Read in input pulse in time domain
df = pd.read_table(pathhead + '/InputEfield_1D.dat', header=0)
eFieldIn = df.values

lowerInd = np.count_nonzero(eFieldIn[:,0] < np.mean(eFieldIn[:,0]) - 3 * taup)
upperInd = np.count_nonzero(eFieldIn[:,0] < np.mean(eFieldIn[:,0]) + 3 * taup)

#%%
plt.figure(2, dpi=200)
plt.clf()
plt.plot(eFieldIn[lowerInd:upperInd,0], eFieldIn[lowerInd:upperInd,1], 'b-o')
plt.title(r'Incident pulse')
plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.xlabel(r'$t$ [s]')
plt.ylabel(r'$E(t)$ [V/m]')
plt.grid(which='major')
plt.margins(x=0)
plt.tight_layout()
plt.savefig(pathhead + '/figs/Input_Et.png')
#plt.show()


#%% Read in input spectrum
df = pd.read_table(pathhead + '/InputSpectrum_1D.dat', header=0)
inSpectrum = df.values


#%%
plt.figure(3, figsize=(8, 6))
plt.clf()
plt.semilogy(inSpectrum[:,0]/omeg0, inSpectrum[:,3]**2, 'r-', label=labstr)
plt.title(r'Input Intensity Spectrum')
#plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
plt.xlabel(r'$\omega/\omega_0$')
plt.ylabel(r'$|E(\omega)|^2$ [V$\cdot$s/m]')
plt.grid(which='major')
plt.legend()
plt.margins(x=0)
plt.minorticks_on()
plt.tight_layout()
plt.savefig(pathhead + '/figs/inSpectrum.png')
#plt.show()
######################################

#%%
""" df = pd.read_table(pathhead + '/kz_plasma.dat', header=0)
kz = df.values

fig = plt.figure(8, dpi=200)
plt.clf()
plt.plot(kz[1:,0]/(2*np.pi)*1e-12, kz[1:,1]*CLIGHT/kz[1:,0], 
    kz[1:,0]/(2*np.pi)*1e-12, kz[1:,2]*CLIGHT/kz[1:,0])
plt.legend(['Re[n]', 'Im[n]'])
plt.xlabel(r'$f$ [THz]')
plt.title('Refractive Index in Plasma')
plt.xlim([20, 100])
plt.ylim([0.0, 3])
plt.axvline(omegPlasma/(2*np.pi)*1e-12, color='k', linestyle='--')
#plt.margins(x=0)
plt.minorticks_on()
#plt.tight_layout()
plt.savefig(pathhead + '/figs/fig_nPlasma.png') """

######################################

#%% Read in spectrum of transmitted pulse
df = pd.read_table(pathhead + '/Spectrum_iteration_' + itnum + '_Transmitted.dat', header=0)
eSpectrumT = df.values

halfT = int(eSpectrumT.shape[0]/2)

dOm = abs(eSpectrumT[1,0] - eSpectrumT[0,0])
print('dOmega = {:.2e}'.format(dOm))

#%%
plt.figure(4, figsize=(8, 6))
plt.clf()
plt.semilogy(eSpectrumT[:halfT,0]/omeg0, eSpectrumT[:halfT,3]**2, label=labstr)
plt.title(r'Transmitted pulse spectrum')
#plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.xlabel(r'$\omega/\omega_0$')
plt.ylabel(r'$|E(\omega)|^2$ [V$\cdot$s/m]')
plt.grid(which='major')
plt.legend()
plt.margins(x=0)
plt.minorticks_on()
#plt.xlim([0.0, np.max(eSpectrumT[:,0])])
plt.tight_layout()
plt.savefig(pathhead + '/figs/EwT_iter{:}.png'.format(itnum))
#plt.show()


#%% Read in spectrum of reflected pulse
df = pd.read_table(pathhead + '/Spectrum_iteration_' + itnum + '_Reflected.dat', header=0)
eSpectrumR = df.values

halfR = int(eSpectrumR.shape[0]/2)

#%%
plt.figure(5, figsize=(8, 6))
plt.clf()
plt.semilogy(eSpectrumR[:halfR,0]/omeg0, eSpectrumR[:halfR,3]**2, label=labstr)
plt.title(r'Reflected pulse spectrum')
#plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.xlabel(r'$\omega/\omega_0$')
plt.ylabel(r'$|E(\omega)|^2$ [V$\cdot$s/m]')
plt.grid(which='major')
plt.legend()
plt.margins(x=0)
plt.minorticks_on()
#plt.xlim([0.0, np.max(eSpectrumR[:,0])])
plt.tight_layout()
plt.savefig(pathhead + '/figs/EwR_iter{:}.png'.format(itnum))
plt.show()
#plt.close('all')
######################################

#%%
#freqUpperCutoff = 750
plt.figure(5, figsize=(8, 6))
plt.clf()
plt.semilogy(
    eSpectrumT[freqLowerCutoff:freqUpperCutoff,0]/omeg0, eSpectrumT[freqLowerCutoff:freqUpperCutoff,3]**2,
    eSpectrumR[freqLowerCutoff:freqUpperCutoff,0]/omeg0, eSpectrumR[freqLowerCutoff:freqUpperCutoff,3]**2
    )
plt.title(r'Both pulse spectra')
plt.legend(['Transmitted', 'Reflected'])
#plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.xlabel(r'$\omega/\omega_0$')
plt.ylabel(r'$|E(\omega)|^2$')
plt.grid(which='major')
plt.margins(x=0)
plt.minorticks_on()
#plt.xlim([0.0, np.max(eSpectrumR[:,0])])
plt.tight_layout()
plt.savefig(pathhead + '/figs/Ew_both_iter{:}.png'.format(itnum))
plt.show()
plt.close('all')
######################################

#%% Overlay transmitted and reflected spectra
fig = plt.figure(6, figsize=(6.4, 4.8), dpi=200)
firstIndex = freqLowerCutoff
plt.clf()
plt.plot(
    inSpectrum[firstIndex:freqUpperCutoff,0], inSpectrum[firstIndex:freqUpperCutoff,3], '-b',
    eSpectrumT[firstIndex:freqUpperCutoff,0], eSpectrumT[firstIndex:freqUpperCutoff,3], '-g', 
    eSpectrumR[firstIndex:freqUpperCutoff,0], eSpectrumR[firstIndex:freqUpperCutoff,3], '-r'
    )
plt.title(r'Transmitted and reflected amplitude spectra')
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.legend(['Incident', 'Transmitted', 'Reflected'])
plt.xlabel(r'$\omega$ [s$^{-1}$]')
plt.ylabel(r'$|E(f)|$ [V$\cdot$s/m]')
plt.grid(which='major')
#plt.xlim([20, 100])
plt.xlim([omeg0 - 7/taup, omeg0 + 7/taup])
plt.text(0.25, 0.75, 
    r'$L=${:.2f} $\mu$m'.format(sampleThickness*1e6), 
    horizontalalignment='center', verticalalignment='center', transform=fig.transFigure,
    fontsize=14)
plt.axvline(omegPlasma, color='k', linestyle='--')
plt.margins(x=0)
plt.minorticks_on()
plt.tight_layout()
plt.savefig(pathhead + '/figs/Ew_Fundamental_iter{:}.png'.format(itnum))
plt.show()
####################################################################
####################################################################

# %%
