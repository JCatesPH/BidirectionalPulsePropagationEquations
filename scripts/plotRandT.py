#%%
import sys
import glob
import math
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

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
    pathhead = '../DATA_DBR'
    itnum = '11'


CLIGHT = 299792458
EPS0 = 8.85418782e-12
CHARGE_E = 1.602176634e-19
MASS_E = 9.10938356e-31

intensityFactor = EPS0*CLIGHT/2 * 1e-4 # Second factor changes it to W/cm^2

mpl.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 16
plt.rcParams['axes.linewidth'] = 2

######################################
#%%
def plotSpectrum(freqArr, spectrumArr, labelArr, filePath, titleStr=None, xlabelStr=r'$\omega/\omega_0$'):
    fig, axs = plt.subplots(figsize=(6.47, 4), dpi=400)

    for freqs, spectrum, label in zip(freqArr, spectrumArr, labelArr):
        axs.semilogy(freqs, spectrum, label=label)

    #ax2.semilogy(df['t [s]'], neVec)
    if titleStr is not None:
        axs.set_title(titleStr)
    #axs.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
    axs.set_xlabel(xlabelStr)
    #axs.set_ylabel(r'Re[$E(t)$] [V/m]')
    #axs.set_xlim([0.0, np.max(freqArr)])
    axs.grid(which='major')
    axs.margins(x=0)
    axs.legend()

    fig.savefig(filePath)
    return fig, axs

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
Z0 = df.loc[df['Variable'] == 'Znaught'].values[0,1]
freqUpperCutoff = int(df.loc[df['Variable'] == 'freqUpperCutoff'].values[0,1])

labstr = r'$I_0={:4.2f}$ [TW/cm$^2$], $\omega_0={:4.2f}$ [THz], $\tau_p={:4.2f}$ [fs]'.format(iIntensity*1e-16, omeg0*1e-12, taup*1e15)

######################################

#%% Plot the structure
df = pd.read_table(pathhead + '/Structure.dat', header=0, index_col=False)
structureDat = df.values

plt.subplots(figsize=(6.47, 4), dpi=400)
plt.plot(structureDat[:,0], structureDat[:,2])
plt.title('Structure Indices of Refraction')
plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
plt.xlabel(r'$z$ [m]')
plt.ylabel(r'$n$')
plt.margins(x=0)
plt.savefig(pathhead + '/figs/Structure.png')
plt.show()

######################################


#%% Read in input spectrum
df = pd.read_table(pathhead + '/InputSpectrum_1D.dat', header=0)
inSpectrum = df.values

omeg = inSpectrum[:,0]
Ez = inSpectrum[:,3]
k0arr = omeg / 3e8

#%% Plot input spectrum
plt.clf()
plt.figure(figsize=(6.47, 4), dpi=400)
plt.semilogy(omeg/omeg0, Ez**2)
plt.title(r'Input Amplitude Spectrum')
#plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
plt.xlabel(r'$\omega/\omega_0$')
plt.grid(which='major')
plt.legend()
#plt.xlim([0.0, np.max(inSpectrum[:,0])])
plt.margins(x=0)
plt.tight_layout()
plt.savefig(pathhead + '/figs/inSpectrum.png')
#plt.show()

#%% Read in the transmitted and reflected spectra
df = pd.read_table(pathhead + '/Spectrum_iteration_' + itnum + '_Transmitted.dat')
omeg = df.values[:,0]
eOmT = df.values[:,1] + 1j * df.values[:,2]

df = pd.read_table(pathhead + '/Spectrum_iteration_' + itnum + '_Reflected.dat')
eOmR = df.values[:,1] + 1j * df.values[:,2]

halflen = int(len(omeg) / 2)

#%%
#plotSpectrum([omeg/omeg0], [np.abs(eOmT)**2], [''], pathhead + '/figs/TransmittedSpectrum.png')
#plotSpectrum([omeg/omeg0], [np.abs(eOmR)**2], [''], pathhead + '/figs/ReflectedSpectrum.png')
plotSpectrum([omeg[:freqUpperCutoff]/omeg0, omeg[:freqUpperCutoff]/omeg0], 
    [np.abs(eOmT[:freqUpperCutoff])**2, np.abs(eOmR[:freqUpperCutoff])**2], 
    ['Transmitted', 'Reflected'], 
    filePath=pathhead + '/figs/BothSpectra.png',
    titleStr='Transmitted and Reflected Spectra')
#plt.show()

#%%
sighat = np.sqrt(4*np.log(2)/taup**2)
lowerInd = np.count_nonzero(omeg[:halflen] < omeg0-15*sighat)
upperInd = np.count_nonzero(omeg[:halflen] < omeg0+15*sighat)
#lowerInd = 1
#upperInd = 1904
tcoef = (np.abs(eOmT[lowerInd:upperInd]) / Ez[lowerInd:upperInd])**2

plt.figure(figsize=(6.47, 4), dpi=400)
plt.plot(omeg[lowerInd:upperInd]/omeg0, tcoef)
plt.xlabel(r'$\omega/\omega_0$')
plt.ylabel(r'$T$')
plt.grid(which='major')
plt.legend()
#plt.xlim([0.8*omeg0, 1.2*omeg0])
#plt.ylim([0, 2])
plt.margins(x=0)
plt.title('Transmission Coefficient')
plt.tight_layout()
plt.savefig(pathhead + '/figs/TransmissionCoef.png')
#plt.show()

rcoef = (np.abs(eOmR[lowerInd:upperInd]) / Ez[lowerInd:upperInd])**2

plt.figure(figsize=(6.47, 4), dpi=400)
plt.plot(omeg[lowerInd:upperInd]/omeg0, rcoef)
plt.xlabel(r'$\omega/\omega_0$')
plt.ylabel(r'$R$')
plt.grid(which='major')
plt.legend()
plt.xlim([0.5, 1.5])
#plt.ylim([0, 2])
plt.margins(x=0)
plt.title('Reflection Coefficient')
plt.tight_layout()
plt.savefig(pathhead + '/figs/ReflectionCoef.png')
plt.show()

#%% Print mean coefficients
print('avg tr = {:}'.format(np.nanmean(np.where(np.isfinite(tcoef), tcoef, 0))))
print('avg rf = {:}'.format(np.nanmean(np.where(np.isfinite(rcoef), rcoef, 0))))

# %% Using Griffiths method naively
n1 = 1.0
n2 = 2.0
n3 = 1.0

R12 = np.abs((n1-n2) / (n1 + n2))**2
R23 = np.abs((n2-n3) / (n2 + n3))**2
T12 = 1 - R12
T23 = 1 - R23  #4*n2*n3/(n2+n3)**2

print("n_d = {:.2f}".format(n2))
print("R12 = {:.4f}".format(R12))
print("T12 = {:.4f}".format(T12))
print("R23 = {:.4f}".format(R23))
print("T23 = {:.4f}".format(T23))

#%%
E0R = 1
E1R = T12 * E0R
E2R = T23 * E1R
E1L = R23 * E1R
E0L = R12 * E0R + T12 * E1L

print('E0R = {:}'.format(E0R))
print('E0L = {:}'.format(E0L))
print('E1R = {:}'.format(E1R))
print('E1L = {:}'.format(E1L))
print('E2R = {:}'.format(E2R))
print('E0L + E2R = {:}'.format(E0L+E2R))

#%%
for i in range(6):
    E1R = T12 * E0R + R12 * E1L
    E2R = T23 * E1R
    E1L = R23 * E1R
    E0L = R12 * E0R + T12 * E1L

print('E0R = {:}'.format(E0R))
print('E0L = {:}'.format(E0L))
print('E1R = {:}'.format(E1R))
print('E1L = {:}'.format(E1L))
print('E2R = {:}'.format(E2R))
print('E0L + E2R = {:}'.format(E0L+E2R))

# %%
