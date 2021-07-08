#%%
import sys
import glob
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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
    pathhead = 'DATA'
    itnum = '5'



#%% Read in input spectrum
df = pd.read_table(pathhead + '/InputSpectrum_1D.dat', header=0)
inSpectrum = df.values


#%%
plt.clf()
plt.figure(figsize=(8, 6))
plt.plot(inSpectrum[:,0], inSpectrum[:,3])
plt.title(r'Input Intensity Spectrum')
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.xlabel(r'$\omega$')
plt.grid(which='major')
plt.tight_layout()
plt.savefig(pathhead + '/figs/inSpectrum.png')
#plt.show()
######################################

#%% Read in spectrum of transmitted pulse
df = pd.read_table(pathhead + '/Spectrum_iteration_' + itnum + '_Transmitted.dat', header=0)
eSpectrumT = df.values


#%%
plt.clf()
plt.figure(figsize=(8, 6))
plt.plot(eSpectrumT[:,0], eSpectrumT[:,3])
plt.title(r'Transmitted pulse spectrum')
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.xlabel(r'$\omega$')
plt.ylabel(r'$|E(\omega)|$ [V$\cdot$s/m]')
plt.grid(which='major')
#plt.show()
plt.tight_layout()
plt.savefig(pathhead + '/figs/Ew_transmitted.png')



#%% Read in spectrum of reflected pulse
df = pd.read_table(pathhead + '/Spectrum_iteration_' + itnum + '_Reflected.dat', header=0)
eSpectrumR = df.values


#%%
plt.clf()
plt.figure(figsize=(8, 6))
plt.plot(eSpectrumR[:,0], eSpectrumR[:,3])
plt.title(r'Reflected pulse spectrum')
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.xlabel(r'$\omega$')
plt.ylabel(r'$|E(\omega)|$ [V$\cdot$s/m]')
plt.grid(which='major')
#plt.show()
plt.tight_layout()
plt.savefig(pathhead + '/figs/Ew_reflected.png')

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
plt.ticklabel_format(axis='x', style='sci', scilimits=(15,15))
plt.legend(['Transmitted', 'Reflected'])
plt.xlabel(r'$\omega$')
plt.ylabel(r'$|E(\omega)|$ [V$\cdot$s/m]')
plt.grid(which='major')
#plt.show()
plt.margins(x=0)
plt.tight_layout()
plt.savefig(pathhead + '/figs/Ew_both.png')
####################################################################
####################################################################

#%% Find point monitor files
pmon_li = [] # List of point monitor files
zmon_li = [] # List of point monitor locations

# Use glob to find all filenames
for name in glob.glob(pathhead + '/PointMon_*'):
    pmon_li.append(name)
    nameli = name.split('_')
    #print(nameli)
    zloc = nameli[-1].split('nm')
    zloc = zloc[0]
    #print('Location of point monitor: {} nm'.format(zloc))
    zmon_li.append(zloc)

print(pmon_li)
print(zmon_li)

#%%
for pointmon in pmon_li:
    zm = zmon_li.pop(0)
    zm_f = float(zm)
    print(zm)
    df = pd.read_table(pointmon, 
        header=None, skiprows=[0,1],
        names=['t [s]', 'Re(Ept) [V/m]', 'Im(Ept) [V/m]', 'Re(Emt) [V/m]', 'Im(Emt) [V/m]', 'Ne', 'Je']
        )

    # Plot the plasma density
    plt.clf()
    plt.figure(figsize=(8, 6))
    plt.semilogy(df['t [s]'], df['Ne'])
    plt.title(r'Electron density at $z={:8.2f}$ $\mu$m'.format(zm_f*1e-3))
    plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
    plt.xlabel(r'$t$ [s]')
    plt.grid(which='major')
    plt.margins(x=0)
    plt.tight_layout()
    plt.savefig(pathhead + '/figs/Ne_' + zm + 'nm.png')
    #plt.show()

    # Plot the forward-propagating pulse
    plt.clf()
    plt.figure(figsize=(8, 6))
    plt.plot(df['t [s]'], df['Re(Ept) [V/m]'])
    plt.title(r'Forward-propagating field at $z={:8.2f}$ $\mu$m'.format(zm_f*1e-3))
    plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
    plt.xlabel(r'$t$ [s]')
    plt.ylabel(r'Re$[E_+(t)]$ [V/m]')
    plt.grid(which='major')
    plt.margins(x=0)
    plt.tight_layout()
    plt.savefig(pathhead + '/figs/Et_forward_' + zm + '.png')
    #plt.show()

    # Plot the backward-propagating pulse
    plt.clf()
    plt.figure(figsize=(8, 6))
    plt.plot(df['t [s]'], df['Re(Emt) [V/m]'])
    plt.title(r'Backward-propagating field at $z={:8.2f}$ $\mu$m'.format(zm_f*1e-3))
    plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
    plt.xlabel(r'$t$ [s]')
    plt.ylabel(r'Re$[E_-(t)]$ [V/m]')
    plt.grid(which='major')
    plt.margins(x=0)
    plt.tight_layout()
    plt.savefig(pathhead + '/figs/Et_backward_' + zm + '.png')
    #plt.show()

    # Find the total field
    eP = df['Re(Ept) [V/m]'].values + 1.0j * df['Im(Ept) [V/m]'].values
    eM = df['Re(Emt) [V/m]'].values + 1.0j * df['Im(Emt) [V/m]'].values

    eTotal = eP + eM

    # Plot the total field
    plt.clf()
    plt.figure(figsize=(8, 6))
    plt.plot(df['t [s]'], np.real(eTotal))
    plt.title(r'Total field at $z={:8.2f}$ $\mu$m'.format(zm_f*1e-3))
    plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
    plt.xlabel(r'$t$ [s]')
    plt.ylabel(r'Re$[E_+(t) + E_-(t)]$ [V/m]')
    plt.grid(which='major')
    plt.margins(x=0)
    plt.tight_layout()
    plt.savefig(pathhead + '/figs/Et_total_' + zm + '.png')
    #plt.show()

    # Find spectra
    eOmP = np.fft.fft(eP)
    eOmM = np.fft.fft(eM)
    eOmTotal = np.fft.fft(eTotal)

    halflen = len(eOmTotal)
    dom_t = float(df.values[-1:,0] - df.values[0,0])
    omeg = [math.pi * i / dom_t for i in range(halflen)]

    plt.clf()
    plt.figure(figsize=(8, 6))
    plt.semilogy(eSpectrumT[:halflen,0], np.abs(eOmTotal[:]))
    plt.title(r'Total spectrum at $z={:8.2f}$ $\mu$m'.format(zm_f*1e-3))
    plt.ticklabel_format(axis='x', style='sci', scilimits=(15,15))
    plt.xlabel(r'$\omega$')
    plt.ylabel(r'$|E(\omega)|$ [V$\cdot$s/m]')
    plt.grid(which='major')
    plt.xlim([0.0, 32e15])
    plt.margins(x=0)
    plt.tight_layout()
    plt.savefig(pathhead + '/figs/Ew_total_' + zm + '.png')
    #plt.show()


#%% Read in forward-propagating pulse
df = pd.read_table(pathhead + '/Efield_iteration_' + itnum + '_Transmitted.dat', header=0)
eFieldT = df.values

#%%
plt.clf()
plt.figure(figsize=(8, 6))
plt.plot(eFieldT[:,0], eFieldT[:,1])
plt.title(r'Transmitted pulse')
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.xlabel(r'$t$ [s]')
plt.ylabel(r'$E(t)$ [V/m]')
plt.grid(which='major')
#plt.show()
plt.tight_layout()
plt.savefig(pathhead + '/figs/Et_transmitted.png')

#%%
idx1 = int(eFieldT.shape[0] / 2 - eFieldT.shape[0] / 20)
idx2 = int(eFieldT.shape[0] / 2 + eFieldT.shape[0] / 20)
plt.clf()
plt.figure(figsize=(8, 6))
plt.plot(eFieldT[idx1:idx2,0], eFieldT[idx1:idx2,1])
plt.title(r'Forward-propagating pulse')
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.xlabel(r'$t$ [s]')
plt.ylabel(r'$E(t)$ [V/m]')
plt.grid(which='major')
#plt.show()
plt.tight_layout()
plt.savefig(pathhead + '/figs/Et_transmitted_zoomed.png')
######################################


#%% Read in backward-propagating pulse
df = pd.read_table(pathhead + '/Efield_iteration_' + itnum + '_Reflected.dat', header=0)
eFieldR = df.values


#%%
plt.clf()
plt.figure(figsize=(8, 6))
plt.plot(eFieldR[:,0], eFieldR[:,1])
plt.title(r'Backward-propagating pulse')
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.xlabel(r'$t$ [s]')
plt.ylabel(r'$E(t)$ [V/m]')
plt.grid(which='major')
plt.tight_layout()
#plt.show()
plt.savefig(pathhead + '/figs/Et_reflected.png')

#%%
idx1 = int(eFieldR.shape[0] / 2 - eFieldR.shape[0] / 20)
idx2 = int(eFieldR.shape[0] / 2 + eFieldR.shape[0] / 20)
plt.clf()
plt.figure(figsize=(8, 6))
plt.plot(eFieldR[idx1:idx2,0], eFieldR[idx1:idx2,1])
plt.title(r'Backward-propagating pulse')
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.xlabel(r'$t$ [s]')
plt.ylabel(r'$E(t)$ [V/m]')
plt.grid(which='major')
#plt.show()
plt.tight_layout()
plt.savefig(pathhead + '/figs/Et_reflected_zoomed.png')
######################################



