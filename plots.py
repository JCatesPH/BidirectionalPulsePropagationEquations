#%%
import sys
import glob
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#%matplotlib widget

if len(sys.argv) > 1:
    pathhead = str(sys.argv[1])
else:
    pathhead = 'DATA'
if len(sys.argv) > 2:
    itnum = str(sys.argv[2])
else:
    itnum = '4'

# Set values if run in interactive mode (VSCode)
if hasattr(sys, 'ps1'):
    print("Interactive mode detected..")
    pathhead = 'DATA'
    itnum = '5'

#plt.ion()
######################################


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
plt.title(r'Structure')
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

#%%
plt.figure(2, figsize=(8, 6))
plt.clf()
plt.plot(eFieldIn[:,0], eFieldIn[:,1], 'b-o')
plt.title(r'Incident pulse')
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.xlabel(r'$t$ [s]')
plt.ylabel(r'$E(t)$ [V/m]')
plt.grid(which='major')

plt.tight_layout()
plt.savefig(pathhead + '/figs/Input_Et.png')
#plt.show()


#%% Read in input spectrum
df = pd.read_table(pathhead + '/InputSpectrum_1D.dat', header=0)
inSpectrum = df.values

halfIn = int(inSpectrum.shape[0]/2)

#%%
plt.figure(3, figsize=(8, 6))
plt.clf()
plt.semilogy(inSpectrum[:halfIn,0]/omeg0, inSpectrum[:halfIn,3]**2, 'r-', label=labstr)
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
plt.savefig(pathhead + '/figs/Ew_transmitted.png')
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
plt.savefig(pathhead + '/figs/Ew_reflected.png')
#plt.show()

######################################


#%% Overlay transmitted and reflected spectra
plt.figure(6, figsize=(8, 6))
plt.clf()
plt.semilogy(eSpectrumT[:halfT,0]/omeg0, eSpectrumT[:halfT,3]**2, '-b', 
    eSpectrumR[:halfR,0]/omeg0, eSpectrumR[:halfR,3]**2, '-r',
    #inSpectrum[:halfR,0], inSpectrum[:halfR,3], '--k',
    )
plt.title(r'Transmitted and reflected pulse spectra')
#plt.ticklabel_format(axis='x', style='sci', scilimits=(15,15))
plt.legend(['Transmitted', 'Reflected'])
plt.xlabel(r'$\omega/\omega_0$')
plt.ylabel(r'$|E(\omega)|^2$ [V$\cdot$s/m]')
plt.grid(which='major')
plt.margins(x=0)
plt.minorticks_on()
plt.tight_layout()
plt.savefig(pathhead + '/figs/Ew_both.png')
#plt.show()

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
    zloc = str(float(zloc[0]) - sourceThickness*1e9)
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

    # Round the plasma density to make cleaner plots
    neVec = df['Ne'].values
    neVec = np.where(neVec < 1e-9, 1e-9, neVec)
    # Plot the plasma density
    plt.figure(7, figsize=(8, 6))
    plt.clf()
    plt.semilogy(df['t [s]'], neVec)
    plt.title(r'Electron density at $z={:8.2f}$ $\mu$m'.format(zm_f*1e-3))
    plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
    plt.xlabel(r'$t$ [s]')
    plt.grid(which='major')
    plt.xlim([np.mean(df['t [s]']) - 2*taup, np.mean(df['t [s]']) + 2*taup])
    plt.margins(x=0)
    plt.minorticks_on()
    plt.tight_layout()
    plt.savefig(pathhead + '/figs/Ne_' + zm + 'nm.png')
    #plt.show()

    # Plot the forward-propagating pulse
    plt.clf()
    plt.plot(df['t [s]'], df['Re(Ept) [V/m]'], label=labstr)
    plt.title(r'Forward-propagating field at $z={:8.2f}$ $\mu$m'.format(zm_f*1e-3))
    plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
    plt.xlabel(r'$t$ [s]')
    plt.ylabel(r'Re$[E_+(t)]$ [V/m]')
    plt.grid(which='major')
    plt.legend()
    plt.margins(x=0)
    plt.tight_layout()
    plt.savefig(pathhead + '/figs/Et_forward_' + zm + '.png')
    #plt.show()

    # Plot the backward-propagating pulse
    plt.clf()
    plt.plot(df['t [s]'], df['Re(Emt) [V/m]'], label=labstr)
    plt.title(r'Backward-propagating field at $z={:8.2f}$ $\mu$m'.format(zm_f*1e-3))
    plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
    plt.xlabel(r'$t$ [s]')
    plt.ylabel(r'Re$[E_-(t)]$ [V/m]')
    plt.grid(which='major')
    plt.legend()
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
    plt.plot(df['t [s]'], np.real(eTotal), label=labstr)
    plt.title(r'Total field at $z={:8.2f}$ $\mu$m'.format(zm_f*1e-3))
    plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
    plt.xlabel(r'$t$ [s]')
    plt.ylabel(r'Re$[E_+(t) + E_-(t)]$ [V/m]')
    plt.grid(which='major')
    plt.legend()
    plt.margins(x=0)
    plt.tight_layout()
    plt.savefig(pathhead + '/figs/Et_total_' + zm + '.png')
    #plt.show()

    # Find spectra
    eOmP = np.fft.fft(eP)
    eOmM = np.fft.fft(eM)
    eOmTotal = np.fft.fft(eTotal)

    halflen = int(len(eOmTotal)/2)
    dom_t = float(df.values[-1:,0] - df.values[0,0])
    omeg = [math.pi * i / dom_t for i in range(halflen)]

    # Plot forward-prop spectrum
    plt.clf()
    plt.semilogy(eSpectrumT[:halflen,0], np.abs(eOmP[:halflen])**2, label=labstr)
    plt.title(r'Forward-propagating spectrum at $z={:6.2f}$ $\mu$m'.format(zm_f*1e-3))
    plt.ticklabel_format(axis='x', style='sci', scilimits=(15,15))
    plt.xlabel(r'$\omega$')
    plt.ylabel(r'$|E(\omega)|^2$ [V$\cdot$s/m]')
    plt.grid(which='major')
    plt.minorticks_on()
    plt.legend()
    plt.xlim([0.0, np.max(omeg)])
    plt.margins(x=0)
    plt.tight_layout()
    plt.savefig(pathhead + '/figs/EwP_' + zm + '.png')
    #plt.show()

    # Plot forward spectrum in THz
    plt.clf()
    plt.semilogy(eSpectrumT[:halflen,0], np.abs(eOmP[:halflen])**2, 'b-o', label=labstr)
    plt.title(r'Total spectrum at $z={:8.2f}$ $\mu$m'.format(zm_f*1e-3))
    plt.ticklabel_format(axis='x', style='sci', scilimits=(12,12))
    plt.xlabel(r'$\omega$')
    plt.ylabel(r'$|E(\omega)|^2$ [V$\cdot$s/m]')
    plt.grid(which='major')
    plt.minorticks_on()
    plt.legend()
    plt.xlim([0.0, 200e12])
    #plt.ylim([1e-6, 1])
    plt.margins(x=0)
    plt.tight_layout()
    plt.savefig(pathhead + '/figs/Ew_fwd_THz_' + zm + '.png')

    # Plot backward-prop spectrum
    plt.clf()
    plt.semilogy(eSpectrumT[:halflen,0], np.abs(eOmM[:halflen])**2, label=labstr)
    plt.title(r'Backward-propagating spectrum at $z={:6.2f}$ $\mu$m'.format(zm_f*1e-3))
    plt.ticklabel_format(axis='x', style='sci', scilimits=(15,15))
    plt.xlabel(r'$\omega$')
    plt.ylabel(r'$|E(\omega)|^2$ [V$\cdot$s/m]')
    plt.grid(which='major')
    plt.minorticks_on()
    plt.legend()
    plt.xlim([0.0, np.max(omeg)])
    plt.margins(x=0)
    plt.tight_layout()
    plt.savefig(pathhead + '/figs/EwM_' + zm + '.png')
    #plt.show()

    # Plot total spectrum
    plt.clf()
    plt.semilogy(eSpectrumT[:halflen,0]/omeg0, np.abs(eOmTotal[:halflen])**2, label=labstr)
    plt.title(r'Total spectrum at $z={:8.2f}$ $\mu$m'.format(zm_f*1e-3))
    #plt.ticklabel_format(axis='x', style='sci', scilimits=(15,15))
    plt.xlabel(r'$\omega/\omega_0$')
    plt.ylabel(r'$|E(\omega)|^2$ [V$\cdot$s/m]')
    plt.grid(which='major')
    plt.minorticks_on()
    plt.legend()
    #plt.xlim([0.0, np.max(omeg)])
    plt.margins(x=0)
    plt.tight_layout()
    plt.savefig(pathhead + '/figs/Ew_total_' + zm + '.png')
    #plt.show()


    # Plot total spectrum in THz
    plt.clf()
    plt.semilogy(eSpectrumT[:halflen,0], np.abs(eOmTotal[:halflen])**2, 'b-o', label=labstr)
    plt.title(r'Total spectrum at $z={:8.2f}$ $\mu$m'.format(zm_f*1e-3))
    plt.ticklabel_format(axis='x', style='sci', scilimits=(12,12))
    plt.xlabel(r'$\omega$')
    plt.minorticks_on()
    plt.ylabel(r'$|E(\omega)|^2$ [V$\cdot$s/m]')
    plt.grid(which='major')
    plt.legend()
    plt.xlim([0.0, 100e12])
    plt.ylim([1e-6, 1])
    plt.margins(x=0)
    plt.tight_layout()
    plt.savefig(pathhead + '/figs/Ew_total_THz_' + zm + '.png')
    #plt.show()

    # Plot the plasma density and EtP
    fig, ax1 = plt.subplots(figsize=(8, 6))

    ax2 = ax1.twinx()

    ax1.plot(df['t [s]'], df['Re(Ept) [V/m]'], '-', color='lightcoral')
    ax1.set_ylabel(r'Re[$E(t)$] [V/m]')

    ax2.semilogy(df['t [s]'], neVec)
    ax1.set_title(r'Electron density and forward-propagating field at $z={:8.2f}$ $\mu$m'.format(zm_f*1e-3))
    ax1.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
    ax1.set_xlabel(r'$t$ [s]')
    ax2.set_ylabel(r'$N_e$')
    ax1.set_xlim([np.mean(df['t [s]']) - 3*taup, np.mean(df['t [s]']) + 3*taup])
    ax1.grid(which='major')
    ax1.margins(x=0)

    plt.savefig(pathhead + '/figs/Ne_EtP_' + zm + 'nm.png')
    #plt.show()

    # Plot total and forward spectrum in THz
    plt.clf()
    plt.semilogy(eSpectrumT[:halflen,0], np.abs(eOmP[:halflen])**2, 'b-o', label='fwd; ' + labstr)
    plt.semilogy(eSpectrumT[:halflen,0], np.abs(eOmTotal[:halflen])**2, 'r-o', label='total; ' + labstr)
    plt.title(r'Both spectra at $z={:8.2f}$ $\mu$m'.format(zm_f*1e-3))
    plt.ticklabel_format(axis='x', style='sci', scilimits=(12,12))
    plt.xlabel(r'$\omega$')
    plt.ylabel(r'$|E(\omega)|^2$ [V$\cdot$s/m]')
    plt.grid(which='major')
    plt.minorticks_on()
    plt.legend()
    plt.xlim([0.0, 100e12])
    plt.ylim([1e-6, 1])
    plt.margins(x=0)
    plt.tight_layout()
    plt.savefig(pathhead + '/figs/Ew_total_THz_' + zm + '.png')
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



