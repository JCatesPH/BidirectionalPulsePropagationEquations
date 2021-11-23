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
Z0 = df.loc[df['Variable'] == 'Znaught'].values[0,1]

labstr = r'$I_0={:4.2f}$ [TW/cm$^2$], $\omega_0={:4.2f}$ [THz], $\tau_p={:4.2f}$ [fs]'.format(iIntensity*1e-16, omeg0*1e-12, taup*1e15)

######################################


#%% Read in input spectrum
df = pd.read_table(pathhead + '/InputSpectrum_1D.dat', header=0)
inSpectrum = df.values

omeg = inSpectrum[:,0]
Ez = inSpectrum[:,3]
k0arr = omeg / 3e8

#%% Plot input spectrum
plt.clf()
plt.figure(figsize=(8, 6))
plt.plot(omeg, Ez, label=labstr)
plt.title(r'Input Amplitude Spectrum')
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.xlabel(r'$\omega$')
plt.grid(which='major')
plt.legend()
plt.xlim([0.0, np.max(inSpectrum[:,0])])
plt.tight_layout()
plt.savefig(pathhead + '/figs/inSpectrum.png')
plt.show()

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

    # Define field variables
    eP = df['Re(Ept) [V/m]'].values + 1.0j * df['Im(Ept) [V/m]'].values
    eM = df['Re(Emt) [V/m]'].values + 1.0j * df['Im(Emt) [V/m]'].values
    # Find spectra
    eOmP = np.fft.fft(eP)
    eOmM = np.fft.fft(eM)

    # Define frequencies
    halflen = int(len(eOmP)/2)
    dom_t = float(df.values[-1:,0] - df.values[0,0])
    omeg = [math.pi * i / dom_t for i in range(halflen)]

    # Plot forward-prop spectrum
    plt.clf()
    plt.semilogy(omeg, np.abs(eOmP[:halflen]), label=r'$E_p$')
    plt.semilogy(omeg, np.abs(eOmM[:halflen]), label=r'$E_m$')
    plt.title(r'Spectra at $z={:6.2f}$ $\mu$m'.format(zm_f*1e-3))
    plt.ticklabel_format(axis='x', style='sci', scilimits=(15,15))
    plt.xlabel(r'$\omega$')
    plt.ylabel(r'$|E(\omega)|$ [V$\cdot$s/m]')
    plt.grid(which='major')
    plt.legend()
    plt.xlim([0.0, np.max(omeg)])
    plt.margins(x=0)
    plt.tight_layout()
    plt.savefig(pathhead + '/figs/Spectra_' + zm + '.png')
    #plt.show()

#%%
tr = np.abs(eOmP[:halflen]) / Ez[:-1]

plt.clf()
plt.plot(omeg, tr)
plt.xlabel(r'$\omega$')
plt.grid(which='major')
plt.legend()
plt.xlim([0.8*omeg0, 1.2*omeg0])
plt.ylim([0, 2])
plt.margins(x=0)
plt.tight_layout()
plt.savefig(pathhead + '/figs/SpectraRatioP_' + zm + '.png')
plt.show()

rf = np.abs(eOmM[:halflen]) / Ez[:-1]

plt.clf()
plt.plot(omeg, rf)
plt.xlabel(r'$\omega$')
plt.grid(which='major')
plt.legend()
plt.xlim([0.8*omeg0, 1.2*omeg0])
plt.ylim([0, 2])
plt.margins(x=0)
plt.tight_layout()
plt.savefig(pathhead + '/figs/SpectraRatioM_' + zm + '.png')
plt.show()

#%% Print mean coefficients
print('avg tr = {:}'.format(np.mean(tr)))
print('avg rf = {:}'.format(np.mean(rf)))

#%%
############################################################################
#    START TRANSFER-MATRIX METHOD
############################################################################

#%% Problem setup
Hz = Ez / Z0 
v = np.array([Ez, Hz])

#%%
k1 = k0arr[1] # Wavenumber in layer 1
kL = k1       # Wavenumber on LHS
kR = k1       # Wavenumber on RHS
L = 1e-6      # Length of layer 1
v0 = v[:,1]   # Initial amplitude
M1 = np.array([[np.cos(k1*L), 1j*Z0*np.sin(k1*L)],[1j/Z0*np.sin(k1*L), np.cos(k1*L)]]) # Transfer matrix for layer 1
EL1 = M1 @ v0 # Compute field at z+L

print('k = {:.2f}'.format(k1))
print('L = {:.2f}'.format(L))
print('v0 =')
print(v0)
print('\nM1 = ')
print(M1)
print('\nEL1 = ')
print(EL1)
print('\nT = {:}'.format(np.abs(EL1[0]/v0[0])))

# %% Using Griffiths method naively
n1 = 1.0
n2 = 1.75
n3 = 1.0
E0R = 1
E0L = np.abs((n1-n2) / (n1 + n2)) * E0R
E1R = 2*n1 / (n1 + n2) * E0R
E1L = np.abs((n2-n3) / (n2 + n3)) * E1R
E2R = 2*n3 / (n3 + n2) * E1R

print('E0R = {:}'.format(E0R))
print('E0L = {:}'.format(E0L))
print('E1R = {:}'.format(E1R))
print('E1L = {:}'.format(E1L))
print('E2R = {:}'.format(E2R))
print('T01 + R01 = {:}'.format(E0L+E1R))
print('T12 + R12 = {:}'.format((E1L+E2R)/E1R))

#%%
for i in range(6):
    E0L = np.abs((n1-n2) / (n1 + n2)) * E0R + 2*n2 / (n1 + n2) * E1L
    E1R = 2*n1 / (n1 + n2) * E0R + np.abs((n1-n2) / (n1 + n2)) * E1L
    E1L = np.abs((n2-n3) / (n2 + n3)) * E1R
    E2R = 2*n3 / (n3 + n2) * E1R

print('E0R = {:}'.format(E0R))
print('E0L = {:}'.format(E0L))
print('E1R = {:}'.format(E1R))
print('E1L = {:}'.format(E1L))
print('E2R = {:}'.format(E2R))
print('T + R = {:}'.format(E0L+E2R))
print('T01 + R01 = {:}'.format(E0L+E1R))
print('T12 + R12 = {:}'.format((E1L+E2R)/E1R))

# %%
