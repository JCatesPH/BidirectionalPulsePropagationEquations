#%%
import sys
import glob
import math
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
#from scipy.fftpack import fft, ifft
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
    #pathhead = '../DATA/silica_indexMatchedKerr_100TWcm2_4umL_050522'
    pathhead = '../DATA/silica_fdtdComp_Kerr_052022'
    itnum = '110'

#%%
mpl.rcParams['font.family'] = 'Tahoma'
plt.rcParams['font.size'] = 14
plt.rcParams['axes.linewidth'] = 2
stdfigsize = (6.66, 5)

#plt.ion()
######################################
#%%
CLIGHT = 299792458
EPS0 = 8.85418782e-12
CHARGE_E = 1.602176634e-19
MASS_E = 9.10938356e-31

def plasmaFreq(nE):
    return np.sqrt(CHARGE_E**2*nE/(EPS0*MASS_E))

#%%
def plotField(timeArr, fieldArr, labelArr, filePath, titleStr=None, xlabelStr=r't [s]', xlims=None):
    fig, axs = plt.subplots(figsize=(11, 9))

    for time, field, label in zip(timeArr, fieldArr, labelArr):
        axs.plot(time, field, label=label)

    #ax2.semilogy(df['t [s]'], neVec)
    if titleStr is not None:
        axs.set_title(titleStr, fontsize=18)
    #axs.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
    axs.set_xlabel(xlabelStr, fontsize=16)
    axs.set_ylabel(r'Re[$E(t)$] [V/m]')
    if xlims is not None:
        axs.set_xlim(xlims)
    axs.tick_params(axis='both', which='major', labelsize=14)
    axs.grid(which='major')
    axs.margins(x=0)
    axs.legend()

    fig.savefig(filePath)
    return fig, axs

def plotSpectrum(freqArr, spectrumArr, labelArr, filePath, titleStr=None, xlabelStr=r'$\omega/\omega_0$', xlims=None):
    fig, axs = plt.subplots(figsize=stdfigsize, dpi=400)

    for freqs, spectrum, label in zip(freqArr, spectrumArr, labelArr):
        axs.semilogy(freqs, spectrum, label=label)

    #ax2.semilogy(df['t [s]'], neVec)
    if titleStr is not None:
        axs.set_title(titleStr)
    #axs.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
    axs.set_xlabel(xlabelStr)
    #axs.set_ylabel(r'Re[$E(t)$] [V/m]')
    if xlims is not None:
        axs.set_xlim(xlims)
    #axs.grid(which='major')
    axs.margins(x=0)
    axs.legend()

    plt.tight_layout()

    fig.savefig(filePath, dpi=400, bbox_inches='tight')
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
nAtoms = df.loc[df['Variable'] == 'num_atoms'].values[0,1]
sourceThickness = df.loc[df['Variable'] == 'LHSsourceLayerThickness'].values[0,1]
freqLowerCutoff = int(df.loc[df['Variable'] == 'freqLowerCutoff'].values[0,1])
freqUpperCutoff = int(df.loc[df['Variable'] == 'freqUpperCutoff'].values[0,1])
#freqUpperCutoff = 8192
numT = df.loc[df['Variable'] == 'num_t'].values[0,1]
domT = df.loc[df['Variable'] == 'domain_t'].values[0,1]

labstr = r'$I_0={:4.2f}$ [TW/cm$^2$], $\omega_0={:4.2f}$ [THz], $\tau_p={:4.2f}$ [fs]'.format(iIntensity*1e-16, omeg0*1e-12/(2*np.pi), taup*1e15)

intensityFactor = EPS0*CLIGHT/2 #/ 32768

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
#%% Load in the structure file
df = pd.read_table(pathhead + '/windowFunc.dat',
        index_col=False
    )
windowfunc = df.values

plt.figure(1, figsize=(8, 6))
plt.clf()
plt.plot(windowfunc[:,0], windowfunc[:,1])
plt.title(r'Window Function')
plt.tight_layout()
plt.savefig(pathhead + '/figs/windowFunc.png')


#%% Read in input pulse in time domain
df = pd.read_table(pathhead + '/InputEfield_1D.dat', header=0)
eFieldIn = df.values

#gaussianEnvelope = np.sqrt(iIntensity/intensityFactor) * np.exp(-2*np.log(2)*(eFieldIn[:,0] - domT/2)**2/taup**2)

#%%
plt.figure(2, figsize=(8, 6))
plt.clf()
plt.plot(eFieldIn[:,0], eFieldIn[:,1], 'b-o')
#plt.semilogy(eFieldIn[:,0], eFieldIn[:,1]**2, 'b-o')
#plt.plot(eFieldIn[:,0], gaussianEnvelope, '-r')
plt.xlim([domT/2 - 3*taup, domT/2 + 3*taup])
plt.title(r'Incident pulse')
plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
plt.xlabel(r'$t$ [s]')
plt.ylabel(r'$E(t)$ [V/m]')
plt.grid(which='major')

plt.tight_layout()
plt.savefig(pathhead + '/figs/Input_Et.png')
plt.show()


#%% Read in input spectrum
df = pd.read_table(pathhead + '/InputSpectrum_1D.dat', header=0)
inSpectrum = df.values

obsMaxIntensity = np.max(inSpectrum[:,3]**2 * intensityFactor)
print("Actual maximum intensity: {:.3e}".format(obsMaxIntensity))
print("  Ratio to input intensity: {:}".format(obsMaxIntensity / iIntensity))

#%% Plot input spectrum
plotSpectrum([inSpectrum[:,0]/omeg0], [inSpectrum[:,3]**2 * intensityFactor], [labstr], 
    filePath=pathhead + '/figs/inSpectrum.png', 
    titleStr='Input Spectrum', 
    xlabelStr=r'$\omega/\omega_0$')
plt.show()


#%% Read in spectrum of transmitted pulse
df = pd.read_table(pathhead + '/Spectrum_iteration_' + itnum + '_Transmitted.dat', header=0)
eSpectrumT = df.values

halfT = int(eSpectrumT.shape[0]/2)

dOm = abs(eSpectrumT[1,0] - eSpectrumT[0,0])
print('dOmega = {:.2e}'.format(dOm))

#%% Plot the transmitted spectrum
plotSpectrum([eSpectrumT[freqLowerCutoff:freqUpperCutoff,0]/omeg0], [eSpectrumT[freqLowerCutoff:freqUpperCutoff,3]**2*intensityFactor], [labstr], 
    filePath=pathhead + '/figs/EwT.png', 
    titleStr='Transmitted Spectrum', 
    xlabelStr=r'$\omega/\omega_0$')



#plt.show()
#%% Read in spectrum of reflected pulse
df = pd.read_table(pathhead + '/Spectrum_iteration_' + itnum + '_Reflected.dat', header=0)
eSpectrumR = df.values

#%% Plot the reflected spectrum
plotSpectrum([eSpectrumR[freqLowerCutoff:freqUpperCutoff,0]/omeg0], [eSpectrumR[freqLowerCutoff:freqUpperCutoff,3]**2*intensityFactor], [labstr], 
    filePath=pathhead + '/figs/EwR.png', 
    titleStr='Reflected Spectrum', 
    xlabelStr=r'$\omega/\omega_0$')
######################################

#%% Plot both reflected and transmitted spectra
plotSpectrum([eSpectrumT[freqLowerCutoff:freqUpperCutoff,0]/omeg0, eSpectrumR[freqLowerCutoff:freqUpperCutoff,0]/omeg0], 
    [eSpectrumT[freqLowerCutoff:freqUpperCutoff,3]**2*intensityFactor, eSpectrumR[freqLowerCutoff:freqUpperCutoff,3]**2*intensityFactor], 
    ['Transmitted', 'Reflected'], 
    filePath=pathhead + '/figs/EwBoth.png', 
    titleStr='Both Spectra', 
    xlabelStr=r'$\omega/\omega_0$')
####################################################################
####################################################################

#%% Find point monitor files
pmon_li = [] # List of point monitor files
zmon_li = [] # List of point monitor locations

# Use glob to find all filenames
for name in glob.glob(pathhead + '/PointMon_iter_{}_*'.format(itnum)):
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
    plt.close('all')
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
    #plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
    plt.xlabel(r'$t$ [s]')
    plt.grid(which='major')
    plt.xlim([np.mean(df['t [s]']) - 2*taup, np.mean(df['t [s]']) + 2*taup])
    plt.margins(x=0)
    plt.minorticks_on()
    plt.tight_layout()
    plt.savefig(pathhead + '/figs/Ne_' + zm + 'nm.png')
    #plt.show()

    plt.clf()
    plt.plot(df['t [s]'], df['Je'])
    plt.title(r'Current density at $z={:8.2f}$ $\mu$m'.format(zm_f*1e-3))
    plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
    plt.xlabel(r'$t$ [s]')
    plt.grid(which='major')
    #plt.xlim([np.mean(df['t [s]']) - 2*taup, np.mean(df['t [s]']) + 2*taup])
    plt.margins(x=0)
    plt.minorticks_on()
    plt.tight_layout()
    plt.savefig(pathhead + '/figs/Je_' + zm + 'nm.png')


    # Use max density to calculate plasma frequency
    omega_pe = plasmaFreq(np.max(neVec))
    #omega_pe = plasmaFreq(1.1e25)
    print('Maximum electron density : {:.1e}\n Plasma frequency: {:.3e}\n wp/w0: {:.3e}'.format(np.max(neVec), omega_pe, omega_pe/omeg0))

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
    #plt.legend()
    plt.xlim([np.mean(df['t [s]']) - 5*taup, np.mean(df['t [s]']) + 5*taup])
    plt.margins(x=0)
    plt.tight_layout()
    plt.savefig(pathhead + '/figs/Et_total_' + zm + '.png')
    #plt.show()

    # Find spectra
    JeOm = np.fft.fft(df['Je'].values) / np.sqrt(len(eP))
    eOmP = np.fft.fft(eP) / np.sqrt(len(eP))
    eOmM = np.fft.fft(eM) / np.sqrt(len(eM))
    eOmTotal = np.fft.fft(eTotal) / np.sqrt(len(eTotal))

    halflen = int(len(eOmTotal)/2) + 1
    dom_t = float(df.values[-1:,0] - df.values[0,0])
    omeg = np.array([2*math.pi * i / dom_t for i in range(halflen)])
    upperFreq_THz = np.max(np.nonzero(omeg / (2*np.pi) < 100e12))

    # Plot forward-prop spectrum
    fig, axs = plotSpectrum([omeg[freqLowerCutoff:freqUpperCutoff]/omeg0], 
        [np.abs(eOmP[freqLowerCutoff:freqUpperCutoff])**2*intensityFactor], 
        [labstr], 
        filePath=pathhead + '/figs/EwP_' + zm + '.png', 
        titleStr=r'Forward-propagating spectrum at $z={:6.2f}$ $\mu$m'.format(zm_f*1e-3), 
        #titleStr='',
        xlabelStr=r'$\omega/\omega_0$')

    axs.axvline(omega_pe / omeg0, color='k', linestyle='--')
    fig.savefig(pathhead + '/figs/EwP_' + zm + '.png')

    """ fig, axs = plotSpectrum([omeg[:upperFreq_THz] / (2*np.pi) * 1e-12], 
        [np.abs(eOmP[:upperFreq_THz])**2*intensityFactor], 
        [labstr], 
        filePath=pathhead + '/figs/EwP_THz_' + zm + '.png', 
        titleStr=r'Forward THz spectrum at $z={:6.2f}$ $\mu$m'.format(zm_f*1e-3), 
        #titleStr='',
        xlabelStr=r'f [THz]')
    #axs.set_ylim([1e-9,1e-6])
    axs.axvline(omega_pe / (2*np.pi) * 1e-12, color='k', linestyle='--')
    fig.savefig(pathhead + '/figs/EwP_THz_' + zm + '.png')
    #plt.show() """

    # Plot backward-prop spectrum
    plotSpectrum([omeg[freqLowerCutoff:freqUpperCutoff]/omeg0], 
        [np.abs(eOmM[freqLowerCutoff:freqUpperCutoff])**2*intensityFactor], 
        [labstr], 
        filePath=pathhead + '/figs/EwM_' + zm + '.png', 
        titleStr=r'Backward-propagating spectrum at $z={:6.2f}$ $\mu$m'.format(zm_f*1e-3), 
        xlabelStr=r'$\omega/\omega_0$')

    # Plot total spectrum
    plotSpectrum([omeg[freqLowerCutoff:freqUpperCutoff]/omeg0], 
        [np.abs(eOmTotal[freqLowerCutoff:freqUpperCutoff])**2*intensityFactor], 
        [labstr], 
        filePath=pathhead + '/figs/Ew_' + zm + '.png', 
        titleStr=r'Spectrum of total field at $z={:6.2f}$ $\mu$m'.format(zm_f*1e-3), 
        xlabelStr=r'$\omega/\omega_0$')

    # Plot both spectra
    fig, axs = plotSpectrum([omeg[freqLowerCutoff:freqUpperCutoff]/omeg0, omeg[freqLowerCutoff:freqUpperCutoff]/omeg0], 
        [np.abs(eOmP[freqLowerCutoff:freqUpperCutoff])**2*intensityFactor, np.abs(eOmM[freqLowerCutoff:freqUpperCutoff])**2*intensityFactor], 
        [r'$A_+$', r'$A_-$'], 
        filePath=pathhead + '/figs/EwBoth_' + zm + '.png', 
        titleStr=r'Both spectra at $z={:6.2f}$ $\mu$m'.format(zm_f*1e-3), 
        xlabelStr=r'$\omega/\omega_0$')

    axs.axvline(omega_pe / omeg0, color='k', linestyle='--')
    fig.savefig(pathhead + '/figs/EwBoth_' + zm + '.png')

    # Plot the plasma density and EtP
    fig, ax1 = plt.subplots(figsize=stdfigsize, dpi=400)

    ax2 = ax1.twinx()

    ax1.plot(df['t [s]'], df['Re(Ept) [V/m]'], '-', color='lightcoral')
    ax1.set_ylabel(r'Re[$E(t)$] [V/m]', fontsize=14)

    ax2.semilogy(df['t [s]'], neVec*1e-6)
    #ax1.set_title(r'Electron density and forward-propagating field at $z={:8.2f}$ $\mu$m'.format(zm_f*1e-3))
    ax1.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
    ax1.tick_params(axis='both', which='major', labelsize=12)
    ax2.tick_params(axis='both', which='major', labelsize=12)
    ax1.set_xlabel(r'$t$ [s]', fontsize=14)
    ax2.set_ylabel(r'$\rho$ [cm$^{-3}$]', fontsize=14)
    ax1.set_xlim([np.mean(df['t [s]']) - 20*taup, np.mean(df['t [s]']) + 20*taup])
    ax1.set_title(r'Field and carrier density at $z={:6.2f}$ $\mu$m'.format(zm_f*1e-3))
    ax2.grid(which='major')
    ax1.margins(x=0)

    plt.savefig(pathhead + '/figs/Ne_EtP_' + zm + 'nm.png', dpi=400, bbox_inches='tight')

    ax1.plot(df['t [s]'], df['Re(Emt) [V/m]'], '-', color='mediumpurple')

    plt.savefig(pathhead + '/figs/Ne_EtP_EtM_' + zm + 'nm.png', dpi=400, bbox_inches='tight')
    #plt.show()

    # Plot transmitted and forward spectrum in THz
    """ plt.clf()
    plt.semilogy(omeg[:upperFreq_THz] / (2*np.pi), np.abs(eOmP[:upperFreq_THz])**2*intensityFactor, 'b-o', label='fwd; ' + labstr)
    #plt.semilogy(omeg[:upperFreq_THz] / (np.pi), np.abs(eOmTotal[:upperFreq_THz])**2, 'r-o', label='total; ' + labstr)
    plt.semilogy(eSpectrumT[:upperFreq_THz,0] / (2*np.pi), np.abs(eSpectrumT[:upperFreq_THz,3])**2*intensityFactor, 'r-o', label='trans; ' + labstr)
    #plt.title(r'Both spectra at $z={:8.2f}$ $\mu$m'.format(zm_f*1e-3))
    plt.ticklabel_format(axis='x', style='sci', scilimits=(12,12))
    plt.xlabel(r'$\omega$')
    plt.ylabel(r'$|E(\omega)|^2$ [V$\cdot$s/m]')
    plt.grid(which='major')
    #plt.ylim([1e-9, 1e-6])
    plt.minorticks_on()
    plt.legend()
    #plt.xlim([0.0, 100e12])
    #plt.ylim([1e-6, 1])
    plt.margins(x=0)
    plt.tight_layout()
    plt.savefig(pathhead + '/figs/Ew-EwT_THz_' + zm + '.png')
    #plt.show() """

    # Plot the FFT of the current
    """ plotSpectrum([omeg[:freqUpperCutoff]/omeg0], 
        [np.abs(JeOm[:freqUpperCutoff])], 
        [r'J_e'], 
        filePath=pathhead + '/figs/JeOm_' + zm + '.png', 
        titleStr=r'Current density at $z={:6.2f}$ $\mu$m'.format(zm_f*1e-3), 
        xlabelStr=r'$\omega/\omega_0$') """

plt.close('all')

#%% Find point monitor files
#omeg = np.array([2*math.pi * i / domT for i in range(int(numT/2))])
#upperFreq_THz = np.max(np.nonzero(omeg / (2*np.pi) < 100e12))

pmon_li = [] # List of point monitor files
zmon_li = [] # List of point monitor locations
time_li = []
fieldP_li = []
fieldM_li = []
omeg_li= []
spectraP_li = []
spectraM_li = []

# Use glob to find all filenames
for name in glob.glob(pathhead + '/PointMon_iter_{}_*'.format(itnum)):
    pmon_li.append(name)
    nameli = name.split('_')
    #print(nameli)
    zloc = nameli[-1].split('nm')
    zloc = str((float(zloc[0]) - sourceThickness*1e9) * 1e-3)
    #print('Location of point monitor: {} nm'.format(zloc))
    zmon_li.append(zloc)
    df = pd.read_table(name, 
        header=None, skiprows=[0,1],
        names=['t [s]', 'Re(Ept) [V/m]', 'Im(Ept) [V/m]', 'Re(Emt) [V/m]', 'Im(Emt) [V/m]', 'Ne', 'Je']
        )

    # Find spectra
    eP = df['Re(Ept) [V/m]'].values + 1.0j * df['Im(Ept) [V/m]'].values
    eOmP = np.fft.fft(eP) / np.sqrt(len(eP))

    eM = df['Re(Emt) [V/m]'].values + 1.0j * df['Im(Emt) [V/m]'].values
    eOmM = np.fft.fft(eM) / np.sqrt(len(eM))

    time_li.append(df['t [s]'].values)
    fieldP_li.append(eP)
    fieldM_li.append(eM)
    omeg_li.append(omeg)
    spectraP_li.append(eOmP)
    spectraM_li.append(eOmM)

print(pmon_li)
print(zmon_li)

#%% Plot the first and last few spectra in domain
redOmega_li = [omeg_li[0], omeg_li[-2], omeg_li[-1]]
redMinus_li = [spectraM_li[0], spectraM_li[-2], spectraM_li[-1]]
redMon_li = [zmon_li[0], zmon_li[-2], zmon_li[-1]]

fig, axs = plotSpectrum([omeg[freqLowerCutoff:freqUpperCutoff] / omeg0 for omeg in redOmega_li], 
    [np.abs(eOmM[freqLowerCutoff:freqUpperCutoff])**2*intensityFactor for eOmM in redMinus_li], 
    redMon_li, 
    filePath=pathhead + '/figs/Ew_FirstLast.png', 
    #titleStr=r'Forward THz spectrum at $z={:6.2f}$ $\mu$m'.format(zm_f*1e-3), 
    titleStr='Backward-propagating Spectra',
    xlabelStr=r'$\omega/\omega_0$')
#axs.set_ylim([1e-9,1e-6])
#axs.set_xlim([0, 60])
#fig.savefig(pathhead + '/figs/EwMz.png')

#%%
plotField(time_li[1:-1:2], 
    fieldP_li[1:-1:2], 
    zmon_li[1:-1:2], 
    filePath=pathhead + '/figs/Etz.png', 
    titleStr=None, 
    xlabelStr=r't [s]', 
    xlims=[domT/2 - 20*taup, domT/2 + 20*taup])
#plt.show()

#%%
#omega_pe = plasmaFreq(4e24)
""" fig, axs = plotSpectrum([omeg[:upperFreq_THz] / (2*np.pi) * 1e-12 for omeg in omeg_li], 
    [np.abs(eOmP[:upperFreq_THz])**2*intensityFactor for eOmP in spectra_li], 
    zmon_li, 
    filePath=pathhead + '/figs/EwPz_THz.png', 
    #titleStr=r'Forward THz spectrum at $z={:6.2f}$ $\mu$m'.format(zm_f*1e-3), 
    titleStr='',
    xlabelStr=r'f [THz]')
#axs.set_ylim([1e-9,1e-6])
#axs.set_xlim([0, 60])
axs.axvline(omega_pe / (2*np.pi) * 1e-12, color='k', linestyle='--')
fig.savefig(pathhead + '/figs/EwPz_THz.png') """

fig, axs = plotSpectrum([omeg[freqLowerCutoff:freqUpperCutoff] / omeg0 for omeg in omeg_li], 
    [np.abs(eOmP[freqLowerCutoff:freqUpperCutoff])**2*intensityFactor for eOmP in spectraP_li], 
    zmon_li, 
    filePath=pathhead + '/figs/EwPz.png', 
    #titleStr=r'Forward THz spectrum at $z={:6.2f}$ $\mu$m'.format(zm_f*1e-3), 
    titleStr='Forward-propagating Spectra',
    xlabelStr=r'$\omega/\omega_0$')
#axs.set_ylim([1e-9,1e-6])
#axs.set_xlim([0, 60])
fig.savefig(pathhead + '/figs/EwPz.png')

fig, axs = plotSpectrum([omeg[freqLowerCutoff:freqUpperCutoff] / omeg0 for omeg in omeg_li], 
    [np.abs(eOmM[freqLowerCutoff:freqUpperCutoff])**2*intensityFactor for eOmM in spectraM_li], 
    zmon_li, 
    filePath=pathhead + '/figs/EwMz.png', 
    #titleStr=r'Forward THz spectrum at $z={:6.2f}$ $\mu$m'.format(zm_f*1e-3), 
    titleStr='Backward-propagating Spectra',
    xlabelStr=r'$\omega/\omega_0$')
#axs.set_ylim([1e-9,1e-6])
#axs.set_xlim([0, 60])
fig.savefig(pathhead + '/figs/EwMz.png')

#plt.show()

#%% Plot the transmitted THz spectrum
""" upperFreq_THz = np.max(np.nonzero(eSpectrumT[:halfT,0] / (2*np.pi) < 100e12))
fig, axs = plotSpectrum([eSpectrumT[:upperFreq_THz,0]/(2*np.pi)*1e-12], [eSpectrumT[:upperFreq_THz,3]**2*intensityFactor], [labstr], 
    filePath=pathhead + '/figs/EwT_THz.png', 
    titleStr='THz Spectrum of Transmitted Field', 
    xlabelStr=r'f [THz]')

axs.axvline(omega_pe / (2*np.pi) * 1e-12, color='k', linestyle='--')
fig.savefig(pathhead + '/figs/EwT_THz.png') """

#%% Read in forward-propagating pulse
df = pd.read_table(pathhead + '/Efield_iteration_' + itnum + '_Transmitted.dat', header=0)
eFieldT = df.values

#%%
plt.clf()
plt.figure(figsize=(8, 6))
plt.plot(eFieldT[:,0], eFieldT[:,1])
#plt.semilogy(eFieldT[:,0], eFieldT[:,1]**2)
plt.title(r'Transmitted pulse')
plt.xlim([domT/2 - 5*taup, domT/2 + 5*taup])
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.xlabel(r'$t$ [s]')
plt.ylabel(r'$E(t)$ [V/m]')
plt.grid(which='major')
#plt.show()
plt.tight_layout()
plt.savefig(pathhead + '/figs/Et_transmitted.png')


#%%
A = np.vstack([eFieldT[:,0].T, np.ones(len(eFieldT[:,0]))]).T
lstSqFit = np.linalg.lstsq(A, eFieldT[:,1], rcond=None)
print("The transmitted field has linear fit with parameters:\n  m={:}\n  c={:}".format(lstSqFit[0][0], lstSqFit[0][1]))

#%%
#pySpectrumT = fft(eFieldT[:,1] - lstSqFit[0][1])
#plt.plot(pySpectrumT)


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
A = np.vstack([eFieldR[:,0].T, np.ones(len(eFieldR[:,0]))]).T
lstSqFit = np.linalg.lstsq(A, eFieldR[:,1], rcond=None)
print("The reflected field has linear fit with parameters:\n  m={:}\n  c={:}".format(lstSqFit[0][0], lstSqFit[0][1]))


#############################################
# %%
# Get data from specific point monitor
index = 1
eForward = fieldP_li[index]
eBackward = fieldM_li[index]
tArr = 1.2e-12 - time_li[index]
zLoc = zmon_li[index]

# %% Plot both with scaled backscattered
rFactor = 1e4 # Scaling factor for backscattered

plt.clf()
plt.figure(figsize=(8, 6))

plt.plot(tArr, eForward)
plt.plot(tArr, rFactor * eBackward)

plt.title(r'Both fields at $z=${:.2f} $\mu$m'.format(float(zLoc)))
plt.xlabel(r'$t$ [s]')
plt.ylabel(r'$E(t)$ [V/m]')
plt.legend([r'$E_+$' , r'$E_-$ (*{:.1e})'.format(rFactor)])

plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.grid(which='major')
plt.tight_layout()
plt.xlim([-2e-13,2e-13])
#plt.show()
plt.savefig(pathhead + '/figs/Et_ScaledBack.png')


# %%
