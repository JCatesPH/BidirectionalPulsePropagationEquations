#%%
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy.fft import fft, ifft

CLIGHT = 299792458
EPS0 = 8.85418782e-12
CHARGE_E = 1.602176634e-19
MASS_E = 9.10938356e-31

#%% Set parameters
num_t = 2**13
dom_t = 1.2e-12
delta_t = dom_t / num_t

print('num_t = {}'.format(num_t))
print('dom_t = {:.3e}'.format(dom_t))
print('delta_t = {:.3e}'.format(delta_t))

#%% Calculate FFT frequencies based on parameters
max_omeg = np.pi / delta_t
delta_om = 2 * np.pi / dom_t
delta_f = 1 / dom_t

print('max_omeg = {:.3e}'.format(max_omeg))
print('delta_omeg = {:.3e}'.format(delta_om))
print('delta_f = {:.3e}'.format(delta_f))
print('min_lambda = {:.3e}'.format(2 * np.pi * CLIGHT / max_omeg))
print('max_lambda = {:.3e}'.format(2 * np.pi * CLIGHT / delta_om))

#%% Use the central frequency for further analysis
lamb0 = 1.9e-6
omeg0 = 2 * np.pi * CLIGHT / lamb0 
maxHarm = int(max_omeg / omeg0)

print('max harmonic = {:d}'.format(maxHarm))
print('points per fund. wavelength = {:d}'.format(int((lamb0/CLIGHT)/delta_t)))
print('points per min. wavelength = {:d}'.format(int(2*np.pi/(12*omeg0)/delta_t)))
#print('points per THz wavelength = {:d}'.format(int((1/10e12)/delta_t)))



#%%
omegArr = np.array([i*delta_om for i in range(int(num_t/2))])
firstHarmonic = 1.5*omeg0
thirdHarmonic = 3.5*omeg0
fifthHarmonic = 5.5*omeg0
sixthHarmonic = 6.5*omeg0
seventhHarmonic = 7.5*omeg0
lastHarmonic = (maxHarm - 0.5)*omeg0

cutoff1st = np.max(np.nonzero(omegArr < firstHarmonic))
cutoff3rd = np.max(np.nonzero(omegArr < thirdHarmonic))
cutoff5th = np.max(np.nonzero(omegArr < fifthHarmonic))
cutoff6th = np.max(np.nonzero(omegArr < sixthHarmonic))
cutoff7th = np.max(np.nonzero(omegArr < seventhHarmonic))
cutoffLast = np.max(np.nonzero(omegArr < lastHarmonic))

print('cutoff for 1st = {:d}'.format(cutoff1st))
print('cutoff for 3rd = {:d}'.format(cutoff3rd))
print('cutoff for 5th = {:d}'.format(cutoff5th))
print('cutoff for 6th = {:d}'.format(cutoff6th))
print('cutoff for 7th = {:d}'.format(cutoff7th))
print('cutoff for last = {:d}'.format(cutoffLast))

customHarmonic = (15 + 0.5)*omeg0
cutoffCustom = np.max(np.nonzero(omegArr < customHarmonic))
print('custom cutoff = {:d}'.format(cutoffCustom))

#%% Find the cutoffs for transparency window of silica
omeg1 = 2 * np.pi * CLIGHT / 0.21e-6
omeg2 = 2 * np.pi * CLIGHT / 6.7e-6
cutoffOm1 = np.max(np.nonzero(omegArr < omeg1))
cutoffOm2 = np.max(np.nonzero(omegArr < omeg2))
print('cutoff for 0.21 = {:d}'.format(cutoffOm1))
print('cutoff for 6.7 = {:d}'.format(cutoffOm2))

# %%
numE = 1.0e27
omegPlasma = np.sqrt(CHARGE_E**2*numE/(EPS0*MASS_E))

print("For {:} electrons, the plasma frequency is: {:.3e} 1/s".format(numE, omegPlasma))
print("Wavelength of {:.2e} [m] gives frequecy {:.2e} 1/s".format(lamb0, omeg0))
print("Approximate plasma skin depth: {:.2e}".format(CLIGHT/omegPlasma))
print("Critial Power: {:.3e} GW".format(17*(omeg0/omegPlasma)**2))
print("Frequecy (not angular) of central wavelength: {:.2e} Hz".format(omeg0/(2*np.pi)))

# %%
numE_li = [1e23, 1e24, 1e25, 1e26, 1e27]
omegp_li = []
lamb_li = []
pcr_li = []
for numE in numE_li:
    omegPlasma = np.sqrt(CHARGE_E**2*numE/(EPS0*MASS_E))
    omegp_li.append(omegPlasma)

    omeg0 = 1.1 * omegPlasma
    lamb0 = 2 * np.pi * CLIGHT / omeg0
    lamb_li.append(lamb0)

    halfCritPower = (17 * (omeg0 / omegPlasma)**2) / 2
    pcr_li.append(halfCritPower)

dat = {'nE [m^-3]' : numE_li, 'omeg_p' : omegp_li, 'lambda_0' : lamb_li, 'Half P_cr [GW]' : pcr_li}
df = pd.DataFrame(data = dat)
df
# %%
lamb0 = 3e-6
omeg0 = 2 * np.pi * CLIGHT / lamb0 
numE_li = [1e23, 1e24, 1e25, 1e26, 1e27]
omegp_li = []
ratio_li = []
pcr_li = []
for numE in numE_li:
    omegPlasma = np.sqrt(CHARGE_E**2*numE/(EPS0*MASS_E))
    omegp_li.append(omegPlasma)

    ratio = omeg0 / omegPlasma
    ratio_li.append(ratio)

    halfCritPower = (17 * (omeg0 / omegPlasma)**2) / 2
    pcr_li.append(halfCritPower)

dat = {'nE [m^-3]' : numE_li, 'omeg_p' : omegp_li, 'omeg0/omegp' : ratio_li, 'Half P_cr [GW]' : pcr_li}
df = pd.DataFrame(data = dat)
df

# %%
