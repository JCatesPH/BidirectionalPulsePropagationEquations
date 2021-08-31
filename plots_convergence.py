#%%
import sys
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors

pathhead = 'DATA'

if len(sys.argv) > 1:
    pathhead = str(sys.argv[1])
else:
    pathhead = 'DATA'

transmitted_li = []
reflected_li = []
#%% Read in spectrum of transmitted pulse
fig = plt.figure(figsize=(8,6))
ax = fig.gca(projection='3d')
z_it = [5, 4, 3, 2, 1, 0]

for itnum in z_it:
    df = pd.read_table(pathhead + '/Spectrum_iteration_' + str(itnum) + '_Transmitted.dat', header=0)
    #df = pd.read_table('DATA/Spectrum_iteration_' + str(itnum) + '_Transmitted.dat', header=0)
    eSpectrumT = df.values
    transmitted_li.append(eSpectrumT)

    idx1 = 0
    idx2 = int(eSpectrumT.shape[0] / 5)
    ax.plot(eSpectrumT[idx1:idx2,0], np.log10(eSpectrumT[idx1:idx2,3]**2), zs=itnum, zdir='y')
    print('itnum = {:}, max(|E|) = {:.2f}'.format(itnum, np.max(eSpectrumT[:,3])))

ax.set_yticks(z_it)
ax.set_xlabel(r'$\omega$')
ax.set_ylabel('Iteration')
ax.set_zlabel('Intensity [log]')
ax.set_title('Transmitted spectrum convergence')

plt.tight_layout()
#plt.show()
plt.savefig(pathhead + '/figs/EwT_convergence.png')
#plt.show()
######################################
#%%
errVec = np.zeros((5,1))
transmitted_li.reverse()
fig, axs = plt.subplots(figsize=(11, 9))

for i in range(1,6):
    errVec[i-1] = np.linalg.norm(transmitted_li[i][idx1:idx2,3]-transmitted_li[i-1][idx1:idx2,3], ord=np.inf) / np.linalg.norm(transmitted_li[i-1][idx1:idx2,3], ord=np.inf)
    axs.semilogy(transmitted_li[i][idx1:idx2,0], (transmitted_li[i][idx1:idx2,3]-transmitted_li[i-1][idx1:idx2,3])**2 + 1e-15, label='{}-{}'.format(i,i-1))

axs.legend()
axs.set_xlabel(r'$\omega$')
axs.set_ylabel(r'$(\Delta A_+)^2$')

plt.savefig(pathhead + '/figs/delta_EwT.png')

#%%
fig, axs = plt.subplots(figsize=(11, 9))
axs.semilogy(np.arange(1,6), errVec)
axs.set_xticks(np.arange(1,6))
axs.set_title('Transmitted Convergence')
plt.savefig(pathhead + '/figs/normErr_EwT.png')

#%%

######################################
#%%
fig = plt.figure(figsize=(8,6))
ax = fig.gca(projection='3d')
z_it = [5, 4, 3, 2, 1, 0]

for itnum in z_it:
    df = pd.read_table(pathhead + '/Spectrum_iteration_' + str(itnum) + '_Reflected.dat', header=0)
    #df = pd.read_table('DATA/Spectrum_iteration_' + str(itnum) + '_Transmitted.dat', header=0)
    eSpectrumR = df.values
    reflected_li.append(eSpectrumR)

    idx1 = 0
    idx2 = int(eSpectrumR.shape[0] / 5)
    ax.plot(eSpectrumR[idx1:idx2,0], np.log10(eSpectrumR[idx1:idx2,3]**2), zs=itnum, zdir='y')
    print('itnum = {:}, max(|E|) = {:.2f}'.format(itnum, np.max(eSpectrumR[:,3])))

ax.set_yticks(z_it)
ax.set_xlabel(r'$\omega$')
ax.set_ylabel('Iteration')
ax.set_zlabel('Intensity [log]')
ax.set_title('Reflected spectrum convergence')

plt.tight_layout()
#plt.show()
plt.savefig(pathhead + '/figs/EwR_convergence.png')


#%%
reflected_li.reverse()
fig, axs = plt.subplots(figsize=(11, 9))

for i in range(1,6):
    errVec[i-1] = np.linalg.norm(reflected_li[i][idx1:idx2,3]-reflected_li[i-1][idx1:idx2,3], ord=np.inf) / np.linalg.norm(reflected_li[i-1][idx1:idx2,3], ord=np.inf)
    axs.semilogy(reflected_li[i][idx1:idx2,0], (reflected_li[i][idx1:idx2,3]-reflected_li[i-1][idx1:idx2,3])**2 + 1e-15, label='{}-{}'.format(i,i-1))

axs.legend()
axs.set_xlabel(r'$\omega$')
axs.set_ylabel(r'$(\Delta A_-)^2$')

plt.savefig(pathhead + '/figs/delta_EwT.png')

#%%
fig, axs = plt.subplots(figsize=(11, 9))
axs.semilogy(np.arange(1,6), errVec)
axs.set_title('Reflected Convergence')
axs.set_xticks(np.arange(1,6))
plt.savefig(pathhead + '/figs/normErr_EwR.png')
plt.show()

#####################################
#%%
z_it = [5, 4, 3, 2, 1, 0]
for itnum in z_it:
    fig = plt.figure(figsize=(8,6))
    ax = fig.gca()
    df = pd.read_table(pathhead + '/Spectrum_iteration_' + str(itnum) + '_Transmitted.dat', header=0)
    #df = pd.read_table('DATA/Spectrum_iteration_' + str(itnum) + '_Transmitted.dat', header=0)
    eSpectrumT = df.values
    idx1 = 0
    idx2 = int(eSpectrumT.shape[0] / 5)
    ax.plot(eSpectrumT[idx1:idx2,0], np.log10(eSpectrumT[idx1:idx2,3]))

    ax.set_title(r'Transmitted $|E(\omega)|$, Iteration {:}'.format(itnum))
    ax.set_xlabel(r'$\omega$')
    print('itnum = {:}, max(|E|) = {:.2f}'.format(itnum, np.max(eSpectrumT[:,3])))

    plt.tight_layout()
    plt.savefig(pathhead + '/figs/EwT_it{:}.png'.format(itnum))

#plt.show()

######################################

#%%
fig = plt.figure(figsize=(8,6))
ax = fig.gca()

itnum1 = 1
infile = pathhead + '/Spectrum_iteration_' + str(itnum1) + '_Reflected.dat'
print(infile)
df1 = pd.read_table(infile, header=0)
edat1 = df1.values

itnum2 = 3
infile = pathhead + '/Spectrum_iteration_' + str(itnum2) + '_Reflected.dat'
print(infile)
df2 = pd.read_table(infile, header=0)
edat2 = df2.values

ax.semilogy(edat1[idx1:idx2,0], edat1[idx1:idx2,3], edat2[idx1:idx2,0], edat2[idx1:idx2,3])

ax.set_title(r'Transmitted $|E(\omega)|$, Iteration {:} and {:}'.format(itnum1, itnum2))
ax.set_xlabel(r'$\omega$')

plt.tight_layout()
#plt.show()
plt.savefig(pathhead + '/figs/Ew_test.png'.format(itnum))


# %%
