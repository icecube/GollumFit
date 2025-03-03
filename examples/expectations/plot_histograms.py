import numpy as np
import matplotlib.pyplot as plt
import h5py
from matplotlib import rc
rc('text', usetex=True)

file = h5py.File('example_label.h5', 'r')

nominal = file['hists'][0]
perturbed = file['hists'][1]

#hist = nominal
hist = perturbed

fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=300)
aspect = (0.4)
vmin=1
vmax=10000

from matplotlib.colors import LogNorm

# Create a log-scaled color map
norm = LogNorm(vmin=vmin, vmax=vmax)

# Plot the first subplot (left)
im1 = axes[0].imshow(hist[0].T, cmap='GnBu', origin='lower',extent=[-1, 0, np.log10(300), np.log10(1e5)],aspect=aspect, norm=norm)
axes[0].set_ylabel(r'$\log_{10}$($E$/GeV)')
axes[0].set_xlabel(r'$\cos(\theta)$')
axes[0].set_title('starting')


# Plot the second subplot (right)
im2 = axes[1].imshow(hist[1].T, cmap='GnBu', origin='lower', extent=[-1, 0, np.log10(300), np.log10(1e5)],aspect=aspect, norm=norm)
axes[1].set_ylabel(r'$\log_{10}$($E$/GeV)')
axes[1].set_xlabel(r'$\cos(\theta)$')
axes[1].set_title('throughgoing')


# Add a colorbar to both subplots #fraction=0.046, pad=0.04
cbar1 = fig.colorbar(im1, ax=axes[0],fraction=0.046, pad=0.04, label=r'$N$')
cbar2 = fig.colorbar(im2, ax=axes[1],fraction=0.046, pad=0.04, label=r'$N$')

xedges = np.linspace(-1, 0, 20 + 1)
yedges = np.linspace(np.log10(300), np.log10(1e5), 24 + 1)

for i in range(len(xedges) - 1):
    for j in range(len(yedges) - 1):
        axes[0].text(xedges[i] + 0.5 * (xedges[i + 1] - xedges[i]),
                 yedges[j] + 0.5 * (yedges[j + 1] - yedges[j]),
                 f'{int(hist[0][i, j])}',
                 color='black',
                 ha='center',
                 va='center', fontsize=4)
        axes[1].text(xedges[i] + 0.5 * (xedges[i + 1] - xedges[i]),
                 yedges[j] + 0.5 * (yedges[j + 1] - yedges[j]),
                 f'{int(hist[1][i, j])}',
                 color='black',
                 ha='center',
                 va='center', fontsize=4)

plt.tight_layout()
plt.savefig('example_hist.pdf')