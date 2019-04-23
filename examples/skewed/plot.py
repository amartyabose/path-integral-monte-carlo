import numpy as np
import matplotlib.pyplot as plt
import h5py
import seaborn as sns
plt.style.use('ggplot')

beta = 3

def hist2d(x, y, wt, title):
    plt.figure()
    hist, xbins, ybins, _ = plt.hist2d(x, y, weights=wt, normed=True, bins=100, range=[[-2,3],[-5,5]])
    plt.xlabel('x')
    plt.ylabel('p')
    plt.title(title)
    plt.colorbar()

f = h5py.File('skewed.h5', 'r')
wigner = f[f'wigner/beta={beta} props=6']
pimc = f[f'pimc/beta={beta} props=5']

_, xbin, _ = plt.hist(pimc[:,2], weights=pimc[:,0], density=True, bins=500, label='PIMC', histtype='step', linewidth=1.5)
plt.hist(wigner[:,2], weights=wigner[:,0], density=True, bins=xbin, label='3X4P', histtype='step', linewidth=1.5)
plt.title('x distribution')
plt.legend()

plt.figure()
plt.hist(wigner[:,3], weights=wigner[:,0], density=True, bins=500, label='3X4P', histtype='step', linewidth=1.5)
plt.title('p distribution')

hist2d(wigner[:,2], wigner[:,3], wigner[:,0], 'Wigner')

plt.show()
