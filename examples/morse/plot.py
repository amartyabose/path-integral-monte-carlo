import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import h5py
import seaborn as sns
plt.style.use('ggplot')

beta = 10

def hist2d(x, y, wt, title):
    plt.figure()
    hist, xbins, ybins, _ = plt.hist2d(x, y, weights=wt, normed=True, bins=100, range=[[-2,3],[-3,3]])
    plt.grid(False)
    plt.xlabel('x')
    plt.ylabel('p')
    plt.title(title)
    plt.colorbar()

f = h5py.File('morse.h5', 'r')
wigner = f[f'wigner/beta={beta} props=16 4X4P']
pimc = f[f'pimc/beta={beta} props=15']

with PdfPages('output.pdf') as pdf:
    hist2d(wigner[:,2], wigner[:,3], wigner[:,0], f'Wigner beta={beta}')
    pdf.savefig()
    plt.close()

    _, xbin, _ = plt.hist(pimc[:,2], weights=pimc[:,0], density=True, bins=300, label='PIMC', histtype='step', linewidth=1.5)
    plt.hist(wigner[:,2], weights=wigner[:,0], density=True, bins=xbin, label='4X4P', histtype='step', linewidth=1.5)
    plt.title(f'x distribution: beta={beta}')
    plt.legend()
    pdf.savefig()
    plt.close()
    
    plt.figure()
    plt.hist(wigner[:,3], weights=wigner[:,0], density=True, bins=300, label='3X4P', histtype='step', linewidth=1.5)
    plt.title(f'p distribution: beta={beta}')
    pdf.savefig()
    plt.close()
