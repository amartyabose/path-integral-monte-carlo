import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')

data_pimc = np.loadtxt('pimc.beta1.G6.props5/phasespace/nd.0.configs.xyz', usecols=(2,))
data_wigner = np.loadtxt('wigner.beta1.G6.props5/phasespace/nd.0.configs.xyz', usecols=(0,2,3))

plt.hist(data_pimc, density=True, bins=100, alpha=0.5, label='PIMC')
plt.hist(data_wigner[:,1], weights=data_wigner[:,0], density=True, bins=100, alpha=0.5, label='Wigner')
plt.xlabel('x')
plt.ylabel('Density')
plt.xlim((-3,3))
plt.legend()
plt.savefig('position_hist.pdf', bbox_inches='tight')
plt.show()

plt.hist2d(data_wigner[:,1], data_wigner[:,2], weights=data_wigner[:,0], density=True, bins=50)
plt.grid(False)
plt.xlabel('x')
plt.ylabel('p')
plt.xlim((-3,3))
plt.ylim((-3,3))
plt.savefig('phasespace.pdf')
plt.show()
