#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('ggplot')

def histogram(data, weight, nbins):
    length = data.shape[0]
    hist, bin_edges = np.histogram(data, weights=weight, bins=nbins, density=True)
    hists = []
    for i in range(10):
        hist_i, _ = np.histogram(data[int(i*length/10) : int((i+1)*length/10)], weights=weight[int(i*length/10) : int((i+1)*length/10)], bins=bin_edges, density=True)
        hists.append(hist_i)
    hists_np = np.array(hists)
    stddevs = np.std(hists_np, axis=0)
    return hist, (bin_edges[:-1]+bin_edges[1:])/2, stddevs

def calculate_hists(filename):
    """Calculates and saves the histograms with errorbars for PIMC data in file."""
    data = pd.read_csv(filename, sep='\t')
    for col in list(data)[1:]:
        print("histogramming "+col)
        hist, bins, stddevs = histogram(data[col].values, data["weight"].values, 100)
        np.savetxt(filename+'_'+col+'_dist', np.array((bins, hist, stddevs)).T)
        print(filename+'_'+col+'_dist written')

def plot(filename, **kwargs):
    data = np.loadtxt(filename)
    plt.errorbar(data[:,0], data[:,1], yerr=data[:,2], capsize=3, **kwargs)
