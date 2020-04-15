# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 09:33:15 2020

@author: austi
"""

import numpy as np
from matplotlib import pyplot as plt
from astroML.correlation import two_point_angular
from astropy.io import fits
import pandas as pd

if "setup_text_plots" not in globals():
    from astroML.plotting import setup_text_plots
setup_text_plots(fontsize=8, usetex=False)

### AllStar file and headers
allStar_dir = 'data/allStar.fits'
star_hdus = fits.open(allStar_dir)
print(star_hdus[1])
star_list = star_hdus[1].data
star_headers = star_hdus[1].header
star_hdus.close()

#%%

if "setup_text_plots" not in globals():
    from astroML.plotting import setup_text_plots
setup_text_plots(fontsize=8, usetex=False)

allstar = star_list
a = allstar['TEFF'] <= 3500
r = allstar[a]
b = r['TEFF'] >= 0
R = r[b]

c = allstar['TEFF'] >= 7500
B = allstar[c]

plt.plot(R['RA'], R['DEC'], 'o', color='r')
plt.plot(B['RA'], B['DEC'], 'o', color='b')
plt.xlabel('RA (degrees)')
plt.ylabel('DEC (degrees)')
plt.title('SDSS Stars')
plt.show()

print("SDSS data size:")
print("  red stars: ", len(R))
print("  blue stars:", len(B))

bins = np.linspace (0.1, 10, 11) # edges for the 10 bins to evaluate
corr_R = two_point_angular(R['RA'], R['DEC'], bins, method='landy-szalay')
corr_B = two_point_angular(B['RA'], B['DEC'], bins, method='landy-szalay')

bin_centers = 0.5 * (bins[1:] + bins[:-1])
plt.plot(bin_centers, corr_R, 'r')
plt.plot(bin_centers, corr_B, 'b')
plt.title('SDSS Stars Correlation')
plt.ylabel(r'$\hat{w}(\theta)$')
plt.xlabel(r'$\theta\ (deg)$')
plt.show() 

print('  Red Star Mean Correlation:', np.nanmean(corr_R))
print('  Blue Star Mean Correlation:', np.nanmean(corr_B))

corr = [corr_R, corr_B]
labels = ['$Red Stars$\n$N=%i$' % len(R),
          '$Blue Stars$\n$N=%i$' % len(B)]

fig = plt.figure(figsize=(5, 2.5))
fig.subplots_adjust(bottom=0.2, top=0.9,
                    left=0.13, right=0.95)

for i in range(2):
    ax = fig.add_subplot(121 + i, xscale='linear', yscale='linear')

    plt.plot(corr[i])

    ax.text(0.95, 0.95, labels[i],
            ha='right', va='top', transform=ax.transAxes)
    ax.set_xlabel(r'$\theta\ (deg)$')
    if i == 0:
        ax.set_ylabel(r'$\hat{w}(\theta)$')

plt.title('SDSS Stars Correlation')
plt.show() 

#%%
#Repeat the above process for Gaia data


R_gaia = pd.read_csv('data/gaia_red2.csv')
B_gaia = pd.read_csv('data/gaia_blue2.csv')

print("Gaia data size:")
print("  red stars: ", len(R_gaia))
print("  blue stars:", len(B_gaia))

plt.plot(R_gaia['ra'], R_gaia['dec'], 'o', color='r')
plt.plot(B_gaia['ra'], B_gaia['dec'], 'o', color='b')
plt.xlabel('RA (degrees)')
plt.ylabel('DEC (degrees)')
plt.title('Gaia Stars')
plt.show()

bins = np.linspace (0.1, 10, 11) # edges for the 10 bins to evaluate
corr_R_gaia = two_point_angular(R_gaia['ra'], R_gaia['dec'], bins, method='landy-szalay')
corr_B_gaia = two_point_angular(B_gaia['ra'], B_gaia['dec'], bins, method='landy-szalay')

bin_centers = 0.5 * (bins[1:] + bins[:-1])
plt.plot(bin_centers, corr_R_gaia, 'r')
plt.plot(bin_centers, corr_B_gaia, 'b')
plt.title('Gaia Stars Correlation')
plt.ylabel(r'$\hat{w}(\theta)$')
plt.xlabel(r'$\theta\ (deg)$')
plt.show() 

print('  Red Star Mean Correlation:', np.nanmean(corr_R_gaia))
print('  Blue Star Mean Correlation:', np.nanmean(corr_B_gaia))

corr = [corr_R_gaia, corr_B_gaia]
labels = ['$Red Stars$\n$N=%i$' % len(R_gaia),
          '$Blue Stars$\n$N=%i$' % len(B_gaia)]

fig = plt.figure(figsize=(5, 2.5))
fig.subplots_adjust(bottom=0.2, top=0.9,
                    left=0.13, right=0.95)

for i in range(2):
    ax = fig.add_subplot(121 + i, xscale='linear', yscale='linear')

    plt.plot(corr[i])

    ax.text(0.95, 0.95, labels[i],
            ha='right', va='top', transform=ax.transAxes)
    ax.set_xlabel(r'$\theta\ (deg)$')
    if i == 0:
        ax.set_ylabel(r'$\hat{w}(\theta)$')

plt.title('Gaia Stars Correlation')
plt.show() 

#%%

#Comparison in different regions using galactic coordinates

plt.plot(R['GLON'], R['GLAT'], 'o', color='r')
plt.plot(B['GLON'], B['GLAT'], 'o', color='b')
plt.xlabel('Galactic Longitude (degrees)')
plt.ylabel('Dalactic Latitude (degrees)')
plt.title('SDSS Stars')
plt.show()

plt.plot(R_gaia['l'], R_gaia['b'], 'o', color='r')
plt.plot(B_gaia['l'], B_gaia['b'], 'o', color='b')
plt.xlabel('Galactic Longitude (degrees)')
plt.ylabel('Galactic Latitude (degrees)')
plt.title('Gaia Stars')
plt.show()

#%%

#Galactic Disk Estimate from SDSS Galactic plot:
    
b_range = [-40,40]

#Getting SDSS disk Stars
R_disk = R[R['GLAT'] > b_range[0]]
R_disk = R_disk[R_disk['GLAT'] < b_range[1]]

B_disk = B[B['GLAT'] > b_range[0]]
B_disk = B_disk[B_disk['GLAT'] < b_range[1]]

#Getting Gaia Disk Stars
R_gaia_disk = R_gaia[R_gaia['b'] > b_range[0]]
R_gaia_disk = R_gaia_disk[R_gaia_disk['b'] < b_range[1]]

B_gaia_disk = B_gaia[B_gaia['b'] > b_range[0]]
B_gaia_disk = B_gaia_disk[B_gaia_disk['b'] < b_range[1]]

#%%

#SDSS disk

bins = np.linspace (0.1, 10, 11) # edges for the 10 bins to evaluate
corr_R_disk = two_point_angular(R_disk['GLON'], R_disk['GLAT'], bins, method='landy-szalay')
corr_B_disk = two_point_angular(B_disk['GLON'], B_disk['GLAT'], bins, method='landy-szalay')

bin_centers = 0.5 * (bins[1:] + bins[:-1])
plt.plot(bin_centers, corr_R_disk, 'r')
plt.plot(bin_centers, corr_B_disk, 'b')
plt.title('SDSS Disk Stars Correlation')
plt.ylabel(r'$\hat{w}(\theta)$')
plt.xlabel(r'$\theta\ (deg)$')
plt.show() 

print('SDSS Disk Stars')
print('  Red Star Mean Correlation:', np.nanmean(corr_R_disk))
print('  Blue Star Mean Correlation:', np.nanmean(corr_B_disk))

corr_disk = [corr_R_disk, corr_B_disk]
labels = ['$Red Stars$\n$N=%i$' % len(R_disk),
          '$Blue Stars$\n$N=%i$' % len(B_disk)]

fig = plt.figure(figsize=(5, 2.5))
fig.subplots_adjust(bottom=0.2, top=0.9,
                    left=0.13, right=0.95)

for i in range(2):
    ax = fig.add_subplot(121 + i, xscale='linear', yscale='linear')

    plt.plot(corr_disk[i])

    ax.text(0.95, 0.95, labels[i],
            ha='right', va='top', transform=ax.transAxes)
    ax.set_xlabel(r'$\theta\ (deg)$')
    if i == 0:
        ax.set_ylabel(r'$\hat{w}(\theta)$')

plt.title('SDSS Disk Stars Correlation')
plt.show() 

#Gaia Disk

bins = np.linspace (0.1, 10, 11) # edges for the 10 bins to evaluate
corr_R_gaia_disk = two_point_angular(R_gaia_disk['l'], R_gaia_disk['b'], bins, method='landy-szalay')
corr_B_gaia_disk = two_point_angular(B_gaia_disk['l'], B_gaia_disk['b'], bins, method='landy-szalay')

bin_centers = 0.5 * (bins[1:] + bins[:-1])
plt.plot(bin_centers, corr_R_gaia_disk, 'r')
plt.plot(bin_centers, corr_B_gaia_disk, 'b')
plt.title('Gaia Disk Stars Correlation')
plt.ylabel(r'$\hat{w}(\theta)$')
plt.xlabel(r'$\theta\ (deg)$')
plt.show() 

print('Gaia Bulge Stars')
print('  Red Star Mean Correlation:', np.nanmean(corr_R_gaia_disk))
print('  Blue Star Mean Correlation:', np.nanmean(corr_B_gaia_disk))

corr_gaia_disk = [corr_R_gaia_disk, corr_B_gaia_disk]
labels = ['$Red Stars$\n$N=%i$' % len(R_gaia_disk),
          '$Blue Stars$\n$N=%i$' % len(B_gaia_disk)]

fig = plt.figure(figsize=(5, 2.5))
fig.subplots_adjust(bottom=0.2, top=0.9,
                    left=0.13, right=0.95)

for i in range(2):
    ax = fig.add_subplot(121 + i, xscale='linear', yscale='linear')

    plt.plot(corr_gaia_disk[i])

    ax.text(0.95, 0.95, labels[i],
            ha='right', va='top', transform=ax.transAxes)
    ax.set_xlabel(r'$\theta\ (deg)$')
    if i == 0:
        ax.set_ylabel(r'$\hat{w}(\theta)$')

plt.title('Gaia Disk Stars Correlation')
plt.show() 

#%%

#Getting SDSS Halo Stars
R_ind = []
for i in range(len(R['GLAT'])):
    if R['GLAT'][i] < b_range[0] or R['GLAT'][i] > b_range[1]:
        R_ind.append(i)
R_halo = R[R_ind]

B_ind = []
for i in range(len(B['GLAT'])):
    if B['GLAT'][i] < b_range[0] or B['GLAT'][i] > b_range[1]:
        B_ind.append(i)
B_halo = B[B_ind]

#Getting Gaia Halo Stars
R_gaia_ind = []
for i in range(len(R_gaia['b'])):
    if R_gaia['b'][i] < b_range[0] or R_gaia['b'][i] > b_range[1]:
        R_gaia_ind.append(i)
R_gaia_halo_l = R_gaia['l'][R_gaia_ind]
R_gaia_halo_b = R_gaia['b'][R_gaia_ind]

B_gaia_ind = []
for i in range(len(B_gaia['b'])):
    if B_gaia['b'][i] < b_range[0] or B_gaia['b'][i] > b_range[1]:
        B_gaia_ind.append(i)
B_gaia_halo_l = B_gaia['l'][B_gaia_ind]
B_gaia_halo_b = B_gaia['b'][B_gaia_ind]

plt.plot(R_gaia_halo_l,R_gaia_halo_b, 'o', color='r')
plt.plot(B_gaia_halo_l,B_gaia_halo_b, 'o', color='b')

#%%

#SDSS Halo

bins = np.linspace (0.1, 10, 11) # edges for the 10 bins to evaluate
corr_R_halo = two_point_angular(R_halo['GLON'], R_halo['GLAT'], bins, method='landy-szalay')
corr_B_halo = two_point_angular(B_halo['GLON'], B_halo['GLAT'], bins, method='landy-szalay')

bin_centers = 0.5 * (bins[1:] + bins[:-1])
plt.plot(bin_centers, corr_R_halo, 'r')
plt.plot(bin_centers, corr_B_halo, 'b')
plt.title('SDSS Halo Stars Correlation')
plt.ylabel(r'$\hat{w}(\theta)$')
plt.xlabel(r'$\theta\ (deg)$')
plt.show() 

print('SDSS Halo Stars')
print('  Red Star Mean Correlation:', np.nanmean(corr_R_halo))
print('  Blue Star Mean Correlation:', np.nanmean(corr_B_halo))

corr_halo = [corr_R_halo, corr_B_halo]
labels = ['$Red Stars$\n$N=%i$' % len(R_halo),
          '$Blue Stars$\n$N=%i$' % len(B_halo)]

fig = plt.figure(figsize=(5, 2.5))
fig.subplots_adjust(bottom=0.2, top=0.9,
                    left=0.13, right=0.95)

for i in range(2):
    ax = fig.add_subplot(121 + i, xscale='linear', yscale='linear')

    plt.plot(corr_halo[i])

    ax.text(0.95, 0.95, labels[i],
            ha='right', va='top', transform=ax.transAxes)
    ax.set_xlabel(r'$\theta\ (deg)$')
    if i == 0:
        ax.set_ylabel(r'$\hat{w}(\theta)$')

plt.title('SDSS Halo Stars Correlation')
plt.show() 

#Gaia Halo

bins = np.linspace (0.1, 10, 11) # edges for the 10 bins to evaluate
corr_R_gaia_halo = two_point_angular(R_gaia_halo_l, R_gaia_halo_b, bins, method='landy-szalay')
corr_B_gaia_halo = two_point_angular(B_gaia_halo_l, B_gaia_halo_b, bins, method='landy-szalay')

bin_centers = 0.5 * (bins[1:] + bins[:-1])
plt.plot(bin_centers, corr_R_gaia_halo, 'r')
plt.plot(bin_centers, corr_B_gaia_halo, 'b')
plt.title('Gaia Halo Stars Correlation')
plt.ylabel(r'$\hat{w}(\theta)$')
plt.xlabel(r'$\theta\ (deg)$')
plt.show() 

print('Gaia Halo Stars')
print('  Red Star Mean Correlation:', np.nanmean(corr_R_gaia_halo))
print('  Blue Star Mean Correlation:', np.nanmean(corr_B_gaia_halo))

corr_gaia_halo = [corr_R_gaia_halo, corr_B_gaia_halo]
labels = ['$Red Stars$\n$N=%i$' % len(R_gaia_halo_l),
          '$Blue Stars$\n$N=%i$' % len(B_gaia_halo_l)]

fig = plt.figure(figsize=(5, 2.5))
fig.subplots_adjust(bottom=0.2, top=0.9,
                    left=0.13, right=0.95)

for i in range(2):
    ax = fig.add_subplot(121 + i, xscale='linear', yscale='linear')

    plt.plot(corr_gaia_halo[i]) 

    ax.text(0.95, 0.95, labels[i],
            ha='right', va='top', transform=ax.transAxes)
    ax.set_xlabel(r'$\theta\ (deg)$')
    if i == 0:
        ax.set_ylabel(r'$\hat{w}(\theta)$')

plt.title('Gaia Halo Stars Correlation')
plt.show() 


