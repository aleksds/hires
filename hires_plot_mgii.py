# Aleks Diamond-Stanic
# 20160517 --> 20160926
#
# the goals of this code include the following:
#
# (1) make plots that visualize the spectra in the wavelength region
#     surrounding Mg II
#
# (2) show the line profiles for the 2796 and 2803 transitions
#     separately for each of 14 galaxies
#
# (3) visualize the velocity structure of the gas traced by MgII using
#     the 2796 lines for the highest velocities and the 2803 for
#     velocities closer to v=0

import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy import units as u
from astropy import constants as const
from astropy.io import ascii
from matplotlib.ticker import AutoMinorLocator

# define the data directory
dir = os.environ['HIRESDIR']

# speed of light in Angstroms per second
c = const.c.to('AA/s')

# wavelengths of relevant absorption lines
mgi2852 = 2852.96328 * u.AA
mgii2803 = 2803.5314853 * u.AA
mgii2796 = 2796.3542699 * u.AA
feii2600 = 2600.1724835 * u.AA
feii2586 = 2586.6495659 * u.AA
feii2382 = 2382.7641781 * u.AA
feii2374 = 2374.4603294 * u.AA
feii2344 = 2344.2129601 * u.AA

# galaxy names
gal = ['J0826', 'J0901', 'J0905', 'J0944', 'J1125', 'J1219', 'J1232', 'J1341', 'J1450', 'J1506', 'J1558', 'J1613', 'J2116', 'J2140']

# galaxy redshifts
zem = np.array([0.603, 0.459, 0.712, 0.514, 0.519, 0.451, 0.400, 0.658, 0.782, 0.608, 0.402, 0.449, 0.728, 0.752])

# approximate centroid velocities 
vcen = np.array([-1250., -1330., -2570., -1240., -2130., -1950., -550., -390., -1790., -1110., -900., -2450., -1480., -770.])

# approximate maximum velocities
vmax = np.array([-1600., -1700., -3000., -2100., -2250., -2240., -750., -1600., -2100., -2090., -1400., -2650., -1800., -1150.])

# sort by maximum outflow velocity
vmax_index = np.flipud(np.argsort(vmax))

# define velocity ranges to plot the profile
vra = np.zeros((len(gal), 3)) * u.km / u.s
# J0826 
vra[0] = [-3500., -1500., 1500.] * u.km / u.s
# J0901 
vra[1] = [-3500., -1700., 1500.] * u.km / u.s
# J0905 
vra[2] = [-3500., -3100., 1500.] * u.km / u.s
# J0944 
vra[3] = [-3500., -1400., 1500.] * u.km / u.s
# J1125
vra[4] = [-3500., -2300., 1500.] * u.km / u.s
# J1219 
vra[5] = [-3500., -2100., 1500.] * u.km / u.s
# J1232 
vra[6] = [-3500., -400., 1500.] * u.km / u.s
# J1341 
vra[7] = [-3500., -400., 1500.] * u.km / u.s
# J1450 
vra[8] = [-3500., -2200., 1500.] * u.km / u.s
# J1506
vra[9] = [-3500., -1600., 1500.] * u.km / u.s
 # J1558
vra[10] = [-3500., -1200., 1500.] * u.km / u.s
# J1613 
vra[11] = [-3500., -2000., 1500.] * u.km / u.s
# J2116 
vra[12] = [-3500., -1600., 1500.] * u.km / u.s
# J2140 
vra[13] = [-3500., -850., 1500.] * u.km / u.s

# make plot showing one galaxy per page
filename = 'all_mgii.pdf'
xls = 7.
minorLocator = AutoMinorLocator()

with PdfPages(filename) as pdf:

    for i in range(0, len(gal)):

        # read in the spectrum
        datafile = dir+gal[i]+'/'+gal[i]+'_stitched_v1.txt'
        print(datafile)
        data = ascii.read(datafile)
        wave = data['wv'] * u.AA
        flux = data['norm']
        xmin = np.array(-3500. * u.km * u.s)
        xmax = np.array(500. * u.km * u.s)

        fig = plt.figure()
        plt.suptitle(gal[i])

        # define the 2803 and 2796 velocity scales
        vel_mgii_2803_aa = (wave - mgii2803*(1+zem[i])) / (mgii2803*(1+zem[i])) * c
        vel_mgii_2803 = vel_mgii_2803_aa.to('km/s')
        vel_mgii_2796_aa = (wave - mgii2796*(1+zem[i])) / (mgii2796*(1+zem[i])) * c
        vel_mgii_2796 = vel_mgii_2796_aa.to('km/s')

        # show the profile in 2796 velocity units
        ax = fig.add_subplot(4,1,1)
        ax.set_xlim(xmin, xmax)
        ax.xaxis.set_minor_locator(minorLocator)
        ax.set_ylim(0., 1.5)
        ax.tick_params(axis='x', labelsize=xls)
        ax.plot(vel_mgii_2796, flux)
        plt.text(xmin+0.03*(xmax-xmin), 0.15, 'MgII 2796', color='blue')
        plt.axvline(x=vcen[i], ymin=0., ymax = 1.5, linewidth=1, color='k', linestyle='dotted')
        plt.axvline(x=vcen[i]+30., ymin=0., ymax = 1.5, linewidth=0.5, color='k')
        plt.axvline(x=vcen[i]-30., ymin=0., ymax = 1.5, linewidth=0.5, color='k')

        # show the profile in 2803 velocity units
        ax = fig.add_subplot(4,1,2)
        ax.set_xlim(xmin, xmax)
        ax.xaxis.set_minor_locator(minorLocator)
        ax.set_ylim(0., 1.5)
        ax.tick_params(axis='x', labelsize=xls)
        ax.plot(vel_mgii_2803, flux, color='red')
        plt.text(xmin+0.03*(xmax-xmin), 0.15, 'MgII 2803', color='black')
        plt.axvline(x=vcen[i], ymin=0., ymax = 1.5, linewidth=1, color='k', linestyle='dotted')
        plt.axvline(x=vcen[i]+30., ymin=0., ymax = 1.5, linewidth=0.5, color='k')
        plt.axvline(x=vcen[i]-30., ymin=0., ymax = 1.5, linewidth=0.5, color='k')

        # show the profile in both 2796 and 2803 velocity units
        ax = fig.add_subplot(4,1,3)
        ax.set_xlim(xmin, xmax)
        ax.xaxis.set_minor_locator(minorLocator)
        ax.set_ylim(0., 1.5)
        ax.tick_params(axis='x', labelsize=xls)
        ax.plot(vel_mgii_2803, flux, color='red')
        ax.plot(vel_mgii_2796, flux, color='blue')
        plt.text(xmin+0.03*(xmax-xmin), 0.15, 'MgII 2796+2803')
        plt.axvline(x=vcen[i], ymin=0., ymax = 1.5, linewidth=1, color='k', linestyle='dotted')
        plt.axvline(x=vcen[i]+30., ymin=0., ymax = 1.5, linewidth=0.5, color='k')
        plt.axvline(x=vcen[i]-30., ymin=0., ymax = 1.5, linewidth=0.5, color='k')

        # show the velocity profile using the 2796 line on the blue
        # side and the 2803 line on the red side
        g2796 = (vel_mgii_2796 > vra[i,0]) & (vel_mgii_2796 < vra[i,1])
        g2803 = (vel_mgii_2803 > vra[i,1]) & (vel_mgii_2803 < vra[i,2])
        ax = fig.add_subplot(4,1,4)
        ax.set_xlim(xmin, xmax)
        ax.xaxis.set_minor_locator(minorLocator)
        ax.set_ylim(0., 1.5)
        ax.tick_params(axis='x', labelsize=xls)
        ax.plot(vel_mgii_2803[g2803], flux[g2803], color='red')
        ax.plot(vel_mgii_2796[g2796], flux[g2796], color='blue')
        plt.text(xmin+0.03*(xmax-xmin), 0.15, 'MgII 2796+2803')
        plt.axvline(x=vcen[i], ymin=0., ymax = 1.5, linewidth=1, color='k', linestyle='dotted')
        plt.axvline(x=vcen[i]+30., ymin=0., ymax = 1.5, linewidth=0.5, color='k')
        plt.axvline(x=vcen[i]-30., ymin=0., ymax = 1.5, linewidth=0.5, color='k')
        

        pdf.savefig()
        plt.close()
    
    os.system("open %s &" % filename)


# make plot showing all galaxies on each page
filename = 'all_mgii_2page.pdf'
xls = 8.
yls = 8.


with PdfPages(filename) as pdf:

    fig = plt.figure()
    
    for i in range(0, len(gal)):

        # read in the spectrum, sorted by outflow velocity
        indx = vmax_index[i]
        datafile = dir+gal[indx]+'/'+gal[indx]+'_stitched_v1.txt'
        print(datafile)
        data = ascii.read(datafile)
        wave = data['wv'] * u.AA
        flux = data['norm']
        xmin = np.array(-3500. * u.km * u.s)
        xmax = np.array(1500. * u.km * u.s)

        # define the 2796 velocity scale
        vel_mgii_2796_aa = (wave - mgii2796*(1+zem[indx])) / (mgii2796*(1+zem[indx])) * c
        vel_mgii_2796 = vel_mgii_2796_aa.to('km/s')

        # plot the profiles using the 2796 velocity scale
        ax = fig.add_subplot(5,3,i+1)
        ax.set_xlim(xmin, xmax)
        ax.xaxis.set_minor_locator(minorLocator)
        ax.set_ylim(0., 1.5)
        ax.tick_params(axis='x', labelsize=xls)
        ax.tick_params(axis='y', labelsize=yls)
        ax.plot(vel_mgii_2796, flux, color='black', linewidth=0.5)
        plt.text(xmin+0.03*(xmax-xmin), 0.15, gal[indx])

    pdf.savefig()
    plt.close()
    
    fig = plt.figure()
        
    for i in range(0, len(gal)):

        # read in the spectrum, sorted by outflow velocity
        indx = vmax_index[i]
        datafile = dir+gal[indx]+'/'+gal[indx]+'_stitched_v1.txt'
        print(datafile)
        data = ascii.read(datafile)
        wave = data['wv'] * u.AA
        flux = data['norm']
        xmin = np.array(-3500. * u.km * u.s)
        xmax = np.array(1500. * u.km * u.s)

        # define the 2803 and 2796 velocity scale
        vel_mgii_2803_aa = (wave - mgii2803*(1+zem[indx])) / (mgii2803*(1+zem[indx])) * c
        vel_mgii_2803 = vel_mgii_2803_aa.to('km/s')
        vel_mgii_2796_aa = (wave - mgii2796*(1+zem[indx])) / (mgii2796*(1+zem[indx])) * c
        vel_mgii_2796 = vel_mgii_2796_aa.to('km/s')

        # define the regions to use the 2796 profile and the regions to use the 2803 profile
        g2796 = (vel_mgii_2796 > vra[indx,0]) & (vel_mgii_2796 < vra[indx,1])
        g2803 = (vel_mgii_2803 > vra[indx,1]) & (vel_mgii_2803 < vra[indx,2])


        # plot the profiles using the 2796 profile on the blue side
        # and the 2803 profile on the red side
        ax = fig.add_subplot(5,3,i+1)
        ax.set_xlim(xmin, xmax)
        ax.xaxis.set_minor_locator(minorLocator)
        ax.set_ylim(0., 1.5)
        ax.tick_params(axis='x', labelsize=xls)
        ax.tick_params(axis='y', labelsize=yls)
        ax.plot(vel_mgii_2803[g2803], flux[g2803], color='red', linewidth=0.5)
        ax.plot(vel_mgii_2796[g2796], flux[g2796], color='blue', linewidth=0.5)
        plt.text(xmin+0.03*(xmax-xmin), 0.15, gal[indx])
        
    pdf.savefig()
    plt.close()
    
    os.system("open %s &" % filename)
