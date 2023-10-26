#%matplotlib inline
import pylab as pl
from astropy.table import Table
from astropy.visualization import simple_norm
import numpy as np
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.colors

basepath = '/orange/adamginsburg/jwst/brick/'

mist = Table.read(f'{basepath}/isochrones/MIST_iso_633a08f2d8bb1.iso.cmd',
                  header_start=12, data_start=13, format='ascii', delimiter=' ', comment='#')
mist['410M405'] = mist['F410M']
mist['405M410'] = mist['F405N']

distance_modulus = 5*np.log10(8500)-5

trilegal = Table.read('/orange/adamginsburg/cmz/sgre/trilegal_nircam.dat', format='ascii.csv', delimiter=' ')

pl.figure(dpi=150)
        
pl.scatter(trilegal['F210M'] - trilegal['F480M'],
          trilegal['F480M'],
          color='red',
          alpha=0.1,
          s=1
          )
        
norm = simple_norm(mist['initial_mass'][mist['log10_isochrone_age_yr'] < 7], stretch='log')
for age in np.unique(mist['log10_isochrone_age_yr']):
    if age in (5,6,7):
    
        agesel = mist['log10_isochrone_age_yr'] == age
        pl.scatter(mist['F210M'][agesel] - mist['F480M'][agesel],
                   (mist['F480M'])[agesel]+distance_modulus,
                   c=mist['initial_mass'][agesel],
                   norm=norm,
                   cmap='inferno',
                   alpha=0.5,
                   s=1,
                  )
        
cb = pl.colorbar()
cb.set_ticks([0.5,1,5,10,50,250])
cb.set_label("Mass")

pl.axhline(0, linestyle='--', color='k')
pl.axvline(0, linestyle='--', color='k')
pl.xlabel("F210M-F480M")
pl.ylabel("F480M")
#pl.axis([-0.2,0.1,-0.1,0.3])
pl.ylim(25,5)
pl.xlim([-0.2, 0.3]);
#agesel.sum()
pl.savefig('isochrone_f210m-f480m.png')
