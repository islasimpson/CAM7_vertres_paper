import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

from CASutils import averaging_utils as avg

pathout="/project/cas/islas/python_savs/L83_paper/DATA_SORT/tape_recorder/"

dat = xr.open_dataset("/project/mojave/observations/SWOOSH/"+\
"swoosh-v02.6-198401-202312-latpress-2.5deg-L31.nc")

# using combined H2O (equivalent latitude filling)
dat = dat.combinedeqfillh2oq

# convert ppmv to kg/kg
dat = dat*18.015280/(1e6*28.964)

# selecting 2005 since there were a lot of gaps in the data prior to that 
# probably best once MLS is in there
# end in 2021 to avoid the Hunga Tonga eruption in January 2022
dat = dat.sel(time=slice("2005-01-01","2021-12-31"))

dat_tr = avg.cosweightlat(dat, -5, 5)
dat_tr = dat_tr.groupby('time.month').mean('time')
dat_tr = dat_tr.rename(level='pre')
dat_tr = dat_tr.rename('Q')

dat_tr.to_netcdf(pathout+'SWOOSH_tropical_q_2005_2021.nc')
