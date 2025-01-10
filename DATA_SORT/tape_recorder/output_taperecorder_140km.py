import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from math import nan
from CASutils import readdata_utils as read
from CASutils import averaging_utils as avg
import sys

exps=['dz1000','dz900','dz800','dz700','dz600','dz500','dz400']
basepath="/project/cas/islas/python_savs/CAM7_vertres_paper/RAW_DATA/140km/mon/"
pathout="/project/cas/islas/python_savs/CAM7_vertres_paper/DATA_SORT/tape_recorder/"

for iexp in exps:
    print(iexp)
    dat = xr.open_mfdataset(basepath+iexp+'/Q_*.nc')
    dat = read.fixcesmtime(dat)
    dat_tr = avg.cosweightlat(dat,-5,5)

    ystart = dat_tr.time.dt.year.isel(time=0).values
    yend = dat_tr.time.dt.year.isel(time=dat_tr.time.size-1).values

    dat_tr = dat_tr.Q.sel(time=slice(str(ystart+5)+'-01-01',str(yend)+'-12-31'))

    dat_tr = dat_tr.groupby('time.month').mean('time')
    dat_tr = dat_tr.mean('lon')

    dat_tr.to_netcdf(pathout+iexp+'_tropical_q.nc')
