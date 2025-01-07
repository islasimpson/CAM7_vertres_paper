import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from math import nan
from CASutils import filter_utils as filt
from CASutils import linfit_utils as linfit
from CASutils import nsq_vortgradients_utils as vortgrad
import sys
import pandas as pd

basepath="/project/cas/islas/python_savs/CAM7_vertres_paper/RAW_DATA/ERA5/U_dayfrom6h/"
pathout="/project/cas/islas/python_savs/CAM7_vertres_paper/DATA_SORT/MRG_filter/"

dat = xr.open_mfdataset(basepath+"/*.nc").U.sel(level=50)
dat = dat.where( ~( (dat.time.dt.month == 2) & (dat.time.dt.day == 29)), drop=True)
dat.time.encoding['calendar'] = 'noleap'
dat = dat.sel(time=slice("1979-01-01","2023-12-31"))

dat = dat.sel(lat=slice(-15,15))

dat_running_100day = dat.rolling(time=100, min_periods=100, center=True)
dat_running_100day = dat_running_100day.construct('segtime')

dat_running_100day_50sep = dat_running_100day.isel(time=slice(0,None,50))

um=[]
mrgfilt=[]
vortgradneg=[]
for isegment in np.arange(0,dat_running_100day_50sep.time.size,1):
    print(isegment)
    datuse = dat_running_100day_50sep.isel(time=isegment)
    datusem = datuse.mean('segtime')
    datusezm = datuse.mean('lon')
    gradvals = vortgrad.barotropicvortgradient(datusezm)    
    neggrad = xr.where(gradvals < 0, 1, 0)
    probneggrad = (neggrad.sum('segtime')/neggrad.segtime.size)*100.
    probneggrad = probneggrad.assign_coords({'time':datuse.time.values}) 
    vortgradneg.append(probneggrad)

    datuseanoms = datuse - datusem
    datuseanoms = datuseanoms.rename(segtime='time')
    datuseanoms['time'] = np.arange(0,datuse.segtime.size,1)
    datuseanoms = datuseanoms.transpose('time','lat','lon')
    mrgfilt_temp = filt.wkfilter(datuseanoms,0.1,-6,3,2.5,5,spd=1)
    mrgfilt_temp = mrgfilt_temp.isel(time=slice(10,mrgfilt_temp.time.size-10))
    mrgfilt_temp = mrgfilt_temp.std('time')
    mrgfilt_temp = mrgfilt_temp.assign_coords({'time':datusem.time.values})    
    mrgfilt.append(mrgfilt_temp)
    um.append(datusem)

mrgfilt = xr.concat(mrgfilt, dim='time').rename('U_mrgfilt')
um = xr.concat(um, dim='time').rename('U')
vortgradneg = xr.concat(vortgradneg, dim='time').rename('prob_vortgradneg')

mrgfilt.to_netcdf(pathout+'U50_MRGfilt_100day_segments_ERA5.nc')
um.to_netcdf(pathout+'U50_MRGfilt_100day_segments_ERA5.nc', mode='a')
vortgradneg.to_netcdf(pathout+'U50_MRGfilt_100day_segments_ERA5.nc', mode='a')
