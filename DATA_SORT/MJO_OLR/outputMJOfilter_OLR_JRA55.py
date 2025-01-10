import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from math import nan
from CASutils import filter_utils as filt
from CASutils import linfit_utils as linfit
from CASutils import averaging_utils as avg
import sys
import pandas as pd

basepath="/project/mojave/observations/JRA55/day/OLR/"
ubasepath="/project/cas/islas/python_savs/CAM7_vertres_paper/RAW_DATA/JRA55/Uzm/"
pathout="/project/cas/islas/python_savs/CAM7_vertres_paper/DATA_SORT/MJO_OLR/"

dat = xr.open_mfdataset(basepath+"*.nc").OLR
dat = dat.where( ~(( dat.time.dt.month == 2) & (dat.time.dt.day == 29)), drop=True)
dat.time.encoding['calendar'] = 'noleap'
dat = dat.sel(time=slice("1979-01-01","2023-12-31"))

udat = xr.open_mfdataset(ubasepath+"*.nc")
monstr = xr.DataArray(udat.indexes['time'].strftime('%Y-%m'), coords = udat.time.coords, name='monstr')
udat = udat.groupby(monstr).mean('time')
udat = udat.rename(monstr='time')
udat['time'] = pd.date_range("1979-01","2024-01", freq='M')
udat = avg.cosweightlat(udat, -5, 5)
udat = udat.interp(pre=50.)


ystart=1979
yend=2023
mjofilt=[]
uqbo=[]
for iyear in np.arange(ystart,yend,1):
    print(iyear)

    udatuse = udat.sel(time=slice(str(iyear)+"-12-01", str(iyear+1)+"-02-28"))
    udatout = udatuse.mean('time').Uzm

    datuse = dat.sel(time=slice(str(iyear)+"-11-01",str(iyear+1)+"-04-30"))
    datclim = datuse.mean('time')
    datanoms = datuse - datclim
    timenew = pd.date_range("1970-11-01","1971-04-30")
    datanoms['time'] = timenew

    mjofilt.append(filt.wkfilter(datanoms,0.15,1,5,20,100, spd=1))
    uqbo.append(udatout)

uqbo = xr.concat(uqbo, dim='year')
uqbo['year'] = np.arange(ystart,yend,1)
uqbo = uqbo.rename('UQBO')

mjofilt = xr.concat(mjofilt, dim='year')
mjofilt['year'] = np.arange(ystart,yend,1)
mjofilt = mjofilt.sel(time=slice("1970-12-01","1971-02-28"))
mjofilt = mjofilt.rename('MJO_OLR')

mjofilt.to_netcdf(pathout+"OLR_mjo_JRA55.nc")
uqbo.to_netcdf(pathout+"OLR_mjo_JRA55.nc", mode="a")
