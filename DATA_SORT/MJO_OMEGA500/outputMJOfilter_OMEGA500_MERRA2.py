import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from math import nan
from CASutils import filter_utils as filt
from CASutils import linfit_utils as linfit
from CASutils import averaging_utils as avg
import pandas as pd

basepath="/project/cas/islas/python_savs/CAM7_vertres_paper/RAW_DATA/MERRA2/OMEGA500/"
pathout="/project/cas/islas/python_savs/CAM7_vertres_paper/DATA_SORT/MJO_OMEGA500/"

w = xr.open_mfdataset(basepath+"omega500_*.nc")
w = w.where( ~( (w.time.dt.month == 2) & (w.time.dt.day == 29)), drop=True)

ystart=w.time.dt.year.isel(time=0) ; yend=w.time.dt.year.isel(time=w.time.size-1)

udat = xr.open_mfdataset("/project/cas/islas/python_savs/L83_paper/RAW_DATA/MERRA2/Uzm/"+
            "Uzm_*.nc")
monstr = xr.DataArray(udat.indexes['time'].strftime('%Y-%m'), coords = udat.time.coords, name='monstr')
udat = udat.groupby(monstr).mean('time')
udat = udat.rename(monstr='time')
udat['time'] = pd.date_range(str(ystart.values)+"-01",str(yend.values+1)+"-01", freq='M')
udat = avg.cosweightlat(udat, -5, 5)
udat = udat.interp(lev=50.)


allmjo=[]
uqbo=[]
for iyear in np.arange(ystart,yend,1):

   udatuse = udat.sel(time=slice(str(iyear)+"-12-01", str(iyear+1)+"-02-28"))
   udatout = udatuse.mean('time').Uzm

   wuse = w.sel(time=slice(str(iyear)+'-11-01',str(iyear+1)+'-03-31')).w
   wuseclim = wuse.mean('time')
   wanoms = wuse - wuseclim
   wanoms['time'] = pd.date_range("1970-11-01","1971-03-31")
   w_mjo = filt.wkfilter(wanoms,0.15,1,5,20,100,spd=1)
   w_mjo = w_mjo.sel(time=slice("1970-12-01","1971-02-28"))
   allmjo.append(w_mjo)
   uqbo.append(udatout)

uqbo = xr.concat(uqbo, dim='year')
uqbo['year'] = np.arange(ystart,yend,1)
uqbo = uqbo.rename('UQBO')

year = xr.DataArray(np.arange(ystart+1,yend+1,1), coords=[np.arange(ystart+1,yend+1,1)],
                    dims=['year'], name='year')
allmjo = xr.concat(allmjo, dim=year)
allmjo = allmjo.rename('MJO_W')
allmjo.to_netcdf(pathout+'W_MJO_MERRA2.nc')
uqbo.to_netcdf(pathout+"W_MJO_MERRA2.nc", mode='a')


