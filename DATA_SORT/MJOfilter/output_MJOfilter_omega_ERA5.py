import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from math import nan
from CASutils import filter_utils as filt
from CASutils import linfit_utils as linfit
import pandas as pd

basepath="/project/cas/islas/python_savs/CAM7_vertres_paper/RAW_DATA/ERA5/W_dayfrom6h/"
pathout="/project/cas/islas/python_savs/CAM7_vertres_paper/DATA_SORT/MJOfilter/"

w = xr.open_mfdataset(basepath+"W_*.nc").W.sel(level=500)
w = w.where( ~( (w.time.dt.month ==2) & (w.time.dt.day == 29)), drop=True)

ystart=w.time.dt.year.isel(time=0) ; yend=w.time.dt.year.isel(time=w.time.size-1)

allmjo=[]
for iyear in np.arange(ystart,yend,1):
   wuse = w.sel(time=slice(str(iyear)+'-11-01',str(iyear+1)+'-03-31'))
   wuseclim = wuse.mean('time')
   wanoms = wuse - wuseclim
   wanoms['time'] = pd.date_range("1970-11-01","1971-03-31")
   w_mjo = filt.wkfilter(wanoms,0.15,1,5,20,100,spd=1)
   w_mjo = w_mjo.sel(time=slice("1970-12-01","1971-02-28"))
   allmjo.append(w_mjo)

year = xr.DataArray(np.arange(ystart+1,yend+1,1), coords=[np.arange(ystart+1,yend+1,1)],
                    dims=['year'], name='year')
allmjo = xr.concat(allmjo, dim=year)
allmjo = allmjo.rename('MJO_W')
allmjo.to_netcdf(pathout+'W_MJO_ERA5.nc')
