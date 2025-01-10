import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from math import nan
from CASutils import filter_utils as filt
from CASutils import linfit_utils as linfit
import pandas as pd

expname=['dz1000','dz900','dz800','dz700','dz600','dz500','dz400']

basepath="/project/cas/islas/python_savs/CAM7_vertres_paper/RAW_DATA/140km/"
pathout="/project/cas/islas/python_savs/CAM7_vertres_paper/DATA_SORT/MJOfilter/"


for iexp in expname:
    print(iexp)
    w = xr.open_mfdataset("/project/cas/islas/python_savs/CAM7_vertres_paper/RAW_DATA/140km/"+\
       "day/plev/day_from_6h/"+iexp+"/OMEGA_*.nc").OMEGA.sel(pre=500)

    ystart = w.time.dt.year.isel(time=0) ; yend = w.time.dt.year.isel(time=w.time.size-1)

    if (iexp == 'dz1000'):
        ystart=1988

    allmjo=[]
    for iyear in np.arange(ystart,yend,1):
        wuse = w.sel(time=slice(str(iyear)+'-11-01',str(iyear+1)+'-03-31'))
        wuseclim = wuse.mean('time')
        wanoms = wuse - wuseclim
        wanoms['time'] = pd.date_range("1970-11-01","1971-03-31")
        w_mjo = filt.wkfilter(wanoms, 0.15,1,5,20,100,spd=1)
        w_mjo = w_mjo.sel(time=slice("1970-12-01","1971-02-28"))
        allmjo.append(w_mjo)

    year = xr.DataArray(np.arange(ystart+1,yend+1,1), coords=[np.arange(ystart+1,yend+1,1)],
                   dims=['year'], name='year')
    allmjo = xr.concat(allmjo, dim=year)
    allmjo = allmjo.rename('MJO_W')
    allmjo.to_netcdf(pathout+'W_mjo_140km_'+iexp+'.nc')
