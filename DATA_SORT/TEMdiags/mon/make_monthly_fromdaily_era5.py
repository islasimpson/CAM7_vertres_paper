import xarray as xr
import numpy as np
import sys

basepath="/project/cas/islas/python_savs/CAM7_vertres_paper/DATA_SORT/TEMdiags/day/"
pathout="/project/cas/islas/python_savs/CAM7_vertres_paper/DATA_SORT/TEMdiags/mon/"

dat = xr.open_dataset(basepath+'ERA5.nc')

monyearstr = xr.DataArray(dat.indexes['time'].strftime('%Y-%m'),
         coords=dat.time.coords, name='monyearstr')

numdays = dat.uzm.isel(lat=0,plev=0).groupby(monyearstr).count()

dat_monthly = dat.groupby(monyearstr).mean('time')
#dat_monthly = dat_monthly.where(numdays > 27, drop=True)

timeout = dat.time.groupby(monyearstr).mean('time')
timeout = timeout.where(numdays > 27, drop=True)
timeout = timeout.rename(monyearstr='time')

dat_monthly = dat_monthly.rename(monyearstr='time')
dat_monthly['time'] = timeout

dat_monthly.to_netcdf(pathout+'ERA5.nc')
