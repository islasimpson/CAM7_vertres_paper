import xarray as xr
import numpy as np
import sys

basepath="/project/cas/islas/python_savs/CAM7_vertres_paper/DATA_SORT/TEMdiags/day/"
pathout="/project/cas/islas/python_savs/CAM7_vertres_paper/DATA_SORT/TEMdiags/mon/"

#mems=['001','002','003']
#mems=['001']
mems=['001','002','003']

for imem in mems:
    dat = xr.open_dataset(basepath+'BWHIST_'+imem+'.nc')
    monyearstr = xr.DataArray(dat.indexes['time'].strftime('%Y-%m'),
             coords=dat.time.coords, name='monyearstr')

    numdays = dat.uzm.isel(lat=0,ilev=0).groupby(monyearstr).count()
    numdays = numdays.drop(['lat','zlon','ilev'])

    dat_monthly = dat.groupby(monyearstr).mean('time')
    dat_monthly = dat_monthly.where(numdays > 27, drop=True)

    timeout = dat.time.groupby(monyearstr).mean('time')
    timeout = timeout.where(numdays > 27, drop=True)
    timeout = timeout.rename(monyearstr='time')

    dat_monthly = dat_monthly.rename(monyearstr='time')
    dat_monthly['time'] = timeout

    dat_monthly.to_netcdf(pathout+'BWHIST_'+imem+'.nc')

