import xarray as xr
import numpy as np
import sys

basepath="/project/cas/islas/python_savs/CAM7_vertres_paper/DATA_SORT/TEMdiags/day/"
pathout="/project/cas/islas/python_savs/CAM7_vertres_paper/DATA_SORT/TEMdiags/mon/"

#fnames=['dz1000_140km','dz900_140km','dz800_140km','dz700_140km','dz600_140km',
#        'dz500_140km','dz400_140km']

#fnames=['dz800_80km','dz700_80km','dz600_80km','dz500_80km']
fnames=['dz500_buggy_80km']



for fname in fnames:
    dat = xr.open_dataset(basepath+fname+'.nc')

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

    dat_monthly.to_netcdf(pathout+fname+'.nc')
