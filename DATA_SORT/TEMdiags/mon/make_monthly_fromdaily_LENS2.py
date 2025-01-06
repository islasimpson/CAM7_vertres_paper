import xarray as xr
import numpy as np
import sys
from CASutils import lensread_utils as lens

basepath="/project/cas/islas/python_savs/L83_paper/DATA_SORT/TEMdiags/day/LENS2/"
pathout="/project/cas/islas/python_savs/L83_paper/DATA_SORT/TEMdiags/mon/LENS2/"

#mems=['001','002','003']
#mems=['001']
#mems=['002','003']
mems = lens.lens2memnamegen(100)
mems = mems[0:30]


for imem in mems:
    dat = xr.open_dataset(basepath+'L83_FHIST_'+imem+'.nc')
    monyearstr = xr.DataArray(dat.indexes['time'].strftime('%Y-%m'),
             coords=dat.time.coords, name='monyearstr')

    numdays = dat.uzm.isel(lat=0,ilev=0).groupby(monyearstr).count()

    dat_monthly = dat.groupby(monyearstr).mean('time')
#    dat_monthly = dat_monthly.where(numdays > 27, drop=True)

    timeout = dat.time.groupby(monyearstr).mean('time')
    timeout = timeout.where(numdays > 27, drop=True)
    timeout = timeout.rename(monyearstr='time')

    dat_monthly = dat_monthly.rename(monyearstr='time')
    dat_monthly['time'] = timeout

    dat_monthly.to_netcdf(pathout+'LENS2_'+imem+'.nc')

