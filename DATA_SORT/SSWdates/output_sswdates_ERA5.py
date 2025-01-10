import xarray as xr
import numpy as np

from CASutils import ssw_utils as ssw
from CASutils import readdata_utils as read

pathout="/project/cas/islas/python_savs/CAM7_vertres_paper/DATA_SORT/SSWdates/"
basepath="/project/cas/islas/python_savs/CAM7_vertres_paper/DATA_SORT/TEMdiags/day/"

dat = xr.open_dataset(basepath+'ERA5.nc')
dat['plev'] = dat.plev/100.
dat = dat.uzm
dat = dat.interp(lat=60.)
dat = dat.interp(plev=10.)
dat = dat.load()

sswdates, nwinters = ssw.ssw_cp(dat)
nssw = len(sswdates)
dateout=[]
for i in np.arange(0,nssw,1):
    dateout.append(sswdates[i].values)

dateout_xr = xr.DataArray(dateout, dims=['ssw'], name='sswdate')
nwinters_xr = xr.DataArray(nwinters, name='nwinters')

dateout_xr.load().to_netcdf(pathout+'sswdates_ERA5.nc')
dat.load().to_netcdf(pathout+'sswdates_ERA5.nc', mode='a')
nwinters_xr.load().to_netcdf(pathout+'sswdates_ERA5.nc', mode='a') 
