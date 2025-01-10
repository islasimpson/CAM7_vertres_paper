import xarray as xr
import numpy as np

from CASutils import ssw_utils as ssw
from CASutils import readdata_utils as read

pathout="/project/cas/islas/python_savs/CAM7_vertres_paper/DATA_SORT/SSWdates/"
basepath="/project/mojave/cesm2/L83/"

for imem in np.arange(3,4,1):
    memstr = str(imem).zfill(3)

    hist = xr.open_mfdataset(basepath+'b.e21.BHISTcmip6.f09_g17.L83_cam6.'+memstr+'/atm/day_1/Uzm/*.nc')
    hist = read.fixcesmtime(hist)
    hist = hist.where( hist.time.dt.hour == 12, drop=True)

    ssp = xr.open_mfdataset(basepath+'b.e21.BSSP370cmip6.f09_g17.L83_cam6.'+memstr+'/atm/day_1/Uzm/*.nc')
    ssp['lat'] = hist.lat.values
    ssp = read.fixcesmtime(ssp)
    ssp = ssp.where( ssp.time.dt.hour == 12, drop=True)

    dat = xr.concat([hist, ssp], dim='time')

    dat = dat.sel(time=slice("1979-01-01","2023-12-31"))

    dat = dat.Uzm
    dat = dat.interp(lat=60.)
    dat = dat.interp(ilev=10.)
    dat = dat.isel(zlon=0)

    sswdates, nwinters = ssw.ssw_cp(dat)
    nssw = len(sswdates)
    dateout=[]
    for i in np.arange(0,nssw,1):
        dateout.append(sswdates[i].values)

    dateout_xr = xr.DataArray(dateout, dims=['ssw'], name='sswdate')
    nwinters_xr = xr.DataArray(nwinters, name='nwinters')

    dateout_xr.to_netcdf(pathout+'sswdates_L83_'+memstr+'_1979_2023.nc')
    dat.to_netcdf(pathout+'sswdates_L83_'+memstr+'_1979_2023.nc', mode='a')
    nwinters_xr.to_netcdf(pathout+'sswdates_L83_'+memstr+'_1979_2023.nc', mode='a')


 