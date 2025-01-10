import xarray as xr
import numpy as np
import pdb
import sys
import glob

from CASutils import ssw_utils as ssw
from CASutils import readdata_utils as read
from CASutils import lensread_utils as lens

pathout="/project/cas/islas/python_savs/CAM7_vertres_paper/DATA_SORT/SSWdates/LENS2/"

memstr = lens.lens2memnamegen(100)

for imem in memstr:
    print(imem)

    hist = xr.open_mfdataset(
     "/project/mojave/cesm2/LENS/atm/day_1/Uzm/"
     +"*.BHIST*"+imem+".*.nc")
    hist = read.fixcesmtime(hist)
    hist = hist.where( hist.time.dt.hour == 12, drop=True)

    ssp = xr.open_mfdataset(
      "/project/mojave/cesm2/LENS/atm/day_1/Uzm/"
 +"*.BSSP370*"+imem+"*.nc")
    ssp = read.fixcesmtime(ssp)
    ssp = ssp.where( ssp.time.dt.hour == 12, drop=True)

    dat = xr.concat([hist, ssp], dim='time')

    dat = dat.sel(time=slice("1979-01-01","2023-12-31"))

    dat = dat.Uzm
    dat = dat.isel(zlon=0)
    dat = dat.interp(lat=60.)
    dat = dat.interp(ilev=10.)
    dat = dat.load()

    sswdates, nwinters = ssw.ssw_cp(dat) 
    nssw = len(sswdates)
    dateout = []
    for i in np.arange(0,nssw,1):
        dateout.append(sswdates[i].values)

    dateout_xr = xr.DataArray(dateout, dims=['ssw'], name='sswdate')
    nwinters_xr = xr.DataArray(nwinters, name='nwinters')

    dateout_xr.load().to_netcdf(pathout+'sswdates_LENS2_'+imem+'_1979_2023.nc')
    dat.load().to_netcdf(pathout+'sswdates_LENS2_'+imem+'_1979_2023.nc', mode='a')
    nwinters_xr.load().to_netcdf(pathout+'sswdates_LENS2_'+imem+'_1979_2023.nc', mode='a')

