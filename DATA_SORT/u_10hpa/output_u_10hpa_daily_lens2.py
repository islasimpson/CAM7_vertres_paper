import xarray as xr
import numpy as np
from CASutils import readdata_utils as read
from CASutils import lensread_utils as lens
import sys

outpath="/project/cas/islas/python_savs/CAM7_vertres_paper/DATA_SORT/u_10hpa/"
basepath="/project/mojave/cesm2/LENS/atm/day_1/Uzm/"

mems = lens.lens2memnamegen_first50(50)

allmems=[]
for imem in mems:
    print(imem)
    dat = xr.open_mfdataset(basepath+"*-"+imem+"*.nc").isel(zlon=0)
    dat = read.fixcesmtime_daily(dat).Uzm
    dat = dat.sel(time=slice("1979-01-01","2023-12-31"))
    dat = dat.where( dat.time.dt.hour == 12, drop=True)
    dat = dat.interp(ilev=10)
    allmems.append(dat)
allmems = xr.concat(allmems, dim='M')
allmems.to_netcdf(outpath+'LENS2.nc')
