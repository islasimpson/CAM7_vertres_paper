import xarray as xr
import numpy as np
from CASutils import readdata_utils as read

outpath="/project/cas/islas/python_savs/CAM7_vertres_paper/DATA_SORT/u_10hpa/"

basepath="/project/cas/islas/python_savs/CAM7_vertres_paper/DATA_SORT/TEMdiags/day/"

mems=['001','002','003']
basepath="/project/cas/islas/python_savs/L83_paper/RAW_DATA/L83_FHIST_BGC/"
for imem in mems:
    dat = xr.open_mfdataset(basepath+imem+'/day_1/*.Uzm.*.nc')
    dat = read.fixcesmtime_daily(dat).Uzm
    dat = dat.sel(time=slice("1979-01-01","2014-12-31"))
    dat = dat.interp(ilev=10)
    dat.to_netcdf(outpath+'L83_Fcase_'+imem+'.nc')

