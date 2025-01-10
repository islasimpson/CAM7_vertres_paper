import xarray as xr
import numpy as np
from CASutils import readdata_utils as read
from functools import partial

outpath="/project/cas/islas/python_savs/CAM7_vertres_paper/DATA_SORT/u_10hpa/"

basepath="/project/cas/islas/python_savs/CAM7_vertres_paper/DATA_SORT/TEMdiags/day/"

def preprocessor(ds):
    ds = ds.u
    ds = ds.sel(pre=10)
    ds = ds.mean('lon')
    return ds

dat = xr.open_mfdataset("/project/cas/islas/python_savs/CAM7_vertres_paper/RAW_DATA/JRA55/dayfrom6h/"+
           "uvtw*.nc", preprocess = partial(preprocessor))

dat = dat.drop_vars(['number','step','pre'])

dat.to_netcdf(outpath+'JRA55.nc')
