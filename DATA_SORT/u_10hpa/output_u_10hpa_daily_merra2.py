import xarray as xr
import numpy as np
from functools import partial

outpath="/project/cas/islas/python_savs/CAM7_vertres_paper/DATA_SORT/u_10hpa/"

basepath="/project/cas/islas/python_savs/CAM7_vertres_paper/RAW_DATA/MERRA2/Uzm/"

def preprocessor(ds):
    ds = ds.Uzm
    ds = ds.sel(lev=10)
    return ds

dat = xr.open_mfdataset(basepath+"*.nc", preprocess=partial(preprocessor))
dat.to_netcdf(outpath+"MERRA2.nc")
