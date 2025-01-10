import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from math import nan
import sys
import pandas as pd
import os

basepath="/project/cas/islas/python_savs/CAM7_vertres_paper/RAW_DATA/MERRA2/OMEGA500/"

ystart=1981 ; yend=2023
for iyear in np.arange(ystart,yend+1,1):
    for imon in np.arange(1,12+1,1):
        filename=basepath+"omega500_merra2_dayfrom6h_"+str(iyear)+str(imon).zfill(2)+".nc"
        dat = xr.open_dataset(filename)
        startdate=str(iyear)+"-"+str(imon).zfill(2)+"-01"
        enddate = pd.Timestamp(startdate) + pd.offsets.MonthEnd(1)
        time = pd.date_range(startdate, enddate, freq='D')

        dat['day'] = time
        dat = dat.rename({"day":"time"})
        os.remove(filename)
        dat.to_netcdf(filename)
