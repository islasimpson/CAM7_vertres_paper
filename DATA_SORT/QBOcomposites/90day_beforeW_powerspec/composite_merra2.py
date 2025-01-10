import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

from CASutils import filter_utils as filt
from CASutils import qbo_utils as qbo
from CASutils import averaging_utils as avg
from CASutils import cospec_utils as cospec
import sys

pre=[50]
latbnd=[5]

pathout="/project/cas/islas/python_savs/CAM7_vertres_paper/DATA_SORT/QBOcomposites/90day_beforeW_powerspec/"
basepath="/project/cas/islas/python_savs/CAM7_vertres_paper/RAW_DATA/MERRA2/dayfrom6h/"

for ip in pre:
    for ilatbnd in latbnd:
        dat = xr.open_mfdataset(basepath+"*.nc").sel(plev=ip, lat=slice(-1.*ilatbnd-2, ilatbnd+2))
        u = dat.u.load()
        w = dat.w.load()

        uzm = u.mean('lon')
        uqbo = avg.cosweightlat(uzm, -5, 5)
    
        # Taking lags -95 to 5, using 10% for tapering, so it's ~90 days before
        lag1=-95 ; lag2=5

        # Take 30 day running mean for QBO
        u30d = filt.runningmean(uqbo, 30, dropna=False)

        # Find the location of E-W transitions
        ewloc = qbo.finde2w(u30d)

        # make sure there's enough room for the lags
        ewloc = ewloc[ (ewloc+lag1 > 0) & (ewloc+lag2 < uzm.time.size) ] 

        lagarr = np.arange(lag1, lag2, 1)

        cospecout_uw = []
        cospecout_uu = []
        cospecout_ww = []

        for icomp in np.arange(0,len(ewloc),1):
            uuse = u.isel(time=slice(int(ewloc[icomp] + lag1), int(ewloc[icomp]) + lag2))
            uuse['time'] = lagarr

            wuse = w.isel(time=slice(int(ewloc[icomp] + lag1), int(ewloc[icomp]) + lag2))
            wuse['time'] = lagarr

            cospec_temp = cospec.cospeccalc(uuse.reset_coords().u, wuse.reset_coords().w,
                            dosymmetries=True, dosegs=False, ftaper=0.1, deseas=False, detrend=False,
                            lat_bnds=(-1.*ilatbnd, ilatbnd))
            cospec_temp = cospec_temp.rename('cospec_uw')
            cospecout_uw.append(cospec_temp)


            cospec_temp = cospec.cospeccalc(uuse.reset_coords().u, uuse.reset_coords().u,
                            dosymmetries=True, dosegs=False, ftaper=0.1, deseas=False, detrend=False,
                            lat_bnds=(-1.*ilatbnd, ilatbnd))
            cospec_temp = cospec_temp.rename('cospec_uu')
            cospecout_uu.append(cospec_temp)

            cospec_temp = cospec.cospeccalc(wuse.reset_coords().w, wuse.reset_coords().w,
                            dosymmetries=True, dosegs=False, ftaper=0.1, deseas=False, detrend=False,
                            lat_bnds=(-1.*ilatbnd, ilatbnd))
            cospec_temp = cospec_temp.rename('cospec_ww')
            cospecout_ww.append(cospec_temp)

        cospecout_uw = xr.concat(cospecout_uw, dim='icomp')
        cospecout_uu = xr.concat(cospecout_uu, dim='icomp')
        cospecout_ww = xr.concat(cospecout_ww, dim='icomp')

        dat = xr.merge([cospecout_uw, cospecout_uu, cospecout_ww])

        dat.to_netcdf(pathout+'MERRA2_'+str(ip)+'hpa_'+str(ilatbnd)+'S'+str(ilatbnd)+'N_90days.nc')















