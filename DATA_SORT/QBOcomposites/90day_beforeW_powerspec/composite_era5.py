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

for ip in pre:
    for ilatbnd in latbnd:
        era5_tem = xr.open_dataset("/project/cas/islas/python_savs/CAM7_vertres_paper/DATA_SORT/TEMdiags/day/"+
                                   "ERA5.nc")
        era5_tem['plev'] = era5_tem.plev/100.
        era5_tem = era5_tem.sel(plev=ip, method='nearest')

        era5_zmfluxes = xr.open_dataset(
            "/project/cas/islas/python_savs/CAM7_vertres_paper_paper/DATA_SORT/E_W_waves_decomp/fluxes/"+
            "fluxes_E_W_stat_era5.nc")
        era5_zmfluxes = era5_zmfluxes.sel(level=ip, method='nearest')

        u = xr.open_mfdataset("/project/cas/islas/python_savs/CAM7_vertres_paper/RAW_DATA/ERA5/"+
            "U_dayfrom6h/U_*.nc").sel(level=ip, lat=slice(-1.*ilatbnd-2, ilatbnd+2)).load()

        w = xr.open_mfdataset("/project/cas/islas/python_savs/CAM7_vertres_paper/RAW_DATA/ERA5/"+
            "W_dayfrom6h/W_*.nc").sel(level=ip, lat=slice(-1.*ilatbnd-2, ilatbnd+2)).load()

        ystart = era5_zmfluxes.isel(time=0).time.dt.year.values
        yend = era5_zmfluxes.isel(time=era5_zmfluxes.time.size-1).time.dt.year.values

        era5_tem = era5_tem.sel(time=slice(str(ystart)+'-01-01',str(yend)+'-12-31'))
        era5_zmfluxes = era5_zmfluxes.sel(time=slice(str(ystart)+'-01-01',str(yend)+'-12-31'))
        u = u.sel(time=slice(str(ystart)+'-01-01',str(yend)+'-12-31'))
        w = w.sel(time=slice(str(ystart)+'-01-01',str(yend)+'-12-31'))

        #---check they have the same length
        nt_tem = era5_tem.time.size
        nt_zmfluxes = era5_zmfluxes.time.size
        nt_u = u.time.size
        nt_w = w.time.size

        if ( (nt_tem != nt_zmfluxes) | (nt_tem != nt_u) | (nt_tem != nt_w) ):
            print(f"Something's wrong, nt_tem={nt_tem},nt_zmfluxes={nt_zmfluxes},nt_u={nt_u},nt_w={nt_w}")

        #---take tropical average (always 5S-5N for QBO winds)
        uqbo = avg.cosweightlat(era5_tem.uzm,-5,5)
        era5_tem_tropics = avg.cosweightlat(era5_tem,-1*ilatbnd,ilatbnd)
        era5_zmfluxes_tropics = avg.cosweightlat(era5_zmfluxes,-1*ilatbnd,ilatbnd)

        # Taking lags -66 to 6, using 10% for tapering, so it's ~2 months before        
        lag1=-95 ; lag2=5

        # 30 day running mean for QBO
        u30d = filt.runningmean(uqbo,30,dropna=False)

        # find the location of E-W transitions
        ewloc = qbo.finde2w(u30d)

        # make sure there's enough room for the lags
        ewloc = ewloc[ (ewloc+lag1 > 0) & (ewloc+lag2 < era5_tem_tropics.time.size) ]

        lagarr = np.arange(lag1,lag2,1)

        comp_uzm = []
        comp_uw_e = []
        cospecout_uw = []
        cospecout_uu = []
        cospecout_ww = []

        for icomp in np.arange(0,len(ewloc),1):
            uzmuse = era5_tem_tropics.uzm.isel(time=slice(int(ewloc[icomp] + lag1),int(ewloc[icomp])+ lag2))
            uzmuse['time'] = lagarr
            comp_uzm.append(uzmuse)

            uweuse = era5_zmfluxes_tropics.UWzm_E.isel(time=slice(int(ewloc[icomp] + lag1),int(ewloc[icomp])+lag2))
            uweuse['time'] = lagarr
            comp_uw_e.append(uweuse)

            uuse = u.U.isel(time=slice(int(ewloc[icomp]+lag1),int(ewloc[icomp]) + lag2))
            uuse['time'] = lagarr

            wuse = w.W.isel(time=slice(int(ewloc[icomp]+lag1),int(ewloc[icomp]) + lag2))
            wuse['time'] = lagarr

            cospec_temp = cospec.cospeccalc(uuse.reset_coords().U,wuse.reset_coords().W,
                            dosymmetries=True,dosegs=False,ftaper=0.1,deseas=False,detrend=False,
                            lat_bnds=(-1*ilatbnd,ilatbnd))
            cospec_temp = cospec_temp.rename('cospec_uw')
            cospecout_uw.append(cospec_temp)

            cospec_temp = cospec.cospeccalc(uuse.reset_coords().U,uuse.reset_coords().U,
                            dosymmetries=True,dosegs=False,ftaper=0.1,deseas=False,detrend=False,
                            lat_bnds=(-1*ilatbnd,ilatbnd))
            cospec_temp = cospec_temp.rename('cospec_uu')
            cospecout_uu.append(cospec_temp)

            cospec_temp = cospec.cospeccalc(wuse.reset_coords().W,wuse.reset_coords().W,
                            dosymmetries=True,dosegs=False,ftaper=0.1,deseas=False,detrend=False,
                            lat_bnds=(-1*ilatbnd,ilatbnd))
            cospec_temp = cospec_temp.rename('cospec_ww')
            cospecout_ww.append(cospec_temp)


        comp_uzm = xr.concat(comp_uzm, dim='icomp')
        comp_uw_e = xr.concat(comp_uw_e, dim='icomp')
        cospecout_uw = xr.concat(cospecout_uw, dim='icomp')
        cospecout_uu = xr.concat(cospecout_uu, dim='icomp')
        cospecout_ww = xr.concat(cospecout_ww, dim='icomp')

        dat = xr.merge([comp_uzm,comp_uw_e,cospecout_uw,cospecout_uu,cospecout_ww])

        dat.to_netcdf(pathout+'ERA5_'+str(ip)+'hpa_'+str(ilatbnd)+'S'+str(ilatbnd)+'N_90days.nc')
