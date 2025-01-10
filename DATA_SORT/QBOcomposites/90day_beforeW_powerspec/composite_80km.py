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

#exp=['dz1000','dz900','dz800','dz700','dz600','dz500','dz400']
#exp=['dz1000']
#exp=['dz800','dz700','dz600','dz500']
exp=['dz800','dz700','dz600','dz500']

for iexp in exp:
    print(iexp)
    for ip in pre:
        for ilatbnd in latbnd:
            dat_tem = xr.open_dataset("/project/cas/islas/python_savs/CAM7_vertres_paper/DATA_SORT/TEMdiags/"+
               "day/"+iexp+"_80km.nc")
            dat_tem = dat_tem.sel(ilev=ip, method='nearest')
            dat_tem = dat_tem.where(dat_tem.time.dt.hour == 12, drop=True)

            dat_zmfluxes = xr.open_dataset("/project/cas/islas/python_savs/CAM7_verters_paper/DATA_SORT/"+
               "E_W_waves_decomp/fluxes/fluxes_E_W_stat_"+iexp+"_80km.nc")
            dat_zmfluxes = dat_zmfluxes.sel(pre=ip, method='nearest')

            u = xr.open_mfdataset("/project/cas/islas/python_savs/CAM7_vertres_paper/RAW_DATA/80km/"+
               "day/plev/day_from_6h/"+iexp+"/U_*.nc")\
               .sel(pre=ip, lat=slice(-1.*ilatbnd-2, ilatbnd+2)).load()

            w = xr.open_mfdataset("/project/cas/islas/python_savs/CAM7_vertres_paper/RAW_DATA/80km/"+
               "day/plev/day_from_6h/"+iexp+"/OMEGA_*.nc")\
               .sel(pre=ip, lat=slice(-1.*ilatbnd-2, ilatbnd+2)).load()

            ystart = dat_zmfluxes.isel(time=0).time.dt.year.values
            yend = dat_zmfluxes.isel(time=dat_zmfluxes.time.size-1).time.dt.year.values

            if (iexp == 'dz1000'):
                ystart=1988
#            else:
#                dat_tem = dat_tem.sel(time=slice(str(ystart)+'-01-01',str(yend)+'-12-31'))
            
            dat_tem = dat_tem.sel(time=slice(str(ystart)+'-01-01',str(yend)+'-12-31'))
            dat_zmfluxes = dat_zmfluxes.sel(time=slice(str(ystart)+'-01-01',str(yend)+'-12-31'))
            u = u.sel(time=slice(str(ystart)+'-01-01',str(yend)+'-12-31'))
            w = w.sel(time=slice(str(ystart)+'-01-01',str(yend)+'-12-31'))

            nt_tem = dat_tem.time.size
            nt_zmfluxes = dat_zmfluxes.time.size
            nt_u = u.time.size
            nt_w = w.time.size


            if ( (nt_tem != nt_zmfluxes) | (nt_tem != nt_u) | (nt_tem != nt_w) ):
                print(f"Something's wrong, nt_tem={nt_tem},nt_zmfluxes={nt_zmfluxes},nt_u={nt_u},nt_w={nt_w}")
                print(f'iexp={iexp}, ip={ip}, ilatbnd={ilatbnd}')
                sys.exit()


            #---take the tropical average (always 5S-5N for QBO winds)
            uqbo = avg.cosweightlat(dat_tem.uzm,-5,5)
            dat_tem_tropics = avg.cosweightlat(dat_tem,-1*ilatbnd,ilatbnd)
            dat_zmfluxes_tropics = avg.cosweightlat(dat_zmfluxes,-1*ilatbnd,ilatbnd)

            # Taking lags -66 to 6, using 10% for tapering, so it's ~ 2 months before
            lag1=-95 ; lag2=5

            # 30 day running mean for QBO
            u30d = filt.runningmean(uqbo,30,dropna=False)

            # find the location of E-W transitions
            ewloc = qbo.finde2w(u30d)

            # make sure there's enough room for the lags
            ewloc = ewloc[ (ewloc+lag1 > 0) & (ewloc+lag2 < dat_tem_tropics.time.size) ]

            lagarr = np.arange(lag1, lag2, 1)

            comp_uzm = []
            comp_uw_e = []
            cospecout_uw = []
            cospecout_uu = []
            cospecout_ww = []
            comp_fz = []


            for icomp in np.arange(0,len(ewloc),1):
                uzmuse = dat_tem_tropics.uzm.isel(time=slice(int(ewloc[icomp] + lag1),int(ewloc[icomp])+ lag2))
                uzmuse['time'] = lagarr
                comp_uzm.append(uzmuse)

                fzuse = dat_tem_tropics.epfz.isel(time=slice(int(ewloc[icomp] + lag1), int(ewloc[icomp])+lag2))
                fzuse['time'] = lagarr
                comp_fz.append(fzuse)


                uweuse = dat_zmfluxes_tropics.UWzm_E.isel(time=slice(int(ewloc[icomp] + lag1),int(ewloc[icomp])+lag2))
                uweuse['time'] = lagarr
                comp_uw_e.append(uweuse)

                uuse = u.U.isel(time=slice(int(ewloc[icomp]+lag1),int(ewloc[icomp]) + lag2))
                uuse['time'] = lagarr

                wuse = w.OMEGA.isel(time=slice(int(ewloc[icomp]+lag1),int(ewloc[icomp]) + lag2))
                wuse['time'] = lagarr

                cospec_temp = cospec.cospeccalc(uuse.reset_coords().U,wuse.reset_coords().OMEGA,
                                dosymmetries=True,dosegs=False,ftaper=0.1,deseas=False,detrend=False,
                                lat_bnds=(-1*ilatbnd,ilatbnd))
                cospec_temp = cospec_temp.rename('cospec_uw')
                cospecout_uw.append(cospec_temp)

                cospec_temp = cospec.cospeccalc(uuse.reset_coords().U,uuse.reset_coords().U,
                                dosymmetries=True,dosegs=False,ftaper=0.1,deseas=False,detrend=False,
                                lat_bnds=(-1*ilatbnd,ilatbnd))
                cospec_temp = cospec_temp.rename('cospec_uu')
                cospecout_uu.append(cospec_temp)

                cospec_temp = cospec.cospeccalc(wuse.reset_coords().OMEGA,wuse.reset_coords().OMEGA,
                                dosymmetries=True,dosegs=False,ftaper=0.1,deseas=False,detrend=False,
                                lat_bnds=(-1*ilatbnd,ilatbnd))
                cospec_temp = cospec_temp.rename('cospec_ww')
                cospecout_ww.append(cospec_temp)


            comp_uzm = xr.concat(comp_uzm, dim='icomp')
            comp_fz = xr.concat(comp_fz, dim='icomp')
            comp_uw_e = xr.concat(comp_uw_e, dim='icomp')
            cospecout_uw = xr.concat(cospecout_uw, dim='icomp')
            cospecout_uu = xr.concat(cospecout_uu, dim='icomp')
            cospecout_ww = xr.concat(cospecout_ww, dim='icomp')

            dat = xr.merge([comp_uzm,comp_fz,comp_uw_e,cospecout_uw,cospecout_uu,cospecout_ww])

            dat.to_netcdf(pathout+iexp+'_'+str(ip)+'hpa_'+str(ilatbnd)+'S'+str(ilatbnd)+'N_80km.nc')

















 
