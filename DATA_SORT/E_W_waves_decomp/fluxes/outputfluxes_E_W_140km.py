import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from math import nan
from CASutils import filter_utils as filt
from CASutils import linfit_utils as linfit
import sys

expname=['dz900','dz800','dz700','dz600','dz500','dz400']

basepath="/project/cas/islas/python_savs/CAM7_vertres_paper/RAW_DATA/140km/day/plev/day_from_6h/"
pathout="/project/cas/islas/python_savs/CAM7_vertres_paper/DATA_SORT/E_W_waves_decomp/fluxes/"

for iexp in expname:
    print(iexp)
    u = xr.open_mfdataset(basepath+iexp+'/U_*.nc').U
    u = u.sel(lat=slice(-30,30))
    print('read in U')
    v = xr.open_mfdataset(basepath+iexp+'/V_*.nc').V
    v = v.sel(lat=slice(-30,30))
    print('read in V')
    w = xr.open_mfdataset(basepath+iexp+'/OMEGA_*.nc').OMEGA
    w = w.sel(lat=slice(-30,30))
    print('read in W')
    t = xr.open_mfdataset(basepath+iexp+'/T_*.nc').T
    t = t.sel(lat=slice(-30,30))
    print('read in T')

    if (iexp == 'dz1000'):
        u = u.sel(time=slice("1988-01-01","2005-12-31"))
        v = v.sel(time=slice("1988-01-01","2005-12-31"))
        w = w.sel(time=slice("1988-01-01","2005-12-31"))
        t = t.sel(time=slice("1988-01-01","2005-12-31"))

    print('before calculating seasonal cycle')
    daystr = xr.DataArray(u.indexes['time'].strftime('%m-%d'), coords = u.time.coords, name='daystr')
    useas = u.groupby(daystr).mean('time')
    vseas = v.groupby(daystr).mean('time')
    wseas = w.groupby(daystr).mean('time')
    tseas = t.groupby(daystr).mean('time')

    u6harm = filt.calc_season_nharm(useas, 6, dimtime=0)
    v6harm = filt.calc_season_nharm(vseas, 6, dimtime=0)
    w6harm = filt.calc_season_nharm(wseas, 6, dimtime=0)
    t6harm = filt.calc_season_nharm(tseas, 6, dimtime=0)
    print('after calculating the seasonal cycle')

    uanoms = u.groupby(daystr) - u6harm
    vanoms = v.groupby(daystr) - v6harm
    wanoms = w.groupby(daystr) - w6harm
    tanoms = t.groupby(daystr) - t6harm

    #---linearly detrending the full time series
    print('before linearly detrending')
    uanoms = linfit.lineardetrend(uanoms, 'time')
    vanoms = linfit.lineardetrend(vanoms, 'time')
    wanoms = linfit.lineardetrend(wanoms, 'time')
    tanoms = linfit.lineardetrend(tanoms, 'time')
    print('after linearly detrending')

    uv = filt.wkfilter_flux(uanoms, vanoms,0.1,spd=1)
    vt = filt.wkfilter_flux(vanoms, tanoms,0.1,spd=1)
    uw = filt.wkfilter_flux(uanoms, wanoms,0.1,spd=1)

    uv_e = uv.eastward.rename('UVzm_E')
    uv_w = uv.westward.rename('UVzm_W')
    uv_stat = ((u6harm - u6harm.mean('lon'))*(v6harm - v6harm.mean('lon'))).mean('daystr')
    uv_stat = uv_stat.rename('UVzm_stat')

    vt_e = vt.eastward.rename('VTzm_E')
    vt_w = vt.westward.rename('VTzm_W')
    vt_stat = ((v6harm - v6harm.mean('lon'))*(t6harm - t6harm.mean('lon'))).mean('daystr')
    vt_stat = vt_stat.rename('VTzm_stat')

    uw_e = uw.eastward.rename('UWzm_E')
    uw_w = uw.westward.rename('UWzm_W')
    uw_stat = ((u6harm - u6harm.mean('lon'))*(w6harm - w6harm.mean('lon'))).mean('daystr')
    uw_stat = uw_stat.rename('UWzm_stat')

    ds = xr.merge([uv_e, uv_w, vt_e, vt_w, uw_e, uw_w, uv_stat, vt_stat, uw_stat])
    ds = ds.mean('lon')

    ds.load().to_netcdf(pathout+'fluxes_E_W_stat_'+iexp+'_140km.nc')




