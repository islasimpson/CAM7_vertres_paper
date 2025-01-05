# Script to calculate TEM diagnostics
# This assumes the data have already been organized into zonal mean fluxes
# Uzm, THzm, VTHzm, Vzm, UVzm, UWzm, Wzm as was done by ctem.F90 for FV dycore
# Following the DynVarMIP protocol (Gerber and Manzini 2016, GMD)
import xarray as xr
import numpy as np
from scipy import integrate
import sys
from CASutils import readdata_utils as read

#expname=[ 'dz1000', 'dz900','dz800','dz700','dz600','dz500','dz400']
#expname=[ 'dz1000' ]
#expname=['dz1000','dz900','dz800','dz700','dz600','dz500','dz400']
expname=['dz500_buggy']

basepath="/project/cas/islas/python_savs/L83_paper/RAW_DATA/80km/ZMfluxes/"

outdir="/project/cas/islas/python_savs/L83_paper/DATA_SORT/TEMdiags/day/"

# set up constants for TEM calculations
p0=101325.
a=6.371e6
om=7.29212e-5
H=7000.
g0=9.80665

for iexp in expname:
    print(iexp)
    filedir=basepath+iexp+"/"
    uzm_xr = xr.open_mfdataset(filedir+"/Uzm*.nc")
    uzm_xr = read.fixcesmtime_daily(uzm_xr)
    thzm_xr = xr.open_mfdataset(filedir+"/THzm*.nc")
    thzm_xr = read.fixcesmtime_daily(thzm_xr)
    vthzm_xr = xr.open_mfdataset(filedir+"/VTHzm*.nc")
    vthzm_xr = read.fixcesmtime_daily(vthzm_xr)
    vzm_xr = xr.open_mfdataset(filedir+"/Vzm*.nc")
    vzm_xr = read.fixcesmtime_daily(vzm_xr)
    uvzm_xr = xr.open_mfdataset(filedir+"/UVzm*.nc")
    uvzm_xr = read.fixcesmtime_daily(uvzm_xr)
    uwzm_xr = xr.open_mfdataset(filedir+"/UWzm*.nc")
    uwzm_xr = read.fixcesmtime_daily(uwzm_xr)
    wzm_xr = xr.open_mfdataset(filedir+"/Wzm*.nc")
    wzm_xr = read.fixcesmtime_daily(wzm_xr)


    lat = uzm_xr.lat
    pre = uzm_xr.ilev 

    uzm_xr = uzm_xr.Uzm.isel(zlon=0)
    thzm_xr = thzm_xr.THzm.isel(zlon=0)
    vthzm_xr = vthzm_xr.VTHzm.isel(zlon=0)
    vzm_xr = vzm_xr.Vzm.isel(zlon=0)
    uvzm_xr = uvzm_xr.UVzm.isel(zlon=0)
    uwzm_xr = uwzm_xr.UWzm.isel(zlon=0)
    wzm_xr = wzm_xr.Wzm.isel(zlon=0)

    npre = pre.size
    nlat = lat.size
    ntime = uzm_xr.time.size

    # convert to numpy arrays
    uzm = np.array(uzm_xr)
    thzm = np.array(thzm_xr)
    vthzm = np.array(vthzm_xr)
    vzm = np.array(vzm_xr)
    uvzm = np.array(uvzm_xr)
    uwzm = np.array(uwzm_xr)
    wzm = np.array(wzm_xr)

    # Compute du/dt
    uzm_xr = uzm_xr.chunk({'time':uzm_xr.time.size, 'ilev':uzm_xr.ilev.size, 'lat':uzm_xr.lat.size})
    dudt = uzm_xr.differentiate('time')
    dudt = dudt.rename('dudt')

    latrad = np.array((lat.values/180.)*np.pi)
    f=2.*om*np.sin(latrad[:])

    # setup pre(ntime,npre,nlat) and latrad(ntime,npre,nlat)
    latradarray = np.tile(latrad,npre*ntime)
    latradarray = np.reshape(latradarray,[ntime,npre,nlat])

    prearray = np.tile(pre,nlat*ntime)
    prearray = np.reshape(prearray,[ntime,nlat,npre])
    prearray = np.moveaxis(prearray,2,1)

    # convert w terms from m/s to Pa/s
    uwzm = -1.*uwzm*prearray*100./H
    wzm = -1.*wzm*prearray*100./H

    # compute the latitudinal gradient of U / cos(phi)
    dudphi = (1./(a*np.cos(latradarray)))*np.gradient(uzm*np.cos(latradarray), latrad, axis=2)

    # compute the vertical gradient of theta and u
    dthdp = np.gradient(thzm, pre*100.,axis=1)
    dudp = np.gradient(uzm, pre*100., axis=1)

    # compute eddy streamfunction and its vertical gradient
    psieddy = vthzm/dthdp
    dpsidp = np.gradient(psieddy,pre*100.,axis=1)

    # (1/acos(phii))**d(psi*cosphi/dphi) for getting w*
    psicos = psieddy*np.cos(latradarray)
    dpsidy = (1./(a*np.cos(latradarray)))*np.gradient(psicos,latrad,axis=2)

    # TEM vertical velocity (Eq A7 of dynvarmip)
    wtem = wzm+dpsidy

    # utendwtem (Eq A10 of dynvarmip)
    utendwtem = -1.*wtem*dudp

    # vtem (Eq A6 of dynvarmip)
    vtem = vzm-dpsidp

    # utendvtem (Eq A9 of dynvarmip)
    farray = np.tile(f,npre*ntime)
    farray = np.reshape(farray,[ntime,npre,nlat])
    utendvtem = vtem*(farray - dudphi)

    # calculate E-P fluxes
    epfy = a*np.cos(latradarray)*(dudp*psieddy - uvzm) # A2
    epfz = a*np.cos(latradarray)*( (farray-dudphi)*psieddy - uwzm) # A3

    # calculate E-P flux divergence and zonal wind tendency due to resolved waves (A5)
    depfydphi = np.gradient(epfy*np.cos(latradarray),latrad, axis=2)*(1./(a*np.cos(latradarray)))
    depfzdp = np.gradient(epfz,pre*100.,axis=1)
    utendepfd = depfydphi + depfzdp
    utendepfd = (1./(a*np.cos(latradarray)))*utendepfd

    # TEM stream function, Eq (A8)
    topvzm = np.zeros([ntime,1,nlat])
    vzmwithzero = np.concatenate((topvzm, vzm), axis=1)
    toppre = np.zeros([1])
    prewithzero = np.concatenate((toppre, pre))
    intv = integrate.cumtrapz(vzmwithzero,prewithzero*100.,axis=1)
    psitem = (2*np.pi*a*np.cos(latradarray)/g0)*(intv - psieddy)

    # final scaling of E-P fluxes and divergence to transform to log-pressure
    epfy = epfy*(prearray*100.)/p0 # A13
    epfz = -1.*(H/p0)*epfz # A14
    wtem = -1.*(H/(prearray*100.))*wtem # A16

    dudt = xr.DataArray(dudt, coords=uzm_xr.coords, name='dudt',
                     attrs={'long_name':'Time rate of change of zonal mean zonal wind',
                            'units':'m/s2'})
    uzm = xr.DataArray(uzm, coords = uzm_xr.coords, name='uzm',
                     attrs={'long_name':'zonal mean zonal wind', 'units':'m/s'})
    epfy = xr.DataArray(epfy, coords = uzm_xr.coords, name='epfy',
                     attrs={'long_name':'northward component of E-P flux', 'units':'m3/s2'})
    epfz = xr.DataArray(epfz, coords = uzm_xr.coords, name='epfz',
                     attrs={'long_name':'upward component of E-P flux', 'units':'m2/s2'})
    vtem = xr.DataArray(vtem, coords = uzm_xr.coords, name='vtem',
                     attrs={'long_name':'Transformed Eulerian mean northward wind', 'units':'m/s'})
    wtem = xr.DataArray(wtem, coords = uzm_xr.coords, name='wtem',
                     attrs={'long_name':'Transformed Eulerian mean upward wind','units':',/s'})
    psitem = xr.DataArray(psitem, coords = uzm_xr.coords, name='psitem',
                     attrs={'long_name':'Transformed Eulerian mean mass stream function','units':'kg/s'})
    utendepfd = xr.DataArray(utendepfd, coords = uzm_xr.coords, name='utendepfd',
                     attrs={'long_name':'tendency of eastward wind due to Eliassen-Palm flux divergence',
                           'units':'m/s2'})
    utendvtem = xr.DataArray(utendvtem, coords = uzm_xr.coords, name='utendvtem',
    attrs={'long_name':'tendency of eastward wind due to TEM northward wind advection and the coriolis term'
              ,'units':'m/s2'})
    utendwtem = xr.DataArray(utendwtem, coords = uzm_xr.coords, name='utendwtem',
    attrs={'long_name':'tendency of eastward wind due to TEM upward wind advection','units':'m/s2'})

    uzm.to_netcdf(outdir+iexp+"_80km.nc")
    epfy.to_netcdf(outdir+iexp+"_80km.nc", mode="a")
    epfz.to_netcdf(outdir+iexp+"_80km.nc", mode="a")
    vtem.to_netcdf(outdir+iexp+"_80km.nc", mode="a")
    wtem.to_netcdf(outdir+iexp+"_80km.nc", mode="a")
    psitem.to_netcdf(outdir+iexp+"_80km.nc", mode="a")
    utendepfd.to_netcdf(outdir+iexp+"_80km.nc", mode="a")
    utendvtem.to_netcdf(outdir+iexp+"_80km.nc", mode="a")
    utendwtem.to_netcdf(outdir+iexp+"_80km.nc", mode="a")
    dudt.to_netcdf(outdir+iexp+"_80km.nc", mode="a")
