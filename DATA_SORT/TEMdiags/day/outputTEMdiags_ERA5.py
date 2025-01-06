# Script to calculate TEM diagnostics
# This assumes the data have already been organized into zonal mean fluxes
# Uzm, THzm, VTHzm, Vzm, UVzm, UWzm, Wzm as was done by ctem.F90 for FV dycore
# Following the DynVarMIP protocol (Gerber and Manzini 2016, GMD)

# Using fluxes that have been computed from the model level data and interpolated
# onto pressure levels on CISL machines at
# /glade/u/home/islas/python/sortera5/grabera5zmfluxes_mlev.ipynb

import xarray as xr
import numpy as np
from scipy import integrate
import sys
from CASutils import readdata_utils as read
import warnings
warnings.filterwarnings('ignore')

expname='ERA5'

basepath="/project/cas/islas/python_savs/CAM7_vertres_paper/RAW_DATA/ERA5/ZMfluxes/plev_from_mlev/"

outdir="/project/cas/islas/python_savs/CAM7_vertres_paper/DATA_SORT/TEMdiags/day/"

# set up constants for TEM calculations
p0=101325.
a=6.371e6
om=7.29212e-5
H=7000.
g0=9.80665

dat_xr = xr.open_mfdataset(basepath+"/zmfluxes*.nc")
dat_xr = dat_xr.squeeze()

lat = dat_xr.lat
pre = dat_xr.plev

pre = pre/100. # converting from Pa to hPa

npre = pre.size
nlat = lat.size
ntime = dat_xr.time.size

# convert to numpy arrays
uzm = np.array(dat_xr.Uzm)
thzm = np.array(dat_xr.THzm)
vthzm = np.array(dat_xr.VTHzm)
vzm = np.array(dat_xr.Vzm)
uvzm = np.array(dat_xr.UVzm)
uwzm = np.array(dat_xr.UWzm)
wzm = np.array(dat_xr.Wzm)

time = np.arange(0,dat_xr.time.size,1)*86400.
uzm_xr = dat_xr.Uzm
uzm_xr['time'] = time

dudt = uzm_xr.differentiate('time')
dudt = dudt.rename('dudt')
dudt['time'] = dat_xr.time

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

# compute the latitudinal gradient of U
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


dudt = xr.DataArray(dudt, coords=dat_xr.coords, name='dudt',
                 attrs={'long_name':'Time rate of change of zonal mean zonal wind',
                        'units':'m/s2'})
uzm = xr.DataArray(uzm, coords = [dat_xr.time, dat_xr.plev, dat_xr.lat],
                        dims=['time','plev','lat'], name='uzm',
                 attrs={'long_name':'zonal mean zonal wind', 'units':'m/s'})
epfy = xr.DataArray(epfy, coords =[dat_xr.time, dat_xr.plev, dat_xr.lat],
                          dims=['time','plev','lat'], name='epfy',
                 attrs={'long_name':'northward component of E-P flux', 'units':'m3/s2'})
epfz = xr.DataArray(epfz, coords = [dat_xr.time, dat_xr.plev, dat_xr.lat],
                          dims=['time','plev','lat'], name='epfz',
                 attrs={'long_name':'upward component of E-P flux', 'units':'m2/s2'})
vtem = xr.DataArray(vtem, coords = [dat_xr.time, dat_xr.plev, dat_xr.lat], 
                          dims=['time','plev','lat'],name='vtem',
                 attrs={'long_name':'Transformed Eulerian mean northward wind', 'units':'m/s'})
wtem = xr.DataArray(wtem, coords =[dat_xr.time, dat_xr.plev, dat_xr.lat],
                          dims=['time','plev','lat'],name='wtem',
                 attrs={'long_name':'Transformed Eulerian mean upward wind','units':',/s'})
psitem = xr.DataArray(psitem, coords = [dat_xr.time, dat_xr.plev, dat_xr.lat],
                              dims=['time','plev','lat'], name='psitem',
                 attrs={'long_name':'Transformed Eulerian mean mass stream function','units':'kg/s'})
utendepfd = xr.DataArray(utendepfd, coords = [dat_xr.time, dat_xr.plev, dat_xr.lat],
                                    dims=['time','plev','lat'], name='utendepfd',
                 attrs={'long_name':'tendency of eastward wind due to Eliassen-Palm flux divergence',
                       'units':'m/s2'})
utendvtem = xr.DataArray(utendvtem, coords = [dat_xr.time, dat_xr.plev, dat_xr.lat],
                                    dims=['time','plev','lat'], name='utendvtem',
attrs={'long_name':'tendency of eastward wind due to TEM northward wind advection and the coriolis term'
          ,'units':'m/s2'})
utendwtem = xr.DataArray(utendwtem, 
                coords = [dat_xr.time,dat_xr.plev,dat_xr.lat],
                dims=['time','plev','lat'], name='utendwtem',
attrs={'long_name':'tendency of eastward wind due to TEM upward wind advection','units':'m/s2'})

print('before output')

uzm.to_netcdf(outdir+expname+".nc")
epfy.to_netcdf(outdir+expname+".nc", mode="a")
epfz.to_netcdf(outdir+expname+".nc", mode="a")
vtem.to_netcdf(outdir+expname+".nc", mode="a")
wtem.to_netcdf(outdir+expname+".nc", mode="a")
psitem.to_netcdf(outdir+expname+".nc", mode="a")
utendepfd.to_netcdf(outdir+expname+".nc", mode="a")
utendvtem.to_netcdf(outdir+expname+".nc", mode="a")
utendwtem.to_netcdf(outdir+expname+".nc", mode="a")
dudt.to_netcdf(outdir+expname+".nc", mode="a")
