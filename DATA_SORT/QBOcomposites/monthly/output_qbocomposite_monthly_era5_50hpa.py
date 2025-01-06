import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from math import nan

from CASutils import averaging_utils as avg
from CASutils import qbo_utils as qbo
from CASutils import plotposition_utils as plotpos
from CASutils import readdata_utils as read
from CASutils import colorbar_utils as cbars

from scipy.ndimage import label
import sys

#expnames=['dz1000']
expnames=['dz900','dz800','dz700','dz600','dz500','dz400']

pathout="/project/cas/islas/python_savs/CAM7_vertres_paper/DATA_SORT/QBOcomposites/monthly/"

def dothecomposite(dat, nlag):
    # compositing based on the transition from easterly to westerly at 60hPa
    try:
#        u = dat.uzm.sel(ilev=60, method='nearest')
        u = dat.uzm.interp(ilev=50)
    except:
#        u = dat.uzm.sel(plev=60, method='nearest')
        u = dat.uzm.interp(plev=50)
    e2w = qbo.finde2w(u)
    e2w = e2w[ (e2w > nlag) & (e2w < dat.time.size - nlag) ]
    
    lagarr = np.arange(-nlag,nlag+1,1)
    for icomp in np.arange(0,len(e2w),1):
        datuse = dat.isel(time=slice( int(e2w[icomp] - nlag), int(e2w[icomp]) + nlag + 1))
        datuse['time'] = lagarr
        if (icomp == 0):
            dat_comp = datuse / len(e2w)
        else:
            dat_comp = dat_comp + datuse / len(e2w)
            
    return dat_comp




#---read in the TEM diagnostics
tem = xr.open_dataset("/project/cas/islas/python_savs/CAM7_vertres_paper/DATA_SORT/TEMdiags/mon/"+
                      "ERA5.nc")
tem['plev'] = tem.plev/100.
#---read in the gravity wave drag
#    basepath="/project/cas/islas/python_savs/L83_paper/DATA_SORT/zonalmean/mon/"
#    var=['BUTGWSPEC','UTGWORO','UTGWSPEC']
#    gwd=[]
#
#    for ivar in var:
#        dat = xr.open_mfdataset(basepath+iexp+"/"+ivar+".nc")
#        dat = read.fixcesmtime(dat)
#        gwd.append(dat[ivar].rename(ivar))
#
#    gwd = xr.merge(gwd)

#---omit the first year
tem = tem.isel(time=slice(12,tem.time.size))
#    gwd = gwd.isel(time=slice(12,gwd.time.size))

#---take the tropical average
w = 5 # 5S to 5N
tem_tr = avg.cosweightlat(tem, -1.*w, w)
#    gwd_tr = avg.cosweightlat(gwd, -1.*w, w)

#    gwd_tr['time'] = tem.time
#    gwd_tr = gwd_tr.rename({'plev':'ilev'})
#    gwd_tr['ilev'] = tem.ilev

#    alldat = xr.merge([tem_tr, gwd_tr])
alldat = tem_tr

nlag=20
comp = dothecomposite(alldat,nlag)

comp.to_netcdf(pathout+"ERA5_composite_50hpa.nc")

