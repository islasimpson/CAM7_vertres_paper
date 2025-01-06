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
#expnames=['dz800','dz700','dz600','dz500']
#expnames=['dz700','dz600','dz500']
expnames=['dz500_buggy']

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



for iexp in expnames:

    #---read in the TEM diagnostics
    tem = xr.open_dataset("/project/cas/islas/python_savs/CAM7_vertres_paper/DATA_SORT/TEMdiags/mon/"+
                          iexp+"_80km.nc")

    #---read in the gravity wave drag
    basepath="/project/cas/islas/python_savs/CAM7_vertres_paper/DATA_SORT/zonalmean/mon/80km/"
    var=['BUTGWSPEC','UTGWORO','UTGWSPEC']
    gwd=[]

    for ivar in var:
        dat = xr.open_mfdataset(basepath+iexp+"/"+ivar+".nc")
        dat = read.fixcesmtime(dat)
        gwd.append(dat[ivar].rename(ivar))

    gwd = xr.merge(gwd)

    #---omit the first year
    tem = tem.isel(time=slice(12,tem.time.size))
    gwd = gwd.isel(time=slice(12,gwd.time.size))

    #---sort out an issue with dz1000, 1987 is missing
    if (iexp == 'dz1000'):
        tem = tem.sel(time=slice("1988-01-01","2005-12-31"))
        gwd = gwd.sel(time=slice("1988-01-01","2005-12-31"))

    #---take the tropical average
    w = 5 # 5S to 5N
    tem_tr = avg.cosweightlat(tem, -1.*w, w)
    gwd_tr = avg.cosweightlat(gwd, -1.*w, w)

    gwd_tr['time'] = tem.time
    gwd_tr = gwd_tr.rename({'plev':'ilev'})
    gwd_tr['ilev'] = tem.ilev

    alldat = xr.merge([tem_tr, gwd_tr])

    nlag=20
    comp = dothecomposite(alldat,nlag)

    comp.to_netcdf(pathout+iexp+"_composite_80km_50hpa.nc")

