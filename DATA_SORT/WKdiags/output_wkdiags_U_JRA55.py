import xarray as xr
import numpy as np
from CASutils import cospec_utils as cospec
from CASutils import calendar_utils as cal
from CASutils import wavenumber_frequency_functions as wf
import sys
from functools import partial

pathout="/project/cas/islas/python_savs/CAM7_vertres_paper/DATA_SORT/WKdiags/"

def preprocessor(ds):
    ds = ds.u
    return ds


#pre=[50,500]
pre=[500]
latbnd=[5]

for ip in pre:
    for ilatbnd in latbnd:
        w = xr.open_mfdataset("/project/cas/islas/python_savs/CAM7_vertres_paper/RAW_DATA/JRA55/dayfrom6h/"+
           "uvtw_*.nc", preprocess = partial(preprocessor)).sel(pre=ip, lat=slice(-1.*ilatbnd-2, ilatbnd+2)).load() 

        wpower = cospec.cospeccalc(w.u,w.u, lat_bnds=(-1.*ilatbnd, ilatbnd),
                   dosymmetries=True, dosegs=True, segsize=100, noverlap=60, ftaper=0.1,
                   spd=1, deseas=True, nharm=4, detrend=True)


        # computing the background
        wpower_avg = wpower.mean(dim='component')
        wpower_avg.loc[{'w':0}] = np.nan

        # reflect the background in frequency
        wpower_avg_r = wpower_avg.reindex(w=list(reversed(wpower_avg.w)))
        wpower_avg_r = wpower_avg_r.where(wpower_avg_r.w != 0, drop=True)
        wpower_avg_r['w'] = -1.*wpower_avg_r.w

        wpower_combined = xr.concat([wpower_avg_r, wpower_avg], dim='w')


        background = wf.smooth_wavefreq(wpower_combined, kern=wf.simple_smooth_kernel(),
                       nsmooth=50, freq_name='w')
        background = background.rename('background')
        background = background.sel(w=slice(np.min(wpower.w), np.max(wpower.w)))

        dout = xr.merge([wpower, background])

        dout.to_netcdf(pathout+'U_JRA55_'+str(ip)+'hpa_'+str(ilatbnd)+'Sto'+str(ilatbnd)+'N.nc') 




