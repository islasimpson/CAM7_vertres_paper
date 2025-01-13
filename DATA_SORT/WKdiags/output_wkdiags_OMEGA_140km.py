import xarray as xr
import numpy as np
from CASutils import cospec_utils as cospec
from CASutils import calendar_utils as cal
from CASutils import wavenumber_frequency_functions as wf
import sys

pathout="/project/cas/islas/python_savs/CAM7_vertres_paper/DATA_SORT/WKdiags/"

pre=[50]
latbnd=[5]

#exps=['dz1000']
exps=['dz1000','dz900','dz800','dz700','dz600','dz500','dz400']

for iexp in exps:
    print(iexp)
    for ip in pre:
        for ilatbnd in latbnd:
            w = xr.open_mfdataset("/project/cas/islas/python_savs/CAM7_vertres_paper/RAW_DATA/140km/"+
               "day/plev/day_from_6h/"+iexp+"/OMEGA_*.nc")\
               .sel(pre=ip, lat=slice(-1.*ilatbnd-2, ilatbnd+2)).load()

            if (iexp == 'dz1000'):
                ystart=1988

            wpower = cospec.cospeccalc(w.OMEGA, w.OMEGA, lat_bnds=(-1.*ilatbnd, ilatbnd),
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

            dout=xr.merge([wpower, background])

            dout.to_netcdf(pathout+iexp+'_'+str(ip)+'hpa_'+str(ilatbnd)+'Sto'+str(ilatbnd)+'N.nc') 

           

