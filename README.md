# CAM7_vertres_paper
Repository for the scripts that are used to process the data and make the figures for Simpson et al (2025) "The path toward the vertical grid options for the Community Atmosphere Model version 7: the impact of vertical resolution on the QBO and tropical waves"

All the data that is used in, and generated from, the scripts below can be downloaded from ???????.  In the following `$DATAROOT` refers to the directory that contains this dataset.

### Pre-computations

* Calculation of daily TEM diagnostics from daily fluxes output directly from CAM.  Computed using scripts in `./DATA_SORT/TEMdiags/day/` and files are saved in `$DATAROOT/DATA_SORT/TEMdiags/day/`
* Calculation of monthly averaged TEM diagnostics from the daily averaged TEM diagnostics.  Computed using scripts in `./DATA_SORT/TEMdiags/mon/` and files are saved in `$DATAROOT/DATA_SORT/TEMdiags/mon`
* Calculation of lagged composites relative to the timing of transition of the QBO from easterly to westerly.  Calculations at `./DATA_SORT/QBOcomposites/monthly` and files are saved in `$DATAROOT/DATA_SORT/QBOcomposites/monthly` and are used in Fig 2.
* Calculation of cospectra for the 90 days prior to the QBO Easterly to Westerly transition at `./DATA_SORT/QBOcomposites/90day_beforeW_powerspec/` and data saved at `$DATAROOT/DATA_SORT/90day_beforeW_powerspec` for use in Fig. 6
* Calculation of dispersion curves for spectra plots at `./DATA_SORT/dispersion_curves/`
* Calculation of wavenumber-frequency power spectra at `./DATA_SORT/WKdiags/` for use in Figs 7 and 8
* Calculation of the power in MRG waves over the lifecycle of the QBO for Fig 9 at `./DATA_SORT/MRG_filter`
* Calculation of MJO filtered vertical velocity for the winter season at `./DATA_SORT/MJOfilter/output_MJOfilter_omega_ERA5.py` for Fig. 10.
* Calculation of the climatological seasonal cycle of 5S to 5N specific humidity at `./DATA_SORT/tape_record/` for Fig 11.
* Obtaining the zonal mean zonal wind at 10 hPa at `./DATA_SORT/u_10hpa/` and then calculating the SSW frequency at `./DATA_SORT/SSWdates/`
* Obtaining the QBO and the MJO filtered vertical velocity variance for Fig 18 at `./DATA_SORT/MJO_OMEGA500/`
