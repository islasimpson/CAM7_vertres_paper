# CAM7_vertres_paper
Repository for the scripts that are used to process the data and make the figures for Simpson et al (2025) "The path toward the vertical grid options for the Community Atmosphere Model version 7: the impact of vertical resolution on the QBO and tropical waves"

All the data that is used in, and generated from, the scripts below can be downloaded from ???????.  In the following `$DATAROOT` refers to the directory that contains this dataset.

### Pre-computations

* Calculation of daily TEM diagnostics from daily fluxes output directly from CAM.  Computed using scripts in `./DATA_SORT/TEMdiags/day/` and files are saved in `$DATAROOT/DATA_SORT/TEMdiags/day/`
* Calculation of monthly averaged TEM diagnostics from the daily averaged TEM diagnostics.  Computed using scripts in `./DATA_SORT/TEMdiags/mon/` and files are saved in `$DATAROOT/DATA_SORT/TEMdiags/mon`
