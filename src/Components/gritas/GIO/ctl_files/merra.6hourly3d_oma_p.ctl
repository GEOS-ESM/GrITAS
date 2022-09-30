DSET     ../../obs/Y%y4/M%m2/D%d2/merra.%e3d_oma_p.%y4%m2%d2_%h2z.hdf
TITLE    6-hourly Gridded Obervation minus Analysis (oma) Values
DTYPE    netcdf
OPTIONS  template
UNDEF    1e+15
XDEF 540 LINEAR  -180 0.6666666666666667
YDEF 361 LINEAR   -90 0.5
ZDEF  42 LEVELS  1000 975 950 925 900 875 850 825 800  775  750  725  700  650
                  600 550 500 450 400 350 300 250 200  150  100   70   50   40
                   30  20  10   7   5   4   3   2   1  0.7  0.5  0.4  0.3  0.1
TDEF   1 LINEAR  00Z01JAN2000 6hr
EDEF   3  NAMES  mean  stdv  nobs
VARS 27
   tv_raob    42  t,z,y,x  RADIOSONDE:    Virtual  temperature [K]
   qv_raob    42  t,z,y,x  RADIOSONDE:    Specific humidity [g/kg]
    u_raob    42  t,z,y,x  RADIOSONDE:    Zonal      wind [m/s]
    v_raob    42  t,z,y,x  RADIOSONDE:    Meridional wind [m/s]
   ps_raob     0  t,y,x    RADIOSONDE:    Surface  pressure  [hPa]
   tv_acraft  42  t,z,y,x  AIRCRAFT:      Virtual  temperature [K]
   qv_acraft  42  t,z,y,x  AIRCRAFT:      Specific humidity [g/kg]
    u_acraft  42  t,z,y,x  AIRCRAFT:      Zonal      wind [m/s]
    v_acraft  42  t,z,y,x  AIRCRAFT:      Meridional wind [m/s]
    u_prof    42  t,z,y,x  PROFILER:      Zonal      wind [m/s]
    v_prof    42  t,z,y,x  PROFILER:      Meridional wind [m/s]
    u_scat     0  t,y,x    SCATTEROMETER: Zonal      wind [m/s]
    v_scat     0  t,y,x    SCATTEROMETER: Meridional wind [m/s]
    w_ssmi     0  t,y,x    SSMI:          Wind      speed [m/s]
    u_amv     42  t,z,y,x  ATMOS MOTION VECTORS: Zonal      wind [m/s]
    v_amv     42  t,z,y,x  ATMOS MOTION VECTORS: Meridional wind [m/s]
  sst_sea      0  t,y,x    SEA SURFACE:   Sea surface temperature [K]
   tv_sea      0  t,y,x    SEA SURFACE:   Virtual     temperature [K]
   ps_sea      0  t,y,x    SEA SURFACE:   Surface      pressure [hPA]
    u_sea      0  t,y,x    SEA SURFACE:   Zonal      wind [m/s]
    v_sea      0  t,y,x    SEA SURFACE:   Meridional wind [m/s]
   qv_sea      0  t,y,x    SEA SURFACE:   Specific humidity [g/kg]
   tv_land     0  t,y,x    LAND SURFACE:  Virtual  temperature [K]
   ps_land     0  t,y,x    LAND SURFACE:  Surface  pressure  [hPA]
    u_land     0  t,y,x    LAND SURFACE:  Zonal      wind [m/s]
    v_land     0  t,y,x    LAND SURFACE:  Meridional wind [m/s]
   qv_land     0  t,y,x    LAND SURFACE:  Specific humidity [g/kg]
ENDVARS
