#!/bin/csh

  setenv RootDir       /discover/nobackup/$user/source_motel/gritas-merra-gio-04/GrITAS
  setenv BinDir        ${RootDir}/Linux/bin
  setenv RunDir        $RootDir/src/Components/gritas/GIO
  setenv ArchRoot      /archive/merra/dao_ops/production/GEOSdas-2_1_4
  setenv n4zip         /discover/nobackup/projects/gmao/share/dasilva/bin/n4zip
# setenv RC_File       ${RunDir}/rc_files/gritas_aircraft_merra.rc
  setenv RC_File       ${RunDir}/rc_files/gritas_upconv_merra2.rc
  #setenv PortArchDir   /discover/nobackup/$user/tai/obs
  setenv PortArchDir   /discover/nobackup/$user/MERRA/obs/conv/1152x721
  setenv Dir           $PortArchDir
  setenv RES           "e"
