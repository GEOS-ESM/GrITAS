#!/bin/csh 



  set Date = 20110601
  set DateE = 20110630
  set YYYY0 = `echo $Date | cut -c1-4`
  set   MM0 = `echo $Date | cut -c5-6`
  set   DD0 = `echo $Date | cut -c7-8`
  set TRUE = 1; set FALSE = 0
  set Comp_Means = $FALSE

  set ExpID   = d5124_m2_jan10
  set ExpID2  = merra2

  set WORK_DIR    = $NOBACKUP/$ExpID/MERRA2
  set STORAGE_DIR = $NOBACKUP/$ExpID/MERRA2/conv
  set OBS_DIR     = /archive/u/dao_ops/MERRA2/gmao_ops/GEOSadas-5_12_4/$ExpID
  set RES     = "d"

  if ( ! -e $WORK_DIR/obs/Y$YYYY0/M$MM0 ) then
   mkdir -p $WORK_DIR/obs/Y$YYYY0/M$MM0
   cd $OBS_DIR/obs/Y$YYYY0/M$MM0
   cp D*/H*/*diag_conv*.ods $WORK_DIR/obs/Y$YYYY0/M$MM0/
  endif

  cd $NOBACKUP/$ExpID/MERRA2

  setenv BinDir  /discover/nobackup/rgovinda/source_motel/GrITAS-MERRA2/GrITAS/Linux/bin
  setenv RunDir  /discover/nobackup/rgovinda/source_motel/GrITAS-MERRA2/GrITAS/src/Components/gritas/GIO
  setenv PortArchDir   $WORK_DIR/conv/$RES
  setenv Dir           $PortArchDir

  set gritas  = ${BinDir}/gritas.x

  source $BinDir/g5_modules


  set RC_File = /discover/nobackup/rgovinda/source_motel/GrITAS-MERRA2/GrITAS/src/Components/gritas/GIO/rc_files/gritas_upconv_merra2.rc
# set Gritas_Core_Opt  = "-nlevs 1 -rc $RC_File -res $RES -ncf -ospl -lb -nopassive"
  set Gritas_Core_Opt  = "-nlevs 15 -rc $RC_File -res $RES -ospl -lb -nopassive"

# unlink obs
# ln -s /gpfsm/dnb02/projects/p53/merra2/intermediate/$ExpID/obs/Y$YYYY/M$MM/D$DD obs

  while ( $Date <= $DateE )
    set YYYY = `echo $Date | cut -c 1-4`
    set   MM = `echo $Date | cut -c 5-6`
    set  Day = `echo $Date | cut -c7-8`

    mkdir -p $STORAGE_DIR/$RES/Y$YYYY/M$MM/D$Day

    foreach Hour ( 00 06 12 18 )
      set DateHr = ${Date}${Hour}
      set DayDir = D${Day}

      set obs_File  = $ExpID.diag_conv.${Date}_${Hour}z.ods

      set out_fileo   = gritaso${Hour}
      /bin/rm -f ${out_fileo}.{bias,stdv,nobs}.nc4

      $gritas -obs -o $out_fileo $Gritas_Core_Opt obs/Y$YYYY/M$MM/${obs_File} &
    
 #    ... for o-f data

      set out_filef   = gritasf${Hour}
      /bin/rm -f ${out_filef}.{bias,stdv,nobs}.nc4

      $gritas -omf -o $out_filef $Gritas_Core_Opt obs/Y$YYYY/M$MM/${obs_File}  &

 #    ... for o-a data

      set out_filea   = gritasa${Hour}
     /bin/rm -f ${out_filea}.{bias,stdv,nobs}.nc4

      $gritas -oma -o $out_filea $Gritas_Core_Opt obs/Y$YYYY/M$MM/${obs_File} &   

#    ... for bias data

#     set out_fileb   = gritasb${Hour}
#    /bin/rm -f ${out_fileb}.{bias,stdv,nobs}.nc4

#     $gritas -obias -o $out_fileb $Gritas_Core_Opt obs/Y$YYYY/M$MM/${obs_File} &   
      wait

      mv ${out_fileo}.bias.nc4 $STORAGE_DIR/$RES/Y$YYYY/M$MM/D$Day/$ExpID2.mean3d_obs_p.${Date}_${Hour}z.nc4
      mv ${out_fileo}.nobs.nc4 $STORAGE_DIR/$RES/Y$YYYY/M$MM/D$Day/$ExpID2.nobs3d_obs_p.${Date}_${Hour}z.nc4
      mv ${out_fileo}.stdv.nc4 $STORAGE_DIR/$RES/Y$YYYY/M$MM/D$Day/$ExpID2.stdv3d_obs_p.${Date}_${Hour}z.nc4

      n4zip $STORAGE_DIR/$RES/Y$YYYY/M$MM/D$Day/$ExpID2.mean3d_obs_p.${Date}_${Hour}z.nc4
      n4zip $STORAGE_DIR/$RES/Y$YYYY/M$MM/D$Day/$ExpID2.nobs3d_obs_p.${Date}_${Hour}z.nc4
      n4zip $STORAGE_DIR/$RES/Y$YYYY/M$MM/D$Day/$ExpID2.stdv3d_obs_p.${Date}_${Hour}z.nc4


      mv ${out_filef}.bias.nc4 $STORAGE_DIR/$RES/Y$YYYY/M$MM/D$Day/$ExpID2.mean3d_omf_p.${Date}_${Hour}z.nc4
      mv ${out_filef}.nobs.nc4 $STORAGE_DIR/$RES/Y$YYYY/M$MM/D$Day/$ExpID2.nobs3d_omf_p.${Date}_${Hour}z.nc4
      mv ${out_filef}.stdv.nc4 $STORAGE_DIR/$RES/Y$YYYY/M$MM/D$Day/$ExpID2.stdv3d_omf_p.${Date}_${Hour}z.nc4

      n4zip $STORAGE_DIR/$RES/Y$YYYY/M$MM/D$Day/$ExpID2.mean3d_omf_p.${Date}_${Hour}z.nc4
      n4zip $STORAGE_DIR/$RES/Y$YYYY/M$MM/D$Day/$ExpID2.nobs3d_omf_p.${Date}_${Hour}z.nc4
      n4zip $STORAGE_DIR/$RES/Y$YYYY/M$MM/D$Day/$ExpID2.stdv3d_omf_p.${Date}_${Hour}z.nc4

      mv ${out_filea}.bias.nc4 $STORAGE_DIR/$RES/Y$YYYY/M$MM/D$Day/$ExpID2.mean3d_oma_p.${Date}_${Hour}z.nc4
      mv ${out_filea}.nobs.nc4 $STORAGE_DIR/$RES/Y$YYYY/M$MM/D$Day/$ExpID2.nobs3d_oma_p.${Date}_${Hour}z.nc4
      mv ${out_filea}.stdv.nc4 $STORAGE_DIR/$RES/Y$YYYY/M$MM/D$Day/$ExpID2.stdv3d_oma_p.${Date}_${Hour}z.nc4

      n4zip $STORAGE_DIR/$RES/Y$YYYY/M$MM/D$Day/$ExpID2.mean3d_oma_p.${Date}_${Hour}z.nc4
      n4zip $STORAGE_DIR/$RES/Y$YYYY/M$MM/D$Day/$ExpID2.nobs3d_oma_p.${Date}_${Hour}z.nc4
      n4zip $STORAGE_DIR/$RES/Y$YYYY/M$MM/D$Day/$ExpID2.stdv3d_oma_p.${Date}_${Hour}z.nc4

#     mv ${out_fileb}.bias.nc4 $STORAGE_DIR/$RES/Y$YYYY/M$MM/D$Day/$ExpID2.mean3d_bias_p.${Date}_${Hour}z.nc4
#     mv ${out_fileb}.nobs.nc4 $STORAGE_DIR/$RES/Y$YYYY/M$MM/D$Day/$ExpID2.nobs3d_bias_p.${Date}_${Hour}z.nc4
#     mv ${out_fileb}.stdv.nc4 $STORAGE_DIR/$RES/Y$YYYY/M$MM/D$Day/$ExpID2.stdv3d_bias_p.${Date}_${Hour}z.nc4

#     n4zip $STORAGE_DIR/$RES/Y$YYYY/M$MM/D$Day/$ExpID2.mean3d_bias_p.${Date}_${Hour}z.nc4
#     n4zip $STORAGE_DIR/$RES/Y$YYYY/M$MM/D$Day/$ExpID2.nobs3d_bias_p.${Date}_${Hour}z.nc4
#     n4zip $STORAGE_DIR/$RES/Y$YYYY/M$MM/D$Day/$ExpID2.stdv3d_bias_p.${Date}_${Hour}z.nc4
  end


     @ Date = $Date + 1
     set Date = `/home/$user/bin/get_date.x $Date`
 end

 # Compute the monthly statistics if desired
# -----------------------------------------
  if ( $Comp_Means == $TRUE ) then
     cd $PortArchDir
     csh -vx ${RunDir}/gritas2means.csh $ExpID2 ${YYYY0}${MM0} -r means
     if ( $status ) then
        echo " Error status returned from gritas2means (-r = means)" 
        echo "  ExpID  YYYY/Mon = ${YYYY0}${MM0}"
        exit 1
     endif
     csh -vx ${RunDir}/gritas2means.csh $ExpID2 ${YYYY0}${MM0} -r rms
     if ( $status ) then
        echo " Error status returned from gritas2means (-r = rms)" 
        echo "   ExpID YYYY/Mon = ${YYYY0}${MM0}"
        exit 1
     endif
     csh -vx ${RunDir}/gritas2means.csh $ExpID2 ${YYYY0}${MM0} -r obrate
     if ( $status ) then
        echo " Error status returned from gritas2means (-r = obrate)" 
        echo "  ExpID  YYYY/Mon = ${YYYY0}${MM0}"
        exit 1
     endif
  endif

