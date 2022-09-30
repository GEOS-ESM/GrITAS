#!/bin/csh 

  set ExpID   = d5124_m2_jan10

  set DAY_TABLE = ( 31 28 31 30 31 30 31 31 30 31 30 31 )
  set Date = 201205
  set DateE = 201205
  set YYYY = `echo $Date | cut -c 1-4`
  set   MM = `echo $Date | cut -c 5-6`
  set DD = 1
  set DAY_MAX = $DAY_TABLE[$MM]
  unsetenv argv
  unset argv
  setenv argv


  set WORK_DIR    = $NOBACKUP/$ExpID/MERRA2/SAT
  set OBS_DIR     = $WORK_DIR
  set STORAGE_DIR = $NOBACKUP/$ExpID/MERRA2
  set ARCHIVE_DIR     = /archive/u/dao_ops/MERRA2/gmao_ops/GEOSadas-5_12_4/$ExpID
  set RC_DIR = /discover/nobackup/rgovinda/source_motel/GrITAS-MERRA2_V2/GrITAS/src/Components/gritas/GIO
  mkdir -p $WORK_DIR

  /bin/rm -rf $OBS_DIR/Y$YYYY/M$MM
  if ( ! -e $OBS_DIR/Y$YYYY/M$MM ) then
    mkdir -p $OBS_DIR/Y$YYYY/M$MM 
    cd $ARCHIVE_DIR/obs/Y$YYYY/M$MM
    while ( $DD <= $DAY_MAX )
     set DDe = $DD
     if ( $DD < 10 ) then
      set DDe = 0$DD
     endif

     if ( $MM == 02 ) then
      set leap  = `/home/rgovinda/bin/get_leap $YYYY $MM`
      if ( $leap == 0 ) then
       set DDe = 29
      endif
     endif
     
    while ( $DD <= $DAY_MAX )
     set DDe = $DD
     if ( $DD < 10 ) then
      set DDe = 0$DD
     endif
     echo "DDe: $DDe" 
     dmget D$DDe/H*/*.ods
     cp D$DDe/H*/*.ods $OBS_DIR/Y$YYYY/M$MM/
     @ DD = $DD + 1
    end
  endif

  cd $WORK_DIR


  set BinDir  = /discover/nobackup/rgovinda/source_motel/GrITAS-MERRA2_V2/GrITAS/Linux/bin
  set gritas  = ${BinDir}/gritas.x

  source $BinDir/g5_modules


  /bin/rm satlist.$Date.txt instrument.txt
  $RC_DIR/get_sat.csh  $OBS_DIR/Y$YYYY/M$MM |& tee satlist.$Date.txt
  set SAT_TABLE = `cat  satlist.$Date.txt`
  $RC_DIR/get_instrument.csh $RC_DIR |& tee instrument.txt
  set INSTRUMENT_TABLE = `cat  instrument.txt`
  set RC_TABLE         = `echo $INSTRUMENT_TABLE`

# set INSTRUMENT_TABLE = ( "amsua_n15" "amsua_n18" "amsua_n19" "amsub_n16" "amsub_n17" "hirs3_n16" "hirs3_n17" "hirs4_metop-a" "hirs4_n18" "hirs4_n19" "iasi_metop-a" "mhs_metop-a" "mhs_n18" "mhs_n19" )
# set INSTRUMENT_TABLE = ( "iasi_metop-a" )
# set RC_TABLE         = `echo $INSTRUMENT_TABLE`

 set kount = 0
 foreach INSTRUMENT ( `echo $INSTRUMENT_TABLE` )
  @ kount = $kount + 1
  set RCI        = $RC_TABLE[$kount]
  set INFILE     = diag_${INSTRUMENT}
  set RC_FILE    = "gritas_${RCI}_merra2.rc"
  echo $RC_FILE

  set RC_File = $RC_DIR/rc_files2/$RC_FILE
  set RES     = "d"
  if ( $INSTRUMENT == "o3lev_aura" ) then
   set Gritas_Core_Opt  = "-nlevs 48 -rc $RC_File -res $RES -ncf -ospl -lb -nopassive"
  else
   set Gritas_Core_Opt  = "-rc $RC_File -res $RES -ospl -lb -nopassive"
  endif

# unlink obs
# ln -s /gpfsm/dnb02/projects/p53/merra2/intermediate/$ExpID/obs/Y$YYYY/M$MM/D$DD obs

    foreach Hour ( all 00 06 12 18 )

      if ( $Hour == all ) then
       set ods_Files = `ls -1 $OBS_DIR/Y$YYYY/M$MM/*${INSTRUMENT}*ods`
       set syn_tag = ""
      else
       set ods_Files = `ls -1 $OBS_DIR/Y$YYYY/M$MM/*${INSTRUMENT}*${Hour}z*ods`
       set syn_tag = "_${Hour}z"
      endif


      set out_fileo   = gritaso${Hour}
      /bin/rm -f ${out_fileo}.{bias,stdv,nobs}.nc4

      $gritas -obs -o $out_fileo $Gritas_Core_Opt ${ods_Files} &
    
 #    ... for o-f data

      set out_filef   = gritasf${Hour}
      /bin/rm -f ${out_filef}.{bias,stdv,nobs}.nc4

      $gritas -omf -o $out_filef $Gritas_Core_Opt ${ods_Files}  &

 #    ... for o-a data

      set out_filea   = gritasa${Hour}
     /bin/rm -f ${out_filea}.{bias,stdv,nobs}.nc4

      $gritas -oma -o $out_filea $Gritas_Core_Opt ${ods_Files} &   

 #    ... for bias data

      set out_fileb   = gritasb${Hour}
     /bin/rm -f ${out_fileb}.{bias,stdv,nobs}.nc4

      $gritas -obias -o $out_fileb $Gritas_Core_Opt ${ods_Files} &   
      wait

      mkdir -p $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM

     if ($INSTRUMENT == "sbuv2_nim07" || $INSTRUMENT == "sbuv2_n07" || $INSTRUMENT == "sbuv2_n11" || $INSTRUMENT == "sbuv2_n14" || $INSTRUMENT == "sbuv2_n16" || $INSTRUMENT == "sbuv2_n17"  ) then
      mv  ${out_fileo}.bias.nc4  $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.mean3d_obs_p.${YYYY}${MM}${syn_tag}.nc4
      mv  ${out_fileo}.nobs.nc4  $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.nobs3d_obs_p.${YYYY}${MM}${syn_tag}.nc4 
      mv  ${out_fileo}.stdv.nc4  $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.stdv3d_obs_p.${YYYY}${MM}${syn_tag}.nc4

      mv  ${out_filef}.bias.nc4  $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.mean3d_omf_p.${YYYY}${MM}${syn_tag}.nc4 
      mv  ${out_filef}.nobs.nc4  $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.nobs3d_omf_p.${YYYY}${MM}${syn_tag}.nc4
      mv  ${out_filef}.stdv.nc4  $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.stdv3d_omf_p.${YYYY}${MM}${syn_tag}.nc4

      mv  ${out_filea}.bias.nc4  $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.mean3d_oma_p.${YYYY}${MM}${syn_tag}.nc4 
      mv  ${out_filea}.nobs.nc4  $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.nobs3d_oma_p.${YYYY}${MM}${syn_tag}.nc4
      mv  ${out_filea}.stdv.nc4  $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.stdv3d_oma_p.${YYYY}${MM}${syn_tag}.nc4

      mv  ${out_fileb}.bias.nc4  $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.mean3d_bias_p.${YYYY}${MM}${syn_tag}.nc4
      mv  ${out_fileb}.nobs.nc4  $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.nobs3d_bias_p.${YYYY}${MM}${syn_tag}.nc4
      mv  ${out_fileb}.stdv.nc4  $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.stdv3d_bias_p.${YYYY}${MM}${syn_tag}.nc4
     else
      lats4d.sh -i ${out_fileo}.bias.nc4 -o $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.mean3d_obs_p.${YYYY}${MM}${syn_tag} -zrev
      lats4d.sh -i ${out_fileo}.nobs.nc4 -o $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.nobs3d_obs_p.${YYYY}${MM}${syn_tag} -zrev
      lats4d.sh -i ${out_fileo}.stdv.nc4 -o $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.stdv3d_obs_p.${YYYY}${MM}${syn_tag} -zrev

      run_ncrcat.csh $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.mean3d_obs_p.${YYYY}${MM}${syn_tag}.nc
      run_ncrcat.csh $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.nobs3d_obs_p.${YYYY}${MM}${syn_tag}.nc
      run_ncrcat.csh $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.stdv3d_obs_p.${YYYY}${MM}${syn_tag}.nc

      lats4d.sh -i ${out_filef}.bias.nc4 -o $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.mean3d_omf_p.${YYYY}${MM}${syn_tag} -zrev
      lats4d.sh -i ${out_filef}.nobs.nc4 -o $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.nobs3d_omf_p.${YYYY}${MM}${syn_tag} -zrev
      lats4d.sh -i ${out_filef}.stdv.nc4 -o $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.stdv3d_omf_p.${YYYY}${MM}${syn_tag} -zrev

      run_ncrcat.csh $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.mean3d_omf_p.${YYYY}${MM}${syn_tag}.nc
      run_ncrcat.csh $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.nobs3d_omf_p.${YYYY}${MM}${syn_tag}.nc
      run_ncrcat.csh $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.stdv3d_omf_p.${YYYY}${MM}${syn_tag}.nc

      lats4d.sh -i ${out_filea}.bias.nc4 -o $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.mean3d_oma_p.${YYYY}${MM}${syn_tag} -zrev
      lats4d.sh -i ${out_filea}.nobs.nc4 -o $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.nobs3d_oma_p.${YYYY}${MM}${syn_tag} -zrev
      lats4d.sh -i ${out_filea}.stdv.nc4 -o $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.stdv3d_oma_p.${YYYY}${MM}${syn_tag} -zrev

      run_ncrcat.csh $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.mean3d_oma_p.${YYYY}${MM}${syn_tag}.nc
      run_ncrcat.csh $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.nobs3d_oma_p.${YYYY}${MM}${syn_tag}.nc
      run_ncrcat.csh $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.stdv3d_oma_p.${YYYY}${MM}${syn_tag}.nc

      lats4d.sh -i ${out_fileb}.bias.nc4 -o $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.mean3d_bias_p.${YYYY}${MM}${syn_tag} -zrev
      lats4d.sh -i ${out_fileb}.nobs.nc4 -o $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.nobs3d_bias_p.${YYYY}${MM}${syn_tag} -zrev
      lats4d.sh -i ${out_fileb}.stdv.nc4 -o $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.stdv3d_bias_p.${YYYY}${MM}${syn_tag} -zrev

      run_ncrcat.csh $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.mean3d_bias_p.${YYYY}${MM}${syn_tag}.nc
      run_ncrcat.csh $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.nobs3d_bias_p.${YYYY}${MM}${syn_tag}.nc
      run_ncrcat.csh $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.stdv3d_bias_p.${YYYY}${MM}${syn_tag}.nc

      if ( -e $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.mean3d_obs_p.${YYYY}${MM}${syn_tag}.nc ) then
       mv  $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.mean3d_obs_p.${YYYY}${MM}${syn_tag}.nc  $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.mean3d_obs_p.${YYYY}${MM}${syn_tag}.nc4
      else
       mv  ${out_fileo}.bias.nc4 $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.mean3d_obs_p.${YYYY}${MM}${syn_tag}.nc4
      endif

      if ( -e $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.nobs3d_obs_p.${YYYY}${MM}${syn_tag}.nc ) then
       mv  $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.nobs3d_obs_p.${YYYY}${MM}${syn_tag}.nc  $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.nobs3d_obs_p.${YYYY}${MM}${syn_tag}.nc4
      else
       mv ${out_fileo}.nobs.nc4 $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.nobs3d_obs_p.${YYYY}${MM}${syn_tag}.nc4
      endif

      if ( -e $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.stdv3d_obs_p.${YYYY}${MM}${syn_tag}.nc ) then
       mv  $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.stdv3d_obs_p.${YYYY}${MM}${syn_tag}.nc  $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.stdv3d_obs_p.${YYYY}${MM}${syn_tag}.nc4
      else
       mv ${out_fileo}.stdv.nc4 $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.stdv3d_obs_p.${YYYY}${MM}${syn_tag}.nc4
      endif


      if ( -e $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.mean3d_omf_p.${YYYY}${MM}${syn_tag}.nc ) then
        mv $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.mean3d_omf_p.${YYYY}${MM}${syn_tag}.nc $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.mean3d_omf_p.${YYYY}${MM}${syn_tag}.nc4
      else
        mv ${out_filef}.bias.nc4  $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.mean3d_omf_p.${YYYY}${MM}${syn_tag}.nc4
      endif

      if ( -e $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.nobs3d_omf_p.${YYYY}${MM}${syn_tag}.nc ) then
       mv $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.nobs3d_omf_p.${YYYY}${MM}${syn_tag}.nc $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.nobs3d_omf_p.${YYYY}${MM}${syn_tag}.nc4
      else
       mv ${out_filef}.nobs.nc4 $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.nobs3d_omf_p.${YYYY}${MM}${syn_tag}.nc4
      endif

      if ( -e $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.stdv3d_omf_p.${YYYY}${MM}${syn_tag}.nc ) then
       mv $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.stdv3d_omf_p.${YYYY}${MM}${syn_tag}.nc $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.stdv3d_omf_p.${YYYY}${MM}${syn_tag}.nc4
      else
       mv ${out_filef}.stdv.nc4 $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.stdv3d_omf_p.${YYYY}${MM}${syn_tag}.nc4
      endif

      if ( -e $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.mean3d_oma_p.${YYYY}${MM}${syn_tag}.nc ) then
       mv $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.mean3d_oma_p.${YYYY}${MM}${syn_tag}.nc $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.mean3d_oma_p.${YYYY}${MM}${syn_tag}.nc4
      else
       mv ${out_filea}.bias.nc4 $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.mean3d_oma_p.${YYYY}${MM}${syn_tag}.nc4
      endif

      if ( -e  $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.nobs3d_oma_p.${YYYY}${MM}${syn_tag}.nc ) then
       mv $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.nobs3d_oma_p.${YYYY}${MM}${syn_tag}.nc $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.nobs3d_oma_p.${YYYY}${MM}${syn_tag}.nc4
      else
       mv ${out_filea}.nobs.nc4 $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.nobs3d_oma_p.${YYYY}${MM}${syn_tag}.nc4
      endif

      if ( -e $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.stdv3d_oma_p.${YYYY}${MM}${syn_tag}.nc ) then
       mv $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.stdv3d_oma_p.${YYYY}${MM}${syn_tag}.nc $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.stdv3d_oma_p.${YYYY}${MM}${syn_tag}.nc4
      else
       mv ${out_filea}.stdv.nc4 $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.stdv3d_oma_p.${YYYY}${MM}${syn_tag}.nc4
      endif

      if ( -e $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.mean3d_bias_p.${YYYY}${MM}${syn_tag}.nc ) then
       mv  $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.mean3d_bias_p.${YYYY}${MM}${syn_tag}.nc  $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.mean3d_bias_p.${YYYY}${MM}${syn_tag}.nc4
      else
       mv ${out_fileb}.bias.nc4 $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.mean3d_bias_p.${YYYY}${MM}${syn_tag}.nc4
      endif

      if ( -e $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.nobs3d_bias_p.${YYYY}${MM}${syn_tag}.nc ) then
       mv  $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.nobs3d_bias_p.${YYYY}${MM}${syn_tag}.nc  $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.nobs3d_bias_p.${YYYY}${MM}${syn_tag}.nc4
      else
       mv ${out_fileb}.nobs.nc4 $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.nobs3d_bias_p.${YYYY}${MM}${syn_tag}.nc4
      endif

      if ( -e $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.stdv3d_bias_p.${YYYY}${MM}${syn_tag}.nc ) then
       mv  $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.stdv3d_bias_p.${YYYY}${MM}${syn_tag}.nc  $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.stdv3d_bias_p.${YYYY}${MM}${syn_tag}.nc4
      else
       mv ${out_fileb}.stdv.nc4 $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.stdv3d_bias_p.${YYYY}${MM}${syn_tag}.nc4
      endif
     endif

      n4zip $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.mean3d_obs_p.${YYYY}${MM}${syn_tag}.nc4
      n4zip $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.nobs3d_obs_p.${YYYY}${MM}${syn_tag}.nc4
      n4zip $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.stdv3d_obs_p.${YYYY}${MM}${syn_tag}.nc4


      n4zip $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.mean3d_omf_p.${YYYY}${MM}${syn_tag}.nc4
      n4zip $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.nobs3d_omf_p.${YYYY}${MM}${syn_tag}.nc4
      n4zip $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.stdv3d_omf_p.${YYYY}${MM}${syn_tag}.nc4


      n4zip $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.mean3d_oma_p.${YYYY}${MM}${syn_tag}.nc4
      n4zip $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.nobs3d_oma_p.${YYYY}${MM}${syn_tag}.nc4
      n4zip $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.stdv3d_oma_p.${YYYY}${MM}${syn_tag}.nc4

      n4zip $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.mean3d_bias_p.${YYYY}${MM}${syn_tag}.nc4
      n4zip $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.nobs3d_bias_p.${YYYY}${MM}${syn_tag}.nc4
      n4zip $STORAGE_DIR/$INSTRUMENT/$RES/Y$YYYY/M$MM/merra2.${INSTRUMENT}.stdv3d_bias_p.${YYYY}${MM}${syn_tag}.nc4
  end
 end
