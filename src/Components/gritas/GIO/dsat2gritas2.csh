#!/bin/csh
  set Usage   = "dsat2gritas2.csh ExpID Instr    YYYYMM [-fetch -cleanup]"
  set Example = "dsat2gritas2.csh merra amsua_am 200401 "

# Set modules
# -----------
# set RootDir  = ${NOBACKUP}/GEOSdas/GrITAS/`uname -s`
# set BinDir   = ${RootDir}/bin

  #set RunDir   = `pwd`
  set RootDir  = $RunDir
  #set BinDir  = $RunDir

  echo " RunDir $RunDir "
  echo " BinDir $BinDir "

  #set ArchRoot = /archive/merra/dao_ops/production/GEOSdas-2_1_4
  echo "ArchRoot $ArchRoot"
#  set ExtraLib = /usr/local/toolworks/totalview.8.6.0-3/lib:/usr/local/intel/mpi/3.2.011/lib64:/usr/local/intel/mkl/9.1.023/lib/em64t:/usr/local/intel/comp/9.1.052/lib:/home/dkokron/play/CSAR_Vis/v1.3.6/lib
######  source ${BinDir}/g5_modules
#  set n4zip   = $SHARE/dasilva/bin/n4zip
#  if ( $?LD_LIBRARY_PATH ) then
#     setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:$ExtraLib
#  else
#     setenv LD_LIBRARY_PATH $ExtraLib
#  endif
  printenv  LD_LIBRARY_PATH
# set n4zip   = ${HOME}/bin/n4zip
  #set n4zip   = /discover/nobackup/projects/gmao/share/dasilva/bin/n4zip
   echo " n4zip $n4zip "

  set gritas  = ${BinDir}/gritas.x
  set grmeans = ${BinDir}/GFIO_mean_r8.x

# Set defaults for options
# ------------------------
  set TRUE = 1; set FALSE = 0 
  set Fetch_ObFiles = $FALSE
  set CleanUp       = $FALSE

# Set options defined by user, if any
# -----------------------------------
  set ReqArgv = ()
  while ( $#argv > 0 )
    switch ( $argv[1] )
       case -fetch:
          set Fetch_ObFiles = $TRUE
          breaksw

       case -cleanup:
          set CleanUp       = $TRUE
          breaksw
       
       default:
          set FirstChar = `echo $argv[1] | awk '{ print substr ($1,1,1)}'`
          if ( "$FirstChar" == "-" ) then
                             # Any other option produces an error
             echo "Illegal option "$argv[1]
             goto err

          else               # ... or is a required argument
             set ReqArgv = ($ReqArgv $argv[1])

          endif
    endsw
    shift
  end

# Get required parameters
# -----------------------
  if ( $#ReqArgv < 3 ) goto err
  set ExpID    = $ReqArgv[1]; shift ReqArgv
  set Instr    = $ReqArgv[1]; shift ReqArgv
  set Date     = $ReqArgv[1]; shift ReqArgv

  set Year     = `echo $Date | awk '{print substr($1,1,4)}'`
  set Month    = `echo $Date | awk '{print substr($1,5,2)}'`

  set LDate    = `echo $Date | awk '{print length($1)}'`
  if ( $LDate < 8 ) then
     set DateFrag = ${Year}${Month}
  else
     set Day      = `echo $Date | awk '{print substr($1,7,2)}'`
     set DateFrag = ${Year}${Month}${Day}
     if ( $LDate  >= 10 ) then
        set Hour     = `echo $Date | awk '{print substr($1,9,2)}'`
        set DateFrag = ${DateFrag}${Hour}
     endif
  endif

# If fetching the data, create the working directory if it does not exists
# ------------------------------------------------------------------------
  set WorkDir     = ${NOBACKUP}/${ExpID}/obs/${Instr}/Y${Year}/M${Month}
  if ( $Fetch_ObFiles == $TRUE ) then
     if ( ! -e $WorkDir ) mkdir -p $WorkDir
  endif

# Set and check directories
# -------------------------
  if ( ! -e $WorkDir ) then
       echo " The directory, $WorkDir, does not exist."
       goto err
  endif
  if ( ! -w $WorkDir ) then
     echo " The directory, $WorkDir, is not set with write permission."
     goto err
  endif
  cd $WorkDir

# Name archived input data files
# ------------------------------
  set name_cmd_basic = "${RunDir}/name_archfiles.csh ${ExpID} ${Instr} ${Date}"
  set arch_anl_Files = ( `${name_cmd_basic} anl` )
  if ( $status ) then
     $name_cmd_basic anl 
     echo "Error status returned from the command, $name_cmd_basic anl"
     goto err
  endif
  set arch_ges_Files = ( `${name_cmd_basic} ges` )
  if ( $status ) then
     $name_cmd_basic ges
     echo "Error status returned from the command, $name_cmd_basic ges"
     goto err
  endif
 
# Nothing to do if there are no input data files
# ----------------------------------------------
  if ( $#arch_anl_Files == 0 ) then
     echo "No input files exist for the instrument, ${Instr} for the"
     echo "  year/month, ${Year}${Month} "
     set OutFiles = ()
     goto setup_portal
  endif

# If desired fetch the data from archive
# --------------------------------------
  if ( $Fetch_ObFiles == $TRUE ) then
     dmget $arch_anl_Files $arch_ges_Files &
     cp $arch_anl_Files $arch_ges_Files . ;set Status = $status
     if ( $status ) then
        echo "Error status (= ${Status}) returned from the command, cp ..."
        goto err
     endif
  endif 

# Set base options
# ----------------
  set SAT_RC_File          = ${RunDir}/rc_files/gritas_${Instr}_merra.rc
#  set Gritas_Core_Opt  = "-nlevs 50 -rc $RC_File -res d -ncf -ospl -lb -nopassive"
  #set Gritas_Core_Opt  = "-rc $SAT_RC_File -res d -ncf -ospl -lb -nopassive"
  set Gritas_Core_Opt  = "-rc $SAT_RC_File -hdf -res m -ncf -ospl -lb -nopassive"

# ... remove any gritas output files
# ----------------------------------
 /bin/rm -f gritas.{bias,stdv,nobs}.hdf

# ... list of output files generated from this script
# ---------------------------------------------------
  set OutFiles = ()

  foreach syn_time ( all 00 06 12 18 )
     if ( $syn_time == all ) then
        set anl_Files = ( `${name_cmd_basic} anl -t `)
        set ges_Files = ( `${name_cmd_basic} ges -t `)
        set syn_tag   = ""

     else
        set anl_Files = ( `${name_cmd_basic} anl -t -syn $syn_time`)
        set ges_Files = ( `${name_cmd_basic} ges -t -syn $syn_time`)
        set syn_tag   = _${syn_time}z
     endif

#    If no files for the given synoptic time ...
#    -------------------------------------------
     if ( $#anl_Files == 0 ) continue  # ... then process the data for
                                       #     the next synoptic time

#    ... run gritas for observation data
#    -----------------------------------
     $gritas -obs $Gritas_Core_Opt ${anl_Files}; set Status = $status
     if ( $Status ) then
        echo "Error status (= ${Status}) returned from ${gritas}"
        goto cleanup
     endif
     set mean_file = merra.mon_mean3d_${Instr}_obs.${Year}${Month}${syn_tag}.hdf
     set stdv_file = merra.mon_stdv3d_${Instr}_obs.${Year}${Month}${syn_tag}.hdf
     set nobs_file = merra.mon_nobs3d_${Instr}_obs.${Year}${Month}${syn_tag}.hdf
     mv gritas.bias.hdf $mean_file
     mv gritas.stdv.hdf $stdv_file
     mv gritas.nobs.hdf $nobs_file
     set OutFiles = ( $OutFiles $mean_file $stdv_file $nobs_file )

#    ... for o-f data
#    ----------------
     $gritas -omf $Gritas_Core_Opt ${ges_Files}; set Status = $status
     if ( $Status ) then
        echo "Error status (= ${Status}) returned from ${gritas}"
        goto cleanup
     endif
     set mean_file = merra.mon_mean3d_${Instr}_omf.${Year}${Month}${syn_tag}.hdf
     set stdv_file = merra.mon_stdv3d_${Instr}_omf.${Year}${Month}${syn_tag}.hdf
     set nobs_file = merra.mon_nobs3d_${Instr}_omf.${Year}${Month}${syn_tag}.hdf
     mv gritas.bias.hdf $mean_file
     mv gritas.stdv.hdf $stdv_file
     mv gritas.nobs.hdf $nobs_file
     set OutFiles = ( $OutFiles $mean_file $stdv_file $nobs_file )

#    ... for o-a data
#    ----------------
    $gritas -omf $Gritas_Core_Opt ${anl_Files}; set Status = $status
     if ( $Status ) then
        echo "Error status (= ${Status}) returned from ${gritas}"
        goto cleanup
     endif
     set mean_file = merra.mon_mean3d_${Instr}_oma.${Year}${Month}${syn_tag}.hdf
     set stdv_file = merra.mon_stdv3d_${Instr}_oma.${Year}${Month}${syn_tag}.hdf
     set nobs_file = merra.mon_nobs3d_${Instr}_oma.${Year}${Month}${syn_tag}.hdf
     mv gritas.bias.hdf $mean_file
     mv gritas.stdv.hdf $stdv_file
     mv gritas.nobs.hdf $nobs_file
     set OutFiles = ( $OutFiles $mean_file $stdv_file $nobs_file )
  end

# ... compress the output data
# ----------------------------
  if ( $#OutFiles > 0 ) then
     $n4zip  ${OutFiles}; set Status = $status
     if ( $Status ) then
        echo "Error status (= ${Status}) returned from ${n4zip}"
       /bin/rm -f ${ExpID}.gritas_*.${Date}.hdf
        goto err
     endif
  endif

# Name the portal directory
# -------------------------
  setup_portal:
#  set PortArchDir = /portal/gmao_ops/gds/merra/obs
# set PortArchDir = /portal/MERRA/obs
   echo $PortArchDir

  if ( ! -e $PortArchDir ) then
       echo " The directory, $PortArchDir, does not exist."
       goto err
  endif
  if ( ! -w $PortArchDir ) then
     echo " The directory, $PortArchDir, is not set with write permission."
     goto err
  endif

  set PortRootDir = ${PortArchDir}/${Instr}/Y${Year}/M${Month}
  if ( ! -e $PortRootDir ) then
     mkdir -p ${PortRootDir}; set Status = $status
     if ( $status ) then
        echo "Error status (= ${Status}) returned from the command, mkdir $PortRootDir"
        goto err
     endif
  endif

# Transfer all output files to the portal directory
# -------------------------------------------------
  if ( $#OutFiles > 0 ) then
     cp $OutFiles ${PortRootDir}; set Status = $status 
     if ( $status ) then
        echo "Error status (= ${Status}) in copying the output files to the directory, ${PortRootDir}"
        goto err
     endif
  endif

# Remove all ob file if desired
# -----------------------------
  set all_anl_Files = ( `${name_cmd_basic} anl -t `)
  set all_ges_Files = ( `${name_cmd_basic} ges -t `)
  if ( $CleanUp == $TRUE ) then
    /bin/rm -fr mon_means $all_anl_Files $all_ges_Files $OutFiles 
  endif

# Greate file flag to denote that the data is processed
# -----------------------------------------------------
  touch data_processed.note

# All is well
# ----------- 
  exit 0

# Usage messages
# --------------
  cleanup:
    /bin/rm -f gritas.{bias,stdv,nobs}.hdf

  err:
     echo "  usage: $Usage"
     echo "example: $Example"
  exit 1
