#!/bin/csh
  set Usage   = "diag2gritas.csh ExpID YYYYMM[DD] [-fetch -restart -cleanup]"
  set Example = "diag2gritas.csh merra 200401 "

# Set modules
# -----------
  set RootDir  = ${NOBACKUP}/GEOSdas/GrITAS/`uname -s`
  set BinDir   = ${RootDir}/bin
  set RunDir   = `pwd`
  set ArchRoot = /archive/merra/dao_ops/production/GEOSdas-2_1_4
  set ExtraLib = /usr/local/toolworks/totalview.8.6.0-3/lib:/usr/local/intel/mpi/3.2.011/lib64:/usr/local/intel/mkl/9.1.023/lib/em64t:/usr/local/intel/comp/9.1.052/lib:/home/dkokron/play/CSAR_Vis/v1.3.6/lib
# source ${BinDir}/g5_modules
#  set n4zip   = $SHARE/dasilva/bin/n4zip
  if ( $?LD_LIBRARY_PATH ) then
     setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:$ExtraLib
  else
     setenv LD_LIBRARY_PATH $ExtraLib
  endif
  printenv  LD_LIBRARY_PATH
# set n4zip   = ${HOME}/bin/n4zip
  set n4zip   = /discover/nobackup/projects/gmao/share/dasilva/bin/n4zip
  set gritas  = ${BinDir}/gritas.x
  set grmeans = ${BinDir}/GFIO_mean_r8.x

# Set defaults for options
# ------------------------
  set TRUE = 1; set FALSE = 0 
  set Fetch_ObFiles = $FALSE
  set Restart       = $FALSE
  set CleanUp       = $FALSE
  set Comp_Means    = $FALSE

# Set options defined by user, if any
# -----------------------------------
  set ReqArgv = ()
  while ( $#argv > 0 )
    switch ( $argv[1] )
       case -means:
          set Comp_Means    = $TRUE
          breaksw

       case -fetch:
          set Fetch_ObFiles = $TRUE
          breaksw

       case -restart:
          set Restart       = $TRUE
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
  if ( $#ReqArgv < 2 ) goto err
  set ExpID    = $ReqArgv[1]; shift ReqArgv
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
  set WorkDir     = ${NOBACKUP}/${ExpID}/obs/Y${Year}/M${Month}
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

# If desired restart the data
# ---------------------------
  if ( $Restart == $TRUE ) then
     if ( -e ob_files && -d ob_files ) then
        set FileList = ( `ls -1 ob_files` )
        if ( $#FileList > 0 ) mv ob_files/* .
     endif
  endif

# If desired fetch the data from archive
# --------------------------------------
  set diag_BName = $ExpID.diag_conv_
  if ( $Fetch_ObFiles == $TRUE ) then
     set ArchDir = ${ArchRoot}/${ExpID}/obs/Y${Year}/M${Month}
     if ( ! -e ob_files ) then
#        dmget ${ArchDir}/${diag_BName}*.${DateFrag} &
        cp ${ArchDir}/${diag_BName}*.${DateFrag}* . ;set Status = $status
        if ( $status ) then
           echo "Error status (= ${Status}) returned from the command, cp ${ArchDir}/..."
           goto err
        endif
     else
        echo "The subdirectory, ob_files, already exists."
        echo "The script continues without fetching archived files from "
        echo "   the directory, $ArchDir"
     endif 
 endif

# Create a directory for processed ob files (if necessary)
# --------------------------------------------------------
  if      ( ! -e ob_files ) then
     mkdir ob_files
  else if ( ! -d ob_files ) then
     echo "The entity, ob_files, exists but it is not a directory"
     goto err
  endif

# Check for existence and permissions of ...
# ------------------------------------------
  set DateHrs    =  ( `ls -1 | egrep $diag_BName | egrep .$DateFrag | awk '{n=split($1,a,"."); print a[n]}' | sort -u`)
#  if ( $#DateHrs == 0 ) then
#     echo "No files with the basename, $BName, exists in the directory, $WorkDir"
#     goto err
#  endif

  foreach DateHr ( $DateHrs )

#    ... analysis files
#    ------------------
     set File = $ExpID.diag_conv_anl.$DateHr
     if ( ! -e $File ) then
        echo " The file, $File, does not exist."
        goto err
     endif
     if ( ! -r $File ) then
        echo " The file, $File, is not set with read permission."
        goto err
     endif

#    ... first guess files
#    ---------------------
     set File = $ExpID.diag_conv_ges.$DateHr
     if ( ! -e $File ) then
        echo " The file, $File, does not exist."
        goto err
     endif
     if ( ! -r $File ) then
        echo " The file, $File, is not set with read permission."
        goto err
     endif
  end

# Set base options
# ----------------
  set RC_File          = ${RunDir}/rc_files/gritas_upconv_merra.rc
  set Gritas_Core_Opt  = "-nlevs 50 -rc $RC_File -res d -ncf -ospl -lb -nopassive"

# For each given synotpic date and hour ...
# -----------------------------------------
  foreach DateHr ( $DateHrs )
     set Date          = `echo $DateHr | awk '{print substr ($1,1,8)}'`
     set Day           = `echo $DateHr | awk '{print substr ($1,7,2)}'`
     set Hour          = `echo $DateHr | awk '{print substr ($1,9,2)}'`

#    ... name the destination directory
#    ----------------------------------
     set DayDir        = D${Day}
     if ( ! -e $DayDir ) mkdir -p ${DayDir}

#    ... input observation files
#    ---------------------------
     set diag_anl_File = ${ExpID}.diag_conv_anl.$DateHr
     set diag_ges_File = ${ExpID}.diag_conv_ges.$DateHr

#    ... remove any gritas output files
#    ----------------------------------
    /bin/rm -f gritas.{bias,stdv,nobs}.hdf

#    ... run gritas for observation data
#    -----------------------------------
     $gritas -obs $Gritas_Core_Opt ${diag_anl_File}; set Status = $status
     if ( $Status ) then
        echo "Error status (= ${Status}) returned from ${gritas}"
        goto cleanup
     endif
     mv gritas.bias.hdf ${DayDir}/merra.mean3d_obs_p.${Date}_${Hour}z.hdf
     mv gritas.stdv.hdf ${DayDir}/merra.stdv3d_obs_p.${Date}_${Hour}z.hdf
     mv gritas.nobs.hdf ${DayDir}/merra.nobs3d_obs_p.${Date}_${Hour}z.hdf

#    ... for o-f data
#    ----------------
     $gritas -omf $Gritas_Core_Opt ${diag_ges_File}; set Status = $status
     if ( $Status ) then
        echo "Error status (= ${Status}) returned from ${gritas}"
        goto cleanup
     endif
     mv gritas.bias.hdf ${DayDir}/merra.mean3d_omf_p.${Date}_${Hour}z.hdf
     mv gritas.stdv.hdf ${DayDir}/merra.stdv3d_omf_p.${Date}_${Hour}z.hdf
     mv gritas.nobs.hdf ${DayDir}/merra.nobs3d_omf_p.${Date}_${Hour}z.hdf

#    ... for o-a data
#    ----------------
     $gritas -omf $Gritas_Core_Opt ${diag_anl_File}; set Status = $status
     if ( $Status ) then
        echo "Error status (= ${Status}) returned from ${gritas}"
        goto cleanup
     endif
     mv gritas.bias.hdf ${DayDir}/merra.mean3d_oma_p.${Date}_${Hour}z.hdf
     mv gritas.stdv.hdf ${DayDir}/merra.stdv3d_oma_p.${Date}_${Hour}z.hdf
     mv gritas.nobs.hdf ${DayDir}/merra.nobs3d_oma_p.${Date}_${Hour}z.hdf

#    ... compress the output data
#    ----------------------------
     $n4zip  ${DayDir}/merra.{mean,stdv,nobs}3d_{obs,omf,oma}_p.${Date}_${Hour}z.hdf; set Status = $status
     if ( $Status ) then
        echo "Error status (= ${Status}) returned from ${n4zip}"
       /bin/rm -f ${ExpID}.gritas_*.${Date}.hdf
        goto err
     endif

#    ... move processed input files to subdirectory
#    ---------------------------------------------- 
     mv $diag_anl_File  $diag_ges_File ob_files
  end

# Name the portal directory
# -------------------------
#  set PortArchDir = /portal/gmao_ops/gds/merra/obs
#  set PortArchDir = /portal/MERRA/obs
  set PortArchDir = /portal/MERRA/obs
  if ( ! -e $PortArchDir ) then
       echo " The directory, $PortArchDir, does not exist."
       goto err
  endif
  if ( ! -w $PortArchDir ) then
     echo " The directory, $PortArchDir, is not set with write permission."
     goto err
  endif

  set PortRootDir = ${PortArchDir}/conv/Y${Year}/M${Month}
  if ( ! -e $PortRootDir ) then
     mkdir -p ${PortRootDir}; set Status = $status
     if ( $status ) then
        echo "Error status (= ${Status}) returned from the command, mkdir $PortRootDir"
        goto err
     endif
  endif

# Transfer all files to the portal directory
# ------------------------------------------
  set DayDir_List  = `ls -1d D??`
  set Mon_OutFiles = `ls -1 | egrep "merra.mon_"`
  cp -r $DayDir_List ${PortRootDir}; set Status = $status 
  if ( $status ) then
     echo "Error status (= ${Status}) in copying the output files to the directory, ${PortRootDir}"
     goto err
  endif
  /bin/rm -fr $DayDir_List

# Remove all ob file if desired
# -----------------------------
  if ( $CleanUp == $TRUE ) then
     /bin/rm -fr ob_files mon_means
  endif

# Compute the monthly statistics if desired
# -----------------------------------------
  if ( $Comp_Means == $TRUE ) then
     cd $PortArchDir
     csh -vx ${RunDir}/gritas2means.csh ${Year}${Month} -r means
     if ( $status ) then
        echo " Error status returned from gritas2means (-r = means)" 
        echo "    Year/Mon = ${Year}${Month}"
        exit 1
     endif
     csh -vx ${RunDir}/gritas2means.csh ${Year}${Month} -r rms
     if ( $status ) then
        echo " Error status returned from gritas2means (-r = rms)" 
        echo "    Year/Mon = ${Year}${Month}"
        exit 1
     endif
     csh -vx ${RunDir}/gritas2means.csh ${Year}${Month} -r obrate
     if ( $status ) then
        echo " Error status returned from gritas2means (-r = obrate)" 
        echo "    Year/Mon = ${Year}${Month}"
        exit 1
     endif
  endif

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
