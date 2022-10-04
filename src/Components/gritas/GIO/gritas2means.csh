#!/bin/csh
#  set Usage   = "gritas2means.csh ExpID YYYYMM [-r result] [-o outbasen]"
#  set Example = "gritas2means.csh merra 200401 "
  set Usage   = "gritas2means.csh YYYYMM [-r result]"
  set Example = "gritas2means.csh 200401 "

# Set modules
# -----------
# set RootDir = ${NOBACKUP}/GEOSdas/GrITAS/`uname -s`
# set BinDir  = ${RootDir}/bin

# set RootDir = $NOBACKUP/MERRA/run
# set BinDir  = $RootDir
  echo "gritas2means: $RootDir" 
  echo "gritas2means: $BinDir" 

  set n4zip = /discover/nobackup/projects/gmao/share/dasilva/bin/n4zip

  echo "gritas2means: $n4zip" 
  set grmeans = ${BinDir}/GFIO_mean_r8.x

# Set defaults for options
# ------------------------
  set TRUE = 1; set FALSE = 0 
  set Result  = means
 
# Set options defined by user, if any
# -----------------------------------
  set ReqArgv = ()
  while ( $#argv > 0 )
    switch ( $argv[1] )
       case -r:
          set Result =  $argv[2]; shift
          breaksw
#       case -o:
#          set OutBaseName = $argv[2]
#          shift
#          breaksw
       
       default:
          set FirstChar = `echo $argv[1] | awk '{ print substr ($1,1,1)}'`
          if ( FirstChar == "-" ) then
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
#  if ( $#ReqArgv < 2 ) goto err
  if ( $#ReqArgv < 1 ) goto err
#  set ExpID  = $ReqArgv[1]; shift ReqArgv
  set ExpID  = merra
  set Date   = $ReqArgv[1]; shift ReqArgv

  set Year   = `echo $Date | awk '{print substr($1,1,4)}'`
  set Month  = `echo $Date | awk '{print substr($1,5,2)}'`

# Check directories
# -----------------
#  set Dir     = ${NOBACKUP}/${ExpID}/obs/Y${Year}/M${Month}
#  set Dir     = /portal/MERRA/obs/Y${Year}/M${Month}
# set Dir     = /portal/MERRA/obs/conv/Y${Year}/M${Month}
# set Dir     = $NOBACKUP/MERRA/obs/conv/Y${Year}/M${Month}
  
  echo "gritas2means.csh: $Dir "

  if ( ! -e $Dir/Y$Year/M$Month ) then
       echo " The directory, $Dir, does not exist."
       goto err
  endif
  if ( ! -r $Dir/Y$Year/M$Month ) then
     echo " The directory, $Dir, is not set with read permission."
     goto err
  endif
  if ( ! -w $Dir/Y$Year/M$Month ) then
     echo " The directory, $Dir, is not set with write permission."
     goto err
  endif

  cd $Dir/Y$Year/M$Month

# Implement result option
# -----------------------
  if      ( $Result == "means" ) then
     set Options       = ""
     set InFile_Result = "mean"
  else if ( $Result == "rms"   ) then
     set Options       = "-rms"
     set InFile_Result = "mean"
  else if ( $Result == "obrate"  ) then
     set Options       = ""
     set InFile_Result = "nobs"
  else
     echo "Invalid assignment for the variable Result (=$Result)"
     goto err
  endif

  set Quants   = ( obs omf oma )
  set SynTimes = ( 00 06 12 18 )
  set OutFiles = ""
  foreach Quant ( $Quants )
     set OutFile  = $ExpID.mon_${Result}_${Quant}.${Year}${Month}.hdf; /bin/rm -f $OutFile
     set InFiles  = D??/$ExpID.${InFile_Result}3d_${Quant}_p.${Year}${Month}??_??z.hdf
     set OutFiles = ( $OutFiles $OutFile )
     $grmeans $Options $InFiles -o $OutFile; set Status = $status
     if ( $Status ) then
        echo "Error status (= ${Status}) returned from $grmeans"
       /bin/rm -f $OutFiles
        goto err
     endif

     foreach SynTime ( $SynTimes )
        set OutFile  = $ExpID.mon_${Result}_${Quant}.${Year}${Month}_${SynTime}z.hdf; /bin/rm -f $OutFile
        set InFiles  = D??/$ExpID.${InFile_Result}3d_${Quant}_p.${Year}${Month}??_${SynTime}z.hdf
        set OutFiles = ( $OutFiles $OutFile )
        $grmeans $Options $InFiles -o $OutFile; set Status = $status
        if ( $Status ) then
           echo "Error status (= ${Status}) returned from $grmeans"
          /bin/rm -f $OutFiles
           goto err
        endif
     end
  end

  $n4zip $OutFiles
  if ( $Status ) then
     echo "Error status (= ${Status}) returned from ${n4zip}"
    /bin/rm -f $OutFiles
     goto err
  endif

# All is well
# ----------- 
  exit 0

# Usage messages
# --------------
  err:
     echo "  usage: $Usage"
     echo "example: $Example"
  exit 1
