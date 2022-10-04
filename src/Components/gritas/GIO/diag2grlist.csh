#!/bin/csh
  set Usage   = "diag2grlist.csh ExpID YYYYMM [-o outbasen]"
  set Example = "diag2grlist.csh merra 200401 "

# Set modules
# -----------
# source ../Linux/bin/g5_modules

# Set defaults for options
# ------------------------
  set TRUE = 1; set FALSE = 0 
 
# Set options defined by user, if any
# -----------------------------------
  set ReqArgv = ()
  while ( $#argv > 0 )
    switch ( $argv[1] )
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
  if ( $#ReqArgv < 2 ) goto err
  set ExpID  = $ReqArgv[1]; shift ReqArgv
  set Date   = $ReqArgv[1]; shift ReqArgv

  set Year   = `echo $Date | awk '{print substr($1,1,4)}'`
  set Month  = `echo $Date | awk '{print substr($1,5,2)}'`

  set RDir   = `pwd`

# Check directories
# -----------------
  set Dir     = ${NOBACKUP}/${ExpID}/obs/Y${Year}/M${Month}
  if ( ! -e $Dir ) then
       echo " The directory, $Dir, does not exist."
       goto err
  endif
  if ( ! -r $Dir ) then
     echo " The directory, $File, is not set with read permission."
     goto err
  endif
  cd $Dir

# Check for existence of files
# ----------------------------
  set BName   = $ExpID.diag_conv_anl.${Year}${Month}
  set InFiles = ( `ls -1 $Dir | egrep $BName` )
  if ( $#InFiles == 0 ) then
     echo "No files with the basename, $BName, exists in the directory, $Dir"
     goto err
  endif
  foreach File ( $InFiles )
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
  set BaseName  = gritas.$$
  set RC_File   = ${RDir}/gritas_upconv_merra.rc
  set Gritas_Core_Opt = "-nlevs 50 -rc $RC_File -res d -ncf -ospl -lb -o $BaseName"

# For each given diag file
# ------------------------
  foreach File ( $InFiles )

#    .... get date tag
#    -----------------
#     set DateHr = `echo $InFile | awk '{n=split($1,a,".");print a[n]}'`

     ${RDir}/gritas.x -obs $Gritas_Core_Opt ${File}; set Status = $status
     if ( $Status ) then
        echo "Error status (= ${Status}) returned from gritas.x"
       /bin/rm -f ${BaseName}*.hdf
        goto err
     endif
    /bin/rm -f ${BaseName}*.hdf
  end

# All is well
# ----------- 
  exit 0

# Usage messages
# --------------
  err:
     echo "  usage: $Usage"
     echo "example: $Example"
  exit 1
