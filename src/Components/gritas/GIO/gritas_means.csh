#!/bin/csh
  set Usage   = "diag2grlist.csh ExpID YYYYMM [-o outbasen]"
  set Example = "diag2grlist.csh merra 200401 "

# Set modules
# -----------
# source ../Linux/bin/g5_modules
  set lats4d = $SHARE/dasilva/opengrads/Contents/lats4d.sh

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

# $SHARE/dasilva/opengrads/Contents/lats4d.sh -i http://opendap:9090/dods/MerraObs/6-hourly/o-f -time 12z15jan2000 12z15jan2000 -func 'mean(@,time=0z1jan2000,time=18z31jan2000)' -o u_scat.omf -vars u_scat -v

# $SHARE/dasilva/opengrads/Contents/lats4d.sh -i http://opendap:9090/dods/MerraObs/6-hourly/o-f -time 12z15jan2000 12z15jan2000 -func 'sum(if(@,==,0,-u,1),time=0z1jan2000,time=18z31jan2000)' -o u_scat.nobs -vars u_scat -v


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
