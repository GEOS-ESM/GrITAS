#!/bin/csh
  set Usage   = "name_archfiles.csh RootDir/ExpID Instr    YYYYMM ftype -t -ods -syn syn_time"
  set Example = "name_archfiles.csh d5_jan98      amsua_am 200501 ges"

  set RootDef = /archive/merra/dao_ops/production/GEOSdas-2_1_4

# Set defaults for options
# ------------------------
  set TRUE = 1; set FALSE = 0 
  set Tail_Names = $FALSE
  set syn_time   = all
  set format     = diag

# Set options defined by user, if any
# -----------------------------------
  set ReqArgv = ()
  while ( $#argv > 0 )
    switch ( $argv[1] )
       case -t:
          set Tail_Names = $TRUE
          breaksw

       case -syn:
          set syn_time   = $argv[2]; shift argv
          breaksw

       case -ods:
          set format     = ods
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
    shift argv
  end

# Get required parameters
# -----------------------
  if ( $#ReqArgv <  3 ) goto err
  set RootDir  = $ReqArgv[1]:h
  set ExpID    = $ReqArgv[1]:t;  shift ReqArgv
  set Instr    = $ReqArgv[1];    shift ReqArgv
  set Date     = $ReqArgv[1];    shift ReqArgv
  if ( $#ReqArgv >= 1 ) then
     set FileType = $ReqArgv[1]; shift ReqArgv
  else
     set FileType = anl
  endif

# For the case when no root directory path is given
# -------------------------------------------------
  if ( $RootDir == $ExpID ) set RootDir = $RootDef # ... set to default

# Name the archive directory
# --------------------------
  set Year     = `echo $Date | awk '{print substr($1,1,4)}'`
  set Month    = `echo $Date | awk '{print substr($1,5,2)}'`
  set ArchDir  = ${RootDir}/${ExpID}/obs/Y${Year}/M${Month}

# Set the names for diag formatted files
# --------------------------------------
  if ( $format == diag ) then

     set sat_band = `echo $Instr | awk '{n=split($1,a,"_"); print a[1]}'`
     set sat_name = `echo $Instr | awk '{n=split($1,a,"_"); print a[n]}'`

#    Set the instrument dependant parameters
#    ---------------------------------------
     set file_frags = ()
     switch ( $sat_band )
        case msu:
           if ( $sat_name == am )     set file_frags = ( msu.006   msu.008   msu.010   msu.012   msu.014   )
           if ( $sat_name == pm )     set file_frags = ( msu.005   msu.007   msu.009   msu.011   )
           if ( $sat_name == tirosn ) set file_frags = ( msu.005 )
           if ( $sat_name == noaa06 ) set file_frags = ( msu.006 )
           if ( $sat_name == noaa07 ) set file_frags = ( msu.007 )
           if ( $sat_name == noaa08 ) set file_frags = ( msu.008 )
           if ( $sat_name == noaa09 ) set file_frags = ( msu.009 )
           if ( $sat_name == noaa10 ) set file_frags = ( msu.010 )
           if ( $sat_name == noaa11 ) set file_frags = ( msu.011 )
           if ( $sat_name == noaa12 ) set file_frags = ( msu.012 )
           if ( $sat_name == noaa14 ) set file_frags = ( msu.014 )
           breaksw
        case ssu:
           if ( $sat_name == am )     set file_frags = ( ssu.006   ssu.008   ssu.014   )
           if ( $sat_name == pm )     set file_frags = ( ssu.005   ssu.007   ssu.009   ssu.011   )
           if ( $sat_name == tirosn ) set file_frags = ( ssu.005 )
           if ( $sat_name == noaa06 ) set file_frags = ( ssu.006 )
           if ( $sat_name == noaa07 ) set file_frags = ( ssu.007 )
           if ( $sat_name == noaa08 ) set file_frags = ( ssu.008 )
           if ( $sat_name == noaa09 ) set file_frags = ( ssu.009 )
           if ( $sat_name == noaa11 ) set file_frags = ( ssu.011 )
           if ( $sat_name == noaa14 ) set file_frags = ( ssu.014 )
           breaksw
        case hirs2:
           if ( $sat_name == am )     set file_frags = ( hirs2.006 hirs2.008 hirs2.010 hirs2.012 hirs2.014 )
           if ( $sat_name == pm )     set file_frags = ( hirs2.005 hirs2.007 hirs2.009 hirs2.011 )
           if ( $sat_name == tirosn ) set file_frags = ( hirs2.005 )
           if ( $sat_name == noaa06 ) set file_frags = ( hirs2.006 )
           if ( $sat_name == noaa07 ) set file_frags = ( hirs2.007 )
           if ( $sat_name == noaa08 ) set file_frags = ( hirs2.008 )
           if ( $sat_name == noaa09 ) set file_frags = ( hirs2.009 )
           if ( $sat_name == noaa10 ) set file_frags = ( hirs2.010 )
           if ( $sat_name == noaa11 ) set file_frags = ( hirs2.011 )
           if ( $sat_name == noaa12 ) set file_frags = ( hirs2.012 )
           if ( $sat_name == noaa14 ) set file_frags = ( hirs2.014 )
           breaksw
        case hirs3:
           if ( $sat_name == am )     set file_frags = ( hirs3.015 hirs3.017 )
           if ( $sat_name == pm )     set file_frags = ( hirs3.016 )
           if ( $sat_name == noaa15 ) set file_frags = ( hirs3.015 )
           if ( $sat_name == noaa16 ) set file_frags = ( hirs3.016 )
           if ( $sat_name == noaa17 ) set file_frags = ( hirs3.017 )
           breaksw
        case hirs4:
           if ( $sat_name == noaa18 ) set file_frags = ( hirs4.018 )
           breaksw
        case amsua:
           if ( $sat_name == am )     set file_frags = ( amsua.015 amsua.017 )
           if ( $sat_name == pm )     set file_frags = ( amsua.016 amsua.018 )
           if ( $sat_name == noaa15 ) set file_frags = ( amsua.015 )
           if ( $sat_name == noaa16 ) set file_frags = ( amsua.016 )
           if ( $sat_name == noaa17 ) set file_frags = ( amsua.017 )
           if ( $sat_name == noaa18 ) set file_frags = ( amsua.018 )
           if ( $sat_name == aqua   ) set file_frags = ( eos_amsua.049 )
           breaksw
        case amsub:
           if ( $sat_name == am )     set file_frags = ( amsub.015 amsub.017 )
           if ( $sat_name == pm )     set file_frags = ( amsub.016 )
           if ( $sat_name == noaa15 ) set file_frags = ( amsub.015 )
           if ( $sat_name == noaa16 ) set file_frags = ( amsub.016 )
           if ( $sat_name == noaa17 ) set file_frags = ( amsub.017 )
           breaksw
        case mhs:
           if ( $sat_name == pm )     set file_frags = ( mhs.018 )
           if ( $sat_name == noaa18 ) set file_frags = ( mhs.018 )
           breaksw
        case ssmi:
           if ( $sat_name == dmsp08 ) set file_frags = ( ssmi.008 )
           if ( $sat_name == dmsp10 ) set file_frags = ( ssmi.010 )
           if ( $sat_name == dmsp11 ) set file_frags = ( ssmi.011 )
           if ( $sat_name == dmsp13 ) set file_frags = ( ssmi.013 )
           if ( $sat_name == dmsp14 ) set file_frags = ( ssmi.014 )
           if ( $sat_name == dmsp15 ) set file_frags = ( ssmi.015 )
           breaksw
        case airs:
           if ( $sat_name == aqua )   set file_frags = ( airs.049 )
           breaksw
        case gsnd:
           if ( $sat_name == goes08 ) set file_frags = ( goes.008 )
           if ( $sat_name == goes10 ) set file_frags = ( goes.010 )
           if ( $sat_name == goes12 ) set file_frags = ( goes.012 )
           breaksw
        default:
     endsw
#     if ( $Instr == airs )      set file_frags = (  airs.049 )
#     if ( $Instr == eos_amsua ) set file_frags = ( amsua.049 )
     if ( $#file_frags == 0 ) then
        echo "Unrecognized instrument, $Instr"
        goto err
     endif

     set files = ""
     foreach name_frag ( $file_frags )
        set frag = diag_${name_frag}.$FileType
        if ( $syn_time == all ) then
           set filesN = ( `ls -1 $ArchDir | egrep $frag` )
        else
           set filesN = ( `ls -1 $ArchDir | egrep $frag | egrep ${syn_time}$` )
        endif
        set files = ( $files $filesN )
     end 

# ... for ods formatted files
# ---------------------------
  else if ( $format == ods ) then
     set frag = ${ExpID}.ana.obs.
     if ( $syn_time == all ) then
        set files = ( `ls -1 $ArchDir | egrep $frag | egrep z.ods$` )
     else
        set files = ( `ls -1 $ArchDir | egrep $frag | egrep _${syn_time}z.ods$` )
     endif

  endif

# Add full path (if desired)
# --------------------------
  if ( $Tail_Names == $FALSE ) then
     set file_list = ""
     foreach file ($files )
        set file_list = ( $file_list ${ArchDir}/$file )
     end

  else
     set file_list = ( $files )

  endif

  echo $file_list

# All is well
# ----------- 
  exit 0

# Usage messages
# --------------
  cleanup:

  err:
     echo "  usage: $Usage"
     echo "example: $Example"
  exit 1
