#!/bin/csh
  set Usage   = "run_dsat.csh2 Date_beg Date_end -ods -noclean -iver tag"
  set Example = "run_dsat.csh2. 200401   200512"

# Set important parameters
# ------------------------
  set RunDir    = `pwd`

# Set defaults for options
# ------------------------
  set TRUE = 1; set FALSE = 0 
  set iVerTag  = "noaa-all"
  set format   =  diag
  set CleanOpt = "-cleanup"

# Set options defined by user, if any
# -----------------------------------
  set ReqArgv = ()
  while ( $#argv > 0 )
    switch ( $argv[1] )
       case -iver:
          set iVerTag   =  $argv[2]; shift argv
          breaksw
       case -ods:
          set format    =   ods
          breaksw
       case -noclean:
          set CleanOpt  = ""
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
  if ( $#ReqArgv < 2 ) goto err
  set Date_beg = $ReqArgv[1]; shift ReqArgv
  set Date_end = $ReqArgv[1]; shift ReqArgv

# Set parameters for while loop
# -----------------------------
  set Year_beg   = `echo $Date_beg | awk '{print substr($1,1,4)}'`
  set Month_beg  = `echo $Date_beg | awk '{print substr($1,5,2)}'`
  set Year_end   = `echo $Date_end | awk '{print substr($1,1,4)}'`
  set Month_end  = `echo $Date_end | awk '{print substr($1,5,2)}'`
    @ JulMon_beg = $Year_beg * 12 + $Month_beg
    @ JulMon_end = $Year_end * 12 + $Month_end
    @ JulMon     = $JulMon_beg
  set ParmLists  = ()

# ... to determine the ExpID's and dates
# --------------------------------------
  while ( $JulMon <= $JulMon_end )
   @ Year  = ( $JulMon - 1 ) / 12
   @ Month = ( $JulMon - 1 ) % 12 + 1
     if ( $Month < 10 ) set Month = 0$Month
     set YYYYMM = $Year$Month 
     if      ( $Year < 1989 ) then
        set ExpID = d5_merra_jan79
     else if ( $Year < 1998 ) then
        set ExpID = d5_merra_jan89
     else
        set ExpID = d5_merra_jan98
     endif
     set ParmLists = ( $ParmLists ${ExpID}:${YYYYMM} )
#     echo $Parm
   @ JulMon++
  end

# ... version of instrument list (as set by the option -ver)
# ----------------------------------------------------------
  switch ( $iVerTag )
     case noaa-all:
        set InstrList =  (   msu_am      msu_pm \
                                         mhs_pm \
                             ssu_am      ssu_pm \
                           hirs2_am    hirs2_pm \
                           hirs3_am    hirs3_pm \
                           amsua_am    amsua_pm \
                           amsub_am    amsub_pm )
        breaksw
     case amsub_pm:
        set InstrList =  ( amsub_pm )
        breaksw
     case noaa_no_amsua_pm:
        set InstrList =  (   msu_am      msu_pm \
                                         mhs_pm \
                             ssu_am      ssu_pm \
                           hirs2_am    hirs2_pm \
                           hirs3_am    hirs3_pm \
                           amsua_am             \
                           amsub_am    amsub_pm )
        breaksw
     case mhs_pm:
        set InstrList =  (               mhs_pm )
        breaksw
     case msu:
        set InstrList =  ( msu_am        msu_pm )
        breaksw
     case amsua:
        set InstrList =  ( amsua_am      amsua_pm )
        breaksw
     case amsua_noaa:
        set InstrList =  ( amsua_noaa15  amsua_noaa16 \
                           amsua_noaa17  amsua_noaa18 )
        breaksw
     case amsub_noaa:
        set InstrList =  ( amsub_noaa15  amsub_noaa16 \
                           amsub_noaa17 )
        breaksw
     case msu_noaa:
        set InstrList =  ( msu_tirosn    msu_noaa06   \
                           msu_noaa07    msu_noaa08   \
                           msu_noaa09    msu_noaa10   \
                           msu_noaa11    msu_noaa12   \
                           msu_noaa14 )
        breaksw
     case ssu_noaa:
#        set InstrList =  ( ssu_tirosn    ssu_noaa06   \
#                           ssu_noaa07    ssu_noaa08   \
#                           ssu_noaa09                 \
#                           ssu_noaa11                 \
#                           ssu_noaa14 )
        set InstrList =  ( ssu_noaa09                 \
                           ssu_noaa11                 \
                           ssu_noaa14 )
        breaksw
     case hirs2_noaa:
        set InstrList =  ( hirs2_tirosn  hirs2_noaa06 \
                           hirs2_noaa07  hirs2_noaa08 \
                           hirs2_noaa09  hirs2_noaa10 \
                           hirs2_noaa11  hirs2_noaa12 \
                           hirs2_noaa14 )
        breaksw
     case hirs3_noaa:
        set InstrList =  ( hirs3_noaa15  hirs3_noaa16 \
                           hirs3_noaa17 )
        breaksw
     case ssmi_dmsp:
        set InstrList =  (  ssmi_dmsp08   ssmi_dmsp10 \
                            ssmi_dmsp11   ssmi_dmsp13 \
                            ssmi_dmsp14   ssmi_dmsp15 ) 
        breaksw
     case gsnd_goes:
        set InstrList =  (  gsnd_goes08   gsnd_goes10 \
                            gsnd_goes12 )
        breaksw
     case airs_aqua:
        set InstrList =  (  airs_aqua )
        breaksw
     default:
        set InstrList =  ( $iVerTag )
#        echo "Unrecognized version of instrument list, $iVerTag"
#        goto err
  endsw

# For each year/month ...
# -----------------------
  foreach ParmList ( $ParmLists )
     set ParmList = ( `echo $ParmList | awk '{n=split($1,a,":"); print a[1],a[2]}'` )
     set ExpID    = $ParmList[1]
     set YYYYMM   = $ParmList[2]

#    ... migrate the archived imput files for all instruments
#    --------------------------------------------------------
     if ( $format == diag ) then
        set ArchList = ()
        foreach Instr ( $InstrList )
           set name_cmd_basic = "${RunDir}/name_archfiles.csh ${ExpID} ${Instr} ${YYYYMM}"
           set anlFiles = ( `${name_cmd_basic} anl` )
           if ( $status ) then
              $name_cmd_basic anl 
              echo "Error status returned from the command, $name_cmd_basic anl"
              goto err
           endif
           set gesFiles = ( `${name_cmd_basic} ges` )
           if ( $status ) then
              $name_cmd_basic ges
              echo "Error status returned from the command, $name_cmd_basic ges"
              goto err
           endif
           set ArchList = ( $ArchList $anlFiles $gesFiles )
        end
        dmget $ArchList &

#       ... run the job for each instrument
#       -----------------------------------
        set Options = ( "-fetch" $CleanOpt )
        foreach Instr ( $InstrList )
           set cmd = "dsat2gritas2.csh $ExpID $Instr $YYYYMM $Options"
           csh -vx $cmd
           if ( $status ) then
              echo "Error status returned from the command, $cmd"
              goto err
           endif
        end 
     else if ( $format == ods ) then
        set First_Instr = $InstrList[1]
        set  Last_Instr = $InstrList[$#InstrList]
        foreach Instr ( $InstrList )
            set Options = ()
           if ( $Instr == $First_Instr ) set Options = ( $Options "-fetch" ) 
           if ( $Instr ==  $Last_Instr ) set Options = ( $Options $CleanOpt )
           set cmd = "dsat2gritas2_ods.csh $ExpID $Instr $YYYYMM $Options"
           csh -vx $cmd
           if ( $status ) then
              echo "Error status returned from the command, $cmd"
              goto err
           endif
        end 
     else
        echo "Invalid input file format $format"
        echo goto err
     endif

  end

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
