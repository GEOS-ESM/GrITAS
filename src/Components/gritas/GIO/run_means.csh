#!/bin/csh


  set ExpID   = d5124_m2_jan79
  set RES     = "d"
  
  set Date = 198102
  set DateE = 198102
  set YYYY = `echo $Date | cut -c1-4`
  set   MM = `echo $Date | cut -c5-6`
  set TRUE = 1; set FALSE = 0
  set Comp_Means = $FALSE

  set    WORK_DIR = $NOBACKUP/$ExpID/MERRA2
  setenv BinDir  /discover/nobackup/rgovinda/source_motel/GrITAS-MERRA2/GrITAS/Linux/bin
  setenv RunDir  /discover/nobackup/rgovinda/source_motel/GrITAS-MERRA2/GrITAS/src/Components/gritas/GIO
  setenv PortArchDir   $WORK_DIR/conv/$RES
  setenv Dir           $PortArchDir


 # Compute the monthly statistics if desired
# -----------------------------------------
# if ( $Comp_Means == $TRUE ) then
     cd $PortArchDir
     csh -vx ${RunDir}/gritas2means.csh $ExpID ${YYYY}${MM} -r means
     if ( $status ) then
        echo " Error status returned from gritas2means (-r = means)" 
        echo "  ExpID  YYYY/MM  = ${YYYY}${MM}"
        exit 1
     endif
     csh -vx ${RunDir}/gritas2means.csh $ExpID ${YYYY}${MM} -r rms
     if ( $status ) then
        echo " Error status returned from gritas2means (-r = rms)" 
        echo "  ExpID  YYYY/MM  = ${YYYY}${MM}"
        exit 1
     endif
     csh -vx ${RunDir}/gritas2means.csh $ExpID ${YYYY}${MM} -r obrate
     if ( $status ) then
        echo " Error status returned from gritas2means (-r = obrate)" 
        echo " ExpID   YYYY/MM  = ${YYYY}${MM}"
        exit 1
     endif
# endif

