#!/bin/bash

# Report/set location of observation system data
#-------------------------------------------------
obsSysLoc()
{
    if [ $# -eq 0 ]; then
	echo "Accepted observing systems: GEOSIT, GEOSFP, MERRA2"
    elif [ $1 == 'GEOSIT' ]; then echo "/home/dao_ops/d5294_geosit_jan18/run/.../archive/obs"
    elif [ $1 == 'GEOSFP' ]; then echo "/home/dao_ops/f516_fp/run/.../archive/obs"
    elif [ $1 == 'MERRA2' ]; then echo "/home/dao_ops/d5124_m2_jan10/run/.../archive/obs"
    else
	echo "Unsupported observing system"
	exit 2
    fi
}

# Function to set list prefix
#------------------------------
getList()
{
    if   [ $1 == 'geosit' ]; then echo GEOSIT
    elif [ $1 == 'geosfp' ]; then echo GEOSFP
    elif [ $1 == 'merra2' ]; then echo MERRA2
    else
	echo "Unmatched instrument = $1"
	exit 1
    fi
}

# Dump ExpID, System, Instrument, and Date - flushed to file otherwise
#------------------------------------------------------------------------
verbose()
{
    echo "Exp ID     = $1"
    echo "System     = $2"
    echo "Instrument = $3"
    echo "DATE       = $4"
    exit 2
}


if [ $# -ne 7 ]; then
    echo "Usage: $0 <YR INI> <MNTH INI> <DAY INI> <YR FIN> <MNTH FIN> <DAY FIN> <OBS SYS>"
    echo `obsSysLoc`
    exit 1
fi


# Parse CL Args
#--------------------------
YR_I=$1; MNTH_I=$2; DAY_I=$3
YR_F=$4; MNTH_F=$5; DAY_F=$6


# Location of observations
#--------------------------
BASE=$(obsSysLoc $7); STATUS=$?
if [ $STATUS -ne 0 ]; then
    echo $BASE
    exit $STATUS
fi

# Directories in use
#---------------------
EXE_DIR=`pwd`

# Build up the list suffix
#---------------------------
SUFFIX=`printf "%d%02d%02d-%d%02d%02d" $YR_I $MNTH_I $DAY_I $YR_F $MNTH_F $DAY_F`


#######################################
# For this YYYYMMDD Date...
# Loop over Synoptic Times and find all ods files
#--------------------------------------
for YR in `seq $YR_I $YR_F`; do
    for MNTH in `seq -f %02g $MNTH_I $MNTH_F`; do
	for DAY in `seq -f %02g $DAY_I $DAY_F`; do

	    # Concat into a single dir to descend into
	    LOC=$BASE/Y$YR/M$MNTH/D$DAY

	    for SYNOPTIC in `seq -f %02g 0 6 18`; do

		# Access all 'diag' type ods files
		for ODS in `ls $LOC/H$SYNOPTIC/*diag*ods`; do #| grep diag.*ods`; do

		    # Split $ODS on '.' delimiter
		    # VARS=(`echo $ODS | awk '{split($0,a,"."); print a[4],a[5],a[6]}'`)y
		    VARS=(`cut -d'/' -f13- <<< $ODS | awk '{split($0,a,"."); print a[1],a[2],a[3]}'`)

		    # Assign
		    EXPID=${VARS[0]}
		    SYS=`cut -d'_' -f2 <<< $EXPID`
		    INSTRUMENT=`cut -d'_' -f2- <<< ${VARS[1]}`
		    DATE=`cut -d'_' -f1 <<< ${VARS[2]}`

		    # # Verbose output
		    # verbose $EXPID $SYS $INSTRUMENT $DATE
		    
		    # Set the list prefix
		    LIST=$(getList $SYS)
		    
		    # Set Header info of list file if non-existent
		    if [ ! -f $EXE_DIR/$LIST.$SUFFIX.H$SYNOPTIC.list ]; then
			echo "Instrument Date ExpID" >> $EXE_DIR/$LIST.$SUFFIX.H$SYNOPTIC.list
		    fi
		    echo "$INSTRUMENT $DATE $EXPID" >> $EXE_DIR/$LIST.$SUFFIX.H$SYNOPTIC.list

		done

	    done
	    #----------------------------------------------------------------------------

	done #DAY
    done #MNTH
done #YR
