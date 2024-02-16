#!/bin/bash

# Return explanation of arguments to user
# ---------------------------------------
explanation()
{
    EXP=$(cat <<EOF

   YYYYMMDD : Initial/final dates in form YYYYMMDD (e.g., 20170104 = 01/04/2017 or Jan. 4, 2017)
   PREFIX   : Prefix of absolute path to experiment
   ABRV     : User-provided shorthand name for experiment
EOF
)
   echo "$EXP"
}

# Check valid command line args
# -----------------------------
if [ $# -ne 4 ]; then
    echo "Usage: $0 <INI YYYYMMDD> <FIN YYYYMMDD> <PREFIX> <ABRV>"
    explanation
    exit 1
fi

# Parse CL Args
#--------------------------
YR_I=${1:0:4}; MNTH_I=${1:4:2}; DAY_I=${1:6:2}
YR_F=${2:0:4}; MNTH_F=${2:4:2}; DAY_F=${2:6:2}
PREFIX=$3

# Directories in use
#---------------------
EXE_DIR=`pwd`

# Build up the list suffix
#---------------------------
SUFFIX="$1-$2"

#######################################
# For this YYYYMMDD Date...
# Loop over Synoptic Times and find all ods files
#--------------------------------------
for YR in `seq $YR_I $YR_F`; do
    for MNTH in `seq -f %02g $MNTH_I $MNTH_F`; do
	for DAY in `seq -f %02g $DAY_I $DAY_F`; do

	    # Concat into a single dir to descend into
	    LOC=$PREFIX/Y$YR/M$MNTH/D$DAY

	    for SYNOPTIC in `seq -f %02g 0 6 18`; do

		# Access all 'diag' type ods files
		for ODS in `ls $LOC/H$SYNOPTIC/*diag*ods`; do

		    # Split $ODS on '.' delimiter
		    VARS=(`cut -d'/' -f13- <<< $ODS | awk '{split($0,a,"."); print a[1],a[2],a[3]}'`)

		    # Assign`
		    INSTRUMENT=`cut -d'_' -f2- <<< ${VARS[1]}`
		    DATE=`cut -d'_' -f1 <<< ${VARS[2]}`

		    # Set the list prefix
		    LIST=$4

		    # Set Header info of list file if non-existent
		    if [ ! -f $EXE_DIR/$LIST.$SUFFIX.H$SYNOPTIC.list ]; then
			echo "Instrument Date ExpID" >> $EXE_DIR/$LIST.$SUFFIX.H$SYNOPTIC.list
		    fi
		    echo "$INSTRUMENT $DATE $4" >> $EXE_DIR/$LIST.$SUFFIX.H$SYNOPTIC.list

		done

	    done
	    #----------------------------------------------------------------------------

	done #DAY
    done #MNTH
done #YR
