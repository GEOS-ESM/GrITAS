#!/bin/csh

   if ( $#argv < 1 ) then
    clear
    echo "USAGE: run_ncrcat.csh file_name"
    exit
   endif

   set FILE = $1
   set NCR_DIR = /usr/local/other/SLES11.1/nco/4.2.3/intel-12.1.0.233/bin
   $NCR_DIR/ncatted -O -h -a units,levels,m,c,"level" $FILE 
   $NCR_DIR/ncatted -O -h -a description,levels,m,c,"satellite channel" $FILE
   $NCR_DIR/ncatted -O -h -a type,levels,m,c,"channels" $FILE
   $NCR_DIR/ncatted -O -h -a long_name,levels,m,c,"satellite channel" $FILE
   $NCR_DIR/ncatted -O -h -a positive,levels,m,c,"up" $FILE
