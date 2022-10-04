#!/bin/csh

      set BIN_DIR = $NOBACKUP/MERRA/run

      set YEAR_TABLE  = ( 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 )
      foreach YYYY ( `echo $YEAR_TABLE` )
       set MONTH_TABLE = ( 01 02 03 04 05 06 07 08 09 10 11 12 )
       touch $BIN_DIR/conv.$YYYY.log
       set stime0 = `time`
       foreach MM ( `echo $MONTH_TABLE` )
        echo $YYYY $MM 
         set stime = `time` 
        $BIN_DIR/set_job $YYYY $MM
         set etime = `time` 
         echo " $YYYY $MM  stime: $stime etime: $etime " >> $BIN_DIR/conv.$YYYY.log
       end
       set etime0 = `time`
       echo " ------------------------------------------------" >> $BIN_DIR/conv.$YYYY.log
       echo " $YYYY stime0: $stime0 etime0: $etime0 " >> $BIN_DIR/conv.$YYYY.log
      end
