#!/bin/csh 

   foreach SAT ( `cat *.txt` )
    set ncount =  `echo $SAT | wc -m`
    set INSTRUMENT = `echo $SAT | cut -c6-$ncount`
    if ( $INSTRUMENT != "conv" ) then
     echo $INSTRUMENT
    endif
   end

