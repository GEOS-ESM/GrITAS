#!/bin/csh

   set OBS_DIR = $1
   foreach FILE ( ` /bin/ls -1 $OBS_DIR/*01_00z.ods` )
    set SAT = $FILE:r:r:e
    echo $SAT
   end
