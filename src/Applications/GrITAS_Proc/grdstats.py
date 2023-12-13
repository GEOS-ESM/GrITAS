#!/usr/bin/env python3

import argparse
import grit_util

parser = argparse.ArgumentParser()
parser.add_argument('ymlfn', type=str)
args = parser.parse_args()


# Instantiate and parse yaml
#-----------------------------
gritas=grit_util.Gritas(args.ymlfn)

# Simple calendar
#------------------
cal=grit_util.calendar(gritas.dateStart, gritas.dateEnd,freq='D')

# Files to access
#------------------
files = [date.strftime(gritas.expDir) for date in cal.dates]

# Conventional obs
#-------------------
convObs=['upconv','upconv2','gps_100lev']

# Vars for creating directory structure for GrITAS output
#---------------------------------------------------------
year_i=cal.years[0];  month_i=cal.months[0];  day_i=cal.days[0]
year_f=cal.years[-1]; month_f=cal.months[-1]; day_f=cal.days[-1]
day_i, day_f, month_i, month_f = map(lambda x: str(x).rjust(2,'0'), [day_i, day_f, month_i, month_f])

# Iterate over all combinations of instruments and residuals, and files in yaml file
#-----------------------------------------------------------------------------------
for res in gritas.resids:
  gritas.wopts(res,year_i,month_i,day_i,year_f,month_f,day_f)
  for instrument in gritas.instruments:
    print("GrITAS reads %s data along %s residual"%(instrument,res))

    # Form Gritas instance's input/output options and call 'out' to execute options at command line
    #-----------------------------------------------------------------------------------------------
    # If instrument == conv, use non-standard named conventional ops rc files
    # Note list a = list b only copies reference to list
    instrs = convObs[:] if instrument == 'conv' else [instrument]

    # Treat instrs list as a stack, calling a separate instance of gritas for each member
    while len(instrs) > 0:
      instr=instrs.pop()
      gritas.iopts(res,instr)
      gritas.nc4Opts(instr)
      gritas.out(files,instrument)
