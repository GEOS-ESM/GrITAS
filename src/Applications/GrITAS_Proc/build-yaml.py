#!/usr/bin/env python

import yaml
import pandas as pd
import optparse, os, sys
import collections

usage = "Usage: %prog [options] "
parser = optparse.OptionParser(usage);

parser.add_option("-c","--confidence",type=float,default=0.95,
                  help='Measurements confidence level (default = 0.95)')
parser.add_option("-d","--date",type=str,default='YYYYMMDD/YYYYMMDD',
                  help='Start/End dates of measurements - delimited via "/" (default = YYYYMMDD/YYYYMMDD)')
parser.add_option("-e","--prefExpDir",type=str,default='./',
                  help='Prefix of directories containing measurements (default = "./")')
parser.add_option("-g","--prefGrITAS",type=str,default='./',
                  help='Prefix of GrITAS src (default = ./)')
parser.add_option("-l","--sysObsListPrefix",type=str,default='foo.list',
                  help='Instrument/Date/ExpID list file prefix (default = foo)')
parser.add_option("-s","--system",type=str,default='',
                  help='Modeling system (default = '')')
parser.add_option("-t","--synopticTimes",type=str,default='00/06/12/18',
                  help='Right-justified, two integer synoptic times - delimited via "/" (default = 00/06/12/18)')
parser.add_option("--dryRun",action="store_true",default=False,
                  help='Build a yaml config for testing')

# Parse the input arguments
(options, args) = parser.parse_args()

conf=options.confidence
pandasDates=pd.date_range(start=options.date.split('/')[0],end=options.date.split('/')[1],freq='D')
startDate=pandasDates[0].strftime('%Y-%m-%d')
endDate=pandasDates[-1].strftime('%Y-%m-%d')
system=options.system.upper()


# List of instruments/dates/expIDs'
# Determined by opening each sysObsList file corresponding to desired synoptic time,
# ..and successively using pandas.merge (def. inner) to determine intersection of all available instrument data
instrDatesIDs=None
for n,t in enumerate(options.synopticTimes.split('/')):
    with open(options.sysObsListPrefix+'.H%s.list'%t,'r') as f:
        fromCSV=pd.read_csv(f,' ')
        instrDatesIDs=fromCSV if n == 0 else pd.merge(instrDatesIDs,fromCSV,how='inner')


# Instruments list
insList=list(instrDatesIDs['Instrument'])

# Get Experiment ID, ensuring it is unique
expID=instrDatesIDs['ExpID'].unique()
if len(expID) > 1:
    raise ValueError('Found multiple unique experiment IDs!')
expID=expID[0]


# Convenience
nSynoptic=len(options.synopticTimes.split('/'))
rcLoc=options.prefGrITAS + '/Components/gritas/etc/'
resids=['omf','oma']
expDirSuffix='/Y%Y/M%m/D%d/H'
if nSynoptic == 1:
    expDirSuffix += options.synopticTimes
elif nSynoptic == 4:
    expDirSuffix += '*'
else:
    expDirSuffix += '{'+','.join(options.synopticTimes.split('/'))+'}'
expDirSuffix += '/%s.diag_'%expID


# Of instruments passed, proceed only with those that have matched *rc & *ods files
#....this is done by ensuring an rc file exists
#....and ods files exist for all dates in range
uniqInstr=[]
for instr in insList:
    thisRCFile='%sgritas_%s.rc'%(rcLoc,instr)
    match=os.path.isfile(thisRCFile)

    # Check all dates in range
    for date in pandasDates:
        for synoptic in options.synopticTimes.split('/'):
            thisODSFile='%s/Y%i/M%s/D%s/H%s/%s.diag_%s.%s_%sz.ods'%(options.prefExpDir,
                                                                    date.year,
                                                                    str(date.month).rjust(2,'0'),
                                                                    str(date.day).rjust(2,'0'),
                                                                    synoptic,
                                                                    expID,
                                                                    instr,
                                                                    date.strftime('%Y%m%d'),
                                                                    synoptic)
            match&=os.path.isfile(thisODSFile)
    
    # If match == True, then append to uniqInstr (if unique) otherwise print Fail message
    if match and instr not in uniqInstr: uniqInstr.append(instr)
    else:
        print("> Failed to match rc file:\n\t\t %s\n  with ods files:\n\t\t %s"\
              %(thisRCFile,'%s/Y*/M*/D*/H*/%s.diag_*.*_*.ods'%(options.prefExpDir,expID)))



# Open a stream for yaml output
OUT=open('%s.gritas.yml'%system,'w')

# # Cleanest ... instantiate 'base' dict and dump immediately to stream
# base = { 'global'     : { 'start date' : startDate,
#                           'end date'   : endDate,
#                           'confidence' : conf },
#          'instruments': uniqInstr,
#          'dryrun'     : options.dryRun,
#          'expid'      : expID,
#          'expdir'     : options.prefExpDir + '/Y%Y/M%m/D%d/H*/*.diag_',
#          'outdir'     : system,
#          'rc location': options.prefGrITAS + '/Components/gritas/etc/',
#          'residuals'  : resids
#      }
# yaml.dump(base,stream=OUT,sort_keys=False,default_flow_style=False)


# To best match examples ... avoid usage of unordered dictionaries
# base = [ ('global', []), ('instruments', uniqInstr), ('dryrun', options.dryRun),\
base = [ ('global', {'start date': startDate, 'end date': endDate, 'confidence': conf}),\
         ('instruments', uniqInstr), ('dryrun', options.dryRun),\
         ('expid', expID),\
         ('expdir', options.prefExpDir + expDirSuffix),\
         ('outdir', system), ('rc location', options.prefGrITAS + '/Components/gritas/etc/'),\
         ('residuals', resids) ]

for k,v in base:
    # if k == 'global':
    #     yaml.dump({'global': [{'start date': startDate}, {'end date': endDate}, {'confidence': conf}]},stream=OUT)
    #     continue
    yaml.dump({k:v},stream=OUT,sort_keys=False,default_flow_style=False)

