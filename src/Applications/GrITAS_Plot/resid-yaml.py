#!/usr/bin/env python
import yaml
import common
import sys, optparse
import pandas as pd
from plot_util import *

usage = "Usage: %prog [options] "
parser = optparse.OptionParser(usage);

parser.add_option("-d","--date",type=str,default='YYYYMMDD/YYYYMMDD',
                  help='Start/End dates of measurements - delimited via "/" (default = YYYYMMDD/YYYYMMDD)')
parser.add_option("-e","--exp",type=str,default='',
                  help='Yaml containing experiment(s) information')
parser.add_option("-r", "--statsInRegions", type=str, default='REG1/REG2/...',
                  help='Access statistics for these regions (default = REG1/REG2/...)')
parser.add_option('--usrDefRegions',type=str,default='',
                  help='CSV file specifying user defined lat/lon regions to consider - ignored if empty (default = '')')
parser.add_option("--dfs", action='store_true', default=False,
                  help='Create yaml for DFS')
parser.add_option("--impact", action='store_true', default=False,
                  help='Create yaml for observation impact')
parser.add_option("--resid", action='store_true', default=False,
                  help='Create yaml for observation residuals')
parser.add_option("--T", action='store_true', default=False,
                  help='Form a time series plot')
parser.add_option("--M", action='store_true', default=False,
                  help='Form a monthly plot of statistics')

# Parse the input arguments
(options, args) = parser.parse_args()

if int(options.dfs)+int(options.impact)+int(options.resid) != 1:
    raise ValueError("Must select either dfs, impact, or resid!")

if options.dfs: yamlOut='dfs.yaml'
elif options.impact: yamlOut='impact.yaml'
else: yamlOut='resid.yaml'

# Output serialization
out=open(yamlOut,'w')
myyam=common.myYML(out)


# Some default regions
#-----------------------------------------------------------------
defRegions={'glo': {'lon': (-180.0,180.0), 'lat': (-90.0,90.0)},\
            'nhe': {'lon': (0.0,360.0),    'lat': (20.0,90.0)},\
            'she': {'lon': (0.0,360.0),    'lat': (-90.0,-20.0)},\
            'tro': {'lon': (0.0,360.0),    'lat': (-20.0,20.0)},\
            'nam': {'lon': (-172.0,-52.0), 'lat': (16.0,72.0)}\
        }

# If user has supplied additional zones of interest, collect those and add to defRegions
#-----------------------------------------------------------------------------------------
if options.usrDefRegions:
    with open(options.usrDefRegions) as f:
        for l in f.readlines():
            name,lon,lat=l.split(' ')
            # lon/lat are strings - need tuple of floats
            lon=tuple(float(x) for x in lon.split(','))
            lat=tuple(float(x) for x in lat.split(','))
            defRegions.update({name: {'lon': lon, 'lat': lat}})
#-----------------------------------------------------------------------------------------

# Actual collection of regions to be considered - selects from defRegions and includes all in usrDefRegions
# ----------------------------------------------------------------------------------------------------------
regions=[Region(k,v['lon'],v['lat']) for k,v in defRegions.items() if k in options.statsInRegions.split('/')]

# Bundle all regions focused on into a super set
Universe=NestedDict('regions',regions)

# Instantiate instruments
instrList = ['amsuan15', 'amsuan19', 'atmsnpp']
Instruments = NestedDict('configure',[instrument(i,-10,10) for i in instrList])


# Declare temporal window to investigate
pandasDates = pd.date_range(start=options.date.split('/')[0],
                            end=options.date.split('/')[1],freq='D') #'MS')

# Get the experiments to consider
# --------------------------------
experiments = NestedDict('experiments');
with open(options.exp, 'r') as f:
    expYaml = yaml.load(f, Loader=yaml.FullLoader)
    experiments.fromYaml(expYaml['experiments'],cls='EXPERIMENT')

# Modify plot params
# -------------------
plotParams = PlotParams(regions=[r.name for r in regions],timeSeries=options.T,monthly=options.M,compVia='ratio')


Global=globalProps(instruments=Instruments,comparator=comparator(True),
                   experiments=experiments,plotParams=plotParams,
                   startDate=pandasDates[0].strftime('%Y-%m-%d'),
                   endDate=pandasDates[-1].strftime('%Y-%m-%d'),
                   obCnt=0,obType='atmsnpp')


Res=Residual(glob=Global,universe=Universe)
Res.serialize(myyam,out)
