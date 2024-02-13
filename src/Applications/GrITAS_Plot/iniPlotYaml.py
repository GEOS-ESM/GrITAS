#!/usr/bin/env python3

import yaml
import common
import os, sys, optparse
import pandas as pd
import defaults
from plot_util import *
from optparse import OptionValueError

usage = "Usage: %prog [options] "
parser = optparse.OptionParser(usage);

parser.add_option("-d","--date",type=str,default='YYYYMMDD/YYYYMMDD',
                  help='Start/End dates of measurements - delimited via "/" (default = YYYYMMDD/YYYYMMDD)')
parser.add_option("-e","--expIDs",type=str,default='/',
                  help='ID(s) for experiment(s) to plot gridded stats for - forward slash delimited (default = /)')
parser.add_option("-i",type=str,dest='instrument',action="callback",
                  default='atmsnpp',callback=defaults.avail_instruments,
                  help='Specify instrument to consider (default = atmsnpp)')
parser.add_option("-r", type=str, dest='statsInRegions',action="callback",
                  default='glo', callback=defaults.avail_regions,
                  help='Geographic regions (delimited via forward slash) over which statistics should be viewed (default = glo)')
parser.add_option("-s", "--statsToView", type=str, default='mean/stdv',
                  help='Statistics to view (default = mean/stdv)')
parser.add_option("--dfs", action='store_true', default=False,
                  help='Create yaml for DFS')
parser.add_option("--impact", action='store_true', default=False,
                  help='Create yaml for observation impact')
parser.add_option("--resid", action='store_true', default=False,
                  help='Create yaml for observation residuals')
parser.add_option("--compVia", type=str, default='ratio',
                  help='If two experiments are provided, compare them according this scheme (default = ratio)')
parser.add_option("--tSeriesVar", type=str, default='',
                  help="Specify stat to view time series for - ignored if options.T is False (default = '')")
parser.add_option('--usrDefRegions',type=str,default='',
                  help='CSV file specifying user defined lat/lon regions to consider - ignored if empty (default = '')')
parser.add_option("--C", action='store_true', default=False,
                  help='Form confidence intervals')
parser.add_option("--L", action='store_true', default=False,
                  help='Plot stats via a line plot; defaults to bars if omitted')
parser.add_option("--M", action='store_true', default=False,
                  help='Form a monthly plot of statistics')
parser.add_option("--T", action='store_true', default=False,
                  help='Form a time series plot')

# Parse the input arguments
(options, args) = parser.parse_args()

if int(options.dfs)+int(options.impact)+int(options.resid) != 1:
    raise OptionValueError("Must select either dfs, impact, or resid!")

if options.dfs:
    yamlOut='dfs.yaml'
    statsFlavor='DFS per Ob'
elif options.impact:
    yamlOut='impact.yaml'
    statsFlavor='Ob count'
else:
    yamlOut='resid.yaml'
    statsFlavor='Standard Deviation'

# Output serialization
out=open(yamlOut,'w')
myyam=common.YML(out)

# If user has supplied additional zones of interest, collect those and add to defaults.regions
#---------------------------------------------------------------------------------------------
if options.usrDefRegions:
    with open(options.usrDefRegions) as f:
        for l in f.readlines():
            name,lon,lat=l.split(' ')
            # lon/lat are strings - need tuple of floats
            lon=tuple(float(x) for x in lon.split(','))
            lat=tuple(float(x) for x in lat.split(','))
            defaults.regions.update({name: {'lon': lon, 'lat': lat}})

# Actual collection of regions to be considered - selects from defaults.regions and includes all in usrDefRegions
# ---------------------------------------------------------------------------------------------------------------
regions=[Region(k,v['lon'],v['lat']) for k,v in defaults.regions.items() if k in options.statsInRegions.split('/')]

# Bundle all regions focused on into a super set
Universe=Collection('regions',regions)

# Declare temporal window to investigate
pandasDates = pd.date_range(start=options.date.split('/')[0],
                            end=options.date.split('/')[1],freq='D')

# Form the experiments to consider
# -------------------------------
experiments = Collection('experiments')
for e in options.expIDs.split('/'):
    # Form an Experiment instance
    exp = Experiment(name=e,nickname=e,pathToFile='%s/%s/oma/%s/%s_gritas.nc4'%
                     (os.getcwd(),e,options.date.replace('/','-'),options.instrument))

    # Determine available instruments
    print(exp.getInstruments())
    _availInstruments = exp.getInstruments()

    # Add instance to collection
    experiments.append(exp)

# Instantiate instruments
# Instruments = Collection('instruments',[defaults.instruments[options.instrument]])
Instruments = Collection('instruments',[Instrument(name,_min=90.0,_max=110.0,vertUnits=unit) for name,unit in _availInstruments])

# Modify plot params
# ------------------
plotParams = PlotParams(regions=[r.name for r in regions],timeSeries=options.T,timeSeriesVar=options.tSeriesVar,
                        monthly=options.M,compVia=options.compVia,linePlot=options.L)

# Construct a GlobalProps instance
# --------------------------------
Global=GlobalProps(instruments=Instruments,experiments=experiments,plotParams=plotParams,
                   stats=Stats(flav=statsFlavor,measures=list(options.statsToView.split('/')),confInterval=options.C),
                   startDate=pandasDates[0].strftime('%Y-%m-%d'),
                   endDate=pandasDates[-1].strftime('%Y-%m-%d'),
                   obCnt=0,obType=options.instrument)

# Construct a StatsViewer instance from a GlobalProps and Collection (Universe) instance
# --------------------------------------------------------------------------------------
SV=StatsViewer(glob=Global,universe=Universe)
SV.serialize(myyam)
