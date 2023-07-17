#!/usr/bin/env python
import common
import sys, optparse
import pandas as pd
from plot_util import *

usage = "Usage: %prog [options] "
parser = optparse.OptionParser(usage);

parser.add_option("-d","--date",type=str,default='YYYYMMDD/YYYYMMDD',
                  help='Start/End dates of measurements - delimited via "/" (default = YYYYMMDD/YYYYMMDD)')
parser.add_option('-r','--usrDefRegions',type=str,default='',
                  help='CSV file specifying user defined lat/lon regions to consider - ignored if empty (default = '')')

# Parse the input arguments
(options, args) = parser.parse_args()

# Output serialization
out=open('regions.yaml','w')
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

regions=[Region(k,v['lon'],v['lat']) for k,v in defRegions.items()]

# Bundle all regions focused on into a super set
Universe=NestedDict('regions',regions)

# Instantiate instruments
instrList = ['amsuan15', 'amsuan19', 'atmsnpp']
Instruments = NestedDict('configure',[instrument(i,-10,10) for i in instrList])


# Declare temporal window to investigate
pandasDates = pd.date_range(start=options.date.split('/')[0],
                            end=options.date.split('/')[1],freq='MS')

Global=globalProps(instruments=Instruments,comparator=monthlyComparator(True),
                   startDate=pandasDates[0].strftime('%Y-%m-%d'),
                   endDate=pandasDates[-1].strftime('%Y-%m-%d'),
                   nicknames=['geosfp', 'geosfpp'],
                   expID=['f5294_fp','f5295_fpp'],
                   fileName='XYZ',obCnt=0,obType='atmsnpp',regions=['glo'])

Res=Residual(glob=Global,universe=Universe)
Res.serialize(myyam,out)
