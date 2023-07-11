#!/usr/bin/env python
import common
import sys, optparse
from plot_util import *

usage = "Usage: %prog [options] "
parser = optparse.OptionParser(usage);

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
    usrRegions=[]
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
Universe=Collection('region',regions)

# Instantiate instruments
instrList = ['amsuan15', 'amsuan19', 'atmsnpp']
Instruments = Collection('configure',[instrument(i,-10,10) for i in instrList])


Res=Residual(Instruments,monthlyComparator(True),Universe)
Res.serialize(myyam,out)
