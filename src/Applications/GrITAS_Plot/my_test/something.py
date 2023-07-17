#!/usr/bin/env python

from helperFns import *
import plot_util
import common
import pandas as pd
import yaml

import optparse

parser = optparse.OptionParser()
parser.add_option("-i","--impacts",action="store_true",default=False,
                  help='Visualize impacts')
parser.add_option("-d","--dfs",action="store_true",default=False,
                  help='Visualize degrees for instrument')
parser.add_option("-r","--residuals",action="store_true",default=False,
                  help='Visualize residuals')
parser.add_option("-f","--ymlConfig",default='',
                  help='Yaml driver (default = "")')
# Parse the arguments
(options, args) = parser.parse_args()


# Only want to do one visualization at a time
if not options.impacts^options.dfs^options.residuals:
    raise ValueError("Scripts accepts only single visualization at a time!")


Global=plot_util.globalProps()
Regions=plot_util.NestedDict('regions')
with open(options.ymlConfig, 'r') as f:
    valsYML = yaml.load(f, Loader=yaml.FullLoader)
    Global.fromYaml(valsYML['global'])
    Regions.fromYaml(valsYML['regions'],cls='REGION') #plot_util.Region()

# Collect all read material into a single Residual instance
Res=plot_util.Residual(glob=Global,universe=Regions)


# Set temporal information
pandasDates = pd.date_range(start=Res.globalProps.startDate, end=Res.globalProps.endDate, freq='MS')
years = pandasDates.year; months = pandasDates.month

# Build up input files from templated file name
inputFiles = ([date.strftime(Res.globalProps.fileName) for date in pandasDates])

fn = inputFiles[0].replace("$expid",Res.globalProps.expID[0])
[nlev,nlat,nlon] = get_dims(fn)
numT = size(years)   # number of times to process
numReg = size(Res.globalProps.regions) # number of regions
numStats = len(Res.globalProps.stats.measures)   # number of statitics to calculate
numExp = size(Res.globalProps.expID)  # number of experiments to handle

total = zeros([nlev,numT,numReg,numStats])
if numExp==2:
   totals = zeros([nlev,numT,numReg,2])

print(fn)
print(totals)
sys.exit(8)

for ii,fname in enumerate(inputFiles):
  for nx,expid in enumerate(Res.globalProps.expID):

    fn = fname.replace("$expid",expid)
    print "this case now: ", fn if numExp>0 else ""

    yyyy = years[ii]
    mm = months[ii]
    mm =  str(mm).rjust(2, '0')

    if Res.globalProps.stats.confidence:
       [lon,lat,lev,mean,stdv,nob,chisqr,chisql,tstud] = get_gritas_data(fn,Res.globalProps.obType)
    else:
       [lon,lat,lev,mean,stdv,nob] = get_gritas_data(fn,Res.globalProps.obType)
    chnindx=lev;


    for nr,region in enumerate(Res.globalProps.regions):
      [reglat,reglon] = reg_def(region);
      [ilat,ilon] = subsetreg(reglat,reglon, lat,lon)


      for ns,stat in enumerate(Res.globalProps.stats.measures):
         if stat == 'sum' or stat == 'mean':
            this = vaccum(stat,Res.globalProps.obCnt,mean[:,ilat,:],nob[:,ilat,:])
            total[:,ii,nr,ns] = Res.globalProps.stats.scale*this
         if stat == 'stdv':
            this = vaccum(stat,Res.globalProps.obCnt,stdv[:,ilat,:],nob[:,ilat,:])
            total[:,ii,nr,ns] = Res.globalProps.stats.scale*this

         if Res.globalProps.stats.confidence:
	    [confl,confr,studt] = getconf(Res.globalProps.obCnt,nob[:,ilat,:],chisql[:,ilat,:],chisqr[:,ilat,:],tstud[:,ilat,:])

      # when only a single experiment is treated ...
      if numExp==0:
        if ns==0:
           if valuesYaml['global']['monthly plot']==True:
              yyyymm = str(yyyy)+mm
              show_bars_monthly(lev,total[:,ii,nr,0],['b'],[Res.globalProps.obType],Res.globalProps.stats.scale,yyyymm)
              pylab.savefig(Res.globalProps.nicknames[nx]+'_'+'monthly_'+Res.globalProps.obType+'_'+region+'_'+yyyymm+'.'+Res.globalProps.figType)
        else:
           if valuesYaml['global']['monthly plot']==True:
              yyyymm = str(yyyy)+mm
              show_bars_monthly(lev,total[:,ii,nr,:],Res.globalProps.colors,Res.globalProps.stats,Res.globalProps.stats.scale,yyyymm,chisql,chisqr,studt)
              pylab.savefig(Res.globalProps.nicknames[nx]+'_'+'monthly_'+Res.globalProps.obType+'_'+region+'_'+yyyymm+'.'+Res.globalProps.figType)
      else:  # one two experiments are treated (compares them) ...
        if nx==0:  # store result
           ii=index([Res.globalProps.stats.measures=='stdv'])
           totals[:,:,:,0] = total[:,:,:,ii]
        if nx==1:  # store result
           yyyymm = str(yyyy)+mm
           if valuesYaml['global']['compare monthly']['doit']==True:
              if   valuesYaml['global']['compare monthly']['type']=='ratio':
                 totals[:,:,:,0] = 100.0*(total[:,:,:,ii]/totals[:,:,:,0])
                 stat=Res.globalProps.stats.measures[0]
                 confl=confl*totals[:,ii,nr,0]
                 confr=confr*totals[:,ii,nr,0]
                 print confl
                 anno=Res.globalProps.expID
                 show_plot_monthly_one(lev,totals[:,ii,nr,0],anno,Res.globalProps.colors[0],'ratio',Res.globalProps.stats.scale,yyyymm,confl,confr,studt)
              elif valuesYaml['global']['compare monthly']['type']=='difference':
                 totals[:,:,:,0] = total[:,:,:,ii] - totals[:,:,:,0]
                 stat=Res.globalProps.stats.measures[0]
                 show_bars_monthly_one(lev,totals[:,ii,nr,0],Res.globalProps.colors[0],'difference',Res.globalProps.stats.scale,yyyymm,chisql,chisqr,studt)
              else:
                 totals[:,:,:,1] = total[:,:,:,ii]
                 show_bars_monthly(lev,totals[:,ii,nr,:],Res.globalProps.colors,Res.globalProps.expID,Res.globalProps.stats.scale,yyyymm,chisql,chisqr,studt)
              pylab.savefig(Res.globalProps.nicknames[1]+'X'+Res.globalProps.nicknames[0]+'_'+'monthly_'+Res.globalProps.obType+'_'+region+'_'+yyyymm+'.'+Res.globalProps.figType)


if numExp==0 and valuesYaml['global']['time series plot']==True:
   for nr,region in enumerate(Res.globalProps.regions):
      for ns,stat in enumerate(Res.globalProps.stats.measures):
         show_bars_time_series(Res.globalProps.obType,years,months,Res.globalProps.stats.scale,total[:,:,nr,ns])
         pylab.savefig(Res.globalProps.nicknames[0]+'_'+'tseries_'+Res.globalProps.obType+'_'+region+'_'+stat+'.'+ext)

