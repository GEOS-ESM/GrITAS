#!/usr/bin/env python

import helperFns
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
pandasDates = pd.date_range(start=Res.globalProps.startDate, end=Res.globalProps.endDate, freq='D') #'MS')
years = pandasDates.year; months = pandasDates.month


nlev, nlat, nlon = (0,0,0)

numT = size(years)   # number of times to process
numReg = size(Res.globalProps.regions) # number of regions
numStats = len(Res.globalProps.stats.measures)   # number of statitics to calculate
numExp = size(Res.globalProps.expID)  # number of experiments to handle


print("Zeros: %i %i %i %i"%(nlev,numT,numReg,numStats))


total = zeros([nlev,numT,numReg,numStats])
if numExp==2:
   totals = zeros([nlev,numT,numReg,2])



for ii,fname in enumerate(Res.globalProps.fileName):
  for nx,expid in enumerate(Res.globalProps.expID):

    yyyy = years[ii]
    mm = months[ii]
    mm =  str(mm).rjust(2, '0')


    foo=files(fname)
    foo.fromGritas(Res.globalProps.obType,Res.globalProps.stats.confidence,'index')# Res.globalProps.instruments.instrument.vertUnits)



    # Iterate over all regions of interest specified in global params
    for nr,region in enumerate(Res.globalProps.regions):

        # Lookup 'region' from universe
        reglat = Res.universe[region].latitude
        reglon = Res.universe[region].longitude      
        
        # Subset of lat/lon defined by region...
        # Will need to adjust these, as they leave out boundary values of region
        ilat = [n for n,l in enumerate(foo.loc['lat']) if reglat._min <= l and l <= reglat._max]
        ilon = [n for n,l in enumerate(foo.loc['lon']) if reglon._min <= l and l <= reglon._max]


        for ns,stat in enumerate(Res.globalProps.stats.measures):

            print(foo.nobs[0,ilat,0]) # no observations at lev[0], for all ilat latitudes and 0th longitude
            print(foo.nobs[1,ilat,0]) # no observations at lev[0], for all ilat latitudes and 0th longitude
            print(foo.nobs[2,ilat,0]) # no observations at lev[0], for all ilat latitudes and 0th longitude
            print(foo.nobs[3,ilat,0]) 
          

            X=foo.getStat(stat,foo.mean[:,ilat,:],thresh=Res.globalProps.obCnt,mask=foo.nobs[:,ilat,:])
            print(X)
            

            this=None
            if stat == 'sum' or stat == 'mean':
                this = vaccum(stat,Res.globalProps.obCnt,foo.mean[:,ilat,:],foo.nobs[:,ilat,:],Res.globalProps.stats.flavor)
            if stat == 'stdv':
                this = vaccum(stat,Res.globalProps.obCnt,foo.stdv[:,ilat,:],foo.nobs[:,ilat,:],Res.globalProps.stats.flavor)

            print(foo.mean[3,ilat,:])

            
            print(foo.mean[3,ilat,:].sum())
            for qq in ilat:
                s=''
                for zz in range(0,foo.getDim('lon')):
                    s+='%s '%str(foo.mean[3,qq,zz])
                print(s+"\n")
            sys.exit()
                
                    

            print(this)

            print(this.sum())
            print(foo.mean[:,ilat,:].sum())
            sys.exit()

            total[:,ii,nr,ns] = Res.globalProps.stats.scale*this


            if Res.globalProps.stats.confidence:
                [confl,confr,studt] = getconf(Res.globalProps.obCnt,nob[:,ilat,:],chisql[:,ilat,:],chisqr[:,ilat,:],tstud[:,ilat,:])

            continue
        sys.exit()

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

