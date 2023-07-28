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


# Assume first and all files have data of same dimensions
#----------------------------------------------------------
tmpFile = Gritas(Res.globalProps.fileName[0])
nlev, nlat, nlon = tmpFile.fromGritas(Res.globalProps.obType,Res.globalProps.stats.confidence,'index').dims()

numT = size(years)   # number of times to process
numReg = size(Res.globalProps.regions) # number of regions
numStats = len(Res.globalProps.stats.measures)   # number of statitics to calculate
numExp = size(Res.globalProps.expID)  # number of experiments to handle


print("Zeros: %i %i %i %i"%(nlev,numT,numReg,numStats))


allStats = zeros([nlev,numT,numReg,numStats])
if numExp==2:
    allStats = zeros([nlev,numT,numReg,2])
# total = zeros([nlev,numT,numReg,numStats])
# if numExp==2:
#    totals = zeros([nlev,numT,numReg,2])


for ii,fname in enumerate(Res.globalProps.fileName):
  for nx,expid in enumerate(filter(None,Res.globalProps.expID)):

    yyyy = years[ii]
    mm = months[ii]
    mm =  str(mm).rjust(2, '0')


    foo=Gritas(fname)
    foo.fromGritas(Res.globalProps.obType,Res.globalProps.stats.confidence,'index')# Res.globalProps.instruments.instrument.vertUnits)

    # Check for matching dimensions
    # exceptions.dims(foo,tmpFile)
    if foo.dims() != tmpFile.dims():
        raise ValueError("Input file [ %s ] of dimensions (%i,%i,%i) does not match dimensions of template file (%i,%i,%i)"%(fname,foo.dims(),tmpFile.dims()))

    # Iterate over all regions of interest specified in global params
    for nr,region in enumerate(Res.globalProps.regions):

        # Lookup 'region' from universe
        reglat = Res.universe[region].latitude
        reglon = Res.universe[region].longitude      
        
        # Subset of lat/lon defined by region...
        # Will need to adjust these, as they leave out boundary values of region
        ilat = [n for n,l in enumerate(foo.loc['lat']) if reglat._min <= l and l <= reglat._max]
        ilon = [n for n,l in enumerate(foo.loc['lon']) if reglon._min <= l and l <= reglon._max]


        # Iterate over desired stats
        for ns,stat in enumerate(Res.globalProps.stats.measures):

            allStats[:,ii,nr,ns]=(Res.globalProps.stats.scale*\
                                  foo.getStat(stat,latSlice=ilat,\
                                              threshold=Res.globalProps.obCnt,mask='nobs'))
            
            # Grab statistics test scores
            if Res.globalProps.stats.confidence:
                foo.getConf(latSlice=ilat,threshold=Res.globalProps.obCnt)

        print(shape(allStats))
        print(ii)
        print(nr)
        #########
        # GET TO WORK W/ FIGURES
        #########
        yyyymm = str(yyyy)+mm

        foo.plotInit(Res.globalProps.nicknames[nx],Res.globalProps.obType,region,Res.globalProps.stats.scale,typ='monthly',yrs=yyyy,mnths=mm)
        foo.monthlyBars(allStats[:,ii,nr,:],stats=Res.globalProps.stats,instruments=Res.globalProps.instruments)

        # Treat a single experiment
        if numExp == 0 and Res.globalProps.monthlyPlot:
            foo.monthlyBars(allStats[:,ii,nr,:],stats=Res.globalProps.stats,instruments=Res.globalProps.instruments)
        else:
            pass


        foo.saveFig()
        plt.show()
        sys.exit()



#         else:  # one two experiments are treated (compares them) ...
#             if nx==0:  # store result
#                 ii=index([Res.globalProps.stats.measures=='stdv'])
#                 totals[:,:,:,0] = total[:,:,:,ii]
#             if nx==1:  # store result
#                 if valuesYaml['global']['compare monthly']['doit']==True:
#                     if   valuesYaml['global']['compare monthly']['type']=='ratio':
#                         totals[:,:,:,0] = 100.0*(total[:,:,:,ii]/totals[:,:,:,0])
#                         stat=Res.globalProps.stats.measures[0]
#                         foo.confl=foo.confl*totals[:,ii,nr,0]
#                         foo.confr=foo.confr*totals[:,ii,nr,0]
#                         print foo.confl
#                         anno=Res.globalProps.expID
#                         show_plot_monthly_one(foo.lev,totals[:,ii,nr,0],anno,Res.globalProps.colors[0],'ratio',Res.globalProps.stats.scale,yyyymm,foo.confl,foo.confr,foo.studt)
#                     elif valuesYaml['global']['compare monthly']['type']=='difference':
#                         totals[:,:,:,0] = total[:,:,:,ii] - totals[:,:,:,0]
#                         stat=Res.globalProps.stats.measures[0]
#                         show_bars_monthly_one(foo.lev,totals[:,ii,nr,0],Res.globalProps.colors[0],'difference',Res.globalProps.stats.scale,yyyymm,chisql,chisqr,foo.studt)
#                     else:
#                         totals[:,:,:,1] = total[:,:,:,ii]
#                         show_bars_monthly(foo.lev,totals[:,ii,nr,:],Res.globalProps.colors,Res.globalProps.expID,Res.globalProps.stats.scale,yyyymm,chisql,chisqr,foo.studt)
#                     pylab.savefig(Res.globalProps.nicknames[1]+'X'+Res.globalProps.nicknames[0]+'_'+'monthly_'+Res.globalProps.obType+'_'+region+'_'+yyyymm+'.'+Res.globalProps.figType)



# print(allStats[:,:,0,0])
# print("***")
# print(allStats[:,0,0,0])
# print("***")
# print(allStats[0,0,0,0])


# Case where we aren't looking at a specific experiment and, rather, desire a time series plot
if numExp==0 and Res.globalProps.tSeriesPlot:
   for nr,region in enumerate(Res.globalProps.regions):
      for ns,stat in enumerate(Res.globalProps.stats.measures):
          tmpFile.getStat(stat,threshold=Res.globalProps.obCnt,mask='nobs')
          print(allStats[:,:,nr,ns])
          print("XXXXXXXXX")
          tmpFile.plotInit(Res.globalProps.nicknames[0],Res.globalProps.obType,region,Res.globalProps.stats.scale,typ='tseries',yrs=years,mnths=months)
          tmpFile.tSeries(allStats[:,:,nr,ns],stats=Res.globalProps.stats,instruments=Res.globalProps.instruments)
          tmpFile.saveFig()

