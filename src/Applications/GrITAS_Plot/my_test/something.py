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
pandasDates = pd.date_range(start=Res.globalProps.startDate, end=Res.globalProps.endDate, freq='MS') #'D') #'MS')
years = pandasDates.year; months = pandasDates.month


# Assume first and all files have data of same dimensions
#----------------------------------------------------------
tmpFile = Gritas(Res.globalProps.fileName[0])
nlev, nlat, nlon = tmpFile.fromGritas(Res.globalProps.obType,Res.globalProps.stats.confidence,'index').dims()

numMnths = size(years)   # number of times to process
numReg = size(Res.globalProps.regions) # number of regions
numStats = len(Res.globalProps.stats.measures)   # number of statitics to calculate
numExp = size(Res.globalProps.expID)  # number of experiments to handle


# Check each experiment has a file per month



# allStats = zeros([nlev,numMnths,numReg,numStats])
# compExp  = zeros([nlev,numMnths,numReg,2]) if numExp == 2 else None

allStats = zeros([nlev,numExp,numReg,numStats])
compExp  = zeros([nlev,numExp,numReg,2]) if numExp == 2 else None

# allStats = zeros([nlev,numMnths,numReg,numStats,numExp])
# compExp  = zeros([nlev,numStats,2]) if numExp == 2 else None

for ii,fname in enumerate(Res.globalProps.fileName):
  for nx,expid in enumerate(filter(None,Res.globalProps.expID)):

    yyyy = years[0] #[ii]
    mm = months[0] #[ii]
    mm =  str(mm).rjust(2, '0')


    datGritas=Gritas(fname)
    datGritas.fromGritas(Res.globalProps.obType,Res.globalProps.stats.confidence,'index')

    # Check for matching dimensions
    if datGritas.dims() != tmpFile.dims():
        raise ValueError("Incompatible dimensions!\n\t - %s\n\t - %s"%(datGritas,tmpFile))


    # Iterate over all regions of interest specified in global params
    for nr,region in enumerate(Res.globalProps.regions):

        # Lookup 'region' from universe
        reglat = Res.universe[region].latitude
        reglon = Res.universe[region].longitude

        # Subset of lat/lon defined by region...
        # Will need to adjust these, as they leave out boundary values of region
        latSubset = [n for n,l in enumerate(datGritas.loc['lat']) if reglat._min <= l and l <= reglat._max]
        lonSubset = [n for n,l in enumerate(datGritas.loc['lon']) if reglon._min <= l and l <= reglon._max]


        # Iterate over desired stats
        for ns,stat in enumerate(Res.globalProps.stats.measures):
            print(shape(allStats[:,ii,nr,ns]))
            print(shape(datGritas.getStat(stat,latSlice=latSubset,threshold=Res.globalProps.obCnt,mask='nobs')))
            allStats[:,ii,nr,ns]=(Res.globalProps.stats.scale*datGritas.getStat(stat,latSlice=latSubset,threshold=Res.globalProps.obCnt,mask='nobs'))

            # Grab statistics test scores
            if Res.globalProps.stats.confidence:
                datGritas.getConfidence(latSlice=latSubset,threshold=Res.globalProps.obCnt)

        print("SHAPE OF STATS: "+str(shape(allStats)))
        print(ii)
        print(nr)


        ##### Treat a single experiment
        # CE: ORIGINAL RT CODE HAD NUMEXP = 0 HERE; SHOULD THIS BE NUMEXP = 1 INSTEAD???
        if numExp == 1:
            datGritas.plotInit(Res.globalProps.nicknames[nx],Res.globalProps.obType,\
                               region,Res.globalProps.stats.scale,typ='monthly',yrs=yyyy,mnths=mm)
            if Res.globalProps.monthlyPlot:
                datGritas.monthlyBars(allStats[:,ii,nr,:],stats=Res.globalProps.stats,instruments=Res.globalProps.instruments)
                datGritas.saveFig()
                plt.show()
            else:
                print("Found numExp = 1, but monthlyPlot is false - skipping monthlyBars plot")

        elif numExp == 2:
            if nx==0:
                # Get index of standard deviations
                stdvIdx = np.argmax([s=='stdv' for s in Res.globalProps.stats.measures])
                compExp[:,:,:,0] = allStats[:,:,:,stdvIdx]
            if nx==1:

                datGritas.plotInit(Res.globalProps.nicknames,Res.globalProps.obType,\
                                   region,Res.globalProps.stats.scale,typ='monthly',yrs=yyyy,mnths=mm)
                print("UGH")
                if Res.globalProps.comparator.monthly:
                    if Res.globalProps.comparator.typ == 'ratio':
                        # print(allStats[:,:,:,stdvIdx])
                        # print(compExp[:,:,:,0])
                        print("stdvIdx = %i"%stdvIdx)
                        print(allStats[:,:,:,stdvIdx]==compExp[:,:,:,0])
                        compExp[:,:,:,0] = 100.0*(allStats[:,:,:,stdvIdx]/compExp[:,:,:,0])
                        print(allStats[:,:,:,stdvIdx]/compExp[:,:,:,0])
                        print(compExp[:,ii,nr,0])
                        print(shape(compExp[:,:,:,0]))
                        datGritas.confl*=compExp[:,ii,nr,0]
                        datGritas.confr*=compExp[:,ii,nr,0]
                        datGritas.monthlyPlot('ratio',compExp[:,ii,nr,0],stats=Res.globalProps.stats,instruments=Res.globalProps.instruments,annotation=Res.globalProps.expID)
                    elif Res.globalProps.comparator.typ == 'difference':
                        compExp[:,:,:,0] = allStats[:,:,:,stdvIdx] - compExp[:,:,:,0]
                        datGritas.monthlyBars(compExp[:,ii,nr,0],stats=Res.globalProps.stats,instruments=Res.globalProps.instruments,flavor='difference')
                    else:
                        compExp[:,:,:,1] = allStats[:,:,:,stdvIdx]
                        datGritas.monthlyBars(compExp[:,ii,nr,:],stats=Res.globalProps.stats,instruments=Res.globalProps.instruments,flavor=Res.globalProps.expID)
                    datGritas.saveFig()
                    plt.show()
        else:
            raise ValueError("Unable to handle more than 2 experiments at once!")



# print(allStats[:,:,0,0])
# print("***")
# print(allStats[:,0,0,0])
# print("***")
# print(allStats[0,0,0,0])

sys.exit(10)

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

