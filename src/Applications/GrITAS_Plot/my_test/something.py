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
                  help='Visualize degrees of freedom for instrument')
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
years = pandasDates.year.to_list(); months = pandasDates.month.to_list()


# Assume first and all files have data of same dimensions
#----------------------------------------------------------
tmpFile = Gritas(Res.globalProps.experiments.members[0].pathToFile.replace('__ID__',Res.globalProps.experiments.members[0].name))
nlev, nlat, nlon = tmpFile.fromGritas(Res.globalProps.obType,Res.globalProps.stats.confidence,'index').dims()

numMnths = size(years)   # number of times to process
numReg = size(Res.globalProps.plotParams.regions) # number of regions
numStats = len(Res.globalProps.supportedStats)  # number of statistics supported
numExp = Res.globalProps.experiments.size()  # number of experiments to handle


print("Handling %i experiments"%numExp)
# Check each experiment has a file per month


allStats = zeros([nlev,numExp,numReg,numStats])
compExp  = zeros([nlev,numExp,numReg,2]) if numExp == 2 else None


# Files to read in
#----------------------
orderedGritKeys=['Control', 'Exp']
datGritas={'Control': None} #, 'Exp': None}

for nx, exp in enumerate(Res.globalProps.experiments.members):
    grit=Gritas(exp.pathToFile.replace('__ID__',exp.name),supportedStats=Res.globalProps.supportedStats)
    grit.fromGritas(Res.globalProps.obType,Res.globalProps.stats.confidence,'index')
    print(grit)
    # Check for matching dimensions
    if grit.dims() != tmpFile.dims():
        raise ValueError("Incompatible dimensions!\n\t - %s\n\t - %s"%(grit,tmpFile))
    # Package
    if nx == 0:
        datGritas['Control']=grit
    elif nx == 1:
        datGritas.update({'Exp': grit})
    else:
        raise ValueError("Unable to handle more than 2 experiments at once!")
#---------------------------------------------------------------------------------------



# Iterate over regions of interest
#---------------------------------------------------------------------------------------
for nr,region in enumerate(Res.globalProps.plotParams.regions):
    # Lookup 'region' from universe
    reglat = Res.universe[region].latitude
    reglon = Res.universe[region].longitude

    # Subset of lat/lon defined by region...
    # Will need to adjust these, as they leave out boundary values of region
    latSubset = [n for n,l in enumerate(datGritas['Control'].loc['lat']) if reglat._min <= l and l <= reglat._max]
    lonSubset = [n for n,l in enumerate(datGritas['Control'].loc['lon']) if reglon._min <= l and l <= reglon._max]

    # Iterate over all supported stats - allows user to toggle which stat(s) to consider
    #---------------------------------------------------------------------------------------
    for ns, stat in enumerate(Res.globalProps.supportedStats):

        # Iterate over control and, potentially, experiment
        # --------------------------------------------------
        for nx, exp in enumerate(orderedGritKeys):
            try:
                allStats[:,nx,nr,ns]=datGritas[exp].getStat(stat,latSlice=latSubset,threshold=Res.globalProps.obCnt,\
                                                            rescale=Res.globalProps.stats.scale,mask='nobs')
                # Optionally grab test statistic scores
                if Res.globalProps.stats.confidence:
                    datGritas[exp].getConfidence(latSlice=latSubset,threshold=Res.globalProps.obCnt)
            except:
                pass



    # Proceed to plot with a single experiment
    # -----------------------------------------
    if numExp == 1:
        if Res.globalProps.plotParams.monthly:
            datGritas['Control'].plotInit(Res.globalProps.experiments.members[0].nickname,Res.globalProps.obType,\
                                          Res.universe[region],Res.globalProps.stats.scale,\
                                          Res.globalProps.plotParams.simpleBars,typ='monthly',\
                                          yrs=years[0],mnths=str(months[0]).rjust(2,'0'))
            datGritas['Control'].monthlyStat(allStats[:,0,nr,:],stats=Res.globalProps.stats,\
                                             instruments=Res.globalProps.instruments)
        elif Res.globalProps.plotParams.timeSeries:
            datGritas['Control'].plotInit(Res.globalProps.experiments.members[0].nickname,Res.globalProps.obType,\
                                          Res.universe[region],Res.globalProps.stats.scale,\
                                          Res.globalProps.plotParams.simpleBars,typ='tseries',yrs=years,mnths=months)
            datGritas['Control'].tSeries(allStats[:,0,nr,:],stats=Res.globalProps.stats,\
                                         instruments=Res.globalProps.instruments)
        else:
            print("Found numExp = 1, but Res.globalProps.plotParams.{monthlyPlot,timeSeries} are False!")

        # Save the figure on hand, if it exists
        if datGritas['Control'].figExist:
            datGritas['Control'].saveFig()
    # -------------------------------------------------------------------------------------

    # Proceed to plot with two experiments
    # -------------------------------------
    if numExp == 2:

        # Grab indicies of mean and stdv - in case user changes order in yaml
        meanIdx = np.argmax([s=='mean' for s in Res.globalProps.stats.measures])
        stdvIdx = np.argmax([s=='stdv' for s in Res.globalProps.stats.measures])

        cntlMean = allStats[:,0,nr,meanIdx]
        cntlStdv = allStats[:,0,nr,stdvIdx]

        # Convenience
        nicknames = [exp.nickname for exp in Res.globalProps.experiments.members]


        print(years)
        print(months)
        print(type(years))
        print(type(months))
        if Res.globalProps.plotParams.monthly:
            datGritas['Exp'].plotInit(nicknames,Res.globalProps.obType,\
                                      Res.universe[region],Res.globalProps.stats.scale,\
                                      Res.globalProps.plotParams.simpleBars,typ='monthly',\
                                      yrs=years,mnths=months)
                                      # yrs=years[0],mnths=str(months[0]).rjust(2,'0'))
            datGritas['Exp'].monthlyComp(Res.globalProps.plotParams.compareVia,
                                         allStats[:,1,nr,stdvIdx],allStats[:,0,nr,stdvIdx],
                                         stats=Res.globalProps.stats,
                                         instruments=Res.globalProps.instruments,annotation=nicknames)
        else:
            print("Found numExp = 2, but plotParams.monthly is False - skipping monthlyComp plot")

        # Save the figure on hand, if it exists
        if datGritas['Exp'].figExist:
            datGritas['Exp'].saveFig()
    # -------------------------------------------------------------------------------------

plt.show()
