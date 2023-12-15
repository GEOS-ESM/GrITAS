#!/usr/bin/env python3

from figutil import *
import plot_util
import common
import pandas as pd
import yaml

import optparse

parser = optparse.OptionParser()
parser.add_option("-d","--dfs",action="store_true",default=False,
                  help='Visualize degrees of freedom for instrument')
parser.add_option("-i","--impacts",action="store_true",default=False,
                  help='Visualize observation impacts')
parser.add_option("-r","--residuals",action="store_true",default=False,
                  help='Visualize residuals')
parser.add_option("-f","--ymlConfig",default='',
                  help='Yaml driver (default = "")')
# Parse the arguments
(options, args) = parser.parse_args()


# Only want to do one visualization at a time
# ---------------------------------------------------------------------------------------
if not options.impacts^options.dfs^options.residuals:
    raise ValueError("Scripts accepts only single visualization at a time!")

# Init and read classes from yaml
# ---------------------------------------------------------------------------------------
Global=plot_util.GlobalProps()
Regions=plot_util.Collection('regions')
with open(options.ymlConfig, 'r') as f:
    valsYML = yaml.load(f, Loader=yaml.FullLoader)
    Global.fromYaml(valsYML['global'])
    Regions.fromYaml(valsYML['regions'],cls='REGION')

# Initialize a StatsViewer instance
# ---------------------------------------------------------------------------------------
SV=plot_util.StatsViewer(glob=Global,universe=Regions)

# Set temporal information
# ---------------------------------------------------------------------------------------
pandasDates = pd.date_range(start=SV.globalProps.startDate, end=SV.globalProps.endDate, freq='MS')
years = pandasDates.year.to_list(); months = pandasDates.month.to_list()

# Convenience
# ---------------------------------------------------------------------------------------
numMnths = size(years)                           # number of times to process
numReg = size(SV.globalProps.plotParams.regions) # number of regions
numStats = len(SV.globalProps.supportedStats)    # number of statistics supported
numExp = SV.globalProps.experiments.size()       # number of experiments to handle

# Numpy arrays to conveniently bundle read statistics
# ---------------------------------------------------------------------------------------
allStats = zeros([])
compExp  = zeros([])

# Read each experiment
# ---------------------------------------------------------------------------------------
gritDict={'Control': None}; orderedGritKeys=['Control', 'Exp']
for nx, exp in enumerate(SV.globalProps.experiments.members):
    # Use current experiment to initialize a Gritas instance
    grit=Gritas(exp.pathToFile.replace('__ID__',exp.name),supportedStats=SV.globalProps.supportedStats)
    grit.fromGritas(SV.globalProps.obType,SV.globalProps.stats.confidence,'index')

    # Use first experiment to set 'allStats' and 'compExp' array dimensions
    if nx == 0:
        allStats.resize([grit.getDim('lev'),numExp,numReg,numStats])
        compExp.resize([grit.getDim('lev'),numExp,numReg,2]) if numExp == 2 else None

    # If a second file is read, assert its dimensions must match those of first file read
    if nx == 1 and grit.dims() != gritDict['Control'].dims():
        raise ValueError("Incompatible dimensions!\n\t - %s\n\t - %s"%(grit,gritDict['Control']))

    # Package
    try:
        gritDict.update({orderedGritKeys[nx]: grit})
    except:
        raise ValueError("Unable to handle more than 2 experiments at once!")

# Iterate over regions of interest
# ---------------------------------------------------------------------------------------
for nr,region in enumerate(SV.globalProps.plotParams.regions):
    # Lookup 'region' from universe
    reglat = SV.universe[region].latitude
    reglon = SV.universe[region].longitude

    # Subset of lat/lon defined by region...
    # Will need to adjust these, as they leave out boundary values of region
    latSubset = [n for n,l in enumerate(gritDict['Control'].loc['lat']) if reglat._min <= l and l <= reglat._max]
    lonSubset = [n for n,l in enumerate(gritDict['Control'].loc['lon']) if reglon._min <= l and l <= reglon._max]

    # Iterate over all supported stats - allows user to toggle which stat(s) to consider
    # ---------------------------------------------------------------------------------------
    for ns, stat in enumerate(SV.globalProps.supportedStats):
        # Iterate over 'Control' and, potentially, 'Experiment'
        # ---------------------------------------------------------------------------------------
        for nx, exp in enumerate(orderedGritKeys):
            try:
                allStats[:,nx,nr,ns]=gritDict[exp].getStat(stat,latSlice=latSubset,threshold=SV.globalProps.obCnt,\
                                                           rescale=SV.globalProps.stats.scale,mask='nobs')
                # Optionally grab test statistic scores
                if SV.globalProps.stats.confidence:
                    gritDict[exp].getConfidence(latSlice=latSubset,threshold=SV.globalProps.obCnt)
            except:
                pass


    # Proceed to plot with a single experiment
    # ---------------------------------------------------------------------------------------
    if numExp == 1:
        gritDict['Control'].plotInit(SV.globalProps.experiments.members[0].nickname,SV.globalProps.obType,\
                                     SV.universe[region],SV.globalProps.stats.scale,\
                                     SV.globalProps.plotParams.simpleBars,yrs=years,mnths=months)
        if SV.globalProps.plotParams.monthly:
            gritDict['Control'].monthlyStat(allStats[:,0,nr,:],stats=SV.globalProps.stats,\
                                            instruments=SV.globalProps.instruments)
        elif SV.globalProps.plotParams.timeSeries:
            gritDict['Control'].tSeries(allStats[:,0,nr,:],stats=SV.globalProps.stats,\
                                        instruments=SV.globalProps.instruments,\
                                        flavor=SV.globalProps.plotParams.timeSeriesVar)
        else:
            print("Found numExp = 1, but SV.globalProps.plotParams.{monthlyPlot,timeSeries} are False!")

        # Save the figure on hand, if it exists
        gritDict['Control'].saveFig()


    # Proceed to plot with two experiments
    # ---------------------------------------------------------------------------------------
    if numExp == 2:
        # Convenience
        nicknames = [exp.nickname for exp in SV.globalProps.experiments.members]

        if SV.globalProps.plotParams.monthly:
            gritDict['Exp'].plotInit(nicknames,SV.globalProps.obType,\
                                     SV.universe[region],SV.globalProps.stats.scale,\
                                     SV.globalProps.plotParams.simpleBars,\
                                     yrs=years,mnths=months)
            gritDict['Exp'].monthlyComp(SV.globalProps.plotParams.compareVia,
                                        allStats[:,1,nr],allStats[:,0,nr],
                                        stats=SV.globalProps.stats,
                                        instruments=SV.globalProps.instruments,annotation=nicknames)
        else:
            print("Found numExp = 2, but plotParams.monthly is False - skipping monthlyComp plot")

        # Save the figure on hand, if it exists
        gritDict['Exp'].saveFig()


plt.show()
