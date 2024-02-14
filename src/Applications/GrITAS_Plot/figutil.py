#!/usr/bin/env python3

import sys
import warnings
import netCDF4 as nc4
import matplotlib as mpl
import matplotlib.pyplot as plt
from pylab import *
from datetime import timedelta
from numpy import nan_to_num as nn

# Global settings for matplotlib
# ------------------------------
mpl.rc('font', family='serif')

def revMaskedArray(arr,action):
    '''
    Reverse an array

    Parameters
    ----------
    arr : numpy.ndarray
       Array of data

    action : bool
       Whether reversal should be performed

    Returns
    -------
    arr, class range
    '''
    if action:
        idx = range(len(arr)-1,-1,-1)
        return arr[idx], idx
    else:
        return arr, range(1,len(arr),1)


def lvlAvg(arr,mask):
    '''
    Compute average of an array, neglecting elements based on a mask

    Parameters
    ----------
    arr : numpy.ndarray

    mask : numpy.ndarray

    Returns
    -------
    float
    '''
    # Access first element of shape(mask) for number lat,lon points for which elements of arr meet mask criterion
    return arr[mask].sum()/shape(mask)[1]

class GritasFig:
    '''
    Manages construction of figure(s) to visualize gridded time-averaged statistics

    Attributes
    ----------
    prefix : str
       Prefix of figure to be saved

    obType : str
       Instrument belonging to experiment(s)

    region : plot_util.Region instance
       A geographic region figure will correspond to

    scale : float
       Rescale statistics

    linePlot : bool
           Plot stats via a line plot; bars used otherwise

    simpleBars : bool
       Embelished confidence intervals

    includeObsCounts : bool
       Whether a side panel should be included to show observation counts per level/channel

    typ : str
       Shorthand for type of figure that will be produced (e.g., 'monthly' or 'tseries')

    figType : str
       Figure format

    yrs : list
       List of years in time window

    mnths : list
       List of months in time window

    yyyymm : str
       Time window as a string, formatted as YYYYMM-YYYYMM

    overlaidAxis : None
       Member variable to hold reference to a twiny() instance, should it exist
    '''
    def __init__(self,prefix,obType,region,scale,linePlot,simpleBars,includeObsCounts,yrs=None,mnths=None):
        '''
        Intialize a GritasFig instance

        Parameters
        ----------
        prefix : str
           Prefix of figure to be saved

        obType : str
           Instrument belonging to experiment(s)

        region : plot_util.Region instance
           A geographic region figure will correspond to

        scale : float
           Rescale statistics

        linePlot : bool
           Plot stats via a line plot; bars used otherwise

        simpleBars : bool
           Embelished confidence intervals

        includeObsCounts : bool
           Whether a side panel should be included to show observation counts per level/channel

        yrs : list
           List of years in time window

        mnths : list
           List of months in time window
        '''
        self.prefix=''.join(prefix) if isinstance(prefix, str) else 'X'.join(prefix)
        self.obType=obType
        self.region=region
        self.scale=scale
        self.linePlot=linePlot
        self.simpleBars=simpleBars
        self.includeObsCounts=includeObsCounts
        self.typ='invalid'
        self.figType='png'
        self.yrs=yrs
        self.mnths=mnths
        self.yyyymm='%s%s'%(self.yrs[0],str(self.mnths[0]).rjust(2,'0'))
        if len(self.yrs) > 1 or len(self.mnths) > 1:
            self.yyyymm += '-%s%s'%(self.yrs[-1],str(self.mnths[-1]).rjust(2,'0'))
        self.overlaidAxis=None

    def monthlyStat(self,allStats,stats=None,instruments=[],annotation=None):
        '''
        Produce a monthly plot of statistic

        Parameters
        ----------
        allStats : numpy.ndarray
           Statistics per level, experiment, region, and stat

        stats : plot_util.Stats instance
           Configure statistics

        instruments : plot_util.Collection instance
           Configure instruments

        annotation : list
           List of strings of annotations (used to label 'CTL' and 'EXP' on figure)

        Returns
        -------
        None
        '''
        self.typ='monthly'
        self.commonFigSetup(allStats,stats.units,instruments)

        # Iterate over all supported stats, skipping if not desired or stat is 'sum' ('sum' handled separately)
        for ns, stat in enumerate(self.supportedStats):
            if stat not in stats.measures or stat == 'sum':
                continue
            statAllLvls=allStats[:,ns]

            # Set center of bars
            barCenter = arange(len(statAllLvls))-1.0/5+ns*self.bar_width

            # Set properties of errorbars
            kwargs={'alpha': self.opacity,'color': stats.colors[ns],'label': stat}#,'error_kw': self.error_config}

            # Confidence on monthlyStat
            # --------------------------
            if stats.confidence:
                # Set xerr to studT score if forming CL for pop. mean, and left/right Chi2 scores if forming
                # CL for pop. stdv.
                xerr=self.studT if stat == 'mean' else [self.leftChi2,self.rightChi2]

                # Assign no error if a stat other than pop. mean/stdv is being estimated
                if stat != 'mean' and stat != 'stdv':
                    xerr=np.ones(np.size(allStats[:,1]))

                # Multiply by sample stdv., thus completing formation of CL
                xerr*=allStats[:,1]

            # Set the left/right error based on stat considered
            if stat == 'mean':
                left = statAllLvls-xerr; right = statAllLvls+xerr
            elif stat == 'stdv':
                left = statAllLvls-xerr[0]; right = statAllLvls+xerr[1]

            # Toggle between plotting stats with bars or line plots
            if self.linePlot:
                self.ax.plot(statAllLvls, barCenter, lw=1, color=kwargs['color'])
                self.ax.fill_betweenx(barCenter, left, right, **kwargs)
            else:
            # Plot bars
                self.ax.barh(barCenter, statAllLvls, self.bar_width, **kwargs)

                # Toggle between simple or embelished bars
                if self.simpleBars:
                    self.ax.errorbar(statAllLvls, barCenter, xerr=xerr, color=stats.colors[ns], capsize=5,ls='')
                else:
                    if stat == 'mean':
                        # Presence of NaNs leads to omissions when barh is passed entire stat/err arrays - iterate instead
                        for _n in range(len(barCenter)):
                            self.ax.barh(barCenter[_n], 2*xerr[_n], self.bar_width, left=statAllLvls[_n]-xerr[_n],\
                                         color=stats.colors[ns],alpha=0.7,hatch='/////',edgecolor=stats.colors[ns])
                    else:
                        # Presence of NaNs leads to omissions when barh is passed entire stat/err arrays - iterate instead
                        for _n in range(len(barCenter)):
                            self.ax.barh(barCenter[_n], xerr[0,_n], self.bar_width, left=left[_n],\
                                         color=stats.colors[ns],alpha=0.7,hatch='/////',edgecolor=stats.colors[ns])
                            self.ax.barh(barCenter[_n], xerr[1,_n], self.bar_width, left=statAllLvls[_n],\
                                         color=stats.colors[ns],alpha=0.7,hatch='/////',edgecolor=stats.colors[ns])

        # Add annotations to figure if not None
        # -------------------------------------
        if annotation:
            self.ax.set_title('CTL: %s\nEXP: %s\n'%tuple(annotation), loc='left', fontsize=self.figTitleSize)

        # Include observation counts if self.includeObsCounts is True
        # -----------------------------------------------------------
        if self.includeObsCounts:
            self.ax_counts.plot(allStats[:,2],arange(len(allStats)),color='g')

        self.ax.legend(loc='lower right')

        # Reset yticks
        self.ax.set_yticks(range(0,self.getDim('lev')))

    def tSeries(self,allStats,stats=None,instruments=[],flavor='None'):
        '''
        Produce a time series plot of statistic

        Parameters
        ----------
        allStats : numpy.ndarray
           Statistics per level, experiment, region, and stat

        stats : plot_util.Stats instance
           Configure statistics

        instruments : plot_util.Collection instance
           Configure instruments

        flavor : str
           Statistic to form time series of (e.g., mean, stdv, or sum)

        Returns
        -------
        None

        Raises
        ------
        ValueError : Variable rval, which is ratio of maxval/minval or its reciprocal, is postive; minval,maxval define the domain of plots produced.
        '''
        self.typ='tseries'
        # Capture minval/maxval
        minval,maxval = self.commonFigSetup(allStats,stats.units,instruments,flavor=flavor)

        # Make copies of allStats for masking
        slices=[allStats.copy(), allStats.copy(), allStats.copy()]
        # Mask
        for n,s in enumerate(slices):
            s[:,(n+1)%len(slices)]=s[:,(n+2)%len(slices)]=np.nan

        cs=None
        cbarLabel=''

        rval = maxval/minval if minval < 0.0 else minval/maxval
        if rval < 0.0:
            maxval= max(abs(minval),abs(maxval))
            minval=-maxval

            if flavor == 'mean':
                cs=self.ax.pcolor(slices[0],cmap=plt.cm.Blues,vmin=np.nanmin(slices[0]),\
                                  vmax=np.nanmax(slices[0]))
                cbarLabel=flavor.capitalize()
            elif flavor == 'stdv':
                cs=self.ax.pcolor(slices[1],cmap=plt.cm.autumn_r,vmin=np.nanmin(slices[1]),\
                                  vmax=np.nanmax(slices[1]))
                cbarLabel=flavor.capitalize()
            elif flavor == 'sum':
                cs=self.ax.pcolor(slices[2],cmap=plt.cm.Greens,vmin=0,vmax=np.nanmax(slices[2]))
                cbarLabel='# Obs.'
        else:
            raise ValueError("rval = maxval/minval is positive! This functionality is not yet established!")
        # Add colorbars to plot
        cbar = plt.colorbar(cs); cbar.ax.set_xlabel(cbarLabel)

        self.ax.axes.xaxis.set_visible(True)
        self.ax.axes.yaxis.set_visible(True)
        self.ax.set_xlabel('Date')
        self.ax.set_xticks(range(len(slices)*len(self.yrs)+1))

        yyyymm = ['/'.join(["{:02d}".format(m),str(y)]) for m,y in zip(self.mnths,self.yrs)]
        yyyymm.insert(0,'') # lazy - major_locator below is causing first entry to be dropped

        self.ax.axes.xaxis.set_major_locator(MultipleLocator(3))
        self.ax.axes.xaxis.set_minor_locator(MultipleLocator(1))
        self.ax.set_xticklabels(yyyymm, minor=False, rotation=45, fontsize=12)
        self.ax.vlines(self.ax.get_xticks(),self.ax.get_ylim()[0],self.ax.get_ylim()[1],color='darkgray',linestyle=':')

    def _baseComp_(self,compareVia,expStats,cntlStats,stats,instruments,annotation):
        allStats=100*(expStats/cntlStats) if compareVia == 'ratio' else expStats - cntlStats
        # Capture minval/maxval
        minval, maxval = self.commonFigSetup(allStats,stats.units,instruments,flavor=annotation)

        # Convenience
        midpnt = minval+0.5*(maxval-minval)
        pos = arange(len(allStats))

        # For the moment, we will only compare difference in means (blue) or ratio of stdv's (red)
        ecolor = 'b' if compareVia == 'difference' else 'r'

        # Include confidence levels
        # --------------------------
        if stats.confidence:
            # xerr based on compareVia
            xerr = self.studT if compareVia == 'difference' else [self.leftChi2,self.rightChi2]
            # Multiply by sample stdv. of Exp, thus completing formation of CL
            xerr*=expStats # this is okay, since this method is only called for ratio

            if compareVia == 'difference':
                self.ax.plot(zeros(len(allStats)), pos, color='k', lw=1, alpha=0.8)
                self.ax.fill_betweenx(pos, allStats-xerr, allStats+xerr, color=ecolor, alpha=0.4)
            if compareVia == 'ratio':
                refline=midpnt*ones(len(allStats))
                self.ax.plot(refline, pos, color='k', lw=1, alpha=0.8)
                self.ax.fill_betweenx(pos, allStats-xerr[0], allStats+xerr[-1], color=ecolor, alpha=0.4)

            # Plot errorbars regardless of compareVia
            self.ax.plot(allStats, pos, lw=2, color=ecolor)
        else:
            self.ax.plot(allStats, pos, lw=2, color=ecolor)

        return minval, maxval

    def _overlayComp_(self,ax,compareVia,expStats,cntlStats,stats,instruments,annotation=''):
        if compareVia != 'difference':
            raise ValueError("Expecting compareVia = 'difference' in overlaid plot!")
        overlayStats=expStats-cntlStats

        ax2_min, ax2_max = np.nanmin(overlayStats), np.nanmax(overlayStats)

        ax.set_xlim([-1.05*max(np.abs([ax2_min,ax2_max])),1.05*max(np.abs([ax2_min,ax2_max]))])
        ax.set_xlabel(r'$<x_{EXP}>-<x_{CTL}>$',fontsize=18)

        for label in ax.get_xticklabels():
            label.set_fontsize(self.maxLabelSize)
            label.set_fontweight('bold')

        if stats.confidence:
            xerr = self.studT
            xerr*=expStats
            ax.plot(overlayStats,arange(len(overlayStats)),'b')
            ax.fill_betweenx(arange(len(overlayStats)), overlayStats-xerr, overlayStats+xerr, color='b', alpha=0.4)
        else:
            ax.plot(overlayStats,arange(len(overlayStats)),'b')

        self.overlaidAxis = ax

    def monthlyComp(self,compareVia,expStats,cntlStats,stats=None,instruments=[],annotation=''):
        '''
        Form a monthly comparison between two experiments according to specified scheme

        Parameters
        ----------
        compareVia : str
           When two experiments are considered, compare them according to scheme (supported: 'ratio' and 'difference')

        expStats : numpy.ndarray
           'Experiment (Exp)' statistics across all vertical levels within a single region. Chosen statistic determined at calling location via 'compIdx'.

        cntlStats : numpy.ndarray
           'Control (Cntl)' statistics across all vertical levels within a single region. Chosen statistic determined at calling location via 'compIdx'.

        stats : plot_util.Stats instance
           Configure statistics

        instruments : plot_util.Collection instance
           Configure instruments

        annotation : list
           List (string) of experiment nicknames

        Returns
        -------
        None
        '''
        self.typ='monthly'

        # Grab indicies of mean, stdv, sum - in case user changes order in yaml
        meanIdx = np.argmax([s=='mean' for s in stats.measures])
        stdvIdx = np.argmax([s=='stdv' for s in stats.measures])
        cntsIdx = np.argmax([s=='sum' for s in stats.measures])

        minval, maxval = 0.0, 0.0
        if compareVia == 'ratio+difference' or compareVia == 'difference+ratio':
            minval, maxval = self._baseComp_('ratio',expStats[:,stdvIdx],
                                             cntlStats[:,stdvIdx],stats,instruments,annotation)
            self._overlayComp_(self.ax.twiny(),'difference',expStats[:,meanIdx],
                               cntlStats[:,meanIdx],stats=stats,instruments=instruments,annotation=annotation)
        elif compareVia == 'ratio':
            minval, maxval = self._baseComp_(compareVia,expStats[:,stdvIdx],
                                             cntlStats[:,stdvIdx],stats,instruments,annotation)
        elif compareVia == 'difference':
            minval, maxval = self._baseComp_(compareVia,expStats[:,meanIdx],
                                             cntlStats[:,meanIdx],stats,instruments,annotation)
        else:
            raise Exception("Error! Should not get here!")

        # Plot deterioration/improvement boxes #midpnt+0.0195*midpnt
        self.ax.annotate('Deterioration', xy=(maxval-0.0195*(maxval-minval),0.5), xycoords='data',
                         xytext=(0,0), textcoords='offset points',
                         size=13, ha='right', va="center",
                         bbox=dict(boxstyle="round", alpha=0.1, color='r'))
        self.ax.annotate('Improvement', xy=(minval+0.0195*(maxval-minval),0.5), xycoords='data',
                         xytext=(0,0), textcoords='offset points',
                         size=13, ha='left', va="center",
                         bbox=dict(boxstyle="round", alpha=0.1, color='g'))

        self.ax.set_title('CTL: %s\nEXP: %s\n'%tuple(annotation), loc='left', fontsize=self.figTitleSize, pad=4)

        # Match self.ax_counts twiny x-axis labels to that of self.overlaidAxis x-axis labels
        if self.includeObsCounts:
            _ax_counts_twin_=self.ax_counts.twiny()
            _ax_counts_twin_.set_xlabel(self.overlaidAxis.get_xlabel(),fontsize=18,fontweight='bold')
            _ax_counts_twin_.set_xticklabels(np.arange(3),fontsize=self.maxLabelSize,fontweight='bold')
            _ax_counts_twin_.set_axis_off()

            ratioCnts = 100*expStats[:,cntsIdx]/cntlStats[:,cntsIdx]
            self.ax_counts.plot(ratioCnts,arange(len(ratioCnts)),color='g')
            self.ax_counts.plot(100*ones(len(ratioCnts)), arange(len(ratioCnts)), color='k', lw=1, alpha=0.8)
            self.ax_counts.xaxis.set_major_locator(MaxNLocator(3))
            self.ax_counts.set_xticks([100-0.001*np.nanmax(ratioCnts),100,100+0.001*np.nanmax(ratioCnts)])
            self.ax_counts.set_xlabel(r'$\%(N^{Obs}_{EXP}/N^{Obs}_{CTL})$',fontsize=18)

            # Reset yticks
            self.ax_counts.set_yticks(range(0,self.getDim('lev')))

    def commonFigSetup(self,allStats,units,instruments,flavor='None'):
        '''
        Common setup for figures

        Parameters
        ----------
        allStats : numpy.ndarray
           Statistic(s) per level, experiment(s), region(s); handles cases where multiple statistics, regions, and experiments are involved.

        units : str
           Units of statistics (used to label horizontal axes of plots)

        instruments : plot_util.Collection instance
           Configure instruments

        flavor : str
           Statistic to form time series of (e.g., mean, stdv, or sum)

        Returns
        -------
        tuple (float)
           min, max range of figure; set by instruments[obType]
        '''
        self.figName='%s_%s_%s_%s_%s.%s'%(self.prefix,self.obType,self.region.name,self.typ,self.yyyymm,self.figType)

        if self.includeObsCounts:
            self.fig, [self.ax, self.ax_counts] = plt.subplots(1,2, sharey=True, figsize=(9,10), gridspec_kw={'width_ratios' :[0.75,0.25]})
        else:
            self.fig, self.ax = plt.subplots(1,1, figsize=(9,10))

        self.figTitleSize=16
        self.minLabelSize=4
        self.maxLabelSize=12
        self.opacity=0.4
        self.bar_width=0.4
        self.error_config = {'ecolor': '0.3'}
        self.ax.margins(y=0)

        for ax in [self.ax, self.ax_counts] if self.includeObsCounts else [self.ax]:
            ax.set_facecolor('#CECECE')
            # Fine tune xtick labels
            for label in ax.get_xticklabels():
                label.set_fontsize(self.maxLabelSize)
                label.set_fontweight('bold')

        # Fine tune ytick labels
        yticks = range(0,self.getDim('lev'))
        for label in self.ax.get_yticklabels():
            label.set_fontsize(self.maxLabelSize) if self.getDim('lev') <= 30 else label.set_fontsize(self.minLabelSize)
            label.set_fontweight('bold')

        vunits = instruments[self.obType].vertUnits
        ylabel='Pressure (%s)'%vunits if vunits == 'hPa' else 'Channel Index'
        xlabels=['(x %s)'%self.scale]
        if units != "1":
            xlabels[0] = units+r'$\left(\sigma_{EXP}/\sigma_{CTL}\right)$' if units == "%" else xlabels[0] + ' (%s' % units + ')'

        self.ax.set_xlabel(xlabels[0],fontsize=18)
        self.ax.set_ylabel(ylabel,fontsize=18)
        self.ax.set_yticks(yticks)
        self.ax.set_yticklabels(int32(self.loc['lev'][yticks]), minor=False, rotation=0)

        prettyTimeWindow = '%i/%i'%(self.mnths[0],self.yrs[0])
        if len(self.mnths) > 1: prettyTimeWindow += ' - %i/%i'%(self.mnths[-1],self.yrs[-1])

        # Always make temporal window centered on figure
        self.fig.suptitle('%s\nTime Window: %s\n'%(self.obType.upper(),prettyTimeWindow), fontsize=self.figTitleSize)

        # Conditional to handle which axis is used to include region information
        if self.includeObsCounts:
            xlabels.append('Obs Count')
            self.ax_counts.set_xlabel(xlabels[1],fontsize=18)
            self.ax_counts.set_title('%s\n'%self.region, loc='right', fontsize=self.figTitleSize)
        else:
            self.ax.set_title('%s\n'%self.region, loc='right', fontsize=self.figTitleSize, pad=4)


        # Include instrument hame and lat/lon coordinates on figure
        _x0,_y0,_width,_height=self.ax.get_position().bounds

        minval, maxval = float(instruments[self.obType]._min), float(instruments[self.obType]._max)
        self.ax.set_xlim([minval,maxval]) if minval != maxval else self.ax.set_xlim([amin(allStats),\
                                                                                     amax(allStats)])
        self.fig.tight_layout()
        return minval, maxval

    def saveFig(self):
        '''
        Saves the current figure should it exist

        Returns
        -------
        None
        '''
        if self.figExist: plt.savefig(self.figName)


class GritasVars:
    '''
    Common variables to all netCDF files produced by GrITAS

    Attributes
    ----------
    var : numpy.ndarray
       Statistics associated with a variable name at all horizontal/vertical gridded cells

    mean : numpy.ndarray
       Sample mean at a vertical level and across all horizontal gridded cells

    stdv : numpy.ndarray
       Sample standard deviation at a vertical level and across all horizontal gridded cells

    nobs : numpy.ndarray
       Number of observations present in each horizontal gridded cell at a set vertical level

    chisqr : numpy.ndarray
       Left Chi-Square value at a set vertical level (same value per horizontal gridded cell)

    chisql : numpy.ndarray
       Right Chi-Square value at a set vertical level (same value per horizontal gridded cell)

    tstud : numpy.ndarray
       Student-t value at a set vertical level (same value per horizontal gridded cell)
    '''
    def read(self,nc4Var,confidence,idx):
        '''
        Set mean, stdv, nobs, chisqr, chisql, tstud

        Parameters
        ----------
        nc4Var : netCDF4.Dataset.variables
           Gridded and time-averaged statistics, associated with a variable name (e.g., instrument name), read from netCDF

        confidence : bool
           Read Chi-Square/Student-t scores to form confidence intervals

        idx : int
           Vertical level where gridded stats should be read

        Returns
        -------
        None
        '''
        self.var = nc4Var[:,:,:,:]
        self.mean = self.var[0,idx,:,:]
        self.stdv = self.var[1,idx,:,:]
        self.nobs = self.var[2,idx,:,:]
        if confidence:
            self.chisqr = self.var[3,idx,:,:]
            self.chisql = self.var[4,idx,:,:]
            self.tstud  = self.var[5,idx,:,:]

    def sliceVar(self,var,levSlice,latSlice,lonSlice):
        '''
        Access a statistic at a subset of vertical levels, longitude, or latitude

        Parameters
        ----------
        var : numpy.ma.core.MaskedArray
           Statistic over all gridded cells from which a subset will be isolated

        levSlice : list
           Indices denoting contiguous collection of vertical levels wherein var will be considered

        latSlice : list
           Indices denoting contiguous collection of latitudes wherein var will be considered

        lonSlice : list
           Indices denoting contigious collection of longitude wherein var will be considered

        Returns
        -------
        var : numpy.ma.core.MaskedArray
           Original var replaced on return by subset defined by either levSlice, latSlice, or lonSlice

        Raises
        ------
        ValueError : levSlice, latSlice, and lonSlice are mutually exclusive
        '''
        if int(bool(levSlice))+int(bool(latSlice))+int(bool(lonSlice)) != 1:
            raise ValueError("Must select levSlice, latSlice, or lonSlice exclusively")

        if levSlice:   var = var[levSlice,:,:]
        elif latSlice: var = var[:,latSlice,:]
        elif lonSlice: var = var[:,:,lonSlice]
        else:          var = var[:,:,:]

        return var

class Gritas(GritasFig,GritasVars):
    '''
    Class to control production of figures from GrITAS-produced netCDF files

    Parameters
    ----------
    fname : str
       Filename containing gridded time-averaged statistics

    loc : dict
       Grid coordinates - 'lev', 'lat', 'lon' keys map to lists of global coordinates

    var : netCDF4.Dataset.variables
       Gridded and time-averaged statistics, associated with a variable name (e.g., instrument name), read from netCDF

    mean : numpy.ndarray
       Sample mean at a vertical level and across all horizontal gridded cells

    stdv : numpy.ndarray
       Sample standard deviation at a vertical level and across all horizontal gridded cells

    nobs : numpy.ndarray
       Number of observations present in each horizontal gridded cell at a set vertical level

    chisqr : numpy.ma.core.MaskedArray
       Unnormalized (i.e., value repeated within each horizontal grid cell) Right Chi-Square score

    chisql : numpy.ma.core.MaskedArray
       Unnormalized (i.e., value repeated within each horizontal grid cell) Left Chi-Square score

    tstud : numpy.ma.core.MaskedArray
       Unnormalized (i.e, value repeated within each horizontal grid cell) Student-t score

    leftChi2 : numpy.ndarray
       Normalized Left Chi-Square score used to compute a confidence interval

    rightChi2 : numpy.ndarray
       Normalized Right Chi-Square score used to compute a confidence interval

    studT : numpy.ndarray
       Normalized Student-t score used to compute a confidence interval

    supportedStats : list
       Raw statistics supported in fname

    figExist : bool
       Whether Gritas instance has a figure associated with it
    '''
    def __init__(self,fname,supportedStats=None):
        '''
        Initialize a Gritas instance

        Parameters
        ----------
        fname : str
           Filename containing gridded time-averaged statistics

        supportedStats : list
           Raw statistics supported in fname (e.g., mean, stdv, etc.)
        '''
        self.fname=fname
        self.loc={}
        self.var=None
        self.mean=None
        self.stdv=None
        self.nobs=None
        self.chisqr=None
        self.chisql=None
        self.tstud=None
        self.leftChi2=None
        self.rightChi2=None
        self.studT=None
        self.supportedStats=supportedStats
        self.figExist=False

    def __repr__(self):
        '''
        String representation of Gritas instance

        Returns
        -------
        str
        '''
        return "NetCDF file "+self.fname+" contains ( levels = %i, nlat = %i, nlon = %i )"%self.dims()

    def getDim(self,var):
        '''
        Report dimension along a single axis of either levels, latitude, or longitude

        Parameters
        ----------
        var : str
           Either 'lev', 'lat', or 'lon'

        Returns
        -------
        int
        '''
        return len(self.loc[var])

    def dims(self):
        '''
        Report dimensions of levels, latitude, and longitude

        Returns
        -------
        tuple (int)

        Raises
        ------
        ValueError : GrITAS file must first be read from
        '''
        try:
            return self.getDim('lev'), self.getDim('lat'), self.getDim('lon')
        except:
            raise ValueError("Unable to report dims! - Must read from GrITAS file first!")

    def fromGritas(self,var,confidence,vunits):
        '''
        Read data from self.fname stored at 'var'

        Parameters
        ----------
        var : str
           Instrument (e.g. variable name) to be read from self.fname

        confidence : bool
           Form confidence intervals after read

        vunits : str
           Vertical units used to label vertical axes of figure(s)

        Returns
        -------
        self
        '''
        f = nc4.Dataset(self.fname,'r', format='NETCDF4')

        # Pick out latitude, longitude and level values stored
        self.loc={k:f.variables[k][:] for k in ['lat','lon','lev']}

        # Reverse order of masked array 'lev', if vertical units are not 'hPa'
        self.loc['lev'], idx = revMaskedArray(self.loc['lev'], ( vunits != 'hPa' ) )

        # Read remaining variables - i.e., the statistics
        try:
            self.read(f.variables[var],confidence,idx)
        except:
            warnings.warn("  > Found instrument name in netCDF file does not match obtype in yaml! Will extract from netCDF.variables keys.")
            varList = list(f.variables.keys())
            self.read(f.variables[varList.pop(-1)],confidence,idx)
        f.close()
        return self

    def getStat(self,stat,levSlice=None,latSlice=None,lonSlice=None,threshold=None,rescale=1.0,mask='nobs'):
        '''
        Get a statistic in a region of globe

        Parameters
        ----------
        stat : str
           Raw statistic under consideration (e.g., mean, stdv, etc.)

        levSlice : list
           Indices denoting contiguous collection of vertical levels wherein stat will be considered

        latSlice : list
           Indices denoting contiguous collection of latitude wherein stat will be considered

        lonSlice : list
           Indices denoting ontiguous collection of longitude wherein stat will be considered

        threshold : int
           Minimum number of observations in a gridded cell to be included in statistics

        rescale : float
           Multiplicatively rescale statistics

        mask : str
           Mask statistics based on 'mask' (defaults to 'nobs'). No mask is applied, if 'mask' != 'nobs'

        Returns
        -------
        numpy.ndarray : statistic (potentially rescaled) within subset of levels, latitudes, or longitudes

        Raises
        ------
        ValueError : stat does not match those supported in netCDF produced by GrITAS
        '''
        if stat not in self.supportedStats:
            raise ValueError("Statistics flavor %s not supported!"%stat)

        # Select the correct member variable based on 'stat'
        var=None
        # if stat == 'sum' or stat == 'mean': var = self.mean
        # if stat == 'stdv':                  var = self.stdv

        if stat == 'sum': var = self.nobs
        if stat == 'mean': var = self.mean
        if stat == 'stdv': var = self.stdv

        # Slice the member variable
        var=self.sliceVar(var,levSlice,latSlice,lonSlice)

        # Convenience - Local dimensions of 'var'
        varNLev, varNLat, varNLon = var.shape

        levelStat=np.zeros(varNLev)
        for n in range(varNLev):
            # Form mask
            ID=np.where(self.nobs[n,latSlice,:]>threshold) if mask == 'nobs' else []

            if stat == 'stdv':
                levelStat[n]=np.sqrt( (var[n,:,:][ID]*var[n,:,:][ID]).sum()/(varNLat*varNLon-1) )
            elif stat == 'mean':
                levelStat[n]=var[n,:,:][ID].sum()/(varNLat*varNLon)
            else:
                levelStat[n]=var[n,:,:][ID].sum()

        return rescale*levelStat

    def getConfidence(self,levSlice=None,latSlice=None,lonSlice=None,threshold=None):
        '''
        Access test statistic scores to form confidence intervals

        Parameters
        ----------
        levSlice : list
           Indices denoting contiguous collection of vertical levels wherein test statistic scores should be accessed

        latSlice : list
           Indices denoting contiguous collection of latitudes wherein test statistic scores should be accessed

        lonSlice : list
           Indices denoting contiguous collection of longitude wherein test statistic scores should be accessed

        threshold : int
           Number of observations within a gridded cell in order for cell to be included

        Returns
        -------
        None
        '''
        # Convenience
        nobs = self.sliceVar(self.nobs,levSlice,latSlice,lonSlice)
        chisql=self.sliceVar(self.chisql,levSlice,latSlice,lonSlice)
        chisqr=self.sliceVar(self.chisqr,levSlice,latSlice,lonSlice)
        tstud= self.sliceVar(self.tstud,levSlice,latSlice,lonSlice)
        nLev, nLat, nLon = nobs.shape

        # Test statistic scores per level
        leftChi2, rightChi2, studT, sampleSize = np.zeros(nLev), np.zeros(nLev), np.zeros(nLev), np.zeros(nLev)

        for n in range(nLev):
            # Form mask
            ID=np.where(nobs[n,:,:]>threshold)
            sampleSize[n] = nobs[n,:,:][ID].sum()
            leftChi2[n] = sqrt((sampleSize[n]-1)/lvlAvg(chisql[n,:,:],ID))
            rightChi2[n] = sqrt((sampleSize[n]-1)/lvlAvg(chisqr[n,:,:],ID))
            studT[n] = lvlAvg(tstud[n,:,:],ID)/sqrt(sampleSize[n])

        self.leftChi2 = leftChi2
        self.rightChi2 = rightChi2
        self.studT = studT

    def plotInit(self,prefix,obType,region,scale,linePlot,simpleBars,includeObsCounts,yrs,mnths):
        '''
        Initialize a figure

        Parameters
        ----------
        prefix : str
           Prefix of figure to be saved

        obType : str
           Instrument name belonging to experiment(s)

        region : plot_util.Region instance
           A geographic region figure will correspond to

        scale : float
           Rescale statistics

        linePlot : bool
           Plot stats via a line plot; bars used otherwise

        simpleBars : bool
           Embelished confidence intervals

        includeObsCounts : bool
           Observation counts are among desired stats to show

        yrs : list
           List of years in time window

        mnths : list
           List of months in time window

        Returns
        -------
        None
        '''
        super().__init__(self,prefix,obType,region,scale,linePlot,simpleBars,includeObsCounts,yrs,mnths)
        self.figExist=True
