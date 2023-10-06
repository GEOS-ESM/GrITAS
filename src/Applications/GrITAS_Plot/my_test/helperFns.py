#!/usr/bin/env python

import sys
from pylab import *
import matplotlib as mpl
import matplotlib.pyplot as plt
from datetime import timedelta

import netCDF4 as nc4


# Global settings for matplotlib
mpl.rc('font', family='serif')
# mpl.rc('text', usetex=True)
# mpl.rcParams.update(mpl.rcParamsDefault)

def revMaskedArray(arr,action):
    if action:
        idx = range(len(arr)-1,-1,-1)
        return arr[idx], idx
    else:
        return arr, range(1,len(arr),1)


def lvlAvg(arr,mask):
    '''
    Return sum(arr)/dim(arr) where elements of arr included in sum and dim are those identified by mask
    -> shape(mask)[1] equates to number of lat x lon points for which elements of arr meet mask criterion
    '''
    return arr[mask].sum()/shape(mask)[1]


class gritasFig:
    def __init__(self,prefix,obType,region,scale,simpleBars,typ='invalid',yrs=None,mnths=None):
        self.prefix=''.join(prefix) if isinstance(prefix, str) else 'X'.join(prefix)
        self.obType=obType
        self.region=region
        self.scale=scale
        self.simpleBars=simpleBars
        self.typ=typ
        self.figType='png'
        # Either a single year/month or (pandas) collection of several
        self.yrs=yrs
        self.mnths=mnths

        self.yyyymm='%s%s'%(yrs[0],str(mnths[0]).rjust(2,'0'))
        if len(yrs) > 1 or len(mnths) > 1:
            self.yyyymm += '-%s%s'%(yrs[-1],str(mnths[-1]).rjust(2,'0'))


        self.figName='%s_%s_%s_%s_%s.%s'%(self.prefix,self.obType,self.region.name,self.typ,\
                                          self.yyyymm,self.figType)\
            if self.typ == 'monthly' else\
               '%s_%s_%s_%s.%s'%(self.prefix,self.typ,self.obType,self.region.name,self.figType)

        self.fig = plt.figure(figsize=(8,10)) if self.typ == 'monthly' else plt.figure(figsize=(10,10))
        self.fig.tight_layout(pad=1.0)
        self.ax=self.fig.gca()
        self.ax.set_facecolor('#CECECE')
        self.minLabelSize=4
        self.maxLabelSize=7

        self.opacity=0.5
        self.bar_width=0.4 #1.0/3 #0.6
        self.error_config = {'ecolor': '0.3'}

    def monthlyStat(self,allStats,stats=None,instruments=[],flavor='None',annotation=None):
        print("IN MONTHLYSTAT (formerly monthlyBars)")
        self.commonFigSetup(allStats,stats.units,instruments,flavor=flavor)


        # Iterate over all supported stats, skipping if not desired
        for ns, stat in enumerate(self.supportedStats):
            if stat not in stats.measures:
                continue
            print("--- with ns, stat = %i, %s"%(ns,stat))
            statAllLvls=allStats[:,ns]

            # Set center of bars
            barCenter = arange(len(statAllLvls))-1.0/5+ns*self.bar_width

            kwargs={'alpha': self.opacity,'color': stats.colors[ns],'label': stat,'error_kw': self.error_config}
            # Confidence on monthlyStat
            # --------------------------
            if stats.confidence:
                '''
                Confidence level assigned to mean:
                    <sample variance> * <test statistic>/sqrt(<# observations - i.e. sample size!>
                The measurements form a random sample {x1,x2,...,xn} drawn from (assumed) a normal distribution
                with (pop.) mean \mu and (pop.) variance \sigma^2. As the population mean and variance of the
                underlying distribution are unknown, we must use a t-score (via t-distribution) to assign a
                (1-\alpha) confidence level for the computed SAMPLE MEAN.
                   CL = [ xbar - t(\alpha/2|n-1)*(s/\sqrt(n)), xbar + t(\alpha/2|n-1)*(s/\sqrt(n)) ]

                where:  - xbar = sample mean
                        - t(\alpha/2|n-1) is area = \alpha/2 of t-distribution w/ n-1 dof (i.e. prob. of test statistic
                          to exceed value st. prob. is \alpha/2)
                        - s^2 = sample variance

                ======

                Confidence level assigned to variance:
                Suppose we want to estimate the variance of (assumed) normal distribution (w/ pop. mean \mu and pop.
                variance \sigma^2) from which random sample {x1,x2,...,xn} are drawn. Assuming pop. mean/variance are
                unknown, the random variable Q = (n-1)s^2/\sigma^2 is chi2 distributed w/ n-1 dof. Thus, we can use
                left/right chi-squared scores to assign a (1-\alpha) confidence level for the computed SAMPLE VARIANCE.

                     CL = [ (n-1)s^2 / chi2(\alpha/2|n-1), (n-1)s^2 / chi2(1-\alpha/2|n-1) ]

                where:  - s^2 = sample variance
                        - n = sample size
                        - chi2(\alpha/2|n-1) is value of Chi2 dist. w/ n-1 dof for which area is \alpha/2
                        - chi2(1-\alpha/2|n-1) is value of Chi2 dist. w/ n-1 dof for which area is 1-\alpha/2
                '''

                # Set xerr to studT score if forming CL for pop. mean, and left/right Chi2 scores if forming
                # CL for pop. stdv.
                xerr=self.studT if stat == 'mean' else [self.leftChi2,self.rightChi2]

                # Assign no error if a stat other than pop. mean/stdv is being estimated
                if stat != 'mean' and stat != 'stdv':
                    xerr=np.ones(np.size(allStats[:,1]))

                # Multiply by sample stdv., thus completing formation of CL
                xerr*=allStats[:,1]


            # Plot bars
            self.ax.barh(barCenter, statAllLvls, self.bar_width, **kwargs)


            # Toggle between simple or prettier bars
            if self.simpleBars:
                self.ax.errorbar(statAllLvls, barCenter, xerr=xerr, color=stats.colors[ns], capsize=5,ls='')
            else:
                self.ax.vlines(statAllLvls,barCenter-0.6*self.bar_width,barCenter+0.6*self.bar_width,color='k')
                if stat == 'mean':
                    self.ax.barh(barCenter, xerr, self.bar_width,left=statAllLvls-xerr,color=stats.colors[ns],\
                                 alpha=1.0,hatch='/////',edgecolor=stats.colors[ns])
                    self.ax.barh(barCenter, xerr, self.bar_width,left=statAllLvls,color=stats.colors[ns],\
                                 alpha=1.0,hatch='/////',edgecolor=stats.colors[ns])
                else:
                    self.ax.barh(barCenter, xerr[0], self.bar_width, left=statAllLvls-xerr[0],\
                                 color=stats.colors[ns],alpha=1.0,hatch='/////',edgecolor=stats.colors[ns])
                    self.ax.barh(barCenter, xerr[1], self.bar_width, left=statAllLvls, color=stats.colors[ns],\
                                 alpha=1.0,hatch='/////',edgecolor=stats.colors[ns])

        # Add annotations if not None
        # Labeling of figure
        # -------------------
        if annotation:
            self.fig.suptitle('CTL: %s'%annotation[0], x=0.125, y=0.93, ha='left', fontsize=14)
            self.ax.set_title('EXP: %s'%annotation[1], loc='left', fontsize=14)

        self.ax.legend(loc='lower right')



    def tSeries(self,allStats,stats=None,instruments=[],flavor='None'):
        # Make copies of allStats for masking
        slices=[allStats.copy(), allStats.copy(), allStats.copy()]
        # Mask
        for n,s in enumerate(slices):
            s[:,(n+1)%len(slices)]=s[:,(n+2)%len(slices)]=np.nan

        # Capture minval/maxval
        minval,maxval = self.commonFigSetup(allStats,stats.units,instruments,flavor=flavor)

        rval = maxval/minval if minval < 0.0 else minval/maxval
        if rval < 0.0:
            maxval= max(abs(minval),abs(maxval))
            minval=-maxval

            cs_mean=self.ax.pcolor(slices[0],cmap=plt.cm.Blues,vmin=np.nanmin(slices[0]),\
                                   vmax=np.nanmax(slices[0]))
            cs_stdv=self.ax.pcolor(slices[1],cmap=plt.cm.autumn_r,vmin=np.nanmin(slices[1]),\
                                   vmax=np.nanmax(slices[1]))
            cs_summ=self.ax.pcolor(slices[2],cmap=plt.cm.Greens,vmin=0,vmax=np.nanmax(slices[2]))

        else:
            cs=plt.pcolor(DUM,cmap=plt.cm.binary,vmin=minval,vmax=maxval)
        # Add colorbars to plot
        cbar_summ=plt.colorbar(cs_summ); cbar_summ.ax.set_xlabel('# Obs.')
        cbar_stdv=plt.colorbar(cs_stdv); cbar_stdv.ax.set_xlabel('Stdv')
        cbar_mean=plt.colorbar(cs_mean); cbar_mean.ax.set_xlabel('Mean')

        self.ax.axes.xaxis.set_visible(True)
        self.ax.axes.yaxis.set_visible(True)
        self.ax.set_xlabel('Date')
        self.ax.set_xticks(range(len(slices)*len(self.yrs)+1))

        yyyymm = ['/'.join(["{:02d}".format(m),str(y)]) for m,y in zip(self.mnths,self.yrs)]
        yyyymm.insert(0,'') # lazy - major_locator below is causing first entry to be dropped

        self.ax.axes.xaxis.set_major_locator(MultipleLocator(3))
        self.ax.axes.xaxis.set_minor_locator(MultipleLocator(1))
        self.ax.set_xticklabels(yyyymm, minor=False, rotation=45)

        self.ax.vlines(self.ax.get_xticks(),self.ax.get_ylim()[0],self.ax.get_ylim()[1],color='darkgray',linestyle=':')


    def monthlyComp(self,compareVia,expStats,cntlStats,stats=None,instruments=[],annotation=''):
        print("IN MONTHLYCOMP (formerly monthlyPlot)")
        allStats=100*(expStats/cntlStats) if compareVia == 'ratio' else expStats - cntlStats

        # Capture minval/maxval
        minval, maxval = self.commonFigSetup(allStats,stats.units,instruments,flavor=annotation)

        # Convenience
        midpnt = minval+0.5*(maxval-minval)
        pos = arange(len(allStats)) #+nf*bar_width

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
            if compareVia == 'ratio':
                refline=midpnt*ones(len(allStats))
                self.ax.plot(refline, pos, color='k', lw=1, alpha=0.8)

            # Plot errorbars regardless of compareVia
            self.ax.errorbar(allStats, pos, xerr=xerr, ecolor=ecolor, lw=2, capsize=5, ls=':', color='gray')
        else:
            self.ax.plot(allStats,pos)

        # Labeling of figure
        # -------------------
        self.fig.suptitle('CTL: %s'%annotation[0], x=0.125, y=0.93, ha='left', fontsize=14)
        self.ax.set_title('EXP: %s'%annotation[1], loc='left', fontsize=14)
        # self.ax.annotate(self.obType.upper(), xy=(maxval-0.06*maxval,amax(pos)-0.5),
        #                  size=12, ha='left', va="bottom")

        # Plot deterioration/improvement boxes #midpnt+0.0195*midpnt
        # -----------------------------------------------------------
        self.ax.annotate('Deterioration', xy=(maxval-0.0195*(maxval-minval),0.5), xycoords='data',
                         xytext=(0,0), textcoords='offset points',
                         size=13, ha='right', va="center",
                         bbox=dict(boxstyle="round", alpha=0.1, color='r'))
        self.ax.annotate('Improvement', xy=(minval+0.0195*(maxval-minval),0.5), xycoords='data',
                         xytext=(0,0), textcoords='offset points',
                         size=13, ha='left', va="center",
                         bbox=dict(boxstyle="round", alpha=0.1, color='g'))



    def commonFigSetup(self,allStats,units,instruments,flavor='None'):
        minval, maxval = float(instruments[self.obType]._min), float(instruments[self.obType]._max)
        self.ax.set_xlim([minval,maxval]) if minval != maxval else self.ax.set_xlim([amin(allStats),\
                                                                                     amax(allStats)])

        self.ax.margins(y=0)
        for label in self.ax.get_xticklabels():
            label.set_fontsize(self.maxLabelSize) #10
            label.set_fontweight('bold')

        yticks = range(0,self.getDim('lev'))
        self.ax.margins(y=0)
        for label in self.ax.get_yticklabels():
            # if vname in ['airs','iasi','cris']:
            #     label.set_fontsize(self.minLabelSize)
            # else:
            #     label.set_fontsize(self.maxLabelSize)

            label.set_fontsize(self.maxLabelSize) if self.getDim('lev') <= 30 else label.set_fontsize(self.minLabelSize)
            label.set_fontweight('bold')


        vunits = instruments[self.obType].vertUnits
        ylabel='Pressure (%s)'%vunits if vunits == 'hPa' else 'Channel Index'
        xlabel= '(x %s)' % str(1/self.scale)
        if units != "1":
            xlabel = units if units == "%" else xlabel + ' (%s' % units + ')'
        self.ax.set_xlabel(xlabel,fontsize=12)
        self.ax.set_ylabel(ylabel,fontsize=12)
        self.ax.set_yticks(yticks)
        self.ax.set_yticklabels(int32(self.loc['lev'][yticks]), minor=False, rotation=0)


        # mytitle = '%s (x %s)'%(flavor,str(1/self.scale)) if self.typ == 'tseries'\
        #           else '%s (x %s)'%(self.obType,str(1/self.scale))
        # mytitle += ' (%s)'%units if units != "1" else ''

        print(self.mnths)
        print(self.yrs)
        prettyTimeWindow = ['/'.join(["{:02d}".format(m),str(y)]) for m,y in zip(self.mnths,self.yrs)]
        print(prettyTimeWindow)
        prettyTimeWindow = '%i/%i'%(self.mnths[0],self.yrs[0])
        if len(self.mnths) > 1: prettyTimeWindow += ' - %i/%i'%(self.mnths[-1],self.yrs[-1])

        self.ax.set_title('%s\nTime Window: %s'%(self.obType.upper(),prettyTimeWindow))

        # Get coordinates of axes
        # ------------------------
        _x0,_y0,_width,_height=self.ax.get_position().bounds
        # Include instrument hame and lat/lon coordinates on figure
        # ----------------------------------------------------------
        # self.fig.text(_x0,_y0+_height+0.005,"%s"%self.obType,ha='left',fontsize=14)
        self.fig.text(_x0+_width,_y0+_height+0.005,'%s\n%s'%(self.region.name,self.region),ha='right',fontsize=14)

        return minval, maxval

    def saveFig(self):
        plt.savefig(self.figName)


class gritasVars:
    def read(self,nc4Var,confidence,idx):
        self.var = nc4Var[:,:,:,:]
        self.mean = self.var[0,idx,:,:]
        self.stdv = self.var[1,idx,:,:]
        self.nobs = self.var[2,idx,:,:]
        if confidence:
            self.chisqr = self.var[3,idx,:,:]
            self.chisql = self.var[4,idx,:,:]
            self.tstud  = self.var[5,idx,:,:]

    def sliceVar(self,var,levSlice,latSlice,lonSlice):
        if not bool(levSlice)^bool(latSlice)^bool(lonSlice):
            raise ValueError("Must select levSlice, latSlice, or lonSlice exclusively")

        if levSlice:   var = var[levSlice,:,:]
        elif latSlice: var = var[:,latSlice,:]
        elif lonSlice: var = var[:,:,lonSlice]
        else:          var = var[:,:,:]

        return var



class Gritas(gritasVars,gritasFig):
    def __init__(self,fname,supportedStats=None):
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
        return "NetCDF file "+self.fname+" contains ( levels = %i, nlat = %i, nlon = %i )"%self.dims()


    # Report dimension along a single axis of either levels, latitude, longitude
    def getDim(self,var):
        return len(self.loc[var])

    # Report tuple of levels, latitude, longitude dimensions
    def dims(self):
        try:
            return self.getDim('lev'), self.getDim('lat'), self.getDim('lon')
        except:
            raise ValueError("Unable to report dims! - Must read from GrITAS file first!")

    def fromGritas(self,var,confidence,vunits):
        '''
        Read data from self.fname stored at 'var'
        - var : instrument name (e.g. atmsnpp)
        - optionally read t-scores and l/r chi-squared scores for assigning confidence intervals
        - vunits : vertical units
        '''
        f = nc4.Dataset(self.fname,'r', format='NETCDF4')

        # Pick out latitude, longitude and level values stored
        self.loc={k:f.variables[k][:] for k in ['lat','lon','lev']}

        # Reverse order of masked array 'lev', if vertical units are not 'hPa'
        self.loc['lev'], idx = revMaskedArray(self.loc['lev'], ( vunits != 'hPa' ) )

        # Read remaining variables - ie. the statistics
        self.read(f.variables[var],confidence,idx)

        f.close()
        return self


    def getStat(self,stat,levSlice=None,latSlice=None,lonSlice=None,threshold=None,rescale=1.0,mask='nobs'):
        # Check for valid statistic
        if stat not in self.supportedStats:
            raise ValueError("Statistics flavor %s not supported!"%stat)

        # Select the correct member variable based on 'stat'
        var=None
        if stat == 'sum' or stat == 'mean': var = self.mean #[:,latSlice,:]
        if stat == 'stdv':                  var = self.stdv #[:,latSlice,:]


        # Slice the member variable
        var=self.sliceVar(var,levSlice,latSlice,lonSlice)

        # Local dimensions of 'var'
        varNLev, varNLat, varNLon = var.shape

        levelStat=np.zeros(varNLev)
        for n in range(varNLev):
            ID=np.where(self.nobs[n,latSlice,:]>threshold) if mask == 'nobs' else []

            if stat == 'stdv':
                levelStat[n]=np.sqrt( (var[n,:,:][ID]*var[n,:,:][ID]).sum()/(varNLat*varNLon-1) )
            elif stat == 'mean':
                levelStat[n]=var[n,:,:][ID].sum()/(varNLat*varNLon)
            else:
                levelStat[n]=var[n,:,:][ID].sum()
        return rescale*levelStat


    def getConfidence(self,levSlice=None,latSlice=None,lonSlice=None,threshold=None):
        # Confidence is repeated for each cell (ie lat x lon) at a given level

        # Slice class member 'nobs'
        nobs = self.sliceVar(self.nobs,levSlice,latSlice,lonSlice)
        # Convenience
        chisql=self.sliceVar(self.chisql,levSlice,latSlice,lonSlice)
        chisqr=self.sliceVar(self.chisqr,levSlice,latSlice,lonSlice)
        tstud= self.sliceVar(self.tstud,levSlice,latSlice,lonSlice)

        # Local dimensions of 'nobs'
        nLev, nLat, nLon = nobs.shape

        # Left/Right Chi2 Scores, Student-T score, and Number of Observations Per Level
        leftChi2, rightChi2, studT, sampleSize = np.zeros(nLev), np.zeros(nLev), np.zeros(nLev), np.zeros(nLev)

        for n in range(nLev):
            ID=np.where(nobs[n,:,:]>threshold)
            sampleSize[n] = nobs[n,:,:][ID].sum()
            # leftChi2[n] = sqrt((sampleSize[n]-1)/(chisql[n,:,:][ID].sum()/(nLat*nLon)))
            # rightChi2[n] = sqrt((sampleSize[n]-1)/(chisqr[n,:,:][ID].sum()/(nLat*nLon)))
            # studT[n] = (tstud[n,:,:][ID].sum()/(nLat*nLon))/sqrt(sampleSize[n])
            leftChi2[n] = sqrt((sampleSize[n]-1)/lvlAvg(chisql[n,:,:],ID))
            rightChi2[n] = sqrt((sampleSize[n]-1)/lvlAvg(chisqr[n,:,:],ID))
            studT[n] = lvlAvg(tstud[n,:,:],ID)/sqrt(sampleSize[n])

        self.leftChi2 = leftChi2
        self.rightChi2 = rightChi2
        self.studT = studT


    def plotInit(self,prefix,obType,region,scale,simpleBars,typ,yrs,mnths):
        gritasFig.__init__(self,prefix,obType,region,scale,simpleBars,typ,yrs,mnths)
        self.figExist=True
