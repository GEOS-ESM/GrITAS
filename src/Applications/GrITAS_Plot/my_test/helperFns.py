#!/usr/bin/env python

import sys
from pylab import *
import matplotlib as mpl
import matplotlib.pyplot as plt
from datetime import timedelta

import netCDF4 as nc4


# Global settings for matplotlib
mpl.rc('font', family='serif')

def revMaskedArray(arr,action):
    if action:
        idx = range(len(arr)-1,-1,-1)
        return arr[idx], idx
    else:
        return arr, range(1,len(arr),1)


def adjustFigAspect(fig,aspect=1):
    '''
    Adjust the subplot parameters so that the figure has the correct
    aspect ratio.
    '''
    xsize,ysize = fig.get_size_inches()
    minsize = min(xsize,ysize)
    xlim = .4*minsize/xsize
    ylim = .4*minsize/ysize
    if aspect < 1:
        xlim *= aspect
    else:
        ylim /= aspect
    fig.subplots_adjust(left=.5-xlim,
                        right=.5+xlim,
                        bottom=.5-ylim,
                        top=.5+ylim)

class gritasFig:
    def __init__(self,prefix,obType,region,scale,typ='invalid',yrs=None,mnths=None):
        self.prefix='X'.join(prefix)
        self.obType=obType
        self.region=region
        self.scale=scale
        self.typ=typ
        self.figType='png'
        # Either a single year/month or (pandas) collection of several
        self.yrs=yrs
        self.mnths=mnths

        self.figName='%s_%s_%s_%s_%s.%s'%(self.prefix,self.typ,self.obType,self.region,\
                                          str(self.yrs)+self.mnths,self.figType)\
            if self.typ == 'monthly' else\
               '%s_%s_%s_%s_%s.%s'%(self.prefix,self.typ,self.obType,self.region,'FIXME-DUMSTAT',self.figType)

        self.fig = plt.figure(figsize=(8,10)) if self.typ == 'monthly' else plt.figure(figsize=(10,10))
        self.ax=self.fig.gca()
        self.ax.set_facecolor('#CECECE')
        self.minLabelSize=4
        self.maxLabelSize=7

        self.opacity=0.5
        self.bar_width=0.6
        self.error_config = {'ecolor': '0.3'}

    def monthlyBars(self,allStats,stats=None,instruments=[],flavor='None'):

        self.commonFigSetup(allStats,stats.units,instruments,flavor=flavor)

        for ns, stat in enumerate(stats.measures):
            omean=allStats[:,ns]

            pos = arange(len(omean))+ns*self.bar_width # hack to center the labels on the center of the bar

            kwargs={'alpha': self.opacity,'color': stats.colors[ns],'label': stat,\
                    'error_kw': self.error_config}
            if stats.confidence:
                xerr=[self.confl,self.confr]
                if stat == 'mean':
                    if stats.measures[1]=='stdv':
                        xerr=self.studt*allStats[:,1]
                    else:
                        raise ValueError("Inconsistent options, aborting...")
                kwargs.update({'xerr': xerr})

            self.ax.barh(pos, omean, self.bar_width, **kwargs)

        self.ax.legend(loc='lower left')
        self.ax.set_title(str(self.yrs)+self.mnths)


    def tSeries(self,allStats,stats=None,instruments=[],flavor='None'):
        # Capture minval/maxval
        minval,maxval = self.commonFigSetup(allStats,stats.units,instruments,flavor=flavor)

        rval = maxval/minval if minval < 0.0 else minval/maxval
        if rval < 0.0:
            maxval= max(abs(minval),abs(maxval))
            minval=-maxval
            cs=plt.pcolor(allStats,cmap=plt.cm.seismic,vmin=minval,vmax=maxval)
        else:
            cs=plt.pcolor(allStats,cmap=plt.cm.binary,vmin=minval,vmax=maxval)
        cbar=plt.colorbar(cs)

        x1 = range(0,len(self.yrs))

        self.ax.axes.xaxis.set_visible(True)
        self.ax.axes.yaxis.set_visible(True)
        self.ax.set_xlabel('Date')
        self.ax.set_xticks(x1)
        yyyymm = 100*self.yrs + self.mnths
        self.ax.set_xticklabels(int32(yyyymm[x1]), minor=False, rotation=45)


    def monthlyPlot(self,case,allStats,stats=None,instruments=[],annotation=''):
        # Capture minval/maxval
        minval, maxval = self.commonFigSetup(allStats,stats.units,instruments,flavor=annotation)

        midpnt = minval+0.5*(maxval-minval)

        pos = arange(len(allStats)) #+nf*bar_width

        # # Set the x-label
        # mylabel = '(x %s)' % str(1/scale)
        # if units != "1":
        #     mylabel = units if units == "%" else mylabel + ' (%s' % units + ')'
        # self.ax.set_xlabel(myLabel,fontsize=12)


        if stats.confidence:
            # xerr based on case
            xerr = self.studt if case == 'difference' or case == 'mean' else [self.confl,self.confr]

            if case == 'difference':
                self.ax.plot(zeros(len(allStats)), pos, color='k', lw=1, alpha=0.8)
            if case == 'ratio':
                refline=midpnt*ones(len(allStats))
                self.ax.plot(refline, pos, color='k', lw=1, alpha=0.8)

                self.fig.suptitle('CTL: %s'%annotation[0], x=0.125, y=0.93, ha='left', fontsize=14)
                self.ax.set_title('EXP: %s'%annotation[1], loc='left', fontsize=14)
                self.ax.annotate(self.obType.upper(), xy=(maxval-0.06*maxval,amax(pos)-0.5),
                                 size=12, ha='left', va="bottom")

            # Plot errorbars regardless of case
            self.ax.errorbar(allStats, pos, xerr=xerr, ecolor='k', lw=2)


            # Plot deterioration/improvement boxes
            self.ax.annotate('Deterioration', xy=(maxval-0.03*maxval,0.5),  xycoords='data',
                             xytext=(0,0), textcoords='offset points',
                             size=13, ha='right', va="center",
                             bbox=dict(boxstyle="round", alpha=0.1, color='r'))
            self.ax.annotate('Improvement',   xy=(midpnt-0.03*midpnt,0.5),  xycoords='data',
                             xytext=(0,0), textcoords='offset points',
                             size=13, ha='right', va="center",
                             bbox=dict(boxstyle="round", alpha=0.1, color='g'))

        else:
            self.ax.plot(allStats,pos)



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
        self.ax.set_ylabel(ylabel, fontsize=12)
        self.ax.set_yticks(yticks)
        self.ax.set_yticklabels(int32(self.loc['lev'][yticks]), minor=False, rotation=0)


        mytitle = '%s (x %s)'%(flavor,str(1/self.scale)) if self.typ == 'tseries'\
                  else '%s (x %s)'%(self.obType,str(1/self.scale))
        mytitle += ' (%s)'%units if units != "1" else ''
        self.ax.set_title(mytitle)

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
    def __init__(self,fname):
        self.fname=fname
        self.loc={}
        self.var=None
        self.mean=None
        self.stdv=None
        self.nobs=None
        self.chisqr=None
        self.chisql=None
        self.tstud=None

        self.confl=None
        self.confr=None
        self.studt=None


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
        f = nc4.Dataset(self.fname,'r', format='NETCDF4')

        self.loc={k:f.variables[k][:] for k in ['lat','lon','lev']}

        # Reverse order of masked array 'lev', if vertical units are not 'hPa'
        self.loc['lev'], idx = revMaskedArray(self.loc['lev'], ( vunits != 'hPa' ) )
        # Read remaining variables - ie. the statistics
        self.read(f.variables[var],confidence,idx)

        f.close()
        return self


    def getStat(self,stat,levSlice=None,latSlice=None,lonSlice=None,threshold=None,mask='nobs'):
        # Check for valid statistic
        if stat != 'sum' and stat != 'mean' and stat != 'stdv':
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
        return levelStat


    def getConfidence(self,levSlice=None,latSlice=None,lonSlice=None,threshold=None):
        # Slice class member 'nobs'
        nobs = self.sliceVar(self.nobs,levSlice,latSlice,lonSlice)
        # Convenience
        chisql=self.sliceVar(self.chisql,levSlice,latSlice,lonSlice)
        chisqr=self.sliceVar(self.chisqr,levSlice,latSlice,lonSlice)
        tstud= self.sliceVar(self.tstud,levSlice,latSlice,lonSlice)

        # Local dimensions of 'nobs'
        nLev, nLat, nLon = nobs.shape

        confl, confr, studt, nccum = np.zeros(nLev), np.zeros(nLev), np.zeros(nLev), np.zeros(nLev)

        for n in range(nLev):
            ID=np.where(nobs[n,:,:]>threshold)
            nccum[n] = nobs[n,:,:][ID].sum()
            confl[n] = sqrt((nccum[n]-1)/(nLat*nLon*chisql[n,1,1]))
            confr[n] = sqrt((nccum[n]-1)/(nLat*nLon*chisqr[n,1,1]))
            studt[n] = tstud[n,1,1]/sqrt(nccum[n]/(nLat*nLon))

        self.confl = confl
        self.confr = confr
        self.studt = studt


    def plotInit(self,prefix,obType,region,scale,typ,yrs,mnths):
        gritasFig.__init__(self,prefix,obType,region,scale,typ,yrs,mnths)
