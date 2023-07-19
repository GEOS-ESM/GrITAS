#!/usr/bin/env python

import sys
from pylab import *
import matplotlib as mpl
import matplotlib.pyplot as plt
from datetime import timedelta

import netCDF4 as nc4

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


class files(gritasVars):
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

    def __repr__(self):
        return "NetCDF file "+self.fname+" contains ( levels = %i, nlat = %i, nlon = %i )"%self.dims()
        
    def dims(self):
        try:
            return len(self.loc['lev']), len(self.loc['lat']), len(self.loc['lon'])
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
        
def get_gritas_data():
    pass

#.................................
def vaccum(stype,tresh,var,nob,valsYaml):
  nlev = size(var[:,0,0])
  nlat = size(var[0,:,0])
  nlon = size(var[0,0,:])
  accum = zeros(nlev)
  nccum = zeros(nlev)
  ndim = nlat*nlon
  for k in range(nlev):
      id=where(nob[k,:,:]>tresh)
      this = var[k,:,:]
      nobs = nob[k,:,:]
      nccum[k] = sum(nobs[id])
      if stype=='sum':
         accum[k] = sum(this[id])
      if stype=='mean':
         accum[k] = sum(this[id])
         accum[k]=accum[k]/(ndim)
      if stype=='stdv':
         accum[k] = sum(this[id]*this[id])
         accum[k]= sqrt(accum[k]/(ndim-1))

  stat=valsYaml['global']['statistics flavor']
  if stat == 'Impact per Ob' or stat == 'DFS per Ob':
     accum = accum/nccum
  if stat == 'Ob count':
     accum = nccum
  
  return accum
#.................................
def getconf(tresh,nob,chisql,chisqr,tstud):
  nlev = size(nob[:,0,0])
  nlat = size(nob[0,:,0])
  nlon = size(nob[0,0,:])
  ndim = nlat*nlon
  confl = zeros(nlev)
  confr = zeros(nlev)
  studt = zeros(nlev)
  nccum = zeros(nlev)
  for k in range(nlev):
      id=where(nob[k,:,:]>tresh)
      nobs = nob[k,:,:]
      nccum[k] = sum(nobs[id])
      confl[k] = sqrt((nccum[k]-1)/(ndim*chisql[k,1,1]))
      confr[k] = sqrt((nccum[k]-1)/(ndim*chisqr[k,1,1]))
      studt[k] = tstud[k,1,1]/sqrt(nccum[k]/ndim)
 
  return confl, confr, studt
      

#.................................
def show_plot_monthly_one(lev,vals,anno,colors,case,scale,this,chisql,chisqr,studt,valsYaml):

    fig = plt.figure(figsize=(8, 10))
    ax = fig.add_subplot(111)
    ax.set_facecolor('#CECECE')

    confidence=valsYaml['global']['confidence']
 
    nf=0
    color = colors
    mpl.rc('font', family='serif')
    opacity = 0.5
    bar_width = 0.6
    label_size = 8
    error_config = {'ecolor': '0.3'}

    obtype=valsYaml['global']['obtype']
    minval = float(valsYaml['global']['configure'][obtype]['min value'])
    maxval = float(valsYaml['global']['configure'][obtype]['max value'])
    if minval == maxval:
       minval = amin(vals)
       maxval = amax(vals)
    midpnt = minval+0.5*(maxval-minval)

    omean = vals
    refline=midpnt*ones(len(omean))
    pos = arange(len(omean))+nf*bar_width # hack to center the labels on the center of the bar
    if confidence==False:
       ax.plot(omean,pos)
    else:
       if case == 'difference' or case == 'mean':
          ax.errorbar(omean, pos, xerr=studt, ecolor='k',linewidth=2)

          if case == 'difference':
             ax.plot(zeros(len(omean)), pos, color='k',linewidth=1, alpha=0.8)

             ax.annotate('Deterioration', xy=(maxval-0.03*maxval,0.5),  xycoords='data',
                xytext=(0,0), textcoords='offset points',
                size=13, ha='right', va="center",
                bbox=dict(boxstyle="round", alpha=0.1, color='r'))
             ax.annotate('Improvement',   xy=(midpnt-0.03*midpnt,0.5),  xycoords='data',
                xytext=(0,0), textcoords='offset points',
                size=13, ha='right', va="center",
                bbox=dict(boxstyle="round", alpha=0.1, color='g'))
       else:
          ax.errorbar(omean, pos, xerr=[confl,confr], ecolor='k', linewidth=2)

          if case == 'ratio':
             ax.plot(refline, pos, color='k',linewidth=1, alpha=0.8)

             suptitle('CTL: '+anno[0], x=0.125, y=0.93, ha='left', fontsize=14)
             title   ('EXP: '+anno[1],                 loc='left', fontsize=14)
             
             ax.annotate('Deterioration', xy=(maxval-0.03*maxval,0.5),  xycoords='data',
                xytext=(0,0), textcoords='offset points',
                size=13, ha='right', va="center",
                bbox=dict(boxstyle="round", alpha=0.1, color='r'))
             ax.annotate('Improvement',   xy=(midpnt-0.03*midpnt,0.5),  xycoords='data',
                xytext=(0,0), textcoords='offset points',
                size=13, ha='right', va="center",
                bbox=dict(boxstyle="round", alpha=0.1, color='g'))
#               arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.1));
             ax.annotate(obtype.upper(), xy=(maxval-0.06*maxval,amax(pos)-0.5), 
                size=12, ha='left', va="bottom")



    xlim([minval,maxval])
 
    y1 = range(0,len(chnindx))
    ax.margins(y=0)
    if len(chnindx) > 30:
       for label in (ax.get_yticklabels()):
          label.set_fontsize(4)
          label.set_fontweight('bold')
    else:
       for label in (ax.get_yticklabels()):
          label.set_fontsize(12)
          label.set_fontweight('bold')
    vunits = valsYaml['global']['configure'][obtype]['vertical units']
    if vunits == "hPa":
       ax.set_ylabel('Pressure ('+vunits+')',fontsize = 12)
    else:
       ax.set_ylabel('Channel index',fontsize = 12)
    ax.set_yticks(y1)
    ax.set_yticklabels(int32(chnindx[y1]), minor=False, rotation=0)

    units=valsYaml['global']['units']
    mylabel = '(x %s)' % str(1/scale)
    if units != "1":
       if units == "%":
          mylabel = units;
       else:
          mylabel = mylabel + ' (%s' % units + ')'
       ax.set_xlabel(mylabel,fontsize = 12)
#   ax.text(0.75,0.90,obtype.upper(),fontsize=14)

#.................................
def show_bars_monthly_one(lev,vals,colors,cases,scale,this,chisql,chisqr,studt,valsYaml):

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)

    confidence=valsYaml['global']['confidence']
 
    nf=0
    case = cases
    color = colors
    mpl.rc('font', family='serif')
    opacity = 0.5
    bar_width = 0.6
    label_size = 7
    error_config = {'ecolor': '0.3'}
    omean = vals
    pos = arange(len(omean))+nf*bar_width # hack to center the labels on the center of the bar
    if confidence==False:
       ax.barh(pos, omean, bar_width,
               alpha=opacity,
               color=color,
               label=case,
               error_kw=error_config)
    else:
       if case == 'mean':
          conf=studt*vals
          ax.barh(pos, omean, bar_width,
                  alpha=opacity,
                  color=color,
                  label=case,
                  xerr=conf,
                  error_kw=error_config)
       else:
          ax.barh(pos, omean, bar_width,
                  alpha=opacity,
                  color=color,
                  label=case,
                  xerr=[confl,confr],
                  error_kw=error_config)

    legend(loc='lower left')
    title(this)

    obtype=valsYaml['global']['obtype']
    minval = float(valsYaml['global']['configure'][obtype]['min value'])
    maxval = float(valsYaml['global']['configure'][obtype]['max value'])
    if minval == maxval:
       minval = amin(vals)
       maxval = amax(vals)
    xlim([minval,maxval])

    y1 = range(0,len(chnindx))
    ax.margins(y=0)
    for label in (ax.get_yticklabels()):
       if case in ['airs','iasi','cris']:
          label.set_fontsize(4)
       else:
          label.set_fontsize(7)
       label.set_fontweight('bold')
    vunits = valsYaml['global']['configure'][obtype]['vertical units']
    if vunits == "hPa":
       ax.set_ylabel('Pressure ('+vunits+')')
    else:
       ax.set_ylabel('Channel index')
    ax.set_yticks(y1)
    ax.set_yticklabels(int32(chnindx[y1]), minor=False, rotation=0)

    units=valsYaml['global']['units']
    mytitle = obtype+' (x %s)' % str(1/scale)
    if units != "1":
       mytitle = mytitle + ' (%s' % units + ')'
    title(mytitle)

#.................................
def show_bars_monthly(lev,vals,colors,cases,scale,this,chisql,chisqr,studt,valsYaml):

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)

    confidence=valsYaml['global']['confidence']
 
    nf=0
    for case in cases:
       color = colors[nf]
       mpl.rc('font', family='serif')
       opacity = 0.5
       bar_width = 0.6
       label_size = 7
       error_config = {'ecolor': '0.3'}
#      mpl.rcParams['ytick.labelsize'] = label_size
       if size(cases)==1:
          omean = vals
       else:
          omean = vals[:,nf]
       pos = arange(len(omean))+nf*bar_width # hack to center the labels on the center of the bar
       if confidence==False:
          ax.barh(pos, omean, bar_width,
                  alpha=opacity,
                  color=color,
                  label=case,
                  error_kw=error_config)
       else:
          if case == 'mean':
             if cases[1]=='stdv':
                conf=studt*vals[:,1]
             else:
                print "inconsistent options, aborting ..."
                sys.exit(1)
             ax.barh(pos, omean, bar_width,
                     alpha=opacity,
                     color=color,
                     label=case,
                     xerr=conf,
                     error_kw=error_config)
          else:
             ax.barh(pos, omean, bar_width,
                     alpha=opacity,
                     color=color,
                     label=case,
                     xerr=[confl,confr],
                     error_kw=error_config)
  
       legend(loc='lower left')
       title(this)
       nf=nf+1

    obtype=valsYaml['global']['obtype']
    minval = float(valsYaml['global']['configure'][obtype]['min value'])
    maxval = float(valsYaml['global']['configure'][obtype]['max value'])
    if minval == maxval:
       minval = amin(vals)
       maxval = amax(vals)
    xlim([minval,maxval])

    y1 = range(0,len(chnindx))
    ax.margins(y=0)
    for label in (ax.get_yticklabels()):
       if case in ['airs','iasi','cris']:
          label.set_fontsize(4)
       else:
          label.set_fontsize(7)
       label.set_fontweight('bold')
    vunits = valsYaml['global']['configure'][obtype]['vertical units']
    if vunits == "hPa":
       ax.set_ylabel('Pressure ('+vunits+')')
    else:
       ax.set_ylabel('Channel index')
    ax.set_yticks(y1)
    ax.set_yticklabels(int32(chnindx[y1]), minor=False, rotation=0)

    units=valsYaml['global']['units']
    mytitle = obtype+' (x %s)' % str(1/scale)
    if units != "1":
       mytitle = mytitle + ' (%s' % units + ')'
    title(mytitle)


#.................................
def get_start_to_end(start_date, end_date):
    date_list = []
    for i in range(0, (end_date - start_date).days + 1):
        date_list.append(  str(start_date + timedelta(days=i))  ) #<-- here
    return date_list
#.................................
def show_bars_time_series(vname,years,months,scale,total,valsYaml):

  fig = plt.figure(figsize=(10, 10))
  minval = float(valsYaml['global']['configure'][vname]['min value'])
  maxval = float(valsYaml['global']['configure'][vname]['max value'])
  if minval == maxval:
     minval = amin(total)
     maxval = amax(total)
  if minval < 0.0:
     rval = maxval/minval
  else:
     rval = minval/maxval
  if rval < 0.0:
     maxval= max(abs(minval),abs(maxval))
     minval=-maxval
     cs=plt.pcolor(total,cmap=plt.cm.seismic,vmin=minval,vmax=maxval)
  else:
     cs=plt.pcolor(total,cmap=plt.cm.binary,vmin=minval,vmax=maxval)
  cbar=plt.colorbar(cs)
  ax = plt.gca()
  ax.axes.xaxis.set_visible(True)
  ax.axes.yaxis.set_visible(True)

  x1 = range(0,len(years))
  ax.margins(y=0)
  for label in (ax.get_xticklabels()):
     label.set_fontsize(10)
     label.set_fontweight('bold')
  ax.set_xlabel('Date')
  ax.set_xticks(x1)
  yyyymm = 100*years + months
  ax.set_xticklabels(int32(yyyymm[x1]), minor=False, rotation=45)
  
  y1 = range(0,len(chnindx))
  ax.margins(y=0)
  for label in (ax.get_yticklabels()):
     if vname in ['airs','iasi','cris']:
        label.set_fontsize(4)
     else:
        label.set_fontsize(7)
     label.set_fontweight('bold')
  vunits = valsYaml['global']['configure'][vname]['vertical units']
  if vunits == "hPa":
     ax.set_ylabel('Pressure ('+vunits+')')
  else:
     ax.set_ylabel('Channel index')
  ax.set_yticks(y1)
  ax.set_yticklabels(int32(chnindx[y1]), minor=False, rotation=0)

  stat=valsYaml['global']['statistics flavor']
  units=valsYaml['global']['units']
  mytitle = stat+' (x %s)' % str(1/scale)
  if units != "1":
     mytitle = mytitle + ' (%s' % units + ')'
  title(mytitle)

#.................................
def index(var):
  for i in range(len(var)):
      if var[i]==True:
         exit
  return i
#.................................
