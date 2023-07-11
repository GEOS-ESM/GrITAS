#!/usr/bin/env python

import sys
import pylab
from pylab import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd
from datetime import date,datetime, timedelta
import yaml

import netCDF4 as nc4


class GrITAS:
  def __init__(self,fname='/dev/null',
                    vname='null',
                    vunits='hPa',
                    confidence=-1.0):
    self.nlon=-1;
    self.nlat=-1;
    self.nlev=-1;
    self.fname=fname;
    self.vname=vname;
    self.vunits=vunits;
    self.confidence=confidence;

#......................................................
  def get_dims(self): 

    f = nc4.Dataset(self.fname,'r', format='NETCDF4')

    lon = f.variables['lon'][:]
    lat = f.variables['lat'][:]
    lev = f.variables['lev'][:]

    self.nlon=len(lon)
    self.nlat=len(lat)
    self.nlev=len(lev)
  
    f.close()
    return self.nlev,self.nlat,self.nlon

#......................................................
  def get_gritas_data(self): 

    f = nc4.Dataset(self.fname,'r', format='NETCDF4')

    #  time,lev,lat,lon
    var = f.variables[self.vname][:,:,:,:]
    lon = f.variables['lon'][:]
    lat = f.variables['lat'][:]
    lev = f.variables['lev'][:]

    f.close()

    if self.vunits != 'hPa':
       nlev= len(lev)
       idx = range(nlev-1,-1,-1)
       lev = lev[idx]
    mean = var[0,idx,:,:]
    stdv = var[1,idx,:,:]
    nobs = var[2,idx,:,:]

    if self.confidence==True:
       chisqr = var[3,idx,:,:]
       chisql = var[4,idx,:,:]
       tstud  = var[5,idx,:,:]
       return lon,lat,lev,mean,stdv,nobs,chisqr,chisql,tstud
    else: 
       return lon,lat,lev,mean,stdv,nobs

#.................................................................
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
#.................................
def vaccum(stype,tresh,var,nob):
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

  stat=valuesYaml['global']['statistics flavor']
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
def subsetreg(reglat, reglon, lat,lon):

   indxlat = []
   for i in range(0, len(lat)) :
     if lat[i] >=reglat[0]:
        if lat[i] <= reglat[1]:
           indxlat.append(i)
   indxlon = []
   for i in range(0, len(lon)) :
     if lon[i] >=reglon[0]:
        if lon[i] <= reglon[1]:
           indxlon.append(i)

   return indxlat, indxlon
#.................................
def show_plot_monthly_one(lev,vals,anno,colors,case,scale,this,chisql,chisqr,studt):

    fig = plt.figure(figsize=(8, 10))
    ax = fig.add_subplot(111)
#   ax.set_facecolor('#CECECE')
    ax.set_facecolor('#E9ECEF')

    confidence=valuesYaml['global']['confidence']
 
    nf=0
    color = colors
    mpl.rc('font', family='serif')
    opacity = 0.5
    bar_width = 0.6
    label_size = 8
    error_config = {'ecolor': '0.3'}

    obtype=valuesYaml['global']['obtype']
    minval = float(valuesYaml['global']['configure'][obtype]['min value'])
    maxval = float(valuesYaml['global']['configure'][obtype]['max value'])
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
    vunits = valuesYaml['global']['configure'][obtype]['vertical units']
    if vunits == "hPa":
       ax.set_ylabel('Pressure ('+vunits+')',fontsize = 12)
    else:
       ax.set_ylabel('Channel index',fontsize = 12)
    ax.set_yticks(y1)
    ax.set_yticklabels(int32(chnindx[y1]), minor=False, rotation=0)

    units=valuesYaml['global']['units']
    mylabel = '(x %s)' % str(1/scale)
    if units != "1":
       if units == "%":
          mylabel = units;
       else:
          mylabel = mylabel + ' (%s' % units + ')'
       ax.set_xlabel(mylabel,fontsize = 12)
#   ax.text(0.75,0.90,obtype.upper(),fontsize=14)

#.................................
def show_bars_monthly_one(lev,vals,colors,cases,scale,this,chisql,chisqr,studt):

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)

    confidence=valuesYaml['global']['confidence']
 
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

    obtype=valuesYaml['global']['obtype']
    minval = float(valuesYaml['global']['configure'][obtype]['min value'])
    maxval = float(valuesYaml['global']['configure'][obtype]['max value'])
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
    vunits = valuesYaml['global']['configure'][obtype]['vertical units']
    if vunits == "hPa":
       ax.set_ylabel('Pressure ('+vunits+')')
    else:
       ax.set_ylabel('Channel index')
    ax.set_yticks(y1)
    ax.set_yticklabels(int32(chnindx[y1]), minor=False, rotation=0)

    units=valuesYaml['global']['units']
    mytitle = obtype+' (x %s)' % str(1/scale)
    if units != "1":
       mytitle = mytitle + ' (%s' % units + ')'
    title(mytitle)

#.................................
def show_bars_monthly(lev,vals,colors,cases,scale,this,chisql,chisqr,studt):

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)

    confidence=valuesYaml['global']['confidence']
 
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

    obtype=valuesYaml['global']['obtype']
    minval = float(valuesYaml['global']['configure'][obtype]['min value'])
    maxval = float(valuesYaml['global']['configure'][obtype]['max value'])
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
    vunits = valuesYaml['global']['configure'][obtype]['vertical units']
    if vunits == "hPa":
       ax.set_ylabel('Pressure ('+vunits+')')
    else:
       ax.set_ylabel('Channel index')
    ax.set_yticks(y1)
    ax.set_yticklabels(int32(chnindx[y1]), minor=False, rotation=0)

    units=valuesYaml['global']['units']
    mytitle = obtype+' (x %s)' % str(1/scale)
    if units != "1":
       mytitle = mytitle + ' (%s' % units + ')'
    title(mytitle)

#.................................
def reg_def(name):
   reglat = zeros(2)
   reglon = zeros(2)
   reglat[0] = float(valuesYaml['region'][name]['lata'])
   reglat[1] = float(valuesYaml['region'][name]['latb'])
   reglon[0] = float(valuesYaml['region'][name]['lona'])
   reglon[1] = float(valuesYaml['region'][name]['lonb'])
   return reglat,reglon
#.................................
def get_start_to_end(start_date, end_date):
    date_list = []
    for i in range(0, (end_date - start_date).days + 1):
        date_list.append(  str(start_date + timedelta(days=i))  ) #<-- here
    return date_list
#.................................
def show_bars_time_series(vname,years,months,scale,total):

  fig = plt.figure(figsize=(10, 10))
  minval = float(valuesYaml['global']['configure'][vname]['min value'])
  maxval = float(valuesYaml['global']['configure'][vname]['max value'])
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
  vunits = valuesYaml['global']['configure'][vname]['vertical units']
  if vunits == "hPa":
     ax.set_ylabel('Pressure ('+vunits+')')
  else:
     ax.set_ylabel('Channel index')
  ax.set_yticks(y1)
  ax.set_yticklabels(int32(chnindx[y1]), minor=False, rotation=0)

  stat=valuesYaml['global']['statistics flavor']
  units=valuesYaml['global']['units']
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

myrc = "impacts.yml";
myrc = "dfs.yml";
myrc = "residuals.yml";
with open(myrc, 'r') as f:
    valuesYaml = yaml.load(f, Loader=yaml.FullLoader)

# Read in relevant global configuration parameters
nicknames = valuesYaml['global']['nickname']
vname = valuesYaml['global']['obtype']
vunits = valuesYaml['global']['configure'][vname]['vertical units']
confidence=valuesYaml['global']['confidence']
scale = float(valuesYaml['global']['scale'])
regions = valuesYaml['global']['regions']
strdate = valuesYaml['global']['start date']
enddate = valuesYaml['global']['end date']
fntmpl= valuesYaml['global']['file name']
ext=valuesYaml['global']['figure type']
stats=valuesYaml['global']['statistics']
colors=valuesYaml['global']['colors']
tresh = valuesYaml['global']['ob count treshold for statistics']
expsid = valuesYaml['global']['experiment identier']


YYYYMMDD = pd.date_range(start=strdate, end=enddate, freq='MS')
years = YYYYMMDD.year
months = YYYYMMDD.month

files = ([date.strftime(fntmpl) for date in YYYYMMDD])

fn = files[0].replace("$expid",expsid[0])
g=GrITAS(fn,vname,vunits,confidence)
[nlev,nlat,nlon] = g.get_dims()
nm = size(years)   # number of times to process
nr = size(regions) # number of regions
ns = size(stats)   # number of statistics to calculate
ne = size(expsid)  # number of experiments to handle

total = zeros([nlev,nm,nr,ns,ne])
if confidence==True:
  confl = zeros([nlev,nm,nr,ne])
  confr = zeros([nlev,nm,nr,ne])
  studt = zeros([nlev,nm,nr,ne])
  cl = zeros([nlev,nm,nr])
  cr = zeros([nlev,nm,nr])
if ne==2:
   totals = zeros([nlev,nm,nr,2])

for ne in range(len(expsid)):
  for ii in range(len(files)):

    expid=expsid[ne];
    fname=files[ii];
    fn = fname.replace("$expid",expid); print "this case now: ", fn
    g = GrITAS(fn,vname,vunits,confidence)

    yyyy = years[ii]
    mm = months[ii]
    mm =  str(mm).rjust(2, '0')

    if confidence==True:
       [lon,lat,lev,mean,stdv,nob,chisqr,chisql,tstud] = g.get_gritas_data()
    else:
       [lon,lat,lev,mean,stdv,nob] = get_gritas_data(fn,vname)
    chnindx=lev;
    for nr in range(len(regions)):
      region=regions[nr]
      [reglat,reglon] = reg_def(region);
      [ilat,ilon] = subsetreg(reglat,reglon, lat,lon)

      for ns in range(len(stats)):
         stat = stats[ns]
         if stat == 'sum' or stat == 'mean':
            this = vaccum(stat,tresh,mean[:,ilat,:],nob[:,ilat,:])
            total[:,ii,nr,ns,ne] = scale*this
         if stat == 'stdv':
            this = vaccum(stat,tresh,stdv[:,ilat,:],nob[:,ilat,:])
            total[:,ii,nr,ns,ne] = scale*this
         if confidence==True:
	    [confl[:,ii,nr,ne],confr[:,ii,nr,ne],studt[:,ii,nr,ne]] = getconf(tresh,nob[:,ilat,:],
                                                                                 chisql[:,ilat,:],
                                                                                 chisqr[:,ilat,:],
                                                                                  tstud[:,ilat,:])
# when only a single experiment is treated ...
if len(expsid)==0:
   if len(stats)==0:
      if valuesYaml['global']['monthly plot']==True:
         yyyymm = str(yyyy)+mm
         show_bars_monthly(lev,total[:,ii,nr,0],['b'],[vname],scale,yyyymm)
         pylab.savefig(nicknames[0]+'_'+'monthly_'+vname+'_'+region+'_'+yyyymm+'.'+ext)
   else:
      if valuesYaml['global']['monthly plot']==True:
         yyyymm = str(yyyy)+mm
         show_bars_monthly(lev,total[:,ii,nr,:],colors,stats,scale,yyyymm,chisql,chisqr,studt)
         pylab.savefig(nicknames[0]+'_'+'monthly_'+vname+'_'+region+'_'+yyyymm+'.'+ext)

else:  # when two experiments are treated (compares them) ...
   for nr in range(len(regions)):
    region = regions[nr]
    for ii in range(len(files)):
     for nx in range(len(expsid)): 
       if nx==0:  # store result
          ii=index([stats=='stdv'])
          totals[:,:,:,0] = total[:,:,:,ii,0]
       if nx==1:  # store result
          yyyymm = str(yyyy)+mm
          if valuesYaml['global']['compare monthly']['doit']==True:
             if   valuesYaml['global']['compare monthly']['type']=='ratio':
                totals[:,:,:,0] = 100.0*(total[:,:,:,ii,1]/totals[:,:,:,ii])
                stat=stats[0]
                cl=confl[:,ii,nr,0] # need care
                cr=confr[:,ii,nr,0] # need care
                cl=cl*totals
                cr=cr*totals
                show_plot_monthly_one(lev,totals[:,ii,nr,0],expsid,colors[0],'ratio',scale,yyyymm,cl,cr,studt)
             elif valuesYaml['global']['compare monthly']['type']=='difference':
                totals[:,:,:,0] = total[:,:,:,ii,1] - totals[:,:,:,0]
                stat=stats[0]
                show_bars_monthly_one(lev,totals[:,ii,nr,0],colors[0],'difference',scale,yyyymm,chisql,chisqr,studt)
             else:
                totals[:,:,:,1] = total[:,:,:,ii,nx]
                show_bars_monthly(lev,totals[:,ii,nr,:],colors,expsid,scale,yyyymm,chisql,chisqr,studt)
             pylab.savefig(nicknames[1]+'X'+nicknames[0]+'_'+'monthly_'+vname+'_'+region+'_'+yyyymm+'.'+ext)

if len(expsid)==0 and valuesYaml['global']['time series plot']==True:
   for nr in range(len(regions)):
       for ns in range(len(stats)):
           show_bars_time_series(vname,years,months,scale,total[:,:,nr,ns])
           pylab.savefig(nicknames[0]+'_'+'tseries_'+vname+'_'+region+'_'+stat+'.'+ext)
