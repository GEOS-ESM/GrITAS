#!/usr/bin/env python

import sys
from datetime import date,datetime, timedelta
import yaml

import netCDF4 as nc4


class GrITAS:
  def __init__(self,fname='/dev/null',
                    ymlfn='null',
                    vname='null',
                    vunits='null',
                    confidence=False):
    self.nlon=-1;
    self.nlat=-1;
    self.nlev=-1;
    self.fname=fname;
    self.vname=vname;
    self.ymlfn=ymlfn;
    self.vunits=vunits;
    self.confidence=confidence;

    self.var=[]

    self.setup()

#......................................................
  def setup(self): 

    print ("setup soon")
#   valuesYaml = yaml.load(self.ymlfn, Loader=yaml.FullLoader)
#   self.confidence=valuesYaml['global']['confidence']
#   self.vunits=valuesYaml['global']['configure'][self.vname]['vertical units']

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
  def get_info(self):

     print ("soon")
#    return self.units, self.confidence

#......................................................
  def get_gritas_data(self): 

    f = nc4.Dataset(self.fname,'r', format='NETCDF4')

    #  time,lev,lat,lon
    self.var = f.variables[self.vname][:,:,:,:]
    lon = f.variables['lon'][:]
    lat = f.variables['lat'][:]
    lev = f.variables['lev'][:]

    f.close()

    nlev= len(lev)
    if self.vunits == 'hPa':
       idx = range(nlev)
    else:
       idx = range(nlev-1,-1,-1)
       lev = lev[idx]
    mean = self.var[0,idx,:,:]
    stdv = self.var[1,idx,:,:]
    nobs = self.var[2,idx,:,:]

    if self.confidence==True:
       chisqr = self.var[3,idx,:,:]
       chisql = self.var[4,idx,:,:]
       tstud  = self.var[5,idx,:,:]
       return lon,lat,lev,mean,stdv,nobs,chisqr,chisql,tstud
    else: 
       return lon,lat,lev,mean,stdv,nobs

