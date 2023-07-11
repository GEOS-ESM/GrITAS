#!/usr/bin/env python

import common

# For debugging
#----------------------------------------
def verboseDict(d,indent=''):
    for k,v in d.items():
        pref='%s%s :'%(indent,k)
        if type(v)==dict:
            print(pref)
            verboseDict(v,indent+'  ')
        else:
            print("%s %s"%(pref,v))



# Domain of latitude/longitude - nothing but a glorified tuple
#---------------------------------------------------
class coordRange:
    def __init__(self,_range):
        self._min, self._max=_range

# Some global properties
#--------------------------
class globalProps:
    def __init__(self,name,*args):
        self.name=name
        self.collect=self.yamlStruc(args)
            
    def yamlStruc(self,*args):
        return {t1:t2 for t1,t2 in args}
        

# A single region to be focused on 
#---------------------------------------------------
class Region(globalProps):
    def __init__(self,name=None,lonRange=(0.0,360.0),latRange=(-90.0,90.0)):
        self.name=name
        self.longitude = coordRange(lonRange)
        self.latitude = coordRange(latRange)
        self.toYaml = self.yamlStruc(('lona',self.longitude._min),('lonb',self.longitude._max),
                                     ('lata',self.latitude._min),('latb',self.latitude._max))

    def __repr__(self):
        return 'Region instance "%s":\n\tLongitude : (%.1f,%.1f)\n\tLatitude  : (%.1f,%.1f)\n'%\
            (self.name,self.longitude._min,self.longitude._max,self.latitude._min,self.latitude._max)

    # Serialize an instance of Region
    def serialize(self,yam):
        yam.writeObj({self.name: self.toYaml})

    def read(self,valsYML):
        self.name=vals
        self.longitude = coordRange(tuple())
        self.latitude = coordRange(tuple())
        

# Collect all regions to be focused on 
#---------------------------------------------------
class Collection:
    def __init__(self,name,regions):
        self.name = name
        self.regions = regions

    
    # Serialize an instance of Collection
    def serialize(self,yam):
        toYaml = {r.name:r.toYaml for r in self.regions}
        yam.writeObj({self.name: toYaml})



# Collect configuration info for different instruments
#-------------------------------------------------------
class instrument(globalProps):
    def __init__(self,name,_min,_max,vertUnits='index'):
        self.name=name
        self._min=_min
        self._max=_max
        self.vertUnits=vertUnits
        self.toYaml=self.yamlStruc(('vertical units',str(self.vertUnits)),
                                   ('min value', self._min),
                                   ('max value', self._max))

    # Serialize instance of instrument
    def serialize(self,yam):
        yam.writeObj({self.name: self.toYaml})



class monthlyComparator(globalProps):
    def __init__(self,doit,typ='ratio'):
        self.name='compare monthly'
        self.doit=doit
        self.typ=typ
        self.toYaml=self.yamlStruc(('doit',self.doit), ('type', self.typ))

    def help(self):
        availTypes=['ratio','difference','trivial']
        print("Available monthly comparator types = \n")
        for n,a in enumerate(availTypes):
            print("\t\t%i) %s"%(n+1,a))

    def serialize(self,yam):
        yam.writeObj({'Comparator': self.toYaml})


class Stats(globalProps):
    def __init__(self,flav,confInterval):
        self.flavor=flav
        self.scale, self.units=self.__set__()
        self.measures=['mean','stdv'] # How is this different from flavor??
        self.colors=[c for c in ['b','r','g','k'][:len(self.measures)]]
        self.confidence=confInterval
        self.toYaml = self.yamlStruc(('flavor',     self.flavor),
                                     ('scale',      self.scale),
                                     ('units',      self.units),
                                     ('measures',   self.measures),
                                     ('colors',     self.colors),
                                     ('confidence', self.confidence))

    def __set__(self):
        if self.flavor   == 'Standard Deviation': return '%.1f'%1,'%'
        elif self.flavor == 'DFS per Ob':         return '%.1e'%10000,'1'
        elif self.flavor == 'Ob count':           return '%.1f'%1,'1'
        else:
            raise ValueError("Unsupported Stats Flavor = %s\n- Supported Flavors: %s, %s, %s"
                             %(self.flavor,'Standard Deviation','DFS per Ob','Ob count'))
        
        
    def serialize(self,yam):
        yam.writeObj(self.toYaml)

# Main class for visualizing residuals from GrITAS
#------------------------------------------
class Residual(Stats):
    '''
    Will incorporate the following structures

    <  GLOBAL  >
    start date : str
    end date : str
    nicknames : list
    experiment identifier : list
    file name : str
    ob count treshold for statistics : int
    obtype : str
    regions : list
    figure type : str (png default)

    \/<Stats>
    \/statistics flavor : str
    \/scale : float
    \/units : str
    \/statistics : list (e.g. ['mean','stdv'])
    \/colors : list (e.g. ['b','r','g'] )
    \/confidence : bool

    monthly plot : bool
    time series plot : bool

    <monthlyComparator>
    doit : bool
    type : str (among 'ratio', 'difference', 'trivial')

    <configure>
      <instrument>
      vertical units : str
      min value : float
      max value : float

    <  REGION  >
    (regions)
    lona : float
    lonb : float
    lata : float
    latb : float
    '''

    def __init__(self,instruments,comparator,universe):
        self.name='global'
        self.startDate='YYYY-MM-DD'
        self.endDate='YYYY-MM-DD'
        self.nicknames=['geosfp', 'geosfpp']
        self.expID=['f5294_fp','f5295_fpp']
        self.fileName='XYZ'
        self.obCnt=0
        self.obType='atmsnpp'
        self.regions=['glo']
        self.figType='png'
        self.monthlyPlot=False
        self.tSeriesPlot=False

        self.stats=Stats('Standard Deviation',True)
        self.instruments=instruments
        self.comparatorMonthly=comparator
        self.Universe=universe
        

    def read(self,valsYML):
        self.startDate=valsYML['start date']
        self.endDate=valsYML['end date']
        self.nicknames=valsYML['nicknames']
        self.obType=valsYML['obType']

    def serialize(self,yam,out):
        toYaml = self.yamlStruc(('start date',                       self.startDate),
                                ('end date',                         self.endDate),
                                ('nicknames',                        self.nicknames),
                                ('experiment identifier',            self.expID),
                                ('file name',                        self.fileName),
                                ('ob count treshold for statistics', self.obCnt),
                                ('obtype',                           self.obType),
                                ('regions',                          self.regions),
                                ('figure type',                      self.figType),
                                ('statistics',                       self.stats.toYaml),
                                ('monthly plot',                     self.monthlyPlot),
                                ('time series plot',                 self.tSeriesPlot),
                                ('Comparator',                       self.comparatorMonthly.toYaml),
                                ('configure',                        self.instruments.toYaml))
        
        yam.writeObj({self.name: toYaml}); out.write("\n")
        self.Universe.serialize(yam)




