#!/usr/bin/env python

import common

# Domain of latitude/longitude - nothing but a glorified tuple
#---------------------------------------------------
class CoordRange:
    '''
    Attributes
    ----------
    _min : int
       Minimum range of coordinate

    _max : int
       Maximum range of coordinate
    '''
    def __init__(self,_range):
        '''
        Initialize a CoordRange instance

        Parameters
        ----------
        _range : tuple
           Float tuple of form (min,max)
        '''
        self._min, self._max=_range

# Some global properties
#--------------------------
class SharedFuncs:
    '''
    Attributes
    ----------
    collect :

    name : str


    Methods
    -------
    yamlStruc(*args)
    '''
    def __init__(self,name,*args):
        '''
        Initialize a SharedFuncs instance

        Parameters...
        '''
        self.name=name
        self.collect=self.yamlStruc(args)

    def yamlStruc(self,*args):
        '''
        Parameters
        ----------
        *args : 
        '''
        return {t1:t2 for t1,t2 in args}


# A single region to be focused on
#---------------------------------------------------
class Region(SharedFuncs):
    def __init__(self,name='globe',lonRange=(0.0,360.0),latRange=(-90.0,90.0)):
        self.name=name
        self.longitude = CoordRange(lonRange)
        self.latitude = CoordRange(latRange)

    def __pretty_lat__(self,l):
        if l < 0:   return '%.1fS'%abs(l)
        elif l > 0: return '%.1fN'%abs(l)
        else:       return '0'

    def __pretty_lon__(self,l):
        if l < 0:   return '%.1fW'%abs(l)
        elif l > 0: return '%.1fE'%abs(l)
        else:       return '0'

    def __repr__(self):
        return '%s - %s\n%s - %s'%(self.__pretty_lat__(self.latitude._min),self.__pretty_lat__(self.latitude._max),\
                                   self.__pretty_lon__(self.longitude._min),self.__pretty_lon__(self.longitude._max))

    def serialize(self,yam):
        yam.writeObj({self.name: self.toYaml})

    def toYaml(self):
        return self.yamlStruc(('lona',self.longitude._min),('lonb',self.longitude._max),
                              ('lata',self.latitude._min),('latb',self.latitude._max))

    def fromYaml(self,yamlTop,**kwargs):
        self.name=kwargs['name']
        self.longitude = CoordRange((yamlTop['lona'],yamlTop['lonb']))
        self.latitude = CoordRange((yamlTop['lata'],yamlTop['latb']))
        return self


# Map a key to a dict of items
#---------------------------------------------------
class NestedDict:
    '''
    Maps a key to list of class instances

    Serializable to yaml
    '''
    def __init__(self,key,members=[]):
        self.key = key
        self.members = members

    def __repr__(self):
        s ='NestedDict instance "%s" with members:\n'%self.key
        try:
            for m in self.members: s += m.__repr__()
        except:
            s += ''
        return s

    def __getitem__(self,key):
        for m in self.members:
            if m.name == key: return m
        raise KeyError('Region "%s" is not a member of Universe regions'%key)

    def size(self):
        return len(self.members)

    def toYaml(self):
        return {mem.name:mem.toYaml() for mem in self.members}

    def fromYaml(self,yamlTop,cls):
        # If reading from yaml, clear contents of members to flush default class instances
        # del self.members[:]
        self.members=[]

        for k in yamlTop.keys():
            if cls=='INSTRUMENT':
                self.append( Instrument().fromYaml(yamlTop[k],name=k) )
            elif cls=='EXPERIMENT':
                self.append( Experiment().fromYaml(yamlTop[k],name=k) )
            else:
                self.append( Region().fromYaml(yamlTop[k],name=k) )

    # Serialize
    def serialize(self,yam):
        yam.writeObj({self.key: self.toYaml() })

    # Append new members
    def append(self,m):
        self.members.append(m)
        # self.toYaml.update({m.name:m.toYaml})


# Collect configuration info for different instruments
#-------------------------------------------------------
class Instrument(SharedFuncs):
    def __init__(self,name='Unspecified',_min=0.0,_max=0.0,vertUnits='index'):
        self.name=name
        self._min=_min
        self._max=_max
        self.vertUnits=vertUnits

    def __repr__(self):
        return 'Instrument instance "%s":\n\t (min,max) = (%.2f,%.2f)\t vertUnits = %s\n'%(self.name,
                                                                                           self._min,self._max,
                                                                                           self.vertUnits)
    def serialize(self,yam):
        yam.writeObj({self.name: self.toYaml()})

    def toYaml(self):
        return self.yamlStruc(('vertical units',str(self.vertUnits)),
                              ('min value', self._min),
                              ('max value', self._max))

    def fromYaml(self,yamlTop,**kwargs):
        self.name=kwargs['name']
        self._min=yamlTop['min value']
        self._max=yamlTop['max value']
        self.vertUnits=yamlTop['vertical units']
        return self

# An experiment class [NB: Experiment & Instrument classes are functionally the same! Consider a factory?!]
class Experiment(SharedFuncs):
    def __init__(self,name='Unspecified',nickname='Unspecified',pathToFile='./'):
        self.name=name
        self.nickname=nickname
        self.pathToFile=pathToFile

    def __repr__(self):
        return 'Experiment instance "%s":\n\t nicknamed = %s\n\t path = %s\n'%(self.name,self.nickname,self.pathToFile)

    def serialize(self,yam):
        yam.writeObj({self.name: self.toYaml()})

    def toYaml(self):
        return self.yamlStruc(('nickname',self.nickname),
                              ('file name', self.pathToFile))

    def fromYaml(self,yamlTop,**kwargs):
        self.name=kwargs['name']
        self.nickname=yamlTop['nickname']
        self.pathToFile=yamlTop['file name']
        return self

class PlotParams(SharedFuncs):
    def __init__(self,regions=[],timeSeries=False,timeSeriesVar='',monthly=True,compVia='ratio'):
        self.regions=regions
        self.timeSeries=timeSeries
        self.timeSeriesVar=timeSeriesVar
        self.monthly=monthly
        self.compareVia=compVia
        self.simpleBars=True
        self.form='png'
        self.__valid__()

    def __valid__(self):
        '''
        Ensure a valid scheme exists to compare two experiments
        '''
        availTypes=['ratio','difference']
        if self.compareVia not in availTypes:
            raise ValueError("%s is an unsupported manner to compare two experiments - please select from %s"%(self.compareVia,availTypes))

    def fromYaml(self,yamlTop):
        self.timeSeries = yamlTop['time series']
        self.timeSeriesVar = yamlTop['time series var']
        self.monthly = yamlTop['monthly']
        self.simpleBars = yamlTop['simple bars']
        self.compareVia = yamlTop['compare via']
        self.regions = yamlTop['regions']
        self.form = yamlTop['format']

    def toYaml(self):
        return self.yamlStruc(('time series', self.timeSeries),
                              ('time series var', self.timeSeriesVar),
                              ('monthly', self.monthly),
                              ('simple bars', self.simpleBars),
                              ('compare via', self.compareVia),
                              ('regions', self.regions),
                              ('format', self.form))

    def serialize(self,yam):
        yam.writeObj(self.toYaml)


class Stats(SharedFuncs):
    def __init__(self,flav='Unspecified',measures=['mean', 'stdv'],confInterval=True):
        self.flavor=flav
        self.scale, self.units=self.__set__()
        self.measures=measures
        self.colors=[c for c in ['b','r','g','k'][:len(self.measures)]]
        self.confidence=confInterval

    def __set__(self):
        if self.flavor   == 'Unspecified': return -1,-1
        if self.flavor   == 'Standard Deviation': return 1.0,'%'
        elif self.flavor == 'DFS per Ob':         return 1e4,'1'
        elif self.flavor == 'Ob count':           return 1.0,'1'
        else:
            raise ValueError("Unsupported Stats Flavor = %s\n- Supported Flavors: %s, %s, %s"
                             %(self.flavor,'Standard Deviation','DFS per Ob','Ob count'))

    def toYaml(self):
        return self.yamlStruc(('flavor',     self.flavor),
                              ('scale',      self.scale),
                              ('units',      self.units),
                              ('measures',   self.measures),
                              ('colors',     self.colors),
                              ('confidence', self.confidence))

    def fromYaml(self,yamlTop):
        self.flavor     = yamlTop['flavor']
        self.scale      = yamlTop['scale']
        self.units      = yamlTop['units']
        self.measures   = yamlTop['measures']
        self.colors     = yamlTop['colors']
        self.confidence = yamlTop['confidence']

    def serialize(self,yam):
        yam.writeObj(self.toYaml)

# Global properties
#------------------------------------------
class GlobalProps(SharedFuncs):
    '''
    Will incorporate the following structures

    <  GLOBAL  >
    start date : str
    end date : str

    supported stats: list

    nicknames : list
    experiment identifier : list
    file name : str
    ob count threshold for statistics : int
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

    <experiments>
      <experiment>
        nickname: str
        file name: str

    <Instruments>
      <Instrument>
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

    def __init__(self,instruments=NestedDict('instruments',[Instrument()]),
                 experiments=NestedDict('experiments',[Experiment()]),
                 plotParams=PlotParams(),stats=Stats(),**kwargs):
        self.name='global'
        self.experiments=experiments
        self.instruments=instruments
        self.plotParams=plotParams
        self.stats=stats
        self.startDate=kwargs.get('startDate')
        self.endDate=kwargs.get('endDate')
        self.obCnt=kwargs.get('obCnt')
        self.obType=kwargs.get('obType')
        self.supportedStats=kwargs.get('supported_stats') if kwargs.get('supported_stats') else ['mean', 'stdv', 'sum']
        self.__parseWildCards__()

    def __parseWildCards__(self):
        try:
            self.fileName=self.fileName.replace("$STRDATE",self.startDate.replace('-',''))
            self.fileName=self.fileName.replace("$ENDDATE",self.endDate.replace('-',''))
            self.fileName=[self.fileName.replace("$TMPNAME",n).replace("$EXPID",e)\
                           for n in self.nicknames for e in self.expID]
        except:
            pass



    def fromYaml(self,yamlTop,**kwargs):
        self.startDate   = yamlTop['start date']
        self.endDate     = yamlTop['end date']

        self.supportedStats = yamlTop['supported stats']

        self.obCnt       = yamlTop['ob count threshold for statistics']
        self.obType      = yamlTop['obtype']


        self.stats.fromYaml(yamlTop['statistics'])
        self.experiments.fromYaml(yamlTop['experiments'],cls='EXPERIMENT')
        self.instruments.fromYaml(yamlTop['instruments'],cls='INSTRUMENT')
        self.plotParams.fromYaml(yamlTop['plot params'])

        self.__parseWildCards__()


    def serialize(self,yam):
        toYaml = self.yamlStruc(('start date',                       self.startDate),
                                ('end date',                         self.endDate),
                                ('supported stats',                  self.supportedStats),
                                ('ob count threshold for statistics', self.obCnt),
                                ('obtype',                           self.obType))

        toYaml.update({'experiments': self.experiments.toYaml()})
        toYaml.update({'plot params': self.plotParams.toYaml()})
        toYaml.update({'statistics': self.stats.toYaml()})
        toYaml.update({'instruments': self.instruments.toYaml()})

        yam.writeObj({self.name: toYaml})


# Main class for visualizing time averaged statistics from GrITAS
#-----------------------------------------------------------------
class StatsViewer:
    def __init__(self,glob=GlobalProps(),universe=NestedDict('regions',[Region()])):
        self.globalProps=glob
        self.universe=universe

    def serialize(self,yam,out):
        self.globalProps.serialize(yam)
        self.universe.serialize(yam)




