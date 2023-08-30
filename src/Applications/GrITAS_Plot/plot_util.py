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
class sharedFuncs:
    def __init__(self,name,*args):
        self.name=name
        self.collect=self.yamlStruc(args)

    def yamlStruc(self,*args):
        return {t1:t2 for t1,t2 in args}


# A single region to be focused on
#---------------------------------------------------
class Region(sharedFuncs):
    def __init__(self,name='globe',lonRange=(0.0,360.0),latRange=(-90.0,90.0)):
        self.name=name
        self.longitude = coordRange(lonRange)
        self.latitude = coordRange(latRange)

    def __repr__(self):
        return 'Region instance "%s":\n\tLongitude : (%.1f,%.1f)\n\tLatitude  : (%.1f,%.1f)\n'%\
            (self.name,self.longitude._min,self.longitude._max,self.latitude._min,self.latitude._max)

    # Serialize an instance of Region
    def serialize(self,yam):
        yam.writeObj({self.name: self.toYaml})

    def toYaml(self):
        return self.yamlStruc(('lona',self.longitude._min),('lonb',self.longitude._max),
                              ('lata',self.latitude._min),('latb',self.latitude._max))

    def fromYaml(self,yamlTop,**kwargs):
        self.name=kwargs['name']
        self.longitude = coordRange((yamlTop['lona'],yamlTop['lonb']))
        self.latitude = coordRange((yamlTop['lata'],yamlTop['latb']))
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
                self.append( instrument().fromYaml(yamlTop[k],name=k) )
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
class instrument(sharedFuncs):
    def __init__(self,name='Unspecified',_min=0,_max=0,vertUnits='index'):
        self.name=name
        self._min=_min
        self._max=_max
        self.vertUnits=vertUnits

    def __repr__(self):
        return 'instrument instance "%s":\n\t (min,max) = (%i,%i)\t vertUnits = %s\n'%(self.name,
                                                                                       self._min,self._max,
                                                                                       self.vertUnits)
    # Serialize instance of instrument
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

# An experiment class [NB: Experiment & instrument classes are functionally the same! Consider a factory?!]
class Experiment:
    def __init__(self,name='Unspecified',nickname='Unspecified',pathToFile='./'):
        self.name=name
        self.nickname=nickname
        self.pathToFile=pathToFile

    def __repr__(self):
        return 'Experiment instance "%s":\n\t nicknamed = %s\n\t path = %s\n'%(self.name,self.nickname,self.pathToFile)

    # Serialize instance of Experiment
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



class comparator(sharedFuncs):
    def __init__(self,monthly,typ='ratio'):
        self.name='compare monthly'
        self.monthly=monthly
        self.typ=typ

    def help(self):
        availTypes=['ratio','difference','trivial']
        print("Available monthly comparator types = \n")
        for n,a in enumerate(availTypes):
            print("\t\t%i) %s"%(n+1,a))

    def toYaml(self):
        return self.yamlStruc(('monthly',self.monthly), ('type', self.typ))

    def fromYaml(self,yamlTop):
        self.monthly=yamlTop['monthly']
        self.typ=yamlTop['type']

    def serialize(self,yam):
        yam.writeObj({'Comparator': self.toYaml})


class PlotParams(sharedFuncs):
    def __init__(self):
        self.timeSeries=False
        self.monthly=False
        self.typ='Unspecified'
        self.regions=[]
        self.form='png'

    def fromYaml(self,yamlTop):
        self.timeSeries = yamlTop['time series']
        self.monthly = yamlTop['monthly']
        self.typ = yamlTop['type']
        self.regions = yamlTop['regions']
        self.form = yamlTop['format']

    def toYaml(self):
        return self.yamlStruc(('time series', self.timeSeries),
                              ('monthly', self.monthly),
                              ('type', self.typ),
                              ('regions', self.regions),
                              ('format', self.form))

    def serialize(self,yam):
        yam.writeObj(self.toYaml)


class Stats(sharedFuncs):
    def __init__(self,flav='Unspecified',confInterval=True):
        self.flavor=flav
        self.scale, self.units=self.__set__()
        self.measures=['mean','stdv'] # How is this different from flavor??
        self.colors=[c for c in ['b','r','g','k'][:len(self.measures)]]
        self.confidence=confInterval

    def __set__(self):
        if self.flavor   == 'Unspecified': return -1,-1
        if self.flavor   == 'Standard Deviation': return '%.1f'%1,'%'
        elif self.flavor == 'DFS per Ob':         return '%.1e'%10000,'1'
        elif self.flavor == 'Ob count':           return '%.1f'%1,'1'
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
class globalProps(sharedFuncs):
    '''
    Will incorporate the following structures

    <  GLOBAL  >
    start date : str
    end date : str
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

    <comparator>
    monthly : bool
    type : str (among 'ratio', 'difference', 'trivial')

    <experiments>
      <experiment>
        nickname: str
        file name: str

    <instruments>
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

    def __init__(self,instruments=NestedDict('instruments',[instrument()]),
                 experiments=NestedDict('experiments',[Experiment()]),
                 comparator=comparator(True),**kwargs):
        self.name='global'
        self.startDate=kwargs.get('startDate')
        self.endDate=kwargs.get('endDate')

        self.obCnt=kwargs.get('obCnt')
        self.obType=kwargs.get('obType')
        # self.regions=kwargs.get('regions')
        # self.figType='png'
        # self.monthlyPlot=False
        # self.tSeriesPlot=False

        self.plotParams=PlotParams()
        self.stats=Stats() #'Standard Deviation',True)
        self.experiments=experiments
        self.instruments=instruments
        self.comparator=comparator

        self.parseWildCards()

    def parseWildCards(self):
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
        self.obCnt       = yamlTop['ob count threshold for statistics']
        self.obType      = yamlTop['obtype']
        # self.regions     = yamlTop['regions']
        # self.figType     = yamlTop['figure type']
        # self.monthlyPlot = yamlTop['monthly plot']
        # self.tSeriesPlot = yamlTop['time series plot']

        print(yamlTop.keys())
        self.stats.fromYaml(yamlTop['statistics'])
        self.experiments.fromYaml(yamlTop['experiments'],cls='EXPERIMENT')
        self.instruments.fromYaml(yamlTop['instruments'],cls='INSTRUMENT') #instrument())
        # self.comparator.fromYaml(yamlTop['Comparator'])
        self.plotParams.fromYaml(yamlTop['plot params'])

        self.parseWildCards()


    def serialize(self,yam):
        toYaml = self.yamlStruc(('start date',                       self.startDate),
                                ('end date',                         self.endDate),
                                ('nicknames',                        self.nicknames),
                                ('experiment identifier',            self.expID),
                                ('file name',                        self.fileName),
                                ('ob count threshold for statistics', self.obCnt),
                                ('obtype',                           self.obType),
                                ('regions',                          self.regions),
                                ('figure type',                      self.figType),
                                ('monthly plot',                     self.monthlyPlot),
                                ('time series plot',                 self.tSeriesPlot))
        toYaml.update({'plot params': self.plotParams.toYaml()})
        toYaml.update({'statistics': self.stats.toYaml()})
        # toYaml.update({'Comparator': self.comparator.toYaml()})
        toYaml.update({'instruments': self.instruments.toYaml()})

        yam.writeObj({self.name: toYaml})


# Main class for visualizing residuals from GrITAS
#------------------------------------------
class Residual:
    def __init__(self,glob=globalProps(),universe=NestedDict('regions',[Region()])):
        self.globalProps=glob
        self.universe=universe

    def serialize(self,yam,out):
        self.globalProps.serialize(yam)
        out.write("\n")
        self.universe.serialize(yam)




