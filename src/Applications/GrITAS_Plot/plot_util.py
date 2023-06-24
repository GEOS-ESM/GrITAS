#!/usr/bin/env python

import common

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
        self.collect=self.yamlStruc(*args)
            
    def yamlStruc(self,*args):
        print(type(args[0]))
        print(args)
        return {self.name : {t1:t2 for t1,t2 in args}}
        # if type(arg) == dict:
        #     yamlStruc(self,*args)
        # else:
        #     return {self.name : {t1:t2 for t1,t2 in args}}
        

# A single region to be focused on 
#---------------------------------------------------
class Region(globalProps):
    def __init__(self,name,lonRange=(0.0,360.0),latRange=(-90.0,90.0)):
        self.name=name
        self.longitude = coordRange(lonRange)
        self.latitude = coordRange(latRange)

    def __repr__(self):
        return 'Region instance "%s":\n\tLongitude : (%.1f,%.1f)\n\tLatitude  : (%.1f,%.1f)\n'%\
            (self.name,self.longitude._min,self.longitude._max,self.latitude._min,self.latitude._max)

    # Serialize an instance of Region
    def serialize(self,yam):
        yam.writeTag(self.yamlStruc(('lona',self.longitude._min),('lonb',self.longitude._max),
                                    ('lata',self.latitude._min),('latb',self.latitude._max)))
        

# Collect all regions to be focused on 
#---------------------------------------------------
class Universe(globalProps):
    def __init__(self,regions):
        self.name = 'region'
        self.regions = regions
        # self.collection = {self.name : {}}

        # for c in regions:
        #     self.collection[self.name].update(c.collection)
    
    # Serialize an instance of Universe
    def serialize(self,yam):
        yam.writeTag(self.yamlStruc(self.regions))



# Collect configuration info for different instruments
#-------------------------------------------------------
class instrument(globalProps):
    def __init__(self,name,_min,_max,vertUnits='index'):
        self.name=name
        self._min=_min
        self._max=_max
        self.vertUnits=vertUnits

    # Serialize instance of instrument
    def serialize(self,yam):
        yam.writeTag(self.yamlStruc(('vertical units',str(self.vertUnits)),
                                    ('min value', self._min),
                                    ('max value', self._max)))



class monthlyComparator(globalProps):
    def __init__(self,doit,typ='ratio'):
        self.name='compare monthly'
        self.doit=doit
        self.typ=typ

    def help(self):
        availTypes=['ratio','difference','trivial']
        print("Available monthly comparator types = \n")
        for n,a in enumerate(availTypes):
            print("\t\t%i) %s"%(n+1,a))

    def serialize(self,yam):
        # collect={'compare monthly': {'doit': self.doit,
        #                              'type': self.typ} }
        yam.writeTag(self.yamlStruc(('doit',self.doit), ('type', self.typ)))


class Stats(globalProps):
    def __init__(self,flav,confInterval):
        self.flavor=flav
        self.scale, self.units=self.__set__()
        self.measures=['mean','stdv'] # How is this different from flavor??
        self.colors=[c for c in ['b','r','g','k'][:len(self.measures)]]
        self.confidence=confInterval
        self.collection = {'flavor': self.flavor,
                           'scale': self.scale,
                           'units': self.units,
                           'measures': self.measures,
                           'colors': self.colors,
                           'confidence': self.confidence}

    def __set__(self):
        if self.flavor   == 'Standard Deviation': return '%.1f'%1,'%'
        elif self.flavor == 'DFS per Ob':         return '%.1e'%10000,'1'
        elif self.flavor == 'Ob count':           return '%.1f'%1,'1'
        else:
            raise ValueError("Unsupported Stats Flavor = %s\n- Supported Flavors: %s, %s, %s"
                             %(self.flavor,'Standard Deviation','DFS per Ob','Ob count'))
        
        
    def serialize(self,yam):
        yam.writeTag(self.collection)

# Main class for visualizing residuals from GrITAS
#------------------------------------------
class Residual(Stats):
    '''
    Will incorporate the following structures

    <  GLOBAL  >
    start date : str
    end date : str
    nickname : list
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

    def __init__(self,instruments,universe,comparator):
        self.name='global'
        self.startDate='YYYY-MM-DD'
        self.endDate='YYYY-MM-DD'
        self.nickname=['geosfp', 'geosfpp']
        self.expID=['f5294_fp','f5295_fpp']
        self.fileName='XYZ'
        self.obCnt=0
        self.obType='atmsnpp'
        self.regions=['glo']
        self.figType='png'
        self.monthlyPlot=False
        self.tSeriesPlot=False

        self.stats=Stats('Standard Deviation',True)
        self.comparatorMonthly=comparator
        self.instruments=instruments
        self.Universe=universe

    def serialize(self,yam):
        collect={self.name: {'start date':                       self.startDate,
                             'end date':                         self.endDate,
                             'nickname':                         self.nickname,
                             'experiment identifier':            self.expID,
                             'file name':                        self.fileName,
                             'ob count treshold for statistics': self.obCnt,
                             'obtype':                           self.obType,
                             'regions':                          self.regions,
                             'figure type':                      self.figType,
                             'statistics':                       self.stats.collection,
                             'configure':                        self.instruments,
                             'monthly plot':                     self.monthlyPlot,
                             'time series plot':                 self.tSeriesPlot,
                             'Comparator':                       self.comparatorMonthly } }
        
        yam.writeTag(collect)
        self.Universe.serialize(yam)

#         self.stats.serialize(yam)
#         self.comparatorMonthly.serialize(yam)
#         self.instruments.serialize(yam) 
#         self.Universe.serialize(yam)
#         # self.stats.serialize(yam)
#         # self.comparatorMonthly.serialize(yam)
# #        self.instruments.serialize(yam)




