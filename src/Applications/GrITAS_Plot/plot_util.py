#!/usr/bin/env python3

import common
from abc import ABC, abstractmethod

class CoordRange:
    '''
    Domain of latitude or longitude - nothing but a glorified tuple

    Attributes
    ----------
    _min : int
       Minimum range of coordinate

    _max : int
       Maximum range of coordinate
    '''
    def __init__(self,_range):
        '''
        Construct a CoordRange instance

        Parameters
        ----------
        _range : tuple
           Float tuple of form (min,max)
        '''
        self._min, self._max=_range

class Serializer(ABC):
    '''
    Manage serialization of all subclasses

    Attributes
    ----------
    name : str
       Key associated with a subclass' properties (held in dictionary)

    toYaml : dict
       Properties of a subclass

    Methods
    -------
    serialize(yam)
       Produces dictionary of form { name: toYaml } and writes to yaml output, if yam is a valid output file, otherwise return to caller.

    fromYaml(self)
      Set class attributes from read-in yaml. Given @abstractmethod decorator, forcing subclasses to implement.
    '''
    def __init__(self,name,*args):
        '''
        Initialize a Serializer instance

        Parameters
        ----------
        name : str
           Key associated with a subclass' properties (held in dictionary)

        *args : tuple(s) as (str,value)
           Tuples collecting subclass' attributes and values
        '''
        self.name=name
        self.toYaml={t1:t2 for t1,t2 in args}

    def serialize(self,yam=None):
        '''
        Serialize class attributes into yaml-formatted output

        Produces dictionary of form { name: toYaml } and writes to yaml output, if yam is a valid output file, otherwise return to caller.

        Parameters
        ----------
        yam : common.YML instance
           Yaml serializer for class attributes

        Returns
        -------
        None if yam is defined, otherwise dict
        '''
        _d={self.name: self.toYaml}
        if yam:
            yam.writeObj(_d)
        else:
            return _d

    @abstractmethod
    def fromYaml(self):
        pass

class Region(Serializer):
    '''
    A geographic region denoted by name and coordinates

    Attributes
    ----------
    longitude : CoordRange instance
       Longitude range of region

    latitude : CoordRange instance
       Latitude range of region

    name : str
       Name given to geographic region

    toYaml : dict
       <Inherited from plot_util.Serializer>

    Methods
    -------
    fromYaml(self,yamlTop,**kwargs)
       Sets Region attributes from yamlTop level in an open yaml input file.
    '''
    def __init__(self,name='globe',lonRange=(0.0,360.0),latRange=(-90.0,90.0)):
        '''
        Initialize a Region instance

        Parameters
        ----------
        name : str
           Name given to geographic region

        lonRange : tuple (float)
           Longitude range of region

        latRange : tuple (float)
           Latitude range of region
        '''
        self.longitude = CoordRange(lonRange)
        self.latitude = CoordRange(latRange)
        Serializer.__init__(self,name,('lona',self.longitude._min),('lonb',self.longitude._max),
                            ('lata',self.latitude._min),('latb',self.latitude._max))

    def __pretty_lat__(self,l):
        '''
        Latitude in degrees North/South

        Returns
        -------
        str
        '''
        if l < 0:   return '%.1fS'%abs(l)
        elif l > 0: return '%.1fN'%abs(l)
        else:       return '0'

    def __pretty_lon__(self,l):
        '''
        Longitude in degrees East/West

        Returns
        -------
        str
        '''
        if l < 0:   return '%.1fW'%abs(l)
        elif l > 0: return '%.1fE'%abs(l)
        else:       return '0'

    def __repr__(self):
        '''
        String representation of Region

        Returns
        -------
        str
        '''
        return '%s\n%s - %s\n%s - %s'%(self.name,self.__pretty_lat__(self.latitude._min),
                                       self.__pretty_lat__(self.latitude._max),
                                       self.__pretty_lon__(self.longitude._min),
                                       self.__pretty_lon__(self.longitude._max))

    def fromYaml(self,yamlTop,**kwargs):
        '''
        Sets Region attributes from yamlTop level in an open yaml input file.

        Parameters
        ----------
        yamlTop : str
           Parse input yaml at level defined by key yamlTop.

        **kwargs

        Returns
        -------
        self
        '''
        self.longitude = CoordRange((yamlTop['lona'],yamlTop['lonb']))
        self.latitude = CoordRange((yamlTop['lata'],yamlTop['latb']))
        Serializer.__init__(self,kwargs['name'],('lona',self.longitude._min),('lonb',self.longitude._max),
                            ('lata',self.latitude._min),('latb',self.latitude._max))
        return self

class Collection:
    '''
    Associate a list of homogeneous class instances to a key

    Attributes
    ----------
    key : str
       Identifier of class instances list

    members : list
       Homogeneous collection of classes

    Methods
    -------
    append(self,m)
       Append a class instance to members list

    fromYaml(self,yamlTop,cls)
       Read class instances from yaml defined by yamlTop.keys()

    serialize(self,yam=None)
       Serialize each class instance in members list

    size(self)
       Number of class instances held in Collection
    '''
    def __init__(self,key,members=[]):
        '''
        Initialize a Collection instance

        Parameters
        ----------
        key : str
           Indentifier of class instances list

        members : list
           Classes to hold in Collection
        '''
        self.key = key
        self.members = members

    def __repr__(self):
        '''
        String representation of Collection

        Returns
        -------
        str
        '''
        s ='Collection instance "%s" with members:\n'%self.key
        try:
            for m in self.members: s += m.__repr__()
        except:
            s += ''
        return s

    def __iter__(self):
        '''
        Return iterator over Collections.members

        Returns
        -------
        List iterator
        '''
        return iter(self.members)

    def __getitem__(self,key):
        '''
        Access class whose name attribute matches key. Class is returned if match is found, otherwise KeyError

        Parameters
        ----------
        key : str
           Name of desired class member

        Returns
        -------
        Class, if class name matches key

        Raises
        ------
        KeyError : class name not found among Collection members
        '''
        for m in self.members:
            if m.name == key:
                return m
        raise KeyError('Member %s not found among Collection members'%m)

    def append(self,m):
        '''
        Append a class instance to members list

        Parameters
        ----------
        m : instance
           Class instance to append to Collection.members

        Returns
        -------
        None
        '''
        self.members.append(m)

    def fromYaml(self,yamlTop,cls):
        '''
        Read class instances from yaml defined by yamlTop.keys()

        Parameters
        ----------
        yamlTop : str
           Parse input yaml at level defined by key yamlTop

        cls : str
           Class instance type used to direct how members are read from yaml

        Returns
        -------
        None
        '''
        # If reading from yaml, clear members contents to flush default class instances
        self.members=[]

        for k in yamlTop.keys():
            if cls.upper()=='INSTRUMENT':
                self.append( Instrument().fromYaml(yamlTop[k],name=k) )
            elif cls.upper()=='EXPERIMENT':
                self.append( Experiment().fromYaml(yamlTop[k],name=k) )
            else:
                self.append( Region().fromYaml(yamlTop[k],name=k) )

    def serialize(self,yam=None):
        '''
        Serialize each class instance in members list

        Parameters
        ----------
        yam : common.YML instance
           Yaml serializer for class attributes

        Returns
        -------
        None if yam is defined, otherwise dict
        '''
        _tmp={}
        for mem in self.members:
            _tmp.update(mem.serialize())
        if yam:
            yam.writeObj({self.key: _tmp})
        else:
            return {self.key: _tmp}

    def size(self):
        '''
        Number of class instances held in Collection

        Returns
        -------
        int
        '''
        return len(self.members)


class Instrument(Serializer):
    '''
    Configuration info for an instrument

    Attributes
    ----------
    _min : float
       Plot range minimum for instrument

    _max : float
       Plot range maximum for instrument

    name : str
       Name given to instrument

    vertUnits : str
       Vertical axis label units

    Methods
    -------
    fromYaml(self,yamlTop,**kwargs)
       Sets Instrument attributes from yamlTop level in an open yaml input file
    '''
    def __init__(self,name='Unspecified',_min=0.0,_max=0.0,vertUnits='index'):
        '''
        Initialize an Instrument instance

        Parameters
        ----------
        name : str
           Name given to instrument

        _min : float
           Plot range minimum for instrument

        _max : float
           Plot range maximum for instrument

        vertUnits : str
           Vertical axis label units (defaults to 'index' - e.g., channel index)
        '''
        self._min=_min
        self._max=_max
        self.vertUnits=vertUnits
        Serializer.__init__(self,name,('vertical units',str(self.vertUnits)),
                            ('min value', self._min),('max value', self._max))

    def __repr__(self):
        '''
        String representation of Instrument

        Returns
        -------
        str
        '''
        return '   %s:\n\t range: (%s,%s)\n\t vert. units: %s'%(self.name,self._min,self._max,self.vertUnits)

    def fromYaml(self,yamlTop,**kwargs):
        '''
        Sets Instrument attributes from yamlTop level in an open yaml input file

        Parameters
        ----------
        yamlTop : str
           Parse input yaml at level defined by key yamlTop.

        **kwargs

        Returns
        -------
        self
        '''
        self._min=yamlTop['min value']
        self._max=yamlTop['max value']
        self.vertUnits=yamlTop['vertical units']
        Serializer.__init__(self,kwargs['name'],('vertical units',str(self.vertUnits)),
                            ('min value', self._min),('max value', self._max))
        return self

class Experiment(Serializer):
    '''
    Experiment information

    Attributes
    ----------
    name : str
       Name given to experiment

    nickname : str
       Shorthand assigned to experiment

    pathToFile : str
       Absolute path of experiment netCDF file

    Methods
    -------
    fromYaml(self,yamlTop,**kwargs)
       Sets Experiment attributes from yamlTop level in an open yaml input file
    '''
    def __init__(self,name='Unspecified',nickname='Unspecified',pathToFile='./'):
        '''
        Initialize an Experiment instance

        Parameters
        ----------
        name : str
           Name given to experiment

        nickname : str
           Shorthand assigned to experiment

        pathToFile : str
           Absolute path of experiment netCDF file
        '''
        self.nickname=nickname
        self.pathToFile=pathToFile
        Serializer.__init__(self,name,('nickname',self.nickname),('file name', self.pathToFile))

    def __repr__(self):
        '''
        String representation of Experiment

        Returns
        -------
        str
        '''
        return 'Experiment instance "%s":\n\t nicknamed = %s\n\t path = %s\n'%(self.name,self.nickname,self.pathToFile)

    def fromYaml(self,yamlTop,**kwargs):
        '''
        Sets Experiment attributes from yamlTop level in an open yaml input file

        Parameters
        ----------
        yamlTop : str
           Parse input yaml at level defined by key yamlTop.

        **kwargs

        Returns
        -------
        self
        '''
        self.nickname=yamlTop['nickname']
        self.pathToFile=yamlTop['file name']
        Serializer.__init__(self,kwargs['name'],('nickname',self.nickname),('file name', self.pathToFile))
        return self

class PlotParams(Serializer):
    '''
    Parameters dictating type of plot to produce, and its properties

    Attributes
    ----------
    regions : list
       List of strings denoting geographic regions within which time-averaged statistics should be plotted

    timeSeries : bool
       Produce a time series plot

    timeSeriesVar : str
       Produce a time series of variable

    monthly : bool
       Produce a monthly plot

    compareVia : str
       When two experiments are considered, compare them according to scheme

    simpleBars : bool
       Embelished confidence intervals

    form : str
       Format of saved figures

    Methods
    -------
    fromYaml(self,yamlTop)
       Sets PlotParams attributes from yamlTop level in an open yaml input file
    '''
    def __init__(self,regions=[],timeSeries=False,timeSeriesVar='',monthly=True,compVia='ratio'):
        '''
        Initialize a PlotParams instance

        Parameters
        ----------
        regions : list
           List of strings denoting geographic regions within which time-average statistics should be plotted

        timeSeries : bool
           Produce a time series plot

        timeSeriesVar : str
           Produce a time series of variable (accepted: mean, stdv, sum)

        monthly : bool
           Produce a monthly plot

        compVia : str
           When two experiments are considered, compare according to scheme (supported: 'ratio' and 'difference')
        '''
        self.regions=regions
        self.timeSeries=timeSeries
        self.timeSeriesVar=timeSeriesVar
        self.monthly=monthly
        self.compareVia=compVia
        self.simpleBars=True
        self.form='png'
        Serializer.__init__(self,'plot params',
                            ('time series', self.timeSeries),
                            ('time series var', self.timeSeriesVar),
                            ('monthly', self.monthly),
                            ('simple bars', self.simpleBars),
                            ('compare via', self.compareVia),
                            ('regions', self.regions),
                            ('format', self.form))
        self.__valid__()

    def __valid__(self):
        '''
        Ensure compareVia scheme is valid/supported

        Returns
        -------
        None

        Raises
        ------
        ValueError : compareVia is neither 'ratio' nor 'difference'
        '''
        valid = ['ratio','difference','difference+ratio','ratio+difference']
        if self.compareVia not in valid:
            raise ValueError("%s is an unsupported manner to compare two experiments - please select from %s"%
                             (self.compareVia,valid))

    def fromYaml(self,yamlTop):
        '''
        Sets PlotParams attributes from yamlTop level in an open yaml input file

        Parameters
        ----------
        yamlTop : str
           Parse input yaml at level defined by key yamlTop.

        Returns
        -------
        None
        '''
        self.timeSeries = yamlTop['time series']
        self.timeSeriesVar = yamlTop['time series var']
        self.monthly = yamlTop['monthly']
        self.simpleBars = yamlTop['simple bars']
        self.compareVia = yamlTop['compare via']
        self.regions = yamlTop['regions']
        self.form = yamlTop['format']
        Serializer.__init__(self,'plot params',
                            ('time series', self.timeSeries),
                            ('time series var', self.timeSeriesVar),
                            ('monthly', self.monthly),
                            ('simple bars', self.simpleBars),
                            ('compare via', self.compareVia),
                            ('regions', self.regions),
                            ('format', self.form))
        self.__valid__()

class Stats(Serializer):
    '''
    Provide control over statistics

    Attributes
    ----------
    flavor : str
       Flavor of statistics

    scale : float
       Multiplicatively rescale statistics by scale

    units : str
       Units of statistics (used to label horizontal axes of plots)

    measures : list
       Time-averaged statistics to access from netCDF

    colors : list
       Colors assigned to measures

    confidence : bool
       Include confidence intervals on plots

    Methods
    -------
    fromYaml(self,yamlTop)
       Sets Stats attributes from yamlTop level in an open yaml input file
    '''
    def __init__(self,flav='Unspecified',measures=['mean', 'stdv'],confInterval=True):
        '''
        Initialize a Stats instance

        Parameters
        ----------
        flav : str
           Flavor of statistics - supported 'Standard Deviation' (residuals), 'DFS per Ob' (DFS), 'Ob count' (obs impact)

        measures : list
           Time-averaged statistics to access from netCDF (among: 'mean', 'stdv', 'sum')

        confInterval : bool
           Include confidence intervals on plots
        '''
        self.flavor=flav
        self.scale, self.units=self.__set__()
        self.measures=measures
        self.colors=[c for c in ['b','r','g','k'][:len(self.measures)]]
        self.confidence=confInterval
        Serializer.__init__(self,'statistics',
                            ('flavor',self.flavor),('scale',self.scale),('units',self.units),
                            ('measures',self.measures),('colors', self.colors),('confidence', self.confidence))

    def __set__(self):
        '''
        Set self.scale and self.units based on self.flavor

        Returns
        -------
        str

        Raises
        ------
        ValueError : Unsupported statistics flavor specified
        '''
        if self.flavor == 'Unspecified': return -1,'-1'
        elif self.flavor == 'Standard Deviation': return 1.0,'%'
        elif self.flavor == 'DFS per Ob':         return 1e4,'1'
        elif self.flavor == 'Ob count':           return 1.0,'1'
        else:
            raise ValueError("Unsupported Stats Flavor = %s\n- Supported Flavors: %s, %s, %s"
                             %(self.flavor,'Standard Deviation','DFS per Ob','Ob count'))

    def fromYaml(self,yamlTop):
        '''
        Sets Stats attributes from yamlTop level in an open yaml input file

        Parameters
        ----------
        yamlTop : str
           Parse input yaml at level defined by key yamlTop

        Returns
        -------
        None
        '''
        self.flavor     = yamlTop['flavor']
        self.scale      = yamlTop['scale']
        self.units      = yamlTop['units']
        self.measures   = yamlTop['measures']
        self.colors     = yamlTop['colors']
        self.confidence = yamlTop['confidence']
        Serializer.__init__(self,'statistics',
                            ('flavor',self.flavor),('scale',self.scale),('units',self.units),
                            ('measures',self.measures),('colors', self.colors),('confidence', self.confidence))

class GlobalProps(Serializer):
    '''
    Global properties of yaml (once produced) that drives plotting; assembled from Collections of instruments and experiment(s), PlotParams, Stats, and several keyword arguments.

    Attributes
    ----------
    name : str
       Name given to instance (defaults to 'global')

    experiments : plot_util.Collection instance
       Experiment(s) targeted in plotting

    instruments : plot_util.Collection instance
       Instrument configuration

    plotParams : plot_util.PlotParams instance
       Properties configuring plotting

    stats : plot_util.Stats instance
       Properties configuring statistics

    startDate : str
       Start date of time-averaged statistics

    endDate : str
       End date of time-averaged statistics

    obCnt : int
       Observation count threshold for presenting statistics

    obType : str
       Instrument to consider in experiment(s)

    supportedStats : list
       Raw time-averaged statistics present in GrITAS-produced netCDF file. Defaults to ['mean', 'stdv', 'sum'], but may be overwritten by user.

    Methods
    -------
    fromYaml(self,yamlTop,**kwargs)
       Sets GlobalProps attributes from yamlTop level in an open yaml input file

    serialize(self,yam)
       Serialize GlobalProps instance

       Overrides Serializer.serialize
    '''
    def __init__(self,instruments=Collection('instruments',[Instrument()]),
                 experiments=Collection('experiments',[Experiment()]),
                 plotParams=PlotParams(),stats=Stats(),**kwargs):
        '''
        Initialize GlobalProps instance

        Parameters
        ----------
        instruments : plot_util.Collection instance
           Instrument configuration

        experiments : plot_util.Collection instance
           Experiment(s) targeted in plotting

        plotParams : plot_util.PlotParams instance
           Properties configuring plotting

        stats : plot_util.Stats instance
           Properties configuring statistics

        **kwargs
           Used to set startDate, endDate, obCnt, obType, supportedStats
        '''
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
        Serializer.__init__(self,self.name,
                            ('start date', self.startDate),
                            ('end date', self.endDate),
                            ('supported stats', self.supportedStats),
                            ('ob count threshold for statistics', self.obCnt),
                            ('obtype', self.obType))

    def fromYaml(self,yamlTop,**kwargs):
        '''
        Sets GlobalProps attributes from yamlTop level in an open yaml input file

        Parameters
        ----------
        yamlTop : str
           Parse input yaml at level defined by key yamlTop.

        **kwargs

        Returns
        -------
        None
        '''
        self.startDate   = yamlTop['start date']
        self.endDate     = yamlTop['end date']
        self.supportedStats = yamlTop['supported stats']
        self.obCnt       = yamlTop['ob count threshold for statistics']
        self.obType      = yamlTop['obtype']
        self.stats.fromYaml(yamlTop['statistics'])
        self.experiments.fromYaml(yamlTop['experiments'],cls='EXPERIMENT')
        self.instruments.fromYaml(yamlTop['instruments'],cls='INSTRUMENT')
        self.plotParams.fromYaml(yamlTop['plot params'])
        Serializer.__init__(self,self.name,
                            ('start date', self.startDate),
                            ('end date', self.endDate),
                            ('supported stats', self.supportedStats),
                            ('ob count threshold for statistics', self.obCnt),
                            ('obtype', self.obType))

    def serialize(self,yam):
        '''
        Serializing each component of GlobalProps instance: experiments, plotParams, stats, instruments, and keyword arguments.

        Overrides Serializer.serialize

        Parameters
        ----------
        yam : common.YML instance
           Yaml serializer for class attributes

        Returns
        -------
        None
        '''
        _tmp=self.toYaml
        _tmp.update(self.experiments.serialize())
        _tmp.update(self.plotParams.serialize())
        _tmp.update(self.stats.serialize())
        _tmp.update(self.instruments.serialize())
        yam.writeObj({self.name: _tmp})

class StatsViewer:
    '''
    Main class for visualizing time-averaged statistics from GrITAS

    Attributes
    ----------
    globalProps : GlobalProps instance
       Properties configuring how time-averaged statistics (via GrITAS) will be viewed

    universe : Collection instance
       Collection of all geographic regions statistics will be visualized within
    '''
    def __init__(self,glob=GlobalProps(),universe=Collection('regions',[Region()])):
        '''
        Construct a StatsViewer instance

        Parameters
        ----------
        glob : GlobalProps instance
           Properties configuring how time-averaged statistics (via GrITAS) will be viewed. Defaults to a generic configuration.

        universe : Collection instance
           Collection of all geographic regions statistics will be visualized within. Defaults to a single region encompassing the entire globe.
        '''
        self.globalProps=glob
        self.universe=universe

    def serialize(self,yam):
        '''
        Serialize StatsViewer attributes into yaml-formatted output

        Parameters
        ----------
        yam : common.YML instance
           Yaml serializer for class attributes

        Returns
        -------
        None
        '''
        self.globalProps.serialize(yam)
        self.universe.serialize(yam)
