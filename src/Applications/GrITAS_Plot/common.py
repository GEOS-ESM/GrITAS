#!/usr/bin/env python

import yaml

class Yaml:
    # def serialize(self):
    # def dump(self):
    #     yaml.dump({k:v},stream=OUT,sort_keys=False,default_flow_style=False)
    pass
        

# Personalized driver of creating yaml input files
# ...because of pyyaml quirkiness, this would be called for each distinct data structure desired in yaml
#-----------------------------------------------------------------------
class myYML:
    def __init__(self,out,sort_keys=False,default_flow_style=False):
        self.out=out
        self.sort_keys=sort_keys
        self.default_flow_style=default_flow_style

    def writeTag(self,obj):
        yaml.dump(obj,stream=self.out,sort_keys=self.sort_keys,\
                  default_flow_style=self.default_flow_style)


# Domain of latitude/longitude - nothing but a glorified tuple
#---------------------------------------------------
class coordRange:
    def __init__(self,_range):
        self._min, self._max=_range


# A single region to be focused on 
#---------------------------------------------------
class Region:
    def __init__(self,name,lonRange=(0.0,360.0),latRange=(-90.0,90.0)):
        self.name=name
        self.longitude = coordRange(lonRange)
        self.latitude = coordRange(latRange)
        self.collection = { self.name : { 'lona': self.longitude._min,
                                          'lonb': self.longitude._max,
                                          'lata': self.latitude._min,
                                          'latb': self.latitude._max } }

    # Serialize an instance of Region
    def serialize(self,yam):
        yam.writeTag(self.collection)
        

# Collect all regions to be focused on 
#---------------------------------------------------
class Universe:
    def __init__(self,regions):
        self.name = 'region'
        self.collection = {self.name : {}}

        for c in regions:
            self.collection[self.name].update(c.collection)
    
    # Serialize an instance of Universe
    def serialize(self,yam):
        yam.writeTag(self.collection)
        
