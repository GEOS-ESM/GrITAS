#!/usr/bin/env python

import yaml

# Personalized driver of creating yaml input files
# ...because of pyyaml quirkiness, this would be called for each distinct data structure desired in yaml
#-----------------------------------------------------------------------
class myYML:
    def __init__(self,out,sort_keys=False,default_flow_style=False):
        self.out=out
        self.indent=3
        self.sort_keys=sort_keys
        self.default_flow_style=default_flow_style

    def writeObj(self,obj):
        yaml.dump(obj,stream=self.out,sort_keys=self.sort_keys,\
                  default_flow_style=self.default_flow_style,indent=self.indent)


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
