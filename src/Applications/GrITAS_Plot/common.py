#!/usr/bin/env python

import yaml

class myYML:
    '''
    Class to parse data structures into a yaml format

    Attributes
    ----------
    default_flow_style : bool
       Organize yaml according to yaml.dump default flow style

    indent : int
       Indentation for producing legible yamls via yaml.dump

    out : file
       Name of file in which formatted output will be rendered

    sort_keys : bool
       Reorder key/values when producing yaml

    Methods
    -------
    writeObj(obj)
       Dump dictionary obj to file using yaml formatting via yaml.dump
    '''
    def __init__(self,out,sort_keys=False,default_flow_style=False):
        self.default_flow_style=default_flow_style
        self.indent=3
        self.out=out
        self.sort_keys=sort_keys

    def writeObj(self,obj):
        '''
        Dump dictionary obj to file using yaml formatting via yaml.dump

        Parameters
        ----------
        obj : dict
           Dictionary dumped to file using yaml formatting
        '''
        yaml.dump(obj,stream=self.out,sort_keys=self.sort_keys,\
                  default_flow_style=self.default_flow_style,indent=self.indent)

def verboseDict(d,indent=''):
    '''
    Recursively print the contents of a dictionary

    Attributes
    ----------
    d : str -or- dictionary
       Key or dictionary to print

    indent : str
       Indendation used to print the values associated with key
    '''
    for k,v in d.items():
        pref='%s%s :'%(indent,k)
        if type(v)==dict:
            print(pref)
            verboseDict(v,indent+'  ')
        else:
            print("%s %s"%(pref,v))
