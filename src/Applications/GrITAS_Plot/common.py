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

    def writeTag(self,obj):
        yaml.dump(obj,stream=self.out,sort_keys=self.sort_keys,\
                  default_flow_style=self.default_flow_style,indent=self.indent)
        # yaml.write_line_break(yaml.Dumper)
        # yaml.dump(u'\n',stream=self.out)
