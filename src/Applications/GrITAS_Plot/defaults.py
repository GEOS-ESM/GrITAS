#!/usr/bin/env python3

import warnings
from optparse import OptionValueError
from plot_util import Instrument

# Default instruments supported and domains for plotting
# ------------------------------------------------------
instruments={ 'amsuan15': Instrument('amsuan15',_min=90.0,_max=110.0),
              'amsuan19': Instrument('amsuan19',_min=-1.5,_max=1.5),
              'atmsnpp': Instrument('atmsnpp',_min=-0.25,_max=1.5),
              'atms': Instrument('atms',_min=90.0,_max=110.0),
              'iasi': Instrument('iasi',_min=90.0,_max=110.0),
              'cris': Instrument('cris',_min=-1.5,_max=1.5) }

# Default regions supported
# -------------------------
regions={'glo': {'lon': (-180.0,180.0), 'lat': (-90.0,90.0)},
         'nhe': {'lon': (0.0,360.0), 'lat': (20.0,90.0)},
         'she': {'lon': (0.0,360.0), 'lat': (-90.0,-20.0)},
         'tro': {'lon': (0.0,360.0), 'lat': (-20.0,20.0)},
         'nam': {'lon': (-172.0,-52.0), 'lat': (16.0,72.0)} }


def avail_regions(option, opt_str, value, parser):
    '''
    Callback function to determine if a provided region specifier is present in default regions dict

    Parameters
    ----------
    option : optparse.Option instance
       Option instance calling this callback function

    opt_str : str
       Option string at command line that triggers this callback

    parser : optparse.OptionParser instance
       Drives parsing of command line options

    value : str
       Argument to option seen on command line

    Returns
    -------
    None

    Raises
    ------
    OptionValueError : region specifier not found among supported regions dictionary
    '''
    # All stats in regions need to match in order to proceed
    valid=True
    for r in value.split('/'):
        if r not in regions:
            valid=False
            break

    # If all region specifiers are matched, set the OptionParser attribute
    if valid:
        setattr(parser.values, option.dest, value)
    else:
        _avail_="Encountered an invalid region specifier - Available default regions: "
        for r in regions.keys():
            _avail_+='%s '%r
        raise OptionValueError(_avail_)

def avail_instruments(option, opt_str, value, parser):
    '''
    Callback function to determine if a provided instrument is present in supported instruments dict

    Parameters
    ----------
    option : optparse.Option instance
       Option instance calling this callback function

    opt_str : str
       Option string at command line that triggers this callback

    parser : optparse.OptionParser instance
       Drives parsing of command line options

    value : str
       Argument to option seen on command line

    Returns
    -------
    None

    Raises
    ------
    Warning : instrument not found among supported instruments dictionary - appended with default domain of (90.0,110.0)
    '''
    if value in instruments:
        setattr(parser.values, option.dest, value)
    else:
        _avail_="\t\tAvailable instruments: "
        for k in instruments.keys():
            _avail_+='%s '%k

        bold='\033[1m'+value+'\033[0m'
        dom='\033[1m'+'(90.0,110.0)'+'\033[0m'
        warnings.warn("\n\tInstrument %s not present in default instruments dictionary!\n%s"%(bold,_avail_) +
                      "\n\tAdding instrument %s to instruments dictionary with default domain = %s"%(bold,dom))
        instruments.update({value: Instrument(value,_min=90.0,_max=110.0)})
        setattr(parser.values, option.dest, value)
