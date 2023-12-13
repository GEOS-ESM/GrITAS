# GrITAS Processing
###### The processing of observation data stream (**ods**) files with GrITAS proceeds through three provided scripts.<br>
---
### Initial Step For Successful Operation
GrITAS is known to fail (segfaults) for `cris_fsr_npp`, `iasi_metop-a`, and `iasi_metop-b` instruments when the host stacksize is not sufficiently large. A quick user workaround is to add the follow line to either of the user's `.bashrc` or `.cshrc` file, and subsequently re-source it.
```bash
# Unlimited stack size in ~/.bashrc
ulimit -s unlimited

# Unlimited stack size in ~./cshrc
unlimit stacksize
```

### Description of Scripts
1. ```driver.sh```
2. ```build-yaml.py```
3. ```grdstats.py```

Unless specified otherwise, a verbose explanation of each script and its input arguments is achieved by calling shell (python) scripts without (with "-h") arguments.

### ```driver.sh```
#### Purpose:
Determine all accessible ***diag***-type **ods** files for a user-specified experiment within two specified calendar dates (inclusive), and across four synoptic times - `00Z, 06Z, 12Z, 18Z`
A description of input arguments is accessed via:
```bash
./driver.sh
```
```
Usage: driver.sh <INI YYYYMMDD> <FIN YYYYMMDD> <PREFIX> <ABRV>

   YYYYMMDD : Initial/final dates in form YYYYMMDD (e.g., 20170104 = 01/04/2017 or Jan. 4, 2017)
   PREFIX   : Prefix of absolute path to experiment
   ABRV     : User-provided shorthand name for experiment
```
| Arguments        | Explanation   | Notes  |
| ---------------- |:-------------:| -----:|
| INI YYYYMMDD   | date of initial observations | Specify as standard four-digit year, two-digit month, two-digit day (ex. 20150110 for January 10, 2015) |
| FIN YYYYMMDD   | date of final observations | [See above] |
| PREFIX         | absolute path prefix to ods files for experiment | For example, experiment x0049 on Discover : /discover/nobackup/projects/gmao/dadev/dao_it/archive/x0049/obs |
| ABRV           | user-defined shorthand for experiment        | - |

Upon error-free termination, four space-delimited `list` files will be generated, one for each synoptic time, cataloging all instrument `*ods` files available within the specified date range for the specified experiment. Each list is of the form:
```
<EXPID>.<YYYYMMDD INI>-<YYYYMMDD FIN>.H<SYNOPTIC TIME>.list
```
Example -- Experiment ABRV = GEOSIT between 2017-01-01 & 2017-01-03:
```bash
          GEOSIT.20170101-20170103.H00.list
          GEOSIT.20170101-20170103.H06.list
          GEOSIT.20170101-20170103.H12.list
          GEOSIT.20170101-20170103.H18.list
```
The contents of each produced list file would then be as follows:
```
Instrument Date ExpID
airs_aqua 20170101 GEOSIT
...
```
Being a space-delimited ASCII file with lines of variable character length, such files may be dense to read, especially for large date ranges. However, the format is easily parsed by the python scripts that follow.

<br>
<br>

### ```build-yaml.py```
#### Purpose:
Seamlessly produce an input yaml file to drive execution of `gritas.x`.<br>
A description of input arguments is accessed via:
```bash
./build-yaml.py -h
```
```bash
Usage: build-yaml.py [options]

Options:
  -h, --help            show this help message and exit
  -c CONFIDENCE, --confidence=CONFIDENCE
                        Measurements confidence level (default = 0.95)
  -d DATE, --date=DATE  Start/End dates of measurements - delimited via "/"
                        (default = YYYYMMDD/YYYYMMDD)
  -e PREFIX, --prefix=PREFIX
                        Absolute path prefix to ods files for experiment
                        (default = "./")
  --rcOverride=RCOVERRIDE
                        Override location of rc files from install/etc (deault
                        = )
  -l EXPODSLISTPREFIX, --expOdsListPrefix=EXPODSLISTPREFIX
                        Instrument/Date/ExpID list file prefix (default =
                        ods.list)
  -s ABRV, --abrv=ABRV  User-defined shorthand for experiment (default = )
  -t SYNOPTICTIMES, --synopticTimes=SYNOPTICTIMES
                        Right-justified, two integer synoptic times -
                        delimited via "/" (default = 00/06/12/18)
  --dryRun              Build a yaml config for testing
```

| Arguments        | Explanation   | Notes  |
| ---------------- |:-------------:| -----:|
| CONFIDENCE       | Confidence level of measurements | Defaults to `0.95` |
| DATE             | Forward slash delimited starting & ending dates | Format `YYYYMMDD/YYYYMMDD`, with `[INI DATE]`/`[FIN DATE]` |
| PREFIX           | Absolute path prefix to ods files for experiment | This should match `PREFIX` passed to `driver.sh` |
| RCOVERRIDE       | Override location of rc files from `<INSTALL DIR>/etc` to specified path | By not passing this flag, rc files will be pulled from `<INSTALL DIR>/etc`. |
| EXPODSLISTPREFIX | Prefix of list files produced by `driver.sh` | Example: `GEOSIT.20170101-20170103` for the GEOSIT experiment identifier between 2017-01-01 & 2017-01-03 |
| ABRV             | User-defined shorthand for experiment | - |
| SYNOPTICTIMES    | Forward slash delimited synoptic times to consider | Any combination within {`00`,`06`,`12`,`18`} may be selected; if option is omitted at command line, all four synoptic times will be considered. |
| dryRun           | Debug | - |

Successful execution will produce a single `yaml` file. Operationally,
* `build-yaml.py` will consider the intersection of each experiment's observation list file (e.g. `$EXPODSLISTPREFIX.H<hh>.list` with `<hh>` among each of chosen `SYNOPTICTIMES`)
* this intersection establishes a set of instruments common across each synoptic time and date range
* for each such instrument, `build-yaml.py` checks that a corresponding resource control (`*rc`) exists
* instruments that remain have an `*rc` file and are included in output `yaml` file only if unique (avoids unnecessary repetition)

**Important**: conventional observations (e.g. `*conv*ods`) require, potentially, several distinct resource control (`*rc`) files. Only select `*rc` files are currently supported
```bash
# Supported conventional obs *rc files
#-------------------
convObs=['upconv','upconv2','gps_100lev']
```

###### To-Do: Improvements
- [x] Set `PREFEXPDIR` internally through an external call to `obsSysLoc()` function within `driver.sh`.
- [x] Set `PREFGRITAS` as a global export statement when environment is sourced.
- [x] Report to user combinations of `*ods` and `*rc` files that do not have matches.
- [ ] Allow user to amend `convObs` list

<br>
<br>

### ```grdstats.py```
#### Purpose: <"gridded statistics">
Drive `gritas.x` for an experiment identifier for user-selected synoptic times within a range of (inclusive) dates. <br>**Most importantly**, `build-yaml.py` and `grdstats.py` are envisioned to work in unison so as to abstract away cumbersome creation of command line input to `gritas.x`. Provided correct input is given to `build-yaml.py`, a user need only do the following:
```bash
./grdstats.py <YAML from build-yaml.py>
```
**Under the hood**, `grdstats.py` forms the following chain of arguments to pass to `gritas.x` for each instrument available for the specified experiment identifier and range of dates:
```bash
gritas.x  -nopassive -<RESID> -rc <INSTALL DIR>/etc/gritas_<INSTRUMENT>.rc \
    -conf $CONFIDENCE -o <OUTDIR> <ODS FILES>
```
where
* `RESID` : hardcoded residuals (currently `omf` & `oma`)
* `INSTALL DIR` : location of GrITAS installation
* `INSTRUMENT` : a given instrument under investigation
* `CONFIDENCE` : confidence level of measurements
* `OUTDIR` : output location of `gritas.x` execution
* `ODS FILES` : templated (where possible) `ods` files on which `gritas.x` will act


Generated NetCDF files from call(s) to `gritas.x` are populated at current working directory according to
```bash
<ABRV>/<RESIDUAL>/<YYYYMMDD INI>-<YYYYMMDD FIN>/<INSTRUMENT>_gritas.nc4
```
Example NetCDF files produced from GEOSIT experiment identifier between 2017-01-01 & 2017-01-03:
```bash
./GEOSIT/oma/20170101-20170103/amsua_metop-a_gritas.nc4
./GEOSIT/oma/20170101-20170103/amsua_metop-b_gritas.nc4
./GEOSIT/oma/20170101-20170103/amsua_n15_gritas.nc4
./GEOSIT/oma/20170101-20170103/amsua_n18_gritas.nc4
...
```

###### To-Do: Improvements
- [ ] Allow for user-specified prefix location of NetCDF files, rather than `./` where `grdstats.py` is executed
- [ ] Allow user to specify residuals
