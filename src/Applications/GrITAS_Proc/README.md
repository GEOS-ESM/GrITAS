# GrITAS Processing
###### The processing of observation data stream (**ods**) files with GrITAS proceeds through three provided scripts.<br>
---
### For Successful Operation
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

Unless specified otherwise, a verbose explanation of each script and its input arguments is achieved by calling ```driver.sh``` without any arguments, and passing **-h** to each ```*.py``` script. The functionality is described below.

### ```driver.sh```
#### Purpose:
Determine all accessible ***diag***-type **ods** files for a user-specified experiment identifier within two specified calendar dates (inclusive), and across four synoptic times - `00Z, 06Z, 12Z, 18Z`

| Arguments        | Explanation   | Notes  |
| ---------------- |:-------------:| -----:|
| YR INI   | year of initial observations | Specify as standard four-digit year (ex. 2015) |
| MNTH INI | month of initial observations| Specify as either one- or two-digit month (ex. 3 or 03 for March) |
| DAY INI  | day of initial observations  | Specify as either one- or two-digit day (ex. 5 or 05 for 5th day of month)|
| YR FIN   | year of final observations   | [See above] |
| MNTH FIN | month of final observations  | [See above] |
| DAY FIN  | day of final observations    | [See above] |
| EXPID    | experiment identifier        | Currently supported: GEOSIT, GEOSFP, MERRA2

Upon error-free termination, four space-delimited `list` files will be generated, one for each synoptic time, cataloging all instrument `*ods` files available within the specified date range. Each list is of the form:
```
<EXPID>.<YYYYMMDD INI>-<YYYYMMDD FIN>.H<SYNOPTIC TIME>.list
```
Example -- GEOSIT experiment identifier between 2017-01-01 & 2017-01-03:
```bash
          GEOSIT.20170101-20170103.H00.list
          GEOSIT.20170101-20170103.H06.list
          GEOSIT.20170101-20170103.H12.list
          GEOSIT.20170101-20170103.H18.list
```
The contents of each produced list file is as follows:
```
Instrument Date ExpID
airs_aqua 20170101 d5294_geosit_jan18
...
```
Being a space-delimited ASCII file with lines of variable character length, such files may be dense to read, especially for large date ranges. However, the format is easily parsed by the python scripts that follow.

__Important__: the precise location of observations corresponding to each supported experiment identifier is hardcoded within the local bash function
```bash
obsSysLoc()
```
Additional supported sources can be added at ```line 12```.

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
  -e PREFEXPDIR, --prefExpDir=PREFEXPDIR
                        Prefix of directories containing measurements (default
                        = "./")
  -g PREFGRITAS, --prefGrITAS=PREFGRITAS
                        Prefix of GrITAS src (default = ./)
  -l EXPOBSLISTPREFIX, --expObsListPrefix=EXPOBSLISTPREFIX
                        Instrument/Date/ExpID list file prefix (default = foo)
  -s EXPID, --exp=EXPID
                        Experiment identifier (default = )
  -t SYNOPTICTIMES, --synopticTimes=SYNOPTICTIMES
                        Right-justified, two integer synoptic times -
                        delimited via "/" (default = 00/06/12/18)
  --dryRun              Build a yaml config for testing
```

| Arguments        | Explanation   | Notes  |
| ---------------- |:-------------:| -----:|
| CONFIDENCE       | Confidence level of measurements (*need to check this one*) | Defaults to `0.95` |
| DATE             | Forward slash delimited starting & ending dates | Format `YYYYMMDD/YYYYMMDD`, with `[INI DATE]`/`[FIN DATE]` |
| PREFEXPDIR       | Directory prefix wherein observations (ods files) are located | User must currently turn to `obsSysLoc()` function of `driver.sh` for precise locations associated with supported experiment identifiers.|
| PREFGRITAS       | Location of GrITAS source code   | User specified |
| EXPOBSLISTPREFIX | Prefix of list files produced by `driver.sh`  | Example: `GEOSIT.20170101-20170103` for the GEOSIT experiment identifier between 2017-01-01 & 2017-01-03 |
| EXPID            | experiment identifier | Currently supported: GEOSIT, GEOSFP, MERRA2 |
| SYNOPTICTIMES    | Forward slash delimited synoptic times to consider | Any combination within {`00`,`06`,`12`,`18`} may be selected; if option is omitted at command line, all four synoptic times will be considered (templated within yaml generated).
| dryRun           | Debug | No additional argument needed.

Successful execution will produce a single `yaml` file. Operationally,
* `build-yaml.py` will consider the intersection of each experiment's observation list file (e.g. `$EXPOBSLISTPREFIX.H<hh>.list` with `<hh>` among each of chosen `SYNOPTICTIMES`)
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
- [ ] Set `PREFEXPDIR` internally through an external call to `obsSysLoc()` function within `driver.sh`.
- [ ] Set `PREFGRITAS` as a global export statement when environment is sourced.
- [ ] Report to user combinations of `*ods` and `*rc` files that do not have matches.
- [ ] Allow user to amend `convObs` list

<br>
<br>

### ```grdstats.py```
#### Purpose: <"gridded statistics">
Drive `gritas.x` for an experiment identifier for user-selected synoptic times within a range of (inclusive) dates. <br>**Most importantly**, `build-yaml.py` and `grdstats.py` are envisioned to work in unison so as to abstract away cumbersome creation of command line input to `gritas.x`. Provided correct input is given to `build-yaml.py`, a user need only do the following:
```bash
./grdstats.py foo.gritas.yml
```
**Under the hood**, `grdstats.py` forms the following chain of arguments to pass to `gritas.x` for each instrument available for the specified experiment identifier and range of dates:
```bash
gritas.x  -nopassive -<RESID> -rc $PREFGRITAS/Components/gritas/etc/gritas_<INSTRUMENT>.rc \
    -conf $CONFIDENCE -o <OUTDIR> <ODS FILES>
```
where
* `RESID` : hardcoded residuals (currently `omf` & `oma`)
* `PREFGRITAS` : location of GrITAS source
* `INSTRUMENT` : a given instrument under investigation
* `CONFIDENCE` : confidence level of measurements
* `OUTDIR` : output location of `gritas.x` execution
* `ODS FILES` : templated (where possible) `ods` files on which `gritas.x` will act


Generated NetCDF files from call(s) to `gritas.x` are populated at current working directory according to
```bash
<EXPID>/<RESIDUAL>/<YYYYMMDD INI>-<YYYYMMDD FIN>/<INSTRUMENT>_gritas.nc4
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
