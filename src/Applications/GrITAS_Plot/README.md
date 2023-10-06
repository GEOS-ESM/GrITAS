# GrITAS Plotting
###### Aids in the visualization of time averaged statistics produced via `gritas.x`.<br>
---

Operation of 'gridded stats view':
- [x] produce yaml for either of residuals, dfs (unsupported), or observation impacts (unsupported).
- [ ] 


```python
./plot-driver-yaml.py

Usage: resid-yaml.py [options]

Options:
  -h, --help            show this help message and exit
  -d DATE, --date=DATE  Start/End dates of measurements - delimited via "/"
                        (default = YYYYMMDD/YYYYMMDD)
  -e EXP, --exp=EXP     Yaml containing experiment(s) information
  -r STATSINREGIONS, --statsInRegions=STATSINREGIONS
                        Access statistics for these regions (default =
                        REG1/REG2/...)
  --usrDefRegions=USRDEFREGIONS
                        CSV file specifying user defined lat/lon regions to
                        consider - ignored if empty (default = )
  --dfs                 Create yaml for DFS
  --impact              Create yaml for observation impact
  --resid               Create yaml for observation residuals
  --T                   Form a time series plot
  --M                   Form a monthly plot of statistics
```
| Arguments        | Explanation   | Notes  |
| ---------------- |:-------------:| -----:|
| DATE   | Time window over which to visualize statistics | Passed at command line as ```YYYYMMDD/YYYYMMDD``` |
| EXP | Standalone yaml config collecting names of experiment(s) to consider, nicknames, and absolute paths | Max of two experiments can be visualized at once |
| STATSINREGIONS | User-selected regions to consider; delimited via forward slash | Default regions available: `glo, nhe, she, tro, nam` denoting global, northern hemisphere, southern hemisphere, tropics, north america, respectively. Additional regions may be specified via `USRDEFREGIONS`. |
| USRDEFREGIONS | CSV file specifying custom lat/lon regions to consider | Specify custom regions (one per line) as ```NAME <MIN LAT>.<MAX LAT> <MIN LON>.<MAX LON>``` |
| DFS | Produces input yaml specialized for DFS | - |
| IMPACT | Produces input yaml specialized for observation impcact | - |
| RESID | Produces input yaml specialized for observation residuals | - |
| T | Form a time series plot | Visualize evolution of an observing system's statistics across an experiment time interval. |
| M | Form a monthly plot of statistics | Visualize statistics for one or two experiments per month, where for two experiments the user should specify how two experiments are to be compared |

Example yamls for each of dfs, impact, resid can be found in `templates`.

To-Do:
- [ ]: address changing plot domains based on different observing systems and type of plot
- [ ]: split time series function call to allow user to toggle whether sum, mean, stdv are shown across time interval; showing all three is acceptable for a single month/year, but becomes busy for larger date ranges


### Running the plotting routine...
- [ ]: describe passing of yaml produced by ```python resid-yaml.py``` to `something.py`
- [ ]: rename `something.py` to a more informative name


