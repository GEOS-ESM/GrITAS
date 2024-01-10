# GrITAS Plotting
#### Aids in the visualization of time-averaged statistics produced via `gritas.x`.

----------------------------

## Basic Operation
1. `iniPlotYaml.py`
  - form yaml configuring how gridded time-averaged statistics should be visualized
2. `gritStatViz.py`
  - produce select plots of time-averaged statistics

### ```iniPlotYaml.py```
```bash
./iniPlotYaml.py -h

Options:
  -h, --help            show this help message and exit
  -d DATE, --date=DATE  Start/End dates of measurements - delimited via "/"
                        (default = YYYYMMDD/YYYYMMDD)
  -e EXPIDS, --expIDs=EXPIDS
                        ID(s) for experiment(s) to plot gridded stats for -
                        forward slash delimited (default = /)
  -i INSTRUMENT         Specify instrument to consider (default = atmsnpp)
  -r STATSINREGIONS     Geographic regions (delimited via forward slash) over
                        which statistics should be viewed (default = glo)
  -s STATSTOVIEW, --statsToView=STATSTOVIEW
                        Statistics to view (default = mean/stdv)
  --dfs                 Create yaml for DFS
  --impact              Create yaml for observation impact
  --resid               Create yaml for observation residuals
  --compVia=COMPVIA     If two experiments are provided, compare them
                        according to this scheme (default = ratio)
  --tSeriesVar=TSERIESVAR
                        Specify stat to view time series for - ignored if
                        options.T is False (default = '')
  --usrDefRegions=USRDEFREGIONS
                        CSV file specifying user defined lat/lon regions to
                        consider - ignored if empty (default = )
  --C                   Form confidence intervals
  --L                   Plot stats via a line plot; defaults to bars if
                        omitted
  --M                   Form a monthly plot of statistics
  --T                   Form a time series plot
```
| Arguments        | Explanation   | Notes  |
| ---------------- |:-------------:| -----:|
| DATE   | Time window over which to visualize statistics. | Passed at command line as ```YYYYMMDD/YYYYMMDD``` |
| EXPIDS | Select one or two experiments to visualize statistics - denoted via user-defined shorthand `ABRV` passsed to `driver.sh`. | Max of two experiments can be visualized at once. If two experiments are provided, the first is treated as the ***control*** (`cntl`) and the second is deemed ***experiment*** (`exp`). See **Additional Notes** below for further information. |
| INSTRUMENT | Single observing instrument among experiment(s) in `EXP` for which statistics should be viewed. | Defaults to *atmsnpp*, if `-i` is not provided. User is alerted if instrument does not match those that are supported. **N.B.**: plotting range associated with each instrument should be adjusted to user's specifications. |
| STATSINREGIONS | User-selected regions to consider. | Default regions available: `glo, nhe, she, tro, nam` denoting Global, Northern Hemisphere, Southern Hemisphere, Tropics, North America, respectively. Custom regions may be specified, provided they are defined via the `USRDEFREGIONS` argument. Multiple regions are forward-slash delimited. |
| STATSTOVIEW | View these statistics, regardless of type of plot desired. | Forward-slash delimited; defaults to (sample) mean and (sample) standard deviation (i.e. 'mean'/'stdv').
| DFS | Produces input yaml specialized for DFS. | - |
| IMPACT | Produces input yaml specialized for observation impact. | - |
| RESID | Produces input yaml specialized for observation residuals. | - |
| COMPVIA | If two experiments are given in `EXP`, compare them according to this scheme. | Supported schemes: **ratio**, **difference**, **ratio+difference**, or **difference+ratio**. When both ratio and difference are chosen, each line plot will be superimposed onto a single figure. See **On Confidence Intervals** for further details. |
| TSERIESVAR | View time series of this statistic. | Acceptable variables: mean, stdv, sum. Option is ignored if `options.T` is not provided at command line. |
| USRDEFREGIONS | CSV file specifying custom lat/lon regions to consider. | Specify custom regions (one per line) as ```NAME <MIN LAT>.<MAX LAT> <MIN LON>.<MAX LON>```. |
| C | Form confidence intervals. | Includes confidence intervals on sample mean and standard deviation for a single experiment, and when two experiments are compared. See **On Confidence Intervals** for further details. |
| L | Plot statistics using a line plot. | When omitted, statistics are visualized using a horizontal bar plot. |
| M | Form a monthly plot of statistics. | Visualize statistics for one or two experiments per month across time window specifed by `DATE`. Should two experiments be passed, a monthly plot of statistics can only proceed by specifying how the two experiments are to be compared (see `COMPVIA`). |
| T | Form a time series plot. | Visualize the time evolution of an observing system's `TSERIESVAR` statistic across time window specified by `DATE`. Only allowed for a single experiment. |

**Additional Notes:**
  The statistics visualization support detailed below will function *with, maximally, only two experiments specified*.

- Example yamls for each of ***DFS***, ***IMPACT***, ***RESID*** can be found in the `templates` directory.

- By default, it is assumed *mean*, *stdv*, *sum* are each supported in the netCDF files produced from `gritas.x`. If this is not the case, modify the instance of `plot_util.GlobalProps` at *line 107* of `iniPlotYaml.py` to explicitly denote the supported stats:
```python
Global=GlobalProps(...)  ---->  Global=GlobalProps(...,supported_stats=[list of supported stats])
```

# On Confidence Intervals
All confidence intervals in `GrITAS` are assigned at a critical value of $`\alpha=0.05`$ or a $`95\%`$ level of confidence. Toggling `options.C` will
include a $`95\%`$ confidence interval on plotted statistics (e.g. mean, standard deviation).

*Confidence Intervals for a Single Experiment*
----------------------------

Assume a random sample of size $`N`$, with sample mean $`\overline{x}`$ and sample variance $`s^2`$, is drawn from an unknown normal distribution $`\mathcal{N}\left(\mu,\sigma^2\right)`$ with
population mean $`\mu`$ and variance $`\sigma^2`$. In the case of `GrITAS`, this would correspond to satellite observations
within a given frequency channel at some latitude/longitude location.
- Assigning a confidence interval for sample mean $`\overline{x}`$:
 - with an unknown population variance, the random variable $`T\equiv\frac{\overline{x}-\mu}{s/\sqrt{N}}`$ is a pivotal quantity and is distributed according to a *Student's t-distribution* with $`N-1`$ degrees of freedom - i.e., $`T\sim t\left(N-1\right)`$
 - $`\left(1-\alpha\right)100\%`$ confidence interval for $`\mu`$ with an unknown population variance $`\sigma^2`$:

 $$
 \mu=\left[\overline{x}-t(N-1\mid\frac{\alpha}{2})\frac{s}{\sqrt{n}},\overline{x}+t(N-1\mid\frac{\alpha}{2})\frac{s}{\sqrt{n}}\right]
 $$

 where $`t(N-1\mid\frac{\alpha}{2})`$ is the real value for which $`P(T>t(N-1))=\frac{\alpha}{2}=0.025`$.

- Assigning a confidence interval for sample variance $`s^2`$:
 - again with an unknown population variance, the random variable $`Q\equiv\frac{\left(N-1\right)s^2}{\sigma^2}`$ is a pivotal quantity and is distributed according to a *Chi-Square Distribution* with $`N-1`$ degrees of freedom - i.e., $`Q\sim\chi^2\left(N-1\right)`$
 - $`\left(1-\alpha\right)100\%`$ confidence interval for $`\sigma^2`$:

 $$
 \sigma^2=\left[s^2-\frac{\left(N-1\right)s^2}{\chi^2\left(N-1\mid\frac{\alpha}{2}\right)},s^2+\frac{\left(N-1\right)s^2}{\chi^2\left(N-1\mid1-\frac{\alpha}{2}\right)}\right]
 $$

 where $`\chi^2\left(N-1\mid\frac{\alpha}{2}\right)`$ is the real value for which $`P\left(Q>\chi^2(N-1)\right)=\frac{\alpha}{2}=0.025`$,
 and $`\chi^2\left(N-1\mid1-\frac{\alpha}{2}\right)`$ is the real value for which $`P\left(Q>\chi^2(N-1)\right)=1-\frac{\alpha}{2}=0.975`$ - or the right/left Chi-Square critical values.

*Confidence Intervals for Two Experiments*
----------------------------

When comparing statistics between two experiments, two options are supported: `ratio` and `difference`. Recall the first experiment listed in `EXP` is deemed the ***control*** (`cntl`), while the second is deemed ***experiment*** (`exp`).
- `ratio`: computes percentage of `cntl` sample standard deviation that is `exp` sample standard deviation:

$$
100\times\left(\frac{s_{\rm exp}}{s_{\rm cntl}}\right)\pm s_{\rm exp}\cdot\left[\chi^2\left(N-1\mid\frac{\alpha}{2}\right),\chi^2\left(N-1\mid1-\frac{\alpha}{2}\right)\right]
$$

 where the right/left Chi-Square scores are taken from `exp`, and the control and experiment datasets are assumed to have the same dimension $`N`$.

- `difference`: computes the difference between the sample means of `exp` and `cntl`:


$$
x_{\rm exp}-x_{\rm cntl}\pm s_{\rm exp}\cdot t(N-1\mid\frac{\alpha}{2})
$$


 where the critical t-value is taken from `exp`, and the control and experiment datasets are assumed to have the same dimension $`N`$.


### ```gritStatViz.py```
```bash
Usage: gritStatViz.py [options]

Options:
  -h, --help            show this help message and exit
  -d, --dfs             Visualize degrees of freedom for instrument
  -i, --impacts         Visualize observation impacts
  -r, --residuals       Visualize residuals
  -f YMLCONFIG, --ymlConfig=YMLCONFIG
                        Yaml driver (default = "")
```
| Arguments        | Explanation   | Notes  |
| ---------------- |:-------------:| -----:|
| DFS | Specifies `YMLCONFIG` is configured for DFS. | No argument needed. |
| IMPACTS | Specifies `YMLCONFIG` is configured for observation impact. | No argument needed. |
| RESIDUALS | Specifies `YMLCONFIG` is configured for observation residuals. | No argument needed. |
| YMLCONFIG | Yaml produced from `iniPlotYaml.py`. |  |

Execution is restricted to one of `DFS`, `IMPACTS`, or `RESIDUALS`. Otherwise, operation should be straightforward.

----------------------------

To-Do:
- [ ] the ratio and difference schemes detailed above are sufficient to visualize how a control experiment is affected by a change in an observing system. However, some changes to the underlying `GrITAS` Fortran code would be needed to adhere to literature:
  - [ ] assignment of a confidence interval to the difference between two sample means requires *Welch's t-interval*, which loosens the assumption of equal dataset sizes and underlying population variances.
  - [ ] assignment of a confidence interval to the ratio of two sample variances requires critical values of the *F-distribution*
- [ ] address changing plot domains based on different observing systems and type of plot
