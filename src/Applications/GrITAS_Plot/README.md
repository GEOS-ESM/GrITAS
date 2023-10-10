# GrITAS Plotting
###### Aids in the visualization of time averaged statistics produced via `gritas.x`.
---

Operation of 'gridded stats view':


```bash
./iniPlotYaml.py -h

Usage: iniPlotYaml.py [options]

Options:
  -h, --help            show this help message and exit
  -d DATE, --date=DATE  Start/End dates of measurements - delimited via "/"
                        (default = YYYYMMDD/YYYYMMDD)
  -e EXP, --exp=EXP     Yaml containing experiment(s) information
  -r STATSINREGIONS, --statsInRegions=STATSINREGIONS
                        Geographic regions over which statistics should be
                        viewed (default = REG1/REG2/...)
  -s STATSTOVIEW, --statsToView=STATSTOVIEW
                        Statistics to view (default = mean/stdv)
  --dfs                 Create yaml for DFS
  --impact              Create yaml for observation impact
  --resid               Create yaml for observation residuals
  --compVia=COMPVIA     If two experiments are provided, compare them
                        according this scheme (default = ratio)
  --tSeriesVar=TSERIESVAR
                        Specify stat to view time series for - ignored if
                        options.T is False (default = '')
  --usrDefRegions=USRDEFREGIONS
                        CSV file specifying user defined lat/lon regions to
                        consider - ignored if empty (default = )
  --C                   Form confidence intervals
  --M                   Form a monthly plot of statistics
  --T                   Form a time series plot
```
| Arguments        | Explanation   | Notes  |
| ---------------- |:-------------:| -----:|
| DATE   | Time window over which to visualize statistics. | Passed at command line as ```YYYYMMDD/YYYYMMDD``` |
| EXP | Standalone yaml config collecting names of experiment(s) to consider, nicknames, and absolute paths. | Max of two experiments can be visualized at once. If two experiments are provided, the first is treated as the ***control*** (`cntl`) and the second is deemed ***experiment*** (`exp`). See **Additional Notes** below for further information. |
| STATSINREGIONS | User-selected regions to consider. | Default regions available: `glo, nhe, she, tro, nam` denoting global, northern hemisphere, southern hemisphere, tropics, north america, respectively. Additional regions may be specified via `USRDEFREGIONS`; Delimited via forward slash. |
| STATSTOVIEW | View these statistics, regardless of type of plot desired. | Forward-slash delimited; defaults to (sample) mean and (sample) standard deviation (i.e. 'mean'/'stdv').
| DFS | Produces input yaml specialized for DFS. | - |
| IMPACT | Produces input yaml specialized for observation impact. | - |
| RESID | Produces input yaml specialized for observation residuals. | - |
| COMPVIA | If two experiments are given in `EXP`, compare them according to this scheme. | Supported schemes: **ratio** or **difference**. See **On Confidence Intervals** for further details. |
| TSERIESVAR | View time series of this statistic. | Acceptable variables: mean, stdv, sum. Options is ignored if `options.T` is not provided at command line. |
| USRDEFREGIONS | CSV file specifying custom lat/lon regions to consider. | Specify custom regions (one per line) as ```NAME <MIN LAT>.<MAX LAT> <MIN LON>.<MAX LON>```. |
| C | Form confidence intervals. | Includes confidence intervals on sample mean and standard deviation for a single experiment, and when two experiments are compared. See **On Confidence Intervals** for further details. |
| M | Form a monthly plot of statistics. | Visualize statistics for one or two experiments per month, where for two experiments the user should specify how two experiments are to be compared. |
| T | Form a time series plot. | Visualize evolution of an observing system's statistics across an experiment time interval. |

**Additional Notes:**
- `EXP` is a standalone YAML denoting experiment(s) to consider. An example, provided below, is for atmsnpp innovations at 0Z for November 2022:
  ```yaml
  experiments:
    f5294_fp:
      nickname: geosfp
      file name: /discover/nobackup/rtodling/Eval/__ID__/obs/NovDec22/__ID__.atms_npp_0h_gritas_omf.nc4
    f5295_fpp:
      nickname: geosfpp
      file name: /discover/nobackup/rtodling/Eval/__ID__/obs/NovDec22/__ID__.atms_npp_0h_gritas_omf.nc4

  ```
  In general, experiments are specified in the following fashion:
  ```yaml
  experiments:
    <some experiment name>:
      nickname: <nickname for experiment>
      file name: <absolute path to netCDF produced as output from gritas.x>
    .
    .
    .
  ```
  NOTE: each instance of `__ID__` in `file name` will be replaced by the corresponding experiment's `nickname`.

  While more than two such experiments may be specified in this yaml file, the statistics visualization support detailed below will function *with, maximally, only two experiments specified*.

- `options.M`: signifies user-selected statistics should be viewed per month across time window specified via `DATE`.
  - NOTE: allowed for either of single or two experiment cases


- `options.T`: signifies the statistic specified by `TSERIESVAR` should be viewed as a time series across the window specified via `DATE`.
  - NOTE: allowed for only the case of a single experiment


- Example yamls for each of *DFS*, *IMPACT*, *RESID* can be found in the `templates` directory.

- By default, it is assumed *mean*, *stdv*, *sum* are each supported in the netCDF files produced from `gritas.x`. If this is not the case, modify the instance of `plot_util.GlobalProps` at *line 107* to explicitly denote the supported stats:
```python
Global=GlobalProps(...)  ---->  Global=GlobalProps(...,supported_stats=[list of supported stats])
```

# On Confidence Intervals
All confidence intervals in `GrITAS` are assigned at a critical value of $`\alpha=0.05`$ or a $`95\%`$ level of confidence. Toggling `options.C` will
therefore include a $`95\%`$ confidence interval on plotted statistics.

**Confidence Intervals for a Single Experiment**
--

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

**Confidence Intervals for Two Experiments**
--

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


To-Do:
- [ ] the ratio and difference schemes detailed above are sufficient to visualize how a control experiment is affected by a change in an observing system. However, some changes to the underlying `GrITAS` Fortran code would be needed to adhere to literature:
  - [ ] assignment of a confidence interval to the difference between two sample means requires *Welch's t-interval*, which loosens the assumption of equal dataset sizes and underlying population variances.
  - [ ] assignment of a confidence interval to the ratio of two sample variances requires critical values of the *F-distribution*
- [ ] address changing plot domains based on different observing systems and type of plot
