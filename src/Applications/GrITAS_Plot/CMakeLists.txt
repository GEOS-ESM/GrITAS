set(PYSCRIPTS common.py defaults.py figutil.py plot_util.py)
# install(FILES ${PYSCRIPTS} DESTINATION lib/Python)
install(PROGRAMS ${PYSCRIPTS} iniPlotYaml.py gritStatViz.py DESTINATION bin)

file(GLOB yml_files templates/*.yaml)
install (
  FILES ${yml_files}
  DESTINATION etc/templates
  )
