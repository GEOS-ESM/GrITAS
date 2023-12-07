# set(PYSCRIPTS build-yaml.py grdstats.py)
# install(PROGRAMS ${PYSCRIPTS} driver.sh DESTINATION bin)
# install(FILES grit_util.py DESTINATION lib/Python)

# set(ENV{PYTHONPATH} "${CMAKE_CURRENT_SOURCE_DIR}/lib/Python")
# export(PYTHONPATH "${CMAKE_INSTALL_PREFIX}/lib/Python")
# message("Setting PYTHONPATH = ${CMAKE_INSTALL_PREFIX}/lib/Python")


set(PYSCRIPTS build-yaml.py grdstats.py grit_util.py)
install(PROGRAMS ${PYSCRIPTS} driver.sh DESTINATION bin)
