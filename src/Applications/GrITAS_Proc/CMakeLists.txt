# set(PYSCRIPTS build-yaml.py grdstats.py)
# install(PROGRAMS ${PYSCRIPTS} driver.sh DESTINATION bin)
# install(FILES grit_util.py DESTINATION lib/Python)

# set(ENV{PYTHONPATH} "${CMAKE_CURRENT_SOURCE_DIR}/lib/Python")
# export(PYTHONPATH "${CMAKE_INSTALL_PREFIX}/lib/Python")
# message("Setting PYTHONPATH = ${CMAKE_INSTALL_PREFIX}/lib/Python")

file(READ build-yaml.py.in FILE_CONTENTS)
string(REPLACE "CMAKE_INSTALL_PREFIX" "${CMAKE_INSTALL_PREFIX}" FILE_CONTENTS "${FILE_CONTENTS}")
file(WRITE build-yaml.py "${FILE_CONTENTS}")

file(READ grit_util.py.in FILE_CONTENTS)
string(REPLACE "CMAKE_INSTALL_PREFIX" "${CMAKE_INSTALL_PREFIX}" FILE_CONTENTS "${FILE_CONTENTS}")
file(WRITE grit_util.py "${FILE_CONTENTS}")

set(PYSCRIPTS build-yaml.py grdstats.py grit_util.py)
install(PROGRAMS ${PYSCRIPTS} driver.sh DESTINATION bin)
