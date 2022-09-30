esma_set_this()

#set (SRCS
#  )

# This is equivalent to FOPT=$(FOPT3) in GNU Make
if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
   set (CMAKE_Fortran_FLAGS_RELEASE "${FOPT3} ${ALIGNCOM}")
endif ()

esma_add_library(${this}
  SRCS ${SRCS}
  DEPENDENCIES GMAO_ods GMAO_mpeu)


ecbuild_add_executable(TARGET grisas.x SOURCES grisas.f LIBS ${this} GMAO_ods GMAO_mpeu) 
ecbuild_add_executable(TARGET gritas.x SOURCES gritas.f LIBS ${this} GMAO_ods GMAO_mpeu)

#set(pythonscripts modis_l2a.py mxd04_l2a.py avhrr_l2a.py patmosx_l2a.py aod_data.py)
#install(PROGRAMS ${pythonscripts} get_aero_obs.csh run_gaas_ana.csh DESTINATION bin)

#install(
#   FILES mxd04_nnr.py avhrr_nnr.py
#   DESTINATION lib/Python
#   )

if (EXTENDED_SOURCE)
  set_target_properties (${this}  PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
  set_target_properties (grisas.x PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
  set_target_properties (gritas.x PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
endif()

#file(GLOB rc_files *.rc)
#file(GLOB pcf_files *.pcf)
#file(GLOB tmpl_files *.tmpl)

#install (
#   FILES ${rc_files} ${pcf_files} ${tmpl_files}
#   DESTINATION etc
#   )
