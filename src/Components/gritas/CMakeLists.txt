esma_set_this()

string (REPLACE " " ";" flags ${FREAL8})
add_compile_options(${flags})
set (CMAKE_Fortran_FLAGS_RELEASE "-O3 ${common_Fortran_flags} ${ALIGNCOM} ${MISMATCH}")


set (srcs
     etc.f
     generate_raob_speed.f90
     gritas_gfio.f
     gritas_grads.f
     gritas.h
     gritas_kxkt.f
     m_gritas_binning.F
     m_gritas_grids.f
     m_gritas_masks.F90
     rdel.f
     stdio.h
  )

target_compile_options(${this} PRIVATE ${flags})

ecbuild_add_executable (TARGET ut_cdfchi.x      SOURCES ut_cdfchi.f        LIBS ${this})

