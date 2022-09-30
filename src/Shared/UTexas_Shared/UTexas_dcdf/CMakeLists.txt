esma_set_this()

string (REPLACE " " ";" flags ${FREAL4})
add_compile_options(${flags})
set (CMAKE_Fortran_FLAGS_RELEASE "-O3 ${common_Fortran_flags} ${ALIGNCOM} ${MISMATCH}")


set (srcs
     algdiv.f
     alngam.f
     alnrel.f
     apser.f
     basym.f
     bcorr.f
     betaln.f
     bfrac.f
     bgrat.f
     bpser.f
     bratio.f
     brcmp1.f
     brcomp.f
     bup.f
     cdfbet.f
     cdfbin.f
     cdfchi.f
     cdfchn.f
     cdff.f
     cdffnc.f
     cdfgam.f
     cdfnbn.f
     cdfnor.f
     cdfpoi.f
     cdft.f
     cumbet.f
     cumbin.f
     cumchi.f
     cumchn.f
     cumf.f
     cumfnc.f
     cumgam.f
     cumnbn.f
     cumnor.f
     cumpoi.f
     cumt.f
     dbetrm.f
     devlpl.f
     dexpm1.f
     dinvnr.f
     dinvr.f
     dlanor.f
     dln1mx.f
     dln1px.f
     dlnbet.f
     dlngam.f
     dstrem.f
     dt1.f
     dzror.f
     erfc1.f
     erf.f
     esum.f
     exparg.f
     fpser.f
     gam1.f
     gaminv.f
     gamln1.f
     gamln.f
     gamma.f
     grat1.f
     gratio.f
     gsumln.f
     ipmpar.f
     psi.f
     rcomp.f
     rexp.f
     rlog1.f
     rlog.f
     spmpar.f
     stvaln.f
  )

target_compile_options(${this} PRIVATE ${flags})

ecbuild_add_executable (TARGET ut_cdfchi.x      SOURCES ut_cdfchi.f        LIBS ${this})

