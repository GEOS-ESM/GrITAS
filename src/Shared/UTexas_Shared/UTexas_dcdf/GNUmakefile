#
# Makefile for ESMA components.
#
# REVISION HISTORY:
#
# 19Jan2004  Zaslavsky  Initial imlementation.
# 20Oct2004  da Silva   Standardization
#

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(PWD)/../../..
endif

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk  # System dependencies

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------

THIS  := $(shell basename `pwd`)
FREAL = $(FREAL4)

LIB  = lib$(THIS).a

BIN_DOTX_SRC = ut_cdfchi.f

BIN_DOTX = $(addsuffix .x, $(basename $(BIN_DOTX_SRC)))
BINS     = $(BIN_DOTX)

esma_install install: $(LIB) $(BINS)
	$(MKDIR) $(ESMALIB) $(ESMAETC) $(ESMAINC)/$(THIS) $(ESMABIN)
	$(CP) -p *.a    $(ESMALIB)
	$(CP) -p $(BINS)        $(ESMABIN)

esma_clean clean: 
	-$(RM) *~ *.[aox] *.[Mm][Oo][Dd] $(BINS)

esma_distclean distclean: 
	-$(RM) *~ *.[aoxd] *.[Mm][Oo][Dd] $(BINS)

esma_doc doc:
	@echo "Target $@ not implemented yet in `pwd`"

esma_help:
	@echo "Standard ESMA targets:"
	@echo "% make esma_install    (builds and install under ESMADIR)"
	@echo "% make esma_clean      (removes deliverables: *.[aox], etc)"
	@echo "% make esma_distclean  (leaves in the same state as cvs co)"
	@echo "% make esma_doc        (generates PDF, installs under ESMADIR)"
	@echo "% make esma_help       (this message)"
	@echo "Environment:"
	@echo "      ESMADIR = $(ESMADIR)"
	@echo "      BASEDIR = $(BASEDIR)"
	@echo "         ARCH = $(ARCH)"
	@echo "         SITE = $(SITE)"
	@echo "         OBJS= $(OBJS) "

help:
	@ echo "Available targets:"
	@ echo "  make lib             (to create libods.a)"
	@ echo "  make lib MODOBJS=    (to omit f90 modules)"
	@ echo "  make test            (to test ODS library)"
	@ echo "  make doc             (to create documentation)"
	@ echo "  make all"  
	@ echo "  make export"  
	@ echo "  make clean"  
	@ echo "  make distclean"
	@ echo "  make depend"  


MODSRC = 

DCDFSRC = \
algdiv.f  betaln.f  brcomp.f  cdff.f    cdft.f    cumfnc.f  dbetrm.f  dln1mx.f  dzror.f   gam1.f    gratio.f  rlog1.f \
alngam.f  bfrac.f   bup.f     cdffnc.f  cumbet.f  cumgam.f  devlpl.f  dln1px.f  erfc1.f   gaminv.f  gsumln.f  rlog.f \
alnrel.f  bgrat.f   cdfbet.f  cdfgam.f  cumbin.f  cumnbn.f  dexpm1.f  dlnbet.f  erf.f     gamln1.f  ipmpar.f  spmpar.f \
apser.f   bpser.f   cdfbin.f  cdfnbn.f  cumchi.f  cumnor.f  dinvnr.f  dlngam.f  esum.f    gamln.f   psi.f     stvaln.f \
basym.f   bratio.f  cdfchi.f  cdfnor.f  cumchn.f  cumpoi.f  dinvr.f   dstrem.f  exparg.f  gamma.f   rcomp.f \
bcorr.f   brcmp1.f  cdfchn.f  cdfpoi.f  cumf.f    cumt.f    dlanor.f  dt1.f     fpser.f   grat1.f   rexp.f \

DOCSRC = 

#                  --------------------
#                  User Defined Targets
#                  --------------------

SRCS = $(MODSRC) $(DCDFSRC)
OBJS := $(addsuffix .o, $(basename $(SRCS)))

ALL_SRCS = $(SRCS) $(BIN_DOTX_SRC)
DEPS := $(addsuffix .d, $(basename $(ALL_SRCS)))

FOPT = $(FOPT3) -traceback
USER_FFLAGS = $(MPFLAG) $(M).

MOD_DIRS = .

vpath % $(MOD_DIRS)

$(LIB) lib : $(DEPS) $(OBJS)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)

%.x : $(LIB) %.o
	$(LD) $(LDFLAGS) -o $@ $*.o $(LIB) $(LIB_SYS)

#                  --------------------
#                      Dependencies
#                  --------------------

#include Makefile.depend

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros
#.
