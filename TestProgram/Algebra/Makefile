TOPDIR = ../..
MKDIR = $(TOPDIR)/Make

EXES = check_algebra_1 
#check_algebra_2 only compiles for SU(3) and REPR=AS
EXES += $(if $(findstring 3,$(NG)),\
		     $(if $(findstring REPR_ANTISYMMETRIC,$(REPR)),\
			      $(if $(findstring GAUGE_SUN,$(GAUGE_GROUP)),\
                        check_algebra_2\
				   ),\
			  ),\
		 )



LIBS += $(TOPDIR)/LibHR/libhr.a


include $(MKDIR)/MkRules

