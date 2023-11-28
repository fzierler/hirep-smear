#REPRESENTATION = REPR_FUNDAMENTAL
#REPRESENTATION = REPR_ADJOINT
#REPRESENTATION = REPR_ANTISYMMETRIC
#REPRESENTATION = REPR_SYMMETRIC

TOPDIR = ../../..
MKDIR = $(TOPDIR)/Make
AUTOSUNDIR = $(TOPDIR)/Make/Utils/autosun


#SRCS = 

spnalgtest: $(MKDIR)/MkFlags spn_vs_sun_algtest.cpp $(AUTOSUNDIR)/adjoint.h \
	$(AUTOSUNDIR)/antisymmetric.h $(AUTOSUNDIR)/complex.h \
	$(AUTOSUNDIR)/fundamental.h $(AUTOSUNDIR)/list.h $(AUTOSUNDIR)/matrix.h \
	$(AUTOSUNDIR)/polynomial.h $(AUTOSUNDIR)/print.h \
	$(AUTOSUNDIR)/representation.h $(AUTOSUNDIR)/sparse.h \
	$(AUTOSUNDIR)/style.h $(AUTOSUNDIR)/sun.h $(AUTOSUNDIR)/symmetric.h
	$(CXX) -I $(AUTOSUNDIR) -o spnalgtest spn_vs_sun_algtest.cpp -g -O0 ${MACRO} &&\
	echo " done." || ( echo " failed." && false )

include $(MKDIR)/MkRules

MACRO = -D_${REPR}_ -D_${GAUGE_GROUP}_ -D_PICA_STYLE_ -DNDEBUG

clean::
	rm -f writeREPR
