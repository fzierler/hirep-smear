#REPRESENTATION = REPR_FUNDAMENTAL
#REPRESENTATION = REPR_ADJOINT
#REPRESENTATION = REPR_ANTISYMMETRIC
#REPRESENTATION = REPR_SYMMETRIC

TOPDIR = ../../..
MKDIR = $(TOPDIR)/Make
AUTOSUNDIR = $(TOPDIR)/Make/Utils/autosun

MACRO = -D_${REPR}_ -D_${GAUGE_GROUP}_ -D_PICA_STYLE_ -DNDEBUG

all: spnalgtest write_spn_sun_algconv

spnalgtest: spn_vs_sun_algtest.cpp lib.h $(MKDIR)/MkFlags \
	$(AUTOSUNDIR)/adjoint.h \
	$(AUTOSUNDIR)/antisymmetric.h $(AUTOSUNDIR)/complex.h \
	$(AUTOSUNDIR)/fundamental.h $(AUTOSUNDIR)/list.h $(AUTOSUNDIR)/matrix.h \
	$(AUTOSUNDIR)/polynomial.h $(AUTOSUNDIR)/print.h \
	$(AUTOSUNDIR)/representation.h $(AUTOSUNDIR)/sparse.h \
	$(AUTOSUNDIR)/style.h $(AUTOSUNDIR)/sun.h $(AUTOSUNDIR)/symmetric.h
	$(CXX) -I $(AUTOSUNDIR) -o spnalgtest spn_vs_sun_algtest.cpp -g -O0 ${MACRO} &&\
	echo "spnalgtest done." || ( echo " failed." && false )

write_spn_sun_algconv: write_spn_sun_algconv.cpp lib.h $(MKDIR)/MkFlags \
	$(AUTOSUNDIR)/adjoint.h \
	$(AUTOSUNDIR)/antisymmetric.h $(AUTOSUNDIR)/complex.h \
	$(AUTOSUNDIR)/fundamental.h $(AUTOSUNDIR)/list.h $(AUTOSUNDIR)/matrix.h \
	$(AUTOSUNDIR)/polynomial.h $(AUTOSUNDIR)/print.h \
	$(AUTOSUNDIR)/representation.h $(AUTOSUNDIR)/sparse.h \
	$(AUTOSUNDIR)/style.h $(AUTOSUNDIR)/sun.h $(AUTOSUNDIR)/symmetric.h
	$(CXX) -I $(AUTOSUNDIR) -o write_spn_sun_algconv write_spn_sun_algconv.cpp -g -O0 ${MACRO} &&\
	echo "write_spn_sun_algconv done." || ( echo " failed." && false )

include $(MKDIR)/MkRules


clean::
	rm -f writeREPR spn_vs_sun_algconv.h spnalgtest lib.o \
		spnalgmacros_adj spnalgmacros_adj.cpp \
		spnalgmacros_asym spnalgmacros_asym.cpp \
		spn_sun_algconv.h write_spn_sun_algconv \
		spN* suN*
