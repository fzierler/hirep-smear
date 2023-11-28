TOPDIR = ../..
MKDIR = $(TOPDIR)/Make

EXES = check_update_1 check_exp

check_update_1_OBJS = ../../HMC/hmc_utils.o

LIBS += $(TOPDIR)/LibHR/libhr.a 



include $(MKDIR)/MkRules

commutator : commutator.cc $(MKDIR)/MkFlags
	echo $(CXX)  _$(REPR)_ _${GAUGE_GROUP}_ &&\
	$(CXX) -o commutator commutator.cc -D_$(REPR)_ -D_${GAUGE_GROUP}_ -D_PICA_STYLE_ -I ../../Make/Utils/autosun

commutator.h: commutator
	./commutator $(NG)





