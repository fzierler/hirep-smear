TOPDIR = ..
MKDIR = $(TOPDIR)/Make

EXES =  measure_spectrum #mk_mesons random_cnfg random_spinor mk_mesons_with_z2semwall mk_mesons_with_z2semwall_new  #measure_formfactor #mk_sfcoupling #trunc_mesons

COMMON_OBJS = 

measure_spectrum_OBJS = meson_measurements.o meson_measurements_ff.o baryon_measurements.o ../GaugeFix/gaugefix.o 
#measure_formfactor_OBJS = meson_measurements.o ../GaugeFix/gaugefix.o

LIBS += $(TOPDIR)/LibHR/libhr.a


include $(MKDIR)/MkRules

