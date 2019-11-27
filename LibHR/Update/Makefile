TOPDIR = ../..
MKDIR = $(TOPDIR)/Make

SRCS = $(wildcard *.c)
LIBNAME = ../libhr.a

#Exclude these source files
EXCLUDE = 

include $(MKDIR)/MkRules

approx_data.db: remez/approx_*
	echo -n "$(_col_write)Writing rational approximation database$(_col_norm)..."
	( cd remez && cat approx_* | ./mkappdata.pl >../$@ && cd .. \
		&& echo " $(_cdone)" ) || ( echo " $(_cfailed)" && false )


ifeq "-O3" "$(findstring -O3,$(CFLAGS))"

CFLAGS_SUB=$(patsubst -O3,-O2,$(CFLAGS))
 
force_fermion_core.o: $(MKDIR)/MkFlags $(RBLD) $(INCDIR)/.autoheaders.h force_fermion_core.c
	if  test "$(call _testdir, $@)" = ""; then \
	  if test "$(DFILES)" = ""; then \
	    $(MAKE) $(DEP) &&\
	    MAKEFILES="$(DEP)" $(MAKE) $@ ;\
	  else \
	    echo -n "$(_col_comp)Compiling$(_col_norm) $(call _cfname,$(call _rname,$(CURDIR)/$@))..." ;\
		echo "$(_col_bold)$(_col_red)INFO: Changing -O3 to -O2 for force_fermion_core.c$(_col_norm)"; \
	    $(CC) -c $(CFLAGS_SUB) $(CPPFLAGS) $(SRC) &&\
	    echo " $(_cdone)" || ( echo " $(_cfailed)" && false) ;\
	  fi ;\
	else \
	  $(MAKE) -C $(dir $@) $(notdir $@) ;\
	fi

endif

# TODO: Fix repeated code 
ifeq "-Ofast" "$(findstring -Ofast,$(CFLAGS))"

CFLAGS_SUB=$(patsubst -Ofast,-O2,$(CFLAGS))
 
force_fermion_core.o: $(MKDIR)/MkFlags $(RBLD) $(INCDIR)/.autoheaders.h force_fermion_core.c
	if  test "$(call _testdir, $@)" = ""; then \
	  if test "$(DFILES)" = ""; then \
	    $(MAKE) $(DEP) &&\
	    MAKEFILES="$(DEP)" $(MAKE) $@ ;\
	  else \
	    echo -n "$(_col_comp)Compiling$(_col_norm) $(call _cfname,$(call _rname,$(CURDIR)/$@))..." ;\
		echo "$(_col_bold)$(_col_red)INFO: Changing -Ofast to -O2 for force_fermion_core.c$(_col_norm)"; \
	    $(CC) -c $(CFLAGS_SUB) $(CPPFLAGS) $(SRC) &&\
	    echo " $(_cdone)" || ( echo " $(_cfailed)" && false) ;\
	  fi ;\
	else \
	  $(MAKE) -C $(dir $@) $(notdir $@) ;\
	fi

endif
