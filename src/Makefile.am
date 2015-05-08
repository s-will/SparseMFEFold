bin_PROGRAMS=SparseMFEFold

# base sources, common to all programs
SparseMFEFold_SOURCES = SparseMFEFold.cc base.hh simple_map.hh	\
	trace_arrow.hh trace_arrow.cc LocARNA/matrix.hh		\
	SparseMFEFold_cmdline.c

############################################################
## gengetopt command line parsers
##
## register parsers and config.h as built sources
BUILT_SOURCES =
#
## register built sources as cleanfiles (otherwise they would not be
## removed by 'make clean')
CLEANFILES = $(BUILT_SOURCES)
#
#%_cmdline.c: $(srcdir)/%.ggo
#	gengetopt < $<

## get a copy of config.h in directory src 
config.h: ../config.h
	if [ -e "$@" ] ; then $(RM) $@ ; fi
	ln -s $< $@


SparseMFEFold_cmdline.c: $(srcdir)/SparseMFEFold.ggo 
	gengetopt < $<
