bin_PROGRAMS=halfzuker

# base sources, common to all programs
halfzuker_SOURCES=halfzuker.cc simple_map.hh trace_arrow.hh

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
