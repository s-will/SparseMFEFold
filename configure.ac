dnl                                                -*- Autoconf -*-
dnl  Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([SparseMFEFold], [0.4], [])

dnl  special dir for aux config files 
AC_CONFIG_AUX_DIR(config)

AM_INIT_AUTOMAKE


dnl  Checks for programs.
AC_PROG_CXX
AC_PROG_CC

dnl TURN ON C++11
AX_CXX_COMPILE_STDCXX_11

AC_LANG_CPLUSPLUS

dnl  Checks for header files.
AC_HEADER_STDC

dnl # for compiling modified ViennaRNA code
if test "$GCC" = yes; then
   AC_DEFINE(UNUSED, __attribute__ ((unused)), avoid warning about unused variables)
else
   AC_DEFINE(UNUSED,, only for gcc)
fi
AC_CHECK_FUNCS(strdup strstr strchr erand48)


dnl # ========================================
dnl # DOXYGEN
dnl 
DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)
DX_INIT_DOXYGEN(HybridPF, Doxyfile, doc)


dnl # c++11
dnl CXXFLAGS="$CXXFLAGS -std=c++11"

dnl  optimization
CXXFLAGS="$CXXFLAGS -O3 -fomit-frame-pointer"

dnl #################
dnl # Debug option
dnl 
AC_MSG_CHECKING([whether to build with debug information])
debuger=no
AC_ARG_ENABLE([debug],
    AC_HELP_STRING(
        [--enable-debug],
        [enable debug data generation (def=no)]
    ),
    debuger="$enableval"
)
AC_MSG_RESULT($debuger)
if test x"$debuger" = x"yes"; then
    AC_DEFINE([DEBUG], [], [Define to 1 for turning on debugging])
    CXXFLAGS="$CXXFLAGS -g -Wall"
else
    AC_DEFINE([NDEBUG], [], [Define to 1 for turning off debugging])
fi


dnl # ========================================
dnl # Vienna RNA
dnl 

dnl  Get Vienna RNA library and include locations
AC_ARG_WITH([vrna],
  [AS_HELP_STRING([--with-vrna],
    [location of the Vienna RNA headers])],
  [VRNA_HOME=$withval;
   VRNA_CFLAGS="-I$withval/include"; 
   VRNA_LDFLAGS="-L$withval/lib"],
  [])

LIBS="-lRNA $LIBS"

LDFLAGS="$VRNA_LDFLAGS $LDFLAGS"
CXXFLAGS="$VRNA_CFLAGS $CXXFLAGS"

AC_CHECK_HEADER("ViennaRNA/fold.h",
	[VRNANOTFOUND=0],
	[VRNANOTFOUND=0],
	[])

dnl  ------------------------------
dnl  handle open mp support
dnl   this is necessary when linking against the vienna lib
AC_OPENMP
AC_SUBST(OPENMP_CXXFLAGS)

CXXFLAGS="$OPENMP_CXXFLAGS $CXXFLAGS"
LDFLAGS="$OPENMP_CXXFLAGS $LDFLAGS"


dnl ###################
dnl  Static linking
dnl 
AC_MSG_CHECKING([whether to link statically])
static_linkage=no
AC_ARG_ENABLE([static],
    AC_HELP_STRING(
        [--enable-static],
        [link statically (def=no)]
    ),
    static_linkage="$enableval"
)
AC_MSG_RESULT($static_linkage)
if test x"$static_linkage" = x"yes"; then
   LDFLAGS="$LDFLAGS -static"
fi


dnl  Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_TYPE_SIZE_T


dnl ##############################################################################
dnl  BEGIN DO NOT COMPILE CHECKS
dnl ##############################################################################

DEPENDENCYNOTFOUND=0;

dnl  error if ViennaRNA not found
if test "$VRNANOTFOUND" = "1"; then
        AC_MSG_NOTICE()
        AC_MSG_NOTICE([The Vienna RNA library is missing.])
        AC_MSG_NOTICE()
        DEPENDENCYNOTFOUND=1;
fi

dnl  error if on of the libraries was not found
if test "$DEPENDENCYNOTFOUND" = "1"; then
        AC_MSG_NOTICE()
        AC_MSG_ERROR([library not found! Try './configure --help'.])
fi


AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile])

AC_OUTPUT([src/Makefile])

