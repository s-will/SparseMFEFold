dnl                                                -*- Autoconf -*-
dnl  Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([SparseMFEFold], [0.4], [])

dnl special dir for aux config files 
AC_CONFIG_AUX_DIR(config)

dnl special dir for m4
AC_CONFIG_MACRO_DIR([m4])

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
dnl rely on pkg-config

dnl check whether pkg-config is available 
AC_CHECK_PROG([PKG_CONFIG_AVAILABLE],pkg-config,[true],[false])

VRNALIB_OK=false

AS_IF([test "$PKG_CONFIG_AVAILABLE" = "true"],[
   dnl call pkg-config (macro PKG_CHECK_MODULES) to set VRNA_LIBS and VRNA_CFLAGS
   PKG_CHECK_MODULES([VRNA],[RNAlib2 >= 2.1.9],[VRNALIB_OK=true],[VRNALIB_OK=false])
   dnl replace ViennaRNA add end of include path -- use perl, since we require it anyway
   dnl This fixes a legacy problem in ViennaRNA <2.2.x pkg-config files
   VRNA_CFLAGS=$(echo $VRNA_CFLAGS | perl -ne 's/(-I.*)\/ViennaRNA(\/\s|\s|$)/\1\2/g;print')

   CXXFLAGS="$VRNA_CFLAGS $CXXFLAGS"
   LDFLAGS="$VRNA_LDFLAGS $LDFLAGS"
   LIBS="$VRNA_LIBS $LIBS"
])

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
if test "$VRNLIB_OK" = "false"; then
        AC_MSG_NOTICE()
        AC_MSG_NOTICE([The Vienna RNA library (or pkg-config) is missing.])
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

