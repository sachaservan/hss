#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([hss], [0.1], [michele.orru@ens.fr])


AC_CONFIG_SRCDIR([src/rms.c])
AC_CONFIG_HEADERS([src/config.h])

# Checks for programs.

AC_CHECK_TARGET_TOOL([latexmk], [latexmk], [no])
#AS_IF([test "x$latexmk" = "xno"], [AC_MSG_ERROR([cannot find pdflatex.])])
AC_PROG_CC([gcc clang cc])
AC_PROG_CC_C99
AC_PROG_RANLIB

# Checks for libraries.
AC_CHECK_LIB(gmp, __gmpz_init, ,
             [AC_MSG_ERROR([GNU MP not found, see https://gmplib.org/])])

# have _GNU_SOURCE defined aroud.
AC_GNU_SOURCE

# Checks for header files.
AC_CHECK_HEADERS([limits.h stdint.h stdlib.h string.h unistd.h])
AC_CHECK_SIZEOF(mp_limb_t, 8, [#include <gmp.h>])
AC_CHECK_SIZEOF(uint32_t, 4, [#include <stdint.h>])
AC_CHECK_SIZEOF(unsigned long long, 8, [])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
#AC_CHECK_FUNCS([dup2 setlocale strdup])


AM_INIT_AUTOMAKE

# Clear out compiler/linker flags
CFLAGS+=" -pedantic -Wall "

# Shut up automake
AM_SILENT_RULES([yes])
AC_SUBST([AM_MAKEFLAGS], [--no-print-directory])

# Adding package options

AC_ARG_ENABLE(debug,
   AS_HELP_STRING([--enable-debug], [enable debugging, default: no]),
   CFLAGS+=" -UNDEBUG -O0 -ggdb -fbounds-check -D_FORTIFY_SOURCE=1",
   CFLAGS+=" -DNDEBUG -O3 -march=native")



AC_DEFINE([SEP], ["\t"], [separator for benchmarks])

AC_ARG_VAR([failure], [failure prob.])
AS_IF([test "x$failure" = x], [failure="17"])
AC_DEFINE_UNQUOTED([FAILURE], [$failure], [log inverse of the failure probability, default: 17.])

AC_ARG_VAR([fb_base], [fixed base])
AS_IF([test "x$fb_base" = x], [fb_base="8"])
AC_DEFINE_UNQUOTED([FB_BASE], [$fb_base], [log of precomputed base for group operation, default: 8.])

AC_ARG_VAR([ss_base], [ss base])
AS_IF([test "x$ss_base" = x], [ss_base="1"])
AC_DEFINE_UNQUOTED([SS_BASE], [$ss_base], [log of secret shares representation, default: 1.])


AC_OUTPUT([Makefile
           src/Makefile
           ])
