#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([ccs], [0.1], [michele.orru@ens.fr])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/rms.c])
AC_CONFIG_HEADERS([src/config.h])

# Checks for programs.
AC_PROG_CC([gcc clang cc])
AC_PROG_CC_C99
AC_PROG_RANLIB

# Checks for libraries.

# have _GNU_SOURCE defined aroud.
AC_GNU_SOURCE

# Checks for header files.
AC_CHECK_HEADERS([limits.h stdint.h stdlib.h string.h unistd.h])
AC_CHECK_LIB(gmp, __gmpz_init, ,
[AC_MSG_ERROR([GNU MP not found, see https://gmplib.org/])])
AC_CHECK_SIZEOF(mp_limb_t, 8, [#include <gmp.h>])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
#AC_CHECK_FUNCS([dup2 setlocale strdup])

# Add compiler/linker flags
CFLAGS="--std=c99 -Wall -march=native"

# Adding package options
AC_ARG_ENABLE(debug,
   AS_HELP_STRING([--enable-debug], [enable debugging, default: no]),
   CFLAGS+=" -DDEBUG -DBN_DEBUG -O0 -ggdb")

AC_OUTPUT([Makefile
           src/Makefile
           ])
