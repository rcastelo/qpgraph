### configure.ac                                        -*- Autoconf -*-
###
### Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)

### analogous to Perl's chop command
define([m4_chop], [m4_substr($1, 0, m4_decr(m4_len($1)))])

### get the version number from the DESCRIPTION file of the package
m4_define([qp_VERSION],
          m4_chop(m4_esyscmd([egrep '^Version:' DESCRIPTION | cut -d ' ' -f 2])))

AC_INIT([qp],[qp_VERSION],[robert.castelo@upf.edu])
AC_CONFIG_AUX_DIR([tools])

AC_PROG_CC

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

### ** Platform. (taken from the configure.ac of R)

AC_CANONICAL_HOST
AC_DEFINE_UNQUOTED(R_PLATFORM, "${host}",
[Define this to be the canonical name (cpu-vendor-os) of your system.])
AC_DEFINE_UNQUOTED(R_CPU, "${host_cpu}",
[Define this to be the name of the CPU of your system.])
AC_DEFINE_UNQUOTED(R_VENDOR, "${host_vendor}",
[Define this to be the name of the vendor of your system.])
AC_DEFINE_UNQUOTED(R_OS, "${host_os}",
[Define this to be the name of the OS of your system.])

case "${host_os}" in
  cygwin*|mingw*|windows*|winnt)
    AC_DEFINE(Win32, 1,
              [Define according to your operating system type.])
    R_OSTYPE="windows"
    ;;
  *)
    AC_DEFINE(Unix, 1,
              [Define according to your operating system type.])
    R_OSTYPE="unix"
   ;;
esac

CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`" -Wall -pedantic"
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
