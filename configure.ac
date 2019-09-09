AC_INIT("DESCRIPTION")

# Check the compiler configured with R
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`

AC_CHECK_LIB([z], [gzeof], , AC_ERROR([zlib not found]))
AC_CHECK_SIZEOF([unsigned long])

AC_OUTPUT(src/Makevars)
