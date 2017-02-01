dnl Copyright (C) 2007-2016  CEA/DEN, EDF R&D, OPEN CASCADE
dnl
dnl Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
dnl CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
dnl
dnl This library is free software; you can redistribute it and/or
dnl modify it under the terms of the GNU Lesser General Public
dnl License as published by the Free Software Foundation; either
dnl version 2.1 of the License, or (at your option) any later version.
dnl
dnl This library is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
dnl Lesser General Public License for more details.
dnl
dnl You should have received a copy of the GNU Lesser General Public
dnl License along with this library; if not, write to the Free Software
dnl Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
dnl
dnl See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
dnl

AC_DEFUN([CHECK_SCOTCH],[

AC_CHECKING(for scotch)

AC_LANG_SAVE
AC_LANG_C

dnl
dnl default values
dnl
SCOTCH_CPPFLAGS="SCOTCH_CPPFLAGS_NOT_DEFINED"
SCOTCH_LIBS="SCOTCH_LIBS_NOT_FOUND"
SCOTCH_LIBSUFFIX="-not-defined"

ENABLE_SCOTCH="no"

AC_CHECKING(for scotch location)

AC_ARG_WITH(scotch,
            [  --with-scotch=DIR      root directory path to SCOTCH library installation ],
            [SCOTCHDIR="$withval"
             AC_MSG_RESULT([Select $withval as path to SCOTCH library])])

if test "x${SCOTCHDIR}" == "x" ; then
  AC_MSG_RESULT(for \${SCOTCH_ROOT_DIR}: ${SCOTCH_ROOT_DIR})
   # --with-scotch option is not used
   if test "x${SCOTCH_ROOT_DIR}" != "x" ; then
      SCOTCHDIR=${SCOTCH_ROOT_DIR}
   fi
fi

if test "x${SCOTCHDIR}" = "x" ; then
  AC_MSG_WARN(SCOTCHDIR is not specified)
  AC_MSG_NOTICE(Trying native scotch...)
  SCOTCHDIR=/usr
fi

CPPFLAGS_old="${CPPFLAGS}"
LIBS_old=$LIBS

scotch_ok=no
scotch_headers_ok=no
scotch_binaries_ok=no

dnl
dnl SCOTCH headers
dnl
AC_CHECKING(for SCOTCH headers)

for d in ${SCOTCHDIR}/include ${SCOTCHDIR}/include/scotch ${SCOTCHDIR}/bin ${SCOTCHDIR} ; do
  dnl
  dnl check scotch.h file availability
  dnl
  AC_CHECK_FILE([${d}/scotch.h],
                [scotch_include_dir_ok=yes],
                [scotch_include_dir_ok=no])

  if test "x${scotch_include_dir_ok}" = "xyes" ; then
    LOCAL_INCLUDES="-DMED_ENABLE_SCOTCH -I${d} ${MPI_INCLUDES}"
    CPPFLAGS="${CPPFLAGS_old} ${LOCAL_INCLUDES} -std=c99"
    AC_TRY_COMPILE([
      #include <stdio.h>
      #include <scotch.h>
      ],
      [SCOTCH_Graph* graph; SCOTCH_graphInit(graph)],
      [scotch_headers_ok=yes],
      [scotch_headers_ok=no])
    if test "x${scotch_headers_ok}" = "xyes" ; then
      break;
    fi
  fi
done

dnl
dnl result for headers
dnl
AC_MSG_RESULT(for scotch headers: $scotch_headers_ok)

dnl
dnl SCOTCH libraries
dnl
if test "x${scotch_headers_ok}" = "xyes" ; then
  AC_CHECKING(for scotch libraries)

  for d in ${SCOTCHDIR}/lib ${SCOTCHDIR}/lib64 ${SCOTCHDIR}/lib/scotch ${SCOTCHDIR}/lib64/scotch ${SCOTCHDIR}/bin ${SCOTCHDIR} ; do
    LOCAL_LIBS="-L${d} -lscotch -lscotcherr"
    LIBS="${LIBS_old} ${LOCAL_LIBS}"
    AC_TRY_LINK([
      #include <stdio.h>
      #include <scotch.h>
      ],
      [SCOTCH_Graph* graph; SCOTCH_graphInit(graph)],
      [scotch_binaries_ok=yes],
      [scotch_binaries_ok=no])
    if test "x${scotch_binaries_ok}" = "xyes" ; then
      break;
    fi
  done
fi

dnl
dnl result for libraries
dnl
AC_MSG_RESULT(for scotch libraries: $scotch_binaries_ok)

dnl
dnl summary
dnl
if test "x${scotch_binaries_ok}" = "xyes" ; then
  SCOTCH_CPPFLAGS=${LOCAL_INCLUDES}
  SCOTCH_LIBS=${LOCAL_LIBS}
  SCOTCH_LIBSUFFIX=""
  ENABLE_SCOTCH="yes"
  scotch_ok=yes
  AC_MSG_RESULT(\$SCOTCH_CPPFLAGS  = ${SCOTCH_CPPFLAGS})
  AC_MSG_RESULT(\$SCOTCH_LIBS      = ${SCOTCH_LIBS})
  AC_MSG_RESULT(\$SCOTCH_LIBSUFFIX = ${SCOTCH_LIBSUFFIX})
fi
AC_MSG_RESULT(for scotch: $scotch_ok)

CPPFLAGS="${CPPFLAGS_old}"
LIBS="${LIBS_old}"

AC_SUBST(SCOTCH_CPPFLAGS)
AC_SUBST(SCOTCH_LIBSUFFIX)
AC_SUBST(SCOTCH_LIBS)
AC_SUBST(ENABLE_SCOTCH)

AC_LANG_RESTORE

])dnl


