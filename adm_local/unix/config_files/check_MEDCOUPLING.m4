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

# Check availability of MEDCOUPLING binary distribution
#
# Author : Roman NIKOLAEV (OPEN CASCADE, 2016)
#

AC_DEFUN([CHECK_MEDCOUPLING],[
AC_REQUIRE([AC_LINKER_OPTIONS])dnl

AC_CHECKING(for MEDCOUPLING)

MEDCOUPLING_ok=no

MEDCOUPLING_LDFLAGS=""
MEDCOUPLING_CXXFLAGS=""

AC_ARG_WITH(medcoupling,
	    [  --with-medcoupling=DIR root directory path of MEDCOUPLING installation ],
	    MEDCOUPLING_DIR="$withval",MEDCOUPLING_DIR="")

if test "x${MEDCOUPLING_DIR}" == "x" ; then
  AC_MSG_RESULT(for \${MEDCOUPLING_ROOT_DIR}: ${MEDCOUPLING_ROOT_DIR})

   # --with-medcoupling option is not used
   if test "x${MEDCOUPLING_ROOT_DIR}" != "x" ; then
    # MEDCOUPLING_ROOT_DIR environment variable defined
      MEDCOUPLING_DIR=${MEDCOUPLING_ROOT_DIR}
   fi

fi

if test -f ${MEDCOUPLING_DIR}/include/InterpKernelValue.hxx ; then
   AC_MSG_RESULT(Using MEDCOUPLING module distribution in ${MEDCOUPLING_DIR})
   MEDCOUPLING_ok=yes

   AC_SUBST(MEDCOUPLING_ROOT_DIR)

   MEDCOUPLING_LDFLAGS=-L${MEDCOUPLING_DIR}/lib${LIB_LOCATION_SUFFIX}
   MEDCOUPLING_CXXFLAGS=-I${MEDCOUPLING_DIR}/include

   AC_SUBST(MEDCOUPLING_LDFLAGS)
   AC_SUBST(MEDCOUPLING_CXXFLAGS)

else
   AC_MSG_WARN("Cannot find MCOUPLING module sources")
fi

AC_MSG_RESULT(for MEDCOUPLING: $MEDCOUPLING_ok)

])dnl
