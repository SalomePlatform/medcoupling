dnl Copyright (C) 2007-2015  CEA/DEN, EDF R&D, OPEN CASCADE
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

dnl define macros :
dnl AC_ENABLE_RENUMBER AC_DISABLE_RENUMBER
dnl author Vincent BERGEAUD
dnl
# AC_ENABLE_RENUMBER
AC_DEFUN([AC_ENABLE_RENUMBER], [dnl
define([AC_ENABLE_RENUMBER_DEFAULT], ifelse($1, yes, yes, no))dnl
AC_ARG_ENABLE(renumber,
changequote(<<, >>)dnl
<<  --enable-renumber[=PKGS]  build with renumber tool [default=>>AC_ENABLE_RENUMBER_DEFAULT],
changequote([, ])dnl
[p=${PACKAGE-default}
case "$enableval" in
yes) enable_renumber=yes ;;
no) enable_renumber=no ;;
*)
  enable_renumber=no
  # Look at the argument we got.  We use all the common list separators.
  IFS="${IFS= 	}"; ac_save_ifs="$IFS"; IFS="${IFS}:,"
  for pkg in $enableval; do
    if test "X$pkg" = "X$p"; then
      enable_renumber=yes
    fi
  done
  IFS="$ac_save_ifs"
  ;;
esac],
enable_renumber=AC_ENABLE_RENUMBER_DEFAULT)dnl

if test "X$enable_renumber" = "Xyes"; then
  ENABLE_RENUMBER="yes"
  AC_SUBST(ENABLE_RENUMBER)
fi
])

# AC_DISABLE_RENUMBER - set the default flag to --disable-renumber
AC_DEFUN([AC_DISABLE_RENUMBER], [AC_ENABLE_RENUMBER(no)])

