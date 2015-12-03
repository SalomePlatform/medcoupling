# Copyright (C) 2015  CEA/DEN, EDF R&D
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#

####################################################################
# _TOHEXA()
# Convert a number (smaller than 16) into hexadecimal representation
# with a leading 0.
MACRO(_TOHEXA num result)
  SET(_hexa_map a b c d e f)
  IF(${num} LESS 10)
    SET(${result} "0${num}")
  ELSE()
    MATH(EXPR _res "${num}-10" )
    LIST(GET _hexa_map ${_res} _out)
    SET(${result} "0${_out}")
  ENDIF()
ENDMACRO(_TOHEXA)

####################################################################
# MEDCOUPLING_XVERSION()
#
# Computes hexadecimal version of MEDCOUPLING package
#
# USAGE: MEDCOUPLING_XVERSION(package)
#
# ARGUMENTS:
#
# package: IN: MEDCOUPLING package name
#
# The macro reads MEDCOUPLING package version from PACKAGE_VERSION variable
# (note package name are uppercase);
# hexadecimal version value in form 0xAABBCC (where AA, BB and CC are
# major, minor and maintenance components of package version in
# hexadecimal form) is put to the PACKAGE_XVERSION variable
MACRO(MEDCOUPLING_XVERSION pkg)
  STRING(TOUPPER ${pkg} _pkg_UC)
  IF(${_pkg_UC}_VERSION)
    SET(_major)
    SET(_minor)
    SET(_patch)
    _TOHEXA(${${_pkg_UC}_MAJOR_VERSION} _major)
    _TOHEXA(${${_pkg_UC}_MINOR_VERSION} _minor)
    _TOHEXA(${${_pkg_UC}_PATCH_VERSION} _patch)
    SET(${_pkg_UC}_XVERSION "0x${_major}${_minor}${_patch}")
  ENDIF()
ENDMACRO(MEDCOUPLING_XVERSION)
