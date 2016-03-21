# Copyright (C) 2013-2016  CEA/DEN, EDF R&D, OPEN CASCADE
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
# See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#
# Author: Adrien Bruneton
#

# Doxygen detection for salome
#
#  !! Please read the generic detection procedure in SalomeMacros.cmake !!
#
# Additional variables:
#
# DOXYGEN_SUPPORT_STL (string) [advanced] : set to YES if doxygen properly manages STL files
#                     or to NO otherwise (version 1.4.4 or older); see description of 
#                     BUILTIN_STL_SUPPORT configuration variable in the doxygen documentation

SALOME_FIND_PACKAGE_AND_DETECT_CONFLICTS(Doxygen DOXYGEN_EXECUTABLE 2)
IF(DOXYGEN_FOUND)
  IF(DOXYGEN_VERSION VERSION_LESS "1.4.5")
    SET(DOXYGEN_SUPPORT_STL NO)
  ELSE()
    SET(DOXYGEN_SUPPORT_STL YES)
  ENDIF()
ENDIF()
MARK_AS_ADVANCED(DOXYGEN_SUPPORT_STL)

IF(DOXYGEN_FOUND)
  SALOME_ACCUMULATE_ENVIRONMENT(PATH ${DOXYGEN_EXECUTABLE})
ENDIF()
