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

# Sphinx detection for Salome
#
#  !! Please read the generic detection procedure in SalomeMacros.cmake !!
#

SALOME_FIND_PACKAGE_AND_DETECT_CONFLICTS(Sphinx SPHINX_EXECUTABLE 2)

# Also retrieve paths to DOCUTILS and SETUPTOOLS:
SET(SETUPTOOLS_ROOT_DIR "$ENV{SETUPTOOLS_ROOT_DIR}" CACHE PATH "Path to the Setuptools installation")
SET(DOCUTILS_ROOT_DIR "$ENV{DOCUTILS_ROOT_DIR}" CACHE PATH "Path to the Docutils installation")

# Ensure the command is run with the given PYTHONPATH
IF(WIN32 AND NOT CYGWIN)
   SET(SPHINX_EXECUTABLE ${SPHINX_EXECUTABLE})
   SET(SPHINX_APIDOC_EXECUTABLE ${SPHINX_APIDOC_EXECUTABLE})
ELSE()
   SET(SPHINX_EXECUTABLE /usr/bin/env PYTHONPATH="${SPHINX_PYTHONPATH}:$$PYTHONPATH" ${SPHINX_EXECUTABLE})
   SET(SPHINX_APIDOC_EXECUTABLE /usr/bin/env PYTHONPATH="${SPHINX_PYTHONPATH}:$$PYTHONPATH" ${SPHINX_APIDOC_EXECUTABLE})
ENDIF()

MARK_AS_ADVANCED(SPHINX_EXECUTABLE)

IF(SPHINX_FOUND)
  SALOME_ACCUMULATE_ENVIRONMENT(PATH ${SPHINX_EXECUTABLE})
  SALOME_ACCUMULATE_ENVIRONMENT(PATH ${SPHINX_APIDOC_EXECUTABLE})
  SALOME_ACCUMULATE_ENVIRONMENT(PYTHONPATH ${SPHINX_PYTHONPATH})
ENDIF()
