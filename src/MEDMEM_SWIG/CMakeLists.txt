# Copyright (C) 2007-2011  CEA/DEN, EDF R&D
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License.
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

FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

SET_SOURCE_FILES_PROPERTIES(libMEDMEM_Swig.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(libMEDMEM_Swig.i PROPERTIES SWIG_FLAGS "-shadow")
SET_SOURCE_FILES_PROPERTIES(libMEDMEM_SwigPYTHON_wrap.cxx PROPERTIES COMPILE_FLAGS "${HDF5_FLAGS} ${MED3_FLAGS} ${XDR_FLAGS}")

INCLUDE_DIRECTORIES(
  ${PYTHON_INCLUDES_DIR}
  ${MED3_INCLUDES_DIR}
  ${HDF5_INCLUDES_DIR}
  ${XDR_INCLUDES_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDMEM
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/GaussPoints
  )

SWIG_ADD_MODULE(libMEDMEM_Swig python libMEDMEM_Swig.i)
SWIG_LINK_LIBRARIES(libMEDMEM_Swig ${PYTHON_LIBRARIES} medmem)

INSTALL(TARGETS _libMEDMEM_Swig DESTINATION ${MED_salomepythondir})
INSTALL(FILES libMEDMEM_Swig.i my_typemap.i MEDMEM_SWIG_Templates.hxx DESTINATION ${MED_salomeinclude_HEADERS})
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/libMEDMEM_Swig.py DESTINATION ${MED_salomescript_PYTHON})
INSTALL(FILES med2sauv sauv2med avs2med DESTINATION ${MED_salomescript_PYTHON}
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
