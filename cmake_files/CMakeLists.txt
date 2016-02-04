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
# See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#

SET(admlocal_cmake_DATA
    FindCppUnit.cmake
    FindGraphviz.cmake
    FindMEDFile.cmake
    FindMetis.cmake
    FindParMetis.cmake
    FindSalomeBoost.cmake
    FindSalomeCppUnit.cmake
    FindSalomeDoxygen.cmake
    FindSalomeGraphviz.cmake
    FindSalomeHDF5.cmake
    FindSalomeLibXml2.cmake
    FindSalomeMEDCoupling.cmake
    FindSalomeMEDFile.cmake
    FindSalomeMetis.cmake
    FindSalomeMPI.cmake
    FindSalomeParMetis.cmake
#    FindSalomePython.cmake
    FindSalomePythonLibs.cmake
    FindSalomePythonInterp.cmake
    FindSalomeScotch.cmake
    FindSalomeSphinx.cmake
    FindSalomeSWIG.cmake
    FindSalomeXDR.cmake
    FindScotch.cmake
    FindSphinx.cmake
    FindXDR.cmake
    SalomeMacros.cmake
    SalomeSetupPlatform.cmake
)

INSTALL(FILES ${admlocal_cmake_DATA} DESTINATION ${MEDCOUPLING_INSTALL_CMAKE_LOCAL})
