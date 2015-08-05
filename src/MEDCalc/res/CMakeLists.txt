# Copyright (C) 2007-2015  CEA/DEN, EDF R&D
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

ADD_SUBDIRECTORY(testfiles)

SET(MED_RESOURCES_FILES
  datasource_add.png
  datasource_changeUnderlyingMesh.png
  datasource_interpolateField.png
  datasource_expandfield.png
  datasource_field.png
  datasource_mesh.png
  datasource.png
  datasource_use.png
  datasource_view.png
  fileimport-32.png
  folder.png
  image_add.png
  MEDCalc.png
  MEDCalc_small.png
  workspace_clean.png
  workspace_save.png
  pqAppIcon16.png
  )
INSTALL(FILES ${MED_RESOURCES_FILES} DESTINATION ${SALOME_MED_INSTALL_RES_DATA})
