# Copyright (C) 2012-2015  CEA/DEN, EDF R&D
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

SET(MED_FILES
  testfield.med
  smallmesh_varfield.med
  smallmesh_constfield.med
  bigmesh.med
  fieldlena.med
  timeseries.med
  parametric_01.med
  parametric_02.med
  enceinte-axis-E500000.med
  enceinte-axis-E600000.med
  entaille.med
  padderoutput.med
  )

INSTALL(FILES ${MED_FILES} DESTINATION ${SALOME_MED_INSTALL_RES_DATA}/medop_testfiles)

SET(MED_IMG_FILES
  irm_test1.png
  irm_test2.png
  # extra
  README.txt
  champ_altitude_MAP.pgm
  )

INSTALL(FILES ${MED_IMG_FILES} DESTINATION ${SALOME_MED_INSTALL_RES_DATA}/medop_testfiles)
