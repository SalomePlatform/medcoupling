# Copyright (C) 2012-2014  CEA/DEN, EDF R&D, OPEN CASCADE
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

SET(dist_admlocalm4_DATA
    ac_check_sizeof_fortran.m4
    check_GUI.m4
    check_Med.m4
    check_Med2.m4
    check_bft.m4
    check_fvm.m4
    check_metis.m4
    check_parmetis.m4
    check_scotch.m4
    med_check_sizeof_medint.m4
    renumber.m4
    splitter.m4
    with_Kernel.m4
)
            
FOREACH(f ${dist_admlocalm4_DATA})
  SET(DEST adm_local/unix/config_files)
  INSTALL(FILES ${f} DESTINATION ${DEST})
ENDFOREACH(f ${dist_admlocalm4_DATA})
