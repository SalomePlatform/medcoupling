#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2007-2016  CEA/DEN, EDF R&D
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
# Author : Anthony GEAY (CEA/DEN/DM2S/STMF/LGLS)

from MEDLoader import *

class CaseIO:
    dictMCTyp={NORM_HEXA8:"hexa8",NORM_POLYHED:"nfaced",NORM_QUAD4:"quad4",NORM_POLYGON:"nsided",NORM_POINT1:"point",NORM_SEG2:"bar2",NORM_SEG3:"bar3",NORM_TRI3:"tria3",NORM_TRI6:"tria6",NORM_QUAD8:"quad8",NORM_TETRA4:"tetra4",NORM_TETRA10:"tetra10",NORM_PYRA5:"pyramid5",NORM_PYRA13:"pyramid13",NORM_PENTA6:"penta6",NORM_PENTA15:"penta15",NORM_HEXA20:"hexa20"}
    discSpatial={ON_CELLS:"element",ON_NODES:"node"}
    dictCompo={1:"scalar",3:"vector",6:"tensor",9:"tensor9"}
    dictMCTyp2 = {v:k for k, v in dictMCTyp.items()}
    discSpatial2 = {v:k for k, v in discSpatial.items()}
    dictCompo2 = {v:k for k, v in dictCompo.items()}
    pass
