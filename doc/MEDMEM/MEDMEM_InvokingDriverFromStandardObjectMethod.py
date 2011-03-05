#  -*- coding: iso-8859-1 -*-
#  Copyright (C) 2007-2010  CEA/DEN, EDF R&D, OPEN CASCADE
#
#  Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
#  CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
#  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#

######################################################################
# This Python script should be executed when the shared library is   #
# generated using SWIG 1.3 (or higher) due to the fact that older    #
# version could not handle the wrapping of several class constructor #
######################################################################
#
from libMEDMEM_Swig import *

medFile = "pointe.med"
medFile2 = "fieldCellDoubleOfpointe.med"
fieldName = "fieldcelldouble"
meshName = "maa1"

try:
    myField = FIELDDOUBLE()

    myDriver1 = myField.addDriver(MED_DRIVER,medFile,fieldName)
    myField.rmDriver()

    myDriver2 = myField.addDriver(MED_DRIVER,medFile2,fieldName)
    myField.rmDriver(myDriver2)

    myMesh = MESH()
    myDriver3 = myMesh.addDriver(MED_DRIVER,medFile,meshName)
    myMesh.read()
    myMesh.rmDriver()

    myMed = MED(MED_DRIVER,medFile)
    myMed.readFileStruct()
    myMed.rmDriver()

except:
    print "There is a problem somewhere !!"
    print "Please consult the error standart output of the python execution !!"
