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
medFile2 = "Field&MeshGeneratedPointe.med"
fieldName = "fieldcelldoublescalar"
meshName = "maa1"

try:
    myField = FIELDDOUBLE()
    myRdOnlyDriver = MED_FIELDDOUBLE_RDONLY_DRIVER(medFile,myField)
    myRdOnlyDriver.setFieldName(fieldName)
    myRdOnlyDriver.open()

    myWrOnlyDriver = MED_FIELDDOUBLE_WRONLY_DRIVER(medFile2,myField)
    myWrOnlyDriver.open()

    myRdOnlyDriver.close()
    myWrOnlyDriver.close()

    print "Invoking field drivers OK"
except :
    print "there is a problem in invoking field drivers !!"
    print "Please consult the error standart output of the python execution !!"

try:
    myMesh = MESH()

    myRdOnlyDriver = MED_MESH_RDONLY_DRIVER(medFile,myMesh)
    myRdOnlyDriver.setMeshName(meshName)
    myRdOnlyDriver.open()
    myRdOnlyDriver.read()
    myRdOnlyDriver.close()

    myWrOnlyDriver = MED_MESH_WRONLY_DRIVER(medFile2,myMesh)
    myWrOnlyDriver.setMeshName(meshName)
    myWrOnlyDriver.open()
    myWrOnlyDriver.write()
    myWrOnlyDriver.close()

    print "Invoking mesh drivers OK"
except :
    print "there is a problem in invoking mesh drivers !!"
    print "Please consult the error standart output of the python execution !!"

