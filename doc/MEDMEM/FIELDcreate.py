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

MedFile = "pointe.med"
meshName = "maa1"

myMesh = MESH(MED_DRIVER,MedFile,meshName)

mySupport = myMesh.getSupportOnAll(MED_CELL)

numberOfComponents = 3
myField = FIELDDOUBLE(mySupport,numberOfComponents)
fieldName = "fieldcelldouble"
myField.setName(fieldName)

for i in range(numberOfComponents):
    if (i == 0):
        name = "Vx"
        desc = "vitesse selon x"
    elif (i == 1):
        name = "Vy"
        desc = "vitesse selon y"
    else:
        name = "Vz"
        desc = "vitesse selon z"
    unit = "m. s-1"
    ip1 = i+1
    myField.setComponentName(ip1,name)
    myField.setComponentDescription(ip1,desc)
    myField.setMEDComponentUnit(ip1,unit)

iterationNumber = 10
myField.setIterationNumber(iterationNumber)

orderNumber = 1
myField.setOrderNumber(orderNumber)

time = 3.435678
myField.setTime(time)

numberOfValue = mySupport.getNumberOfElements(MEDMEM_ALL_ELEMENTS)

for i in range(numberOfValue):
    ip1 = i+1
    for j in range(numberOfComponents):
        jp1 = j+1
        value = (ip1+jp1)*0.1
        myField.setValueIJ(ip1,jp1,value)

id = myField.addDriver(MED_DRIVER)
