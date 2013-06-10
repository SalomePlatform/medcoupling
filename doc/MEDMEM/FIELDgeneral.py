#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2007-2013  CEA/DEN, EDF R&D, OPEN CASCADE
#
# Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
# CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
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

######################################################################
# This Python script should be executed when the shared library is   #
# generated using SWIG 1.3 (or higher) due to the fact that older    #
# version could not handle the wrapping of several class constructor #
######################################################################
#
from libMEDMEM_Swig import *

MedFile = "pointe.med"
meshName = "maa1"
fieldName = "fieldcelldoublescalar"

myMesh = MESH(MED_DRIVER,MedFile,meshName)

mySupport = myMesh.getSupportOnAll(MED_CELL)

myField = FIELDDOUBLE(mySupport,MED_DRIVER,MedFile,fieldName,-1,-1)

numberOfComponents = myField.getNumberOfComponents()

for i in range(numberOfComponents):
    ip1 = i+1
    name = myField.getComponentName(ip1)
    desc = myField.getComponentDescription(ip1)
    unit = myField.getMEDComponentUnit(ip1)

    print "Component ",ip1
    print "  - name       : ",name
    print "  - decription : ",desc
    print "  - unit       : ", unit

iterationNumber = myField.getIterationNumber()
orderNumber = myField.getOrderNumber()
time = myField.getTime()
print "Iteration ",iterationNumber,"  at time ",time,\
      " (and order number ",orderNumber,")"

numberOfValue = mySupport.getNumberOfElements(MED_ALL_ELEMENTS)
value = myField.getValue()

for i in range(numberOfValue):
    print "  * ",value[i*numberOfComponents:(i+1)*numberOfComponents]
