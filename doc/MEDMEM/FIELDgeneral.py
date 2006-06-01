# Copyright (C) 2005  OPEN CASCADE, CEA, EDF R&D, LEG
#           PRINCIPIA R&D, EADS CCR, Lip6, BV, CEDRAT
# 
######################################################################
#                                                                    #
# This Python script should be executed when the shared library is   #
# generated using SWIG 1.3 (or higher) due to the fact that older    #
# version could not handle the wrapping of several class constructor #
#                                                                    #
######################################################################
from libMEDMEM_Swig import *

MedFile = "pointe.med"
meshName = "maa1"
fieldName = "fieldcelldouble"

myMesh = MESH(MED_DRIVER,MedFile,meshName)

mySupport = SUPPORT(myMesh,"Support on CELLs",MED_CELL)

myField = FIELDDOUBLE(mySupport,MED_DRIVER,MedFile,fieldName)

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
value = myField.getValue(MED_FULL_INTERLACE)

for i in range(numberOfValue):
    print "  * ",value[i*numberOfComponents:(i+1)*numberOfComponents]
