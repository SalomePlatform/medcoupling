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

medFile = "pointe.med"
medFile2 = "fieldCellDoubleOfpointe.me"
fieldName = "fieldcelldouble"
meshName = "maa1"

try:
    myField = FIEDLDOUBLE()

    myDriver1 = myField->addDriver(MED_DRIVER,medFile,fieldName)
    myField.rmDriver()

    myDriver2 = myField->addDriver(MED_DRIVER,medFile2,fieldName)
    myField.rmDriver(myDriver2)

    myMesh = MESH()
    myDriver3 = myMesh->addDriver(MED_DRIVER,medFile,meshName)
    myMesh.read()
    myMesh.rmDriver()

    myMed = MED()
    myMed.readFileStruct()
    myMed.rmDriver()

except:
    print "There is a problem somewhere !!"
    print "Please consult the error standart output of the python execution !!"
