######################################################################
#                                                                    #
# This Python script should be executed when the shared library is   #
# generated using SWIG 1.3 (or higher) due to the fact that older    #
# version could not handle the wrapping of several class constructor #
#                                                                    #
######################################################################
from libMEDMEM_Swig import *

medFile = "pointe.med"
medFile2 = "Field&MeshGeneratedPointe.med"
fieldName = "fieldcelldouble"
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

    myWrOnlyDriver = MED_MESH_WRONLY_DRIVER(medFile,myMesh)
    myWrOnlyDriver.setMeshName(meshName)
    myWrOnlyDriver.open()
    myWrOnlyDriver.write()

    myWrOnlyDriver.close()

    print "Invoking mesh drivers OK"
except :
    print "there is a problem in invoking mesh drivers !!"
    print "Please consult the error standart output of the python execution !!"

try:
    myMed = MED()
    myRdOnlyDriver = MED_MED_RDONLY_DRIVER(medFile,myMed)
    myRdOnlyDriver.open() 
    myRdOnlyDriver.readFileStruct()
    myRdOnlyDriver.close()
    myMed.updateSupport()

    print "Invoking Med drivers OK"
except :
    print "There is a problem in invoking MED drivers !!"
    print "Please consult the error standart output of the python execution !!"
