######################################################################
#                                                                    #
# This Python script should be executed when the shared library is   #
# generated using SWIG 1.3 (or higher) due to the fact that older    #
# version could not handle the wrapping of several class constructor #
#                                                                    #
######################################################################
from libMEDMEM_Swig import *

medFile = "pointe.med"
fieldName = "fieldcelldouble"
meshName = "maa1"

try:
    print "Creation of MESH object"
    myMesh = MESH(MED_DRIVER,medFile,meshName)

    print "Creation of MED object"
    myMed = MED(MED_DRIVER,medFile)

    print "Test the driver removal dor MESH"
    myMesh.rmDriver()

    print "Test the driver removal dor MED"
    myMed.rmDriver()

    print "End of Python script"
    
except:
    print "There is a problem somewhere !!"
    print "Consult the error standart output of the python execution !!"
