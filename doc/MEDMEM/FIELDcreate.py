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

myMesh = MESH(MED_DRIVER,MedFile,meshName)

mySupport = SUPPORT(myMesh,"Support on all CELLs",MED_CELL)

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

numberOfValue = mySupport.getNumberOfElements(MED_ALL_ELEMENTS)

for i in range(numberOfValue):
    ip1 = i+1
    for j in range(numberOfComponents):
        jp1 = j+1
        value = (ip1+jp1)*0.1
        myField.setValueIJ(ip1,jp1,value)

id = myField.addDriver(MED_DRIVER)
