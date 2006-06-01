# Copyright (C) 2005  OPEN CASCADE, CEA, EDF R&D, LEG
#           PRINCIPIA R&D, EADS CCR, Lip6, BV, CEDRAT
# 
from libMEDMEM_Swig import *

MedFile = "pointe.med"
meshName = "maa1"

myMesh = MESH(MED_DRIVER,MedFile,meshName)

name = myMesh.getName()

if (name != meshName) :
    print "Error when reading mesh name : We ask for mesh #",meshName,"#"
    print "and we get mesh #",name
else :
    print "Mesh name : ",name
    spaceDimension = myMesh.getSpaceDimension()
    meshDimension = myMesh.getMeshDimension()
    print "Space Dimension : ",spaceDimension
    print "Mesh Dimension : ",meshDimension
