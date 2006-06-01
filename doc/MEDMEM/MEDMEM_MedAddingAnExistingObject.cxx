// Copyright (C) 2005  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
// CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
// 
#include "MEDMEM_Exception.hxx"
#include "MEDMEM_define.hxx"

#include "MEDMEM_Field.hxx"
#include "MEDMEM_Mesh.hxx"
#include "MEDMEM_Med.hxx"

main () {
  
  const char * fileName   = "pointe.med";
  const char * fileName2  = "Field&MeshGeneratedPointe.med";
  const char * fileName3  = "MedGeneratedPointe.med";
  const char * fieldName1  = "fieldcelldouble";
  const char * fieldName2  = "fieldcelldoublebis";
  const char * meshName1   = "maa1";
  const char * meshName2   = "maa1bis";

  try {

    // FAIRE LE TEST AVEC LES CHAMPS AUSSI !.

    MESH myMesh(MED_DRIVER,fileName,meshName1);
    myMesh.setName(meshName2);
    myMesh.rmDriver();

    MED  myMed(MED_DRIVER,fileName);
    myMed.read();
    myMed.addMesh(&myMesh);
    int myMedDriver = myMed.addDriver(MED_DRIVER,fileName3);
    myMed.write(myMedDriver);

    // FAIRE LE TEST AVEC LES CHAMPS AUSSI !.

  } catch (MEDEXCEPTION& ex){
    MESSAGE(ex.what()) ;
  }
}
