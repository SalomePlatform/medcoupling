#include "MEDMEM_Exception.hxx"
#include "MEDMEM_define.hxx"

#include "MEDMEM_Field.hxx"
#include "MEDMEM_Mesh.hxx"
#include "MEDMEM_Med.hxx"

main () {
  
  const char * fileName   = "pointe.med";
  const char * fileName2  = "fieldCellDoubleOfpointe.med";
  const char * fieldName  = "fieldcelldouble";
  const char * meshName   = "maa1";
    
  try {
    // Test creation of drivers from the standard driver method of an object
    FIELD<double> * myField = new FIELD<double>();
    int myDriver1 = myField->addDriver(MED_DRIVER, fileName, fieldName);
    //myField->read();
    //This test failed due to inadequate Support implementation
    myField->rmDriver();  // TESTER LA VALIDITE DE myDriver2 !!!!

    int myDriver2 = myField->addDriver(MED_DRIVER, fileName2, fieldName);
    //myField->write(myDriver2);
    //This test failed due to inadequate Support implementation
    myField->rmDriver(myDriver2);

    MESH * myMesh  = new MESH();
    int myDriver3  = myMesh->addDriver(MED_DRIVER, fileName, meshName);
    myMesh->read();
    myMesh->rmDriver();

    MED  *  myMed  = new MED();
    int myDriver4  = myMed->addDriver(MED_DRIVER, fileName);
    myMed->readFileStruct();
    myMed->rmDriver();

    delete myField;
    delete myMesh;
    delete myMed;

  } catch (MEDEXCEPTION& ex){
    MESSAGE(ex.what()) ;
  }
}
