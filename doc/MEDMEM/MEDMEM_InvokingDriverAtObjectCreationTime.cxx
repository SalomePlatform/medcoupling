#include "MEDMEM_Exception.hxx"
#include "MEDMEM_define.hxx"

#include "MEDMEM_Field.hxx"
#include "MEDMEM_Mesh.hxx"
#include "MEDMEM_Med.hxx"

using namespace MEDMEM ;
using namespace MED_EN ;

main () {
  
  const char * fileName  = "pointe.med";
  const char * fieldName = "fieldcelldouble";
  const char * meshName  = "maa1";
  
  try {
    
    // Test creation of drivers at object Creation time

    //This test failed due to inadequate Support implementation   
    // FIELD<double> myField (MED_DRIVER,fileName,fieldName); 
    MESH          myMesh  (MED_DRIVER,fileName,meshName);
    MED           myMed   (MED_DRIVER,fileName);

    // Test removal of drivers
    //myField.rmDriver();
    myMesh.rmDriver ();
    myMed.rmDriver  ();

  } catch (MEDEXCEPTION& ex){
    MESSAGE(ex.what()) ;
  }
}
