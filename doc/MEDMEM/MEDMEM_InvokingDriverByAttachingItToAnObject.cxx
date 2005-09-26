#include "MEDMEM_Exception.hxx"
#include "MEDMEM_define.hxx"

#include "MEDMEM_Field.hxx"
#include "MEDMEM_Mesh.hxx"
#include "MEDMEM_Med.hxx"
#include "MEDMEM_MedMedDriver.hxx"
#include "MEDMEM_MedMeshDriver.hxx"

using namespace MEDMEM ;
using namespace MED_EN ;

main () {
  
  const char * fileName   = "pointe.med";
  const char * fileName2  = "Field&MeshGeneratedPointe.med";
  const char * fileName3  = "MedGeneratedPointe.med";
  const char * fieldName  = "fieldcelldouble";
  const char * meshName   = "maa1";

  try {
    // Test creation of drivers from the standard driver method of an object
    {  
      FIELD<double> * myField = new FIELD<double>();
      MED_FIELD_RDONLY_DRIVER<double> myRdOnlyDriver(fileName,myField);
      myRdOnlyDriver.setFieldName(fieldName);
      myRdOnlyDriver.open(); 
      //This test failed due to inadequate Support implementation   
      // myRdOnlyDriver.read();
      // try  { myRdOnlyDriver.write(); } catch  (MEDEXCEPTION& ex) 
      //   { MESSAGE(ex.what()); }
      MED_FIELD_WRONLY_DRIVER<double> myWrOnlyDriver(fileName2,myField);
      myWrOnlyDriver.open(); 
      //This test failed due to inadequate Support implementation   
      // myWrOnlyDriver.write(); 
      // try  myWrOnlyDriver.read(); catch  (MEDEXCEPTION& ex) 
      //   { MESSAGE(ex.what()); }
      myRdOnlyDriver.close();
      myWrOnlyDriver.close();
      delete myField;
    }

    {
      MESH * myMesh = new MESH();
      MED_MESH_RDONLY_DRIVER myRdOnlyDriver(fileName,myMesh);
      myRdOnlyDriver.setMeshName(meshName);
      myRdOnlyDriver.open(); 
      myRdOnlyDriver.read();
      myRdOnlyDriver.close(); 
      // try  { myRdOnlyDriver.write(); } catch  (MEDEXCEPTION& ex)
      //   { MESSAGE(ex.what()); }
      MED_MESH_WRONLY_DRIVER myWrOnlyDriver(fileName2,myMesh);
      myWrOnlyDriver.setMeshName(meshName);
      myWrOnlyDriver.open(); 
      myWrOnlyDriver.write(); 
      // try  myWrOnlyDriver.read(); catch  (MEDEXCEPTION& ex)
      //   { MESSAGE(ex.what()); }
      // myRdOnlyDriver.close(); 
      //While we use H5close() in the MESH/FIELD drivers, the next
      //line will fail, because all files are previously closed !
      myWrOnlyDriver.close();
      delete myMesh;
    }

    {
      MED * myMed = new MED();
      MED_MED_RDONLY_DRIVER myRdOnlyDriver(fileName,myMed);
      myRdOnlyDriver.open(); 
      myRdOnlyDriver.readFileStruct();
      myRdOnlyDriver.close(); 
      myMed->updateSupport(); // DOIT ETRE SUPPRIMEE
      //      myRdOnlyDriver.read();
      // try { myRdOnlyDriver.write(); } catch  (MEDEXCEPTION& ex) 
      //   { MESSAGE(ex.what()); }
      //MED_MED_WRONLY_DRIVER myWrOnlyDriver(fileName3,myMed);
      //myWrOnlyDriver.open(); 
      //myWrOnlyDriver.write(); // Not implemented yet.
      //myWrOnlyDriver.close();
      delete myMed;
    }

  } catch (MEDEXCEPTION& ex){
    cout << "MAIN BLOCK EXCEPTION" << endl;
    MESSAGE(ex.what()) ;
  }
}
