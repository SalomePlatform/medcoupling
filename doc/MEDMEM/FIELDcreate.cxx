using namespace std;
#include "MEDMEM_Mesh.hxx"
#include "MEDMEM_Field.hxx"

using namespace MEDMEM;
using namespace MED_EN ;

int main (int argc, char ** argv) {

  const string MedFile = "pointe.med" ;
  const string MeshName = "maa1" ;

  /* read MESH */
  MESH * myMesh = new MESH(MED_DRIVER,MedFile,MeshName) ;
  //  myMesh->read() ;

  // we need a support :
  SUPPORT * mySupport = new SUPPORT(myMesh,"Support on all CELLs",MED_CELL);

  /* create FIELD on mySupport, with 3 components */
  int NumberOfCompoennts = 3 ;
  FIELD<double> myField(mySupport,NumberOfCompoennts) ;
  const string FieldName = "fieldcelldouble" ;
  myField.setName(FieldName) ;

  // Components information
  string * ComponentsNames = new string[NumberOfCompoennts] ;
  ComponentsNames[0] = "Vx" ;
  ComponentsNames[1] = "Vy" ;
  ComponentsNames[2] = "Vz" ;
  myField.setComponentsNames(ComponentsNames) ;

  string * ComponentsDescriptions = new string[NumberOfCompoennts] ;
  ComponentsDescriptions[0] = "vitesse selon x" ;
  ComponentsDescriptions[1] = "vitesse selon y" ;
  ComponentsDescriptions[2] = "vitesse selon z" ;
  myField.setComponentsDescriptions(ComponentsDescriptions) ;

  string * ComponentsUnits = new string[NumberOfCompoennts] ;
  ComponentsUnits[0] = "m.s-1" ;
  ComponentsUnits[1] = "m.s-1" ;
  ComponentsUnits[2] = "m.s-1" ;
  myField.setMEDComponentsUnits(ComponentsUnits) ;
  
  // Iteration information :
  int IterationNumber = 10 ; // set value to MED_NOPDT if undefined (default)
  myField.setIterationNumber(IterationNumber) ;

  int OrderNumber = 1 ; // set value to MED_NONOR if undefined (default)
  myField.setOrderNumber(OrderNumber) ;

  double Time = 3.435678 ; // in second
  myField.setTime(Time) ;

  // Value :
  int NumberOfValue = mySupport->getNumberOfElements(MED_ALL_ELEMENTS);
  for(int i=1; i<=NumberOfValue; i++) // i^th element
    for (int j=1; j<=NumberOfCompoennts; j++) { // j^th component
      double myValue = (i+j) * 0.1 ;
      myField.setValueIJ(i,j,myValue);
    }
  
  // save this new field
  int id = myField.addDriver(MED_DRIVER) ;

  return 0 ;
}
