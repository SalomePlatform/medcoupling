using namespace std;
#include "MEDMEM_Mesh.hxx"
#include "MEDMEM_Field.hxx"

using namespace MEDMEM;
using namespace MED_EN ;

int main (int argc, char ** argv) {

  const string MedFile = "pointe.med" ;
  const string MeshName = "maa1" ;
  const string FieldName = "fieldcelldoublevector" ;

  /* read MESH */
  MESH * myMesh = new MESH(MED_DRIVER,MedFile,MeshName) ;
  //  myMesh->read() ;

  /* read FIELD */
  // we need a support :
  SUPPORT * mySupport = new SUPPORT(myMesh,"Support on all Cells",MED_CELL);
  FIELD<double> myField(mySupport,MED_DRIVER,MedFile,FieldName) ;
  //  myField.read() ;

  /* what in Field ? */
  // How many components
  int NumberOfCompoennts = myField.getNumberOfComponents() ;

  const string * ComponentsNames = myField.getComponentsNames();
  const string * ComponentsDescriptions = myField.getComponentsDescriptions();
  const string * ComponentsUnits = myField.getMEDComponentsUnits();

  for(int i=0;i<NumberOfCompoennts; i++) {
    cout << "Component " << i << " :" <<endl ;
    cout << "  - name        : " << ComponentsNames[i] << endl ;
    cout << "  - description : " << ComponentsDescriptions[i] << endl ;
    cout << "  - unit        : " << ComponentsUnits[i] << endl ;
  }

  // Which iteration :
  int IterationNumber = myField.getIterationNumber() ; // negative mean undefined
  int OrderNumber = myField.getOrderNumber() ;
  // internal iteration at this time iteration, negative mean undefined
  double Time = myField.getTime() ;

  cout << "Iteration " << IterationNumber << " at time " << Time <<
    " (and order number " << OrderNumber << ")" << endl ;

  // How many Value :
  int NumberOfValue = mySupport->getNumberOfElements(MED_ALL_ELEMENTS);
  // Value
  const double * Value = myField.getValue(MED_FULL_INTERLACE);
  for(int i=0; i<NumberOfValue; i++) {
    for(int j=0; j<NumberOfCompoennts; j++)
      cout << Value[i*NumberOfCompoennts+j] << " " ;
    cout << endl ;
  }

  delete mySupport;
  delete myMesh;

  return 0 ;
}
