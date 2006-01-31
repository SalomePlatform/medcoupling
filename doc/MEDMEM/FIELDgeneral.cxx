// Copyright (C) 2005  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
// CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either 
// version 2.1 of the License.
// 
// This library is distributed in the hope that it will be useful 
// but WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public  
// License along with this library; if not, write to the Free Software 
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
// See http://www.salome-platform.org/
//
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
  const double * Value = myField.getValue();
  for(int i=0; i<NumberOfValue; i++) {
    for(int j=0; j<NumberOfCompoennts; j++)
      cout << Value[i*NumberOfCompoennts+j] << " " ;
    cout << endl ;
  }

  delete mySupport;
  delete myMesh;

  return 0 ;
}
