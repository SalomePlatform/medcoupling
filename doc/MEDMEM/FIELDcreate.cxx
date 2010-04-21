//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D, OPEN CASCADE
//
//  Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
//  CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

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
