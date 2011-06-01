// Copyright (C) 2007-2011  CEA/DEN, EDF R&D, OPEN CASCADE
//
// Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
// CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

#include "MEDMEM_Mesh.hxx"

using namespace MEDMEM ;
using namespace MED_EN ;

int main (int argc, char ** argv) {

  const string MedFile = "pointe.med" ;
  const string MeshName = "maa1" ;
  MESH myMesh(MED_DRIVER,MedFile,MeshName) ;

  cout << "Mesh name : " << myMesh.getName()  << endl << endl ; 

  int SpaceDimension = myMesh.getSpaceDimension() ;
  int NumberOfNodes = myMesh.getNumberOfNodes() ;
  cout << "Space dimension  : " << SpaceDimension << endl << endl ; 
  cout << "Number of nodes  : " << NumberOfNodes  << endl << endl ; 

  cout << "Show Nodes Coordinates : " << endl ;

  // coordinates names :
  cout << "Name :" << endl ;
  const string * CoordinatesNames = myMesh.getCoordinatesNames() ;
  for(int i=0; i<SpaceDimension ; i++) {
    cout << " - " << CoordinatesNames[i] << endl ;
  }
  // coordinates unit :
  cout << "Unit :" << endl ;
  const string * CoordinatesUnits = myMesh.getCoordinatesUnits() ;
  for(int i=0; i<SpaceDimension ; i++) {
    cout << " - " << CoordinatesUnits[i] << endl ;
  }
  // coordinates value
  const double * Coordinates = 
    myMesh.getCoordinates(MED_FULL_INTERLACE) ;
  for(int i=0; i<NumberOfNodes ; i++) {
    cout << "Nodes " << i+1 << " : " ;
    for (int j=0; j<SpaceDimension ; j++)
      cout << Coordinates[i*SpaceDimension+j] << " " ;
    cout << endl ;
  }
  
  return 0 ;
}
