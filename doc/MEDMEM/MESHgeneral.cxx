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

using namespace std;
#include "MEDMEM_Mesh.hxx"

using namespace MEDMEM ;

int main (int argc, char ** argv) {

  const string MedFile = "pointe.med" ;
  const string MeshName = "maa1" ;
  
  // create a MESH object by reading it on file :
  MESH myMesh(MED_DRIVER,MedFile,MeshName) ;

  string Name = myMesh.getName() ;
  if (Name != MeshName) {
    cout << "Error when reading mesh name : We ask for mesh #"
         << MeshName <<"# and we get mesh #"<< Name <<"#"<< endl << endl ;
    return -1;
  }

  cout << "Mesh name : " << Name  << endl << endl ; 

  int SpaceDimension = myMesh.getSpaceDimension() ;
  int MeshDimension = myMesh.getMeshDimension() ;
  
  cout << "Space Dimension : " << SpaceDimension << endl << endl ; 
  cout << "Mesh Dimension : " << MeshDimension << endl << endl ; 

  return 0 ;
}
