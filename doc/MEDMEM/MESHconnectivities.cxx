//  Copyright (C) 2007-2011  CEA/DEN, EDF R&D, OPEN CASCADE
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

#include "MEDMEM_Mesh.hxx"

using namespace MEDMEM ;
using namespace MED_EN ;

int main (int argc, char ** argv) {

//   const string MedFile = "polyedres.med" ;
//   const string MeshName = "Erreur orientation" ;
//   const string MedFile = "polygones.med" ;
//   const string MeshName = "Bord" ;
  const string MedFile = "pointe.med" ;
  const string MeshName = "maa1" ;
  MESH myMesh(MED_DRIVER,MedFile,MeshName) ;
  myMesh.read() ;

  cout << "Mesh name : " << myMesh.getName()  << endl << endl ; 

  // we get all type for cell entity :
  int NumberOfTypes = myMesh.getNumberOfTypes(MED_CELL) ;
  cout << "Show Connectivity (Nodal) :" << endl ;
  // this example use access with a specified medGeometryElement array
  const medGeometryElement * Types = myMesh.getTypes(MED_CELL);
  string * cellTypeNames =  myMesh.getCellTypeNames(MED_CELL);
  for (int i=0; i<NumberOfTypes; i++) {
    cout << "For type " << cellTypeNames[i] << " : " << endl ;
    medGeometryElement myType = Types[i] ;
    int NumberOfElements = myMesh.getNumberOfElements(MED_CELL,myType);
    int NomberOfNodesPerCell = Types[i]%100 ;
    const int * Connectivity = 
      myMesh.getConnectivity(MED_FULL_INTERLACE,
			     MED_NODAL,
			     MED_CELL,
			     myType);
    for (int j=0; j<NumberOfElements; j++){
      cout << "Element "<< j+1 <<" : " ;
      for (int k=0; k<NomberOfNodesPerCell; k++)
	cout << Connectivity[j*NomberOfNodesPerCell+k]<<" ";
      cout << endl ;
    }
  }
  cout << "Show Reverse Nodal Connectivity :" << endl ;
  // this example use global access with index array
  int NumberOfNodes = myMesh.getNumberOfNodes() ;
  const int * ReverseNodalConnectivity = 
    myMesh.getReverseConnectivity(MED_NODAL) ;
  const int * ReverseNodalConnectivityIndex = 
    myMesh.getReverseConnectivityIndex(MED_NODAL) ;
  for (int i=0; i<NumberOfNodes; i++) {
    cout << "Node "<<i+1<<" : " ;
    int IndexBegin = ReverseNodalConnectivityIndex[i] ;
    int IndexEnd = ReverseNodalConnectivityIndex[i+1] ;
    for (int j=IndexBegin; j<IndexEnd; j++)
      // Index value begin at 1 so use j-1
      cout << ReverseNodalConnectivity[j-1] << " " ; 
    cout << endl ;
  }

  cout << "Show Connectivity (Descending) :" << endl ;
  // this example use global access with index array
  int NumberOfElements = myMesh.getNumberOfElements(MED_CELL,MED_ALL_ELEMENTS);
  const int * DescendingConnectivity =  
    myMesh.getConnectivity(MED_FULL_INTERLACE,
			   MED_DESCENDING,
			   MED_CELL,
			   MED_ALL_ELEMENTS);
  const int * DescendingConnectivityIndex =
    myMesh.getConnectivityIndex(MED_DESCENDING,MED_CELL);
  for (int i=0; i<NumberOfElements; i++) {
    cout << "Element "<<i+1<<" : " ;
    int IndexBegin = DescendingConnectivityIndex[i] ;
    int IndexEnd = DescendingConnectivityIndex[i+1] ;
    for (int j=IndexBegin; j<IndexEnd; j++)
      // Index value begin at 1 so use j-1
      cout << DescendingConnectivity[j-1] << " " ;
    cout << endl ;
  }

  cout << "Show Reverse Descending Connectivity :" << endl ;
  // this example use global access with Index array
  const int * ReverseDescendingConnectivity = 
    myMesh.getReverseConnectivity(MED_DESCENDING) ;
  const int * ReverseDescendingConnectivityIndex = 
    myMesh.getReverseConnectivityIndex(MED_DESCENDING) ;

  int MeshDimension = myMesh.getMeshDimension() ;
  int NumberOfConstituents  = 0;
  string Constituent ;
  medEntityMesh ConstituentEntity ;
  // test if we have face (3D) or edge (2D)
  if (MeshDimension==3) {
    Constituent = "Face" ;
    ConstituentEntity = MED_FACE ;
  }
  if (MeshDimension==2) {
    Constituent = "Edge" ;
    ConstituentEntity = MED_EDGE ;
  }

  NumberOfConstituents = 
    myMesh.getNumberOfElements(ConstituentEntity,MED_ALL_ELEMENTS);
  
  if (MeshDimension==1) {
    MESSAGE("ERROR : MeshDimension = 1 !");
    MESSAGE("We could not see Reverse Descending Connectivity.") ;
  } else {
    for (int i=0; i<NumberOfConstituents; i++) {
      cout << Constituent << " " << i+1 << " : " ;
      int IndexBegin = ReverseDescendingConnectivityIndex[i] ;
      int IndexEnd = ReverseDescendingConnectivityIndex[i+1] ;
      for (int j=IndexBegin;j<IndexEnd;j++)
      // Index value begin at 1 so use j-1
	cout << ReverseDescendingConnectivity[j-1] << " " ;
      cout << endl ;
    }
  }
  cout << "Show "<< Constituent <<" Connectivity (Nodal) :" << endl ;
  // this example use global access with index array
  const int * ConstituentConnectivity =  
    myMesh.getConnectivity(MED_FULL_INTERLACE,
			   MED_NODAL,
			   ConstituentEntity,
			   MED_ALL_ELEMENTS);
  const int * ConstituentConnectivityIndex =
    myMesh.getConnectivityIndex(MED_NODAL,ConstituentEntity);
  for (int i=0; i<NumberOfConstituents; i++) {
    cout << Constituent << " " << i+1 << " : " ;
    int IndexBegin = ConstituentConnectivityIndex[i] ;
    int IndexEnd = ConstituentConnectivityIndex[i+1] ;
    for (int j=IndexBegin; j<IndexEnd; j++)
      // Index value begin at 1 so use j-1
      cout << ConstituentConnectivity[j-1]<<" ";
    cout << endl ;
  }

  int nbPolygons = myMesh.getNumberOfPolygons();
  if ( nbPolygons > 0 )
  {
    cout << "Show Connectivity (Nodal) of POLYGONS:" << endl ;
    const int* Connectivity = myMesh.getPolygonsConnectivity(MED_NODAL,MED_CELL);
    const int* ConnectivityIndex = myMesh.getPolygonsConnectivityIndex(MED_NODAL,MED_CELL);
    for (int j=0; j<nbPolygons; j++){
      cout << "Polygon "<< j+1 <<" : " ;
    int IndexBegin = ConnectivityIndex[j];
    int IndexEnd   = ConnectivityIndex[j+1];
      for (int k=IndexBegin; k<IndexEnd; k++)
	cout << Connectivity[k-1]<<" ";
      cout << endl ;
    }
  }

  int nbPolyhedrons = myMesh.getNumberOfPolyhedron();
  if ( nbPolyhedrons > 0 )
  {
    cout << "Show Connectivity (Nodal) of POLYHEDRONS:" << endl ;
    const int* Connectivity = myMesh.getPolyhedronConnectivity(MED_NODAL);
    const int* FaceIndex    = myMesh.getPolyhedronFacesIndex();
    const int* Index        = myMesh.getPolyhedronIndex(MED_NODAL);
    for (int j=0; j<nbPolyhedrons; j++){
      cout << "Polyhedron "<< j+1 <<" : " << endl;
      int FaceIndexBegin = Index[j];
      int FaceIndexEnd = Index[j+1];
      for (int k=FaceIndexBegin; k<FaceIndexEnd; k++) {
        cout << "  Face " << k - FaceIndexBegin + 1 << " : ";
        int IndexBegin = FaceIndex[k-1];
        int IndexEnd   = FaceIndex[k];
        for (int i=IndexBegin; i<IndexEnd; i++)
          cout << Connectivity[i-1]<<" ";
        cout << endl ;
      }
    }
  }

  return 0 ;
}
