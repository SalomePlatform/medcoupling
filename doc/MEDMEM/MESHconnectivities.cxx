#include "MEDMEM_Mesh.hxx"
#include "MEDMEM_CellModel.hxx"

int main (int argc, char ** argv) {

  const string MedFile = "pointe.med" ;
  const string MeshName = "maa1" ;
  MESH myMesh(MED_DRIVER,MedFile,MeshName) ;
  myMesh.read() ;

  cout << "Mesh name : " << myMesh.getName()  << endl << endl ; 

  // we get all type for cell entity :
  int NumberOfTypes = myMesh.getNumberOfTypes(MED_CELL) ;
  const CELLMODEL * Types = myMesh.getCellsTypes(MED_CELL) ;

  cout << "Show Connectivity (Nodal) :" << endl ;
  // this example use access with a specified medGeometryElement through
  // CELLMODEL class
  for (int i=0; i<NumberOfTypes; i++) {
    cout << "For type " << Types[i].getName() << " : " << endl ;
    medGeometryElement myType = Types[i].getType() ;
    int NumberOfElements = myMesh.getNumberOfElements(MED_CELL,myType);
    int NomberOfNodesPerCell = Types[i].getNumberOfNodes() ;
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

  return 0 ;
}
