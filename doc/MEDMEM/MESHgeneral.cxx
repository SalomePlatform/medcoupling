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
