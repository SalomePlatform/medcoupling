#include "MEDMEM_Mesh.hxx"

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
