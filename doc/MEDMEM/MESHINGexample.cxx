#include "MEDMEM_Meshing.hxx"
#include "MEDMEM_Group.hxx"

using namespace std;

int main (int argc, char ** argv) {

  // filename to save the generated MESH
  string filename = "meshing.med" ;

  MESHING myMeshing ;
  myMeshing.setName("meshing") ;

  // define coordinates

  int SpaceDimension = 3 ;
  int NumberOfNodes = 19 ;
  double Coordinates[57] = {
    0.0, 0.0, 0.0, 
    0.0, 0.0, 1.0, 
    2.0, 0.0, 1.0, 
    0.0, 2.0, 1.0, 
    -2.0, 0.0, 1.0, 
    0.0, -2.0, 1.0, 
    1.0, 1.0, 2.0, 
    -1.0, 1.0, 2.0, 
    -1.0, -1.0, 2.0, 
    1.0, -1.0, 2.0, 
    1.0, 1.0, 3.0, 
    -1.0, 1.0, 3.0, 
    -1.0, -1.0, 3.0, 
    1.0, -1.0, 3.0, 
    1.0, 1.0, 4.0, 
    -1.0, 1.0, 4.0, 
    -1.0, -1.0, 4.0, 
    1.0, -1.0, 4.0,
    0.0, 0.0, 5.0
  };

  myMeshing.setCoordinates(SpaceDimension,NumberOfNodes,Coordinates,"CARTESIAN",MED_FULL_INTERLACE);

  string Names[3] = { "X","Y","Z" } ;
  myMeshing.setCoordinatesNames(Names);

  string Units[3] = { "cm","cm","cm" } ;
  myMeshing.setCoordinatesUnits(Units) ;

  // define conectivities

  // cell part
  
  const int NumberOfTypes = 3 ;
  medGeometryElement Types[NumberOfTypes] = {MED_TETRA4,MED_PYRA5,MED_HEXA8} ;
  const int NumberOfElements[NumberOfTypes] = {12,2,2} ;

  myMeshing.setNumberOfTypes(NumberOfTypes,MED_CELL);
  myMeshing.setTypes(Types,MED_CELL);
  myMeshing.setNumberOfElements(NumberOfElements,MED_CELL);

  const int sizeTetra = 12*4 ;
  int ConnectivityTetra[sizeTetra]=
  {
    1,2,3,6,
    1,2,4,3,
    1,2,5,4,
    1,2,6,5,
    2,7,4,3,
    2,8,5,4,
    2,9,6,5,
    2,10,3,6,
    2,7,3,10,
    2,8,4,7,
    2,9,5,8,
    2,10,6,9
  };
  
  myMeshing.setConnectivity(ConnectivityTetra,MED_CELL,MED_TETRA4);

  int ConnectivityPyra[2*5]=
  {
    7,8,9,10,2,
    15,18,17,16,19
  };

  myMeshing.setConnectivity(ConnectivityPyra,MED_CELL,MED_PYRA5);

  int ConnectivityHexa[2*8]=
  {
    11,12,13,14,7,8,9,10,
    15,16,17,18,11,12,13,14
  };

  myMeshing.setConnectivity(ConnectivityHexa,MED_CELL,MED_HEXA8);

  // face part

  const int NumberOfFacesTypes = 2 ;
  medGeometryElement FacesTypes[NumberOfFacesTypes] = {MED_TRIA3,MED_QUAD4} ;
  const int NumberOfFacesElements[NumberOfFacesTypes] = {4,4} ;

  myMeshing.setNumberOfTypes(NumberOfFacesTypes,MED_FACE);
  myMeshing.setTypes(FacesTypes,MED_FACE);
  myMeshing.setNumberOfElements(NumberOfFacesElements,MED_FACE);

  const int sizeTria = 3*4 ;
  int ConnectivityTria[sizeTria]=
  {
    1,4,3,
    1,5,4,
    1,6,5,
    1,3,6
  };
  
  myMeshing.setConnectivity(ConnectivityTria,MED_FACE,MED_TRIA3);

  int ConnectivityQua[4*4]=
  {
    7,8,9,10,
    11,12,13,14,
    11,7,8,12,
    12,8,9,13
  };

  myMeshing.setConnectivity(ConnectivityQua,MED_FACE,MED_QUAD4);

  // edge part

  // not yet implemented : if set, results are unpredictable.

  // Some groups :

  // Node :
  {
    GROUP myGroup ;
    myGroup.setName("SomeNodes");
    myGroup.setMesh(&myMeshing);
    myGroup.setEntity(MED_NODE);
    myGroup.setNumberOfGeometricType(1);
    medGeometryElement myTypes[1] = {MED_NONE};
    myGroup.setGeometricType(myTypes);
    const int myNumberOfElements[1] = {4} ;
    myGroup.setNumberOfElements(myNumberOfElements);
    const int index[1+1] = {1,5} ;
    const int value[4]= { 1,4,5,7} ;
    myGroup.setNumber(index,value);
    
    myMeshing.addGroup(myGroup);
  }
  {
    GROUP myGroup ;
    myGroup.setName("OtherNodes");
    myGroup.setMesh(&myMeshing);
    myGroup.setEntity(MED_NODE);
    myGroup.setNumberOfGeometricType(1);
    medGeometryElement myTypes[1] = {MED_NONE};
    myGroup.setGeometricType(myTypes);
    const int myNumberOfElements[1] = {3} ;
    myGroup.setNumberOfElements(myNumberOfElements);
    const int index[1+1] = {1,4} ;
    const int value[3]= { 2,3,6} ;
    myGroup.setNumber(index,value);
    
    myMeshing.addGroup(myGroup);
  }

  // Cell :
  {
    GROUP myGroup ;
    myGroup.setName("SomeCells");
    myGroup.setMesh(&myMeshing);
    myGroup.setEntity(MED_CELL);
    myGroup.setNumberOfGeometricType(3);
    medGeometryElement myTypes[3] = {MED_TETRA4,MED_PYRA5,MED_HEXA8};
    myGroup.setGeometricType(myTypes);
    const int myNumberOfElements[3] = {4,1,2} ;
    myGroup.setNumberOfElements(myNumberOfElements);
    const int index[3+1] = {1,5,6,8} ;
    const int value[4+1+2]=
    {
      2,7,8,12,
      13,
      15,16
    };
    myGroup.setNumber(index,value);
    
    myMeshing.addGroup(myGroup);
  }
  {
    GROUP myGroup ;
    myGroup.setName("OtherCells");
    myGroup.setMesh(&myMeshing);
    myGroup.setEntity(MED_CELL);
    myGroup.setNumberOfGeometricType(2);
    medGeometryElement myTypes[] = {MED_TETRA4,MED_PYRA5};
    myGroup.setGeometricType(myTypes);
    const int myNumberOfElements[] = {4,1} ;
    myGroup.setNumberOfElements(myNumberOfElements);
    const int index[3+1] = {1,5,6} ;
    const int value[4+1]=
    {
      3,4,5,9,
      14
    };
    myGroup.setNumber(index,value);
    
    myMeshing.addGroup(myGroup);
  }

  // Face :
  {
    GROUP myGroup ;
    myGroup.setName("SomeFaces");
    myGroup.setMesh(&myMeshing);
    myGroup.setEntity(MED_FACE);
    myGroup.setNumberOfGeometricType(2);
    medGeometryElement myTypes[2] = {MED_TRIA3,MED_QUAD4};
    myGroup.setGeometricType(myTypes);
    const int myNumberOfElements[2] = {2,3} ;
    myGroup.setNumberOfElements(myNumberOfElements);
    const int index[2+1] = {1,3,6} ;
    const int value[2+3]=
    {
      2,4,
      5,6,8
    } ;
    myGroup.setNumber(index,value);
    
    myMeshing.addGroup(myGroup);
  }
  {
    GROUP myGroup ;
    myGroup.setName("OtherFaces");
    myGroup.setMesh(&myMeshing);
    myGroup.setEntity(MED_FACE);
    myGroup.setNumberOfGeometricType(1);
    medGeometryElement myTypes[1] = {MED_TRIA3};
    myGroup.setGeometricType(myTypes);
    const int myNumberOfElements[1] = {2} ;
    myGroup.setNumberOfElements(myNumberOfElements);
    const int index[1+1] = {1,3} ;
    const int value[2]=
    {
      1,3
    } ;
    myGroup.setNumber(index,value);
    
    myMeshing.addGroup(myGroup);
  }

  // all rigtht, we save it !

  int id = myMeshing.addDriver(MED_DRIVER,filename,myMeshing.getName());
  myMeshing.write(id) ;

}
