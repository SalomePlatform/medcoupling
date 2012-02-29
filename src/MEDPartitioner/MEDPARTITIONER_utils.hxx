//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
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
#ifndef MEDPARTITIONER_UTILS_HXX_
#define MEDPARTITIONER_UTILS_HXX_
#include "MEDCouplingUMesh.hxx"

#include <string>
#include <vector>
#include <map>

#ifdef LOCALIZED
#undef LOCALIZED
#endif

#if defined(_DEBUG_) || defined(_DEBUG)
//# define LOCALIZED(message) #message , __FILE__ , __FUNCTION__ , __LINE__
# define LOCALIZED(message) #message , __FUNCTION__ , __LINE__
#else
# define LOCALIZED(message) #message
#endif

namespace MEDPARTITIONER
{
  using namespace std;
  using namespace ParaMEDMEM;
  
  string trim(string& s,const string& drop);
  string intToStr(int i);
  string doubleToStr(double i);
  int strToInt(string s);
  double strToDouble(string s);
  bool testArg(const char *arg, const char *argExpected, string& argValue);
  vector<int> createRandomSize(int size);
  void randomizeAdj(int* xadj, int* adjncy, vector<int>& ran, vector<int>& vx, vector<int>& va);
  void testRandomize();
  
  string reprVectorOfString(const vector<string>& vec);
  string reprVectorOfString(const vector<string>& vec, string sep);
  string reprMapOfStringInt(const map<string,int>& mymap);
  string reprMapOfStringVectorOfString(const map< string,vector<string> >& mymap);
  string reprFieldDescriptions(const vector<string>& vec, string sep);
  
  string serializeFromString(const string& s);
  string serializeFromVectorOfString(const vector<string>& vec);
  vector<string> deserializeToVectorOfString(const string& str);
  string eraseTagSerialized(string fromStr, string tag);
  
  vector<string> vectorizeFromMapOfStringInt(const map<string,int>& mymap);
  map<string,int> devectorizeToMapOfStringInt(const vector<string>& vec);
  
  vector<string> vectorizeFromMapOfStringVectorOfString(const map< string,vector<string> >& mymap);
  map< string,vector<string> > devectorizeToMapOfStringVectorOfString(const vector<string>& vec);
  
  vector<string> selectTagsInVectorOfString(const vector<string>& vec, string tag);
  vector<string> deleteDuplicatesInVectorOfString(const vector<string>& vec);
  map< string,vector<string> > deleteDuplicatesInMapOfStringVectorOfString(const map< string,vector<string> >& mymap);
  
  string cle1ToStr(string s, int inew);
  void cle1ToData(string cle, string& s, int& inew);
  
  string cle2ToStr(string s, int inew, int iold);
  void cle2ToData(string cle, string& s, int& inew, int& iold);
  
  string extractFromDescription(string description, string tag);
  void fieldDescriptionToData(string description,
    int& idomain, string& fileName, string& meshName, string& fieldName,
    int& typeField, int& DT, int& IT);
  void fieldShortDescriptionToData(string description,
    string& fieldName, int& typeField, int& entity, int& DT, int& IT);
  
  ParaMEDMEM::DataArrayInt* createDataArrayIntFromVector(vector<int>& v);
  ParaMEDMEM::DataArrayInt* createDataArrayIntFromVector(vector<int>& v, int nbComponents);
  ParaMEDMEM::DataArrayDouble* createDataArrayDoubleFromVector(vector<double>& v);

  void sendVectorOfString(const vector<string>& vec, const int target);
  vector<string> recvVectorOfString(const int source);
  //TODO void sendrecvVectorOfString(const vector<string>& vec, const int source, const int target);
  vector<string> sendAndReceiveVectorOfString(const vector<string>& vec, const int source, const int target);
  vector<string> allgathervVectorOfString(const vector<std::string>& vec);
  
  vector<string> browseFieldDouble(const MEDCouplingFieldDouble* fd);
  vector<string> browseAllFields(const string& myfile);
  vector<string> browseAllFieldsOnMesh(const string& myfile, const string& mymesh, int idomain);
  vector<string> GetInfosOfField(const char *fileName, const char *meshName, int idomain );

  void sendDoubleVec(const std::vector<double>& vec, int target);
  std::vector<double>* recvDoubleVec(int source);
  void recvDoubleVec(std::vector<double>& vec, int source);
    
  void sendIntVec(const std::vector<int>& vec, int target);
  std::vector<int>* recvIntVec(int source);
  void recvIntVec(std::vector<int>& vec, int source);
  
  void sendDataArrayInt(ParaMEDMEM::DataArrayInt* da, int target);
  ParaMEDMEM::DataArrayInt* recvDataArrayInt(int source);
  void sendDataArrayDouble(ParaMEDMEM::DataArrayDouble* da, int target);
  ParaMEDMEM::DataArrayDouble* recvDataArrayDouble(int source);
  
  ParaMEDMEM::MEDCouplingUMesh* createEmptyMEDCouplingUMesh();
  
  void testVectorOfStringMPI();
  void testMapOfStringIntMPI();
  void testMapOfStringVectorOfStringMPI();
  void testDataArrayMPI();
  void testPersistantMpi0To1(int taille, int nb);
  void testPersistantMpiRing(int taille, int nb);
  void testPersistantMpiRingOnCommSplit(int taille, int nb);

  class MyGlobals
  {
    public : static int _verbose;  //0 to 1000 over 200 is debug
    public : static int _rank;
    public : static int _world_size;
    public : static int _randomize;
    public : static int _atomize;
    public : static int _creates_boundary_faces;
    public : static int _is0verbose; //cout if rank 0 and verbose
    public : static vector<string> _fileNames;    //on [iold]
    public : static vector<string> _meshNames;    //on [iold]
    public : static vector<string> _fieldDescriptions;
    //used for descriptions of components of fields for example...
    public : static vector<string> _generalInformations;
    
  };
  
  /*int MyGlobals::_verbose=0;
  int MyGlobals::_is0verbose=0;
  int MyGlobals::_rank=-1;
  int MyGlobals::_world_size=-1;*/
}
#endif /*MEDPARTITIONER_UTILS_HXX_*/
