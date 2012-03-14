// Copyright (C) 2007-2011  CEA/DEN, EDF R&D
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

#ifndef __MEDPARTITIONER_UTILS_HXX__
#define __MEDPARTITIONER_UTILS_HXX__

#include "MEDCouplingUMesh.hxx"

#include <string>
#include <vector>
#include <map>

//# define LOCALIZED(message) #message , __FILE__ , __FUNCTION__ , __LINE__

namespace MEDPARTITIONER
{
  std::string Trim(const std::string& s,const std::string& drop);
  std::string IntToStr(const int i);
  std::string DoubleToStr(const double i);
  int StrToInt(const std::string& s);
  double StrToDouble(const std::string& s);
  bool TestArg(const char *arg, const char *argExpected, std::string& argValue);
  std::vector<int> CreateRandomSize(const int size);
  void RandomizeAdj(int* xadj, int* adjncy, std::vector<int>& ran, std::vector<int>& vx, std::vector<int>& va);
  void TestRandomize();
  
  std::string ReprVectorOfString(const std::vector<std::string>& vec);
  std::string ReprVectorOfString(const std::vector<std::string>& vec, const std::string separator);
  std::string ReprMapOfStringInt(const std::map<std::string,int>& mymap);
  std::string ReprMapOfStringVectorOfString(const std::map< std::string,std::vector<std::string> >& mymap);
  std::string ReprFieldDescriptions(const std::vector<std::string>& vec,const  std::string separator);
  
  std::string SerializeFromString(const std::string& s);
  std::string SerializeFromVectorOfString(const std::vector<std::string>& vec);
  std::vector<std::string> DeserializeToVectorOfString(const std::string& str);
  std::string EraseTagSerialized(const std::string& fromStr, const std::string& tag);
  
  std::vector<std::string> VectorizeFromMapOfStringInt(const std::map<std::string,int>& mymap);
  std::map<std::string,int> DevectorizeToMapOfStringInt(const std::vector<std::string>& vec);
  
  std::vector<std::string> VectorizeFromMapOfStringVectorOfString(const std::map< std::string,std::vector<std::string> >& mymap);
  std::map< std::string,std::vector<std::string> > DevectorizeToMapOfStringVectorOfString(const std::vector<std::string>& vec);
  
  std::vector<std::string> SelectTagsInVectorOfString(const std::vector<std::string>& vec, const std::string tag);
  std::vector<std::string> DeleteDuplicatesInVectorOfString(const std::vector<std::string>& vec);
  std::map< std::string,std::vector<std::string> > DeleteDuplicatesInMapOfStringVectorOfString(const std::map< std::string,std::vector<std::string> >& mymap);
  
  std::string Cle1ToStr(const std::string& s, const int inew);
  void Cle1ToData(const std::string& cle, std::string& s, int& inew);
  
  std::string Cle2ToStr(const std::string& s,const int inew,const int iold);
  void Cle2ToData(const std::string& cle, std::string& s, int& inew, int& iold);
  
  std::string ExtractFromDescription(const std::string& description,const  std::string& tag);
  void FieldDescriptionToData(const std::string& description,
                              int& idomain, std::string& fileName, std::string& meshName, std::string& fieldName,
                              int& typeField, int& DT, int& IT);
  void FieldShortDescriptionToData(const std::string& description,
                                   std::string& fieldName, int& typeField, int& entity, int& DT, int& IT);
  
  ParaMEDMEM::DataArrayInt *CreateDataArrayIntFromVector(const std::vector<int>& v);
  ParaMEDMEM::DataArrayInt *CreateDataArrayIntFromVector(const std::vector<int>& v, const int nbComponents);
  ParaMEDMEM::DataArrayDouble *CreateDataArrayDoubleFromVector(const std::vector<double>& v);

  //not adviced, interblocking, use sendAndReceive
  //void SendVectorOfString(const std::vector<std::string>& vec, const int target);
  //std::vector<std::string> RecvVectorOfString(const int source);
  //TODO void sendRecvVectorOfString(const std::vector<std::string>& vec, const int source, const int target);
  std::vector<std::string> SendAndReceiveVectorOfString(const std::vector<std::string>& vec, const int source, const int target);
  std::vector<std::string> AllgathervVectorOfString(const std::vector<std::string>& vec);
  
  std::vector<std::string> BrowseFieldDouble(const ParaMEDMEM::MEDCouplingFieldDouble* fd);
  std::vector<std::string> BrowseAllFields(const std::string& myfile);
  std::vector<std::string> BrowseAllFieldsOnMesh(const std::string& myfile, const std::string& mymesh, const int idomain);
  std::vector<std::string> GetInfosOfField(const char *fileName, const char *meshName, const int idomain );

  void SendDoubleVec(const std::vector<double>& vec, const int target);
  std::vector<double> *RecvDoubleVec(const int source);
  void RecvDoubleVec(std::vector<double>& vec, const int source);
    
  void SendIntVec(const std::vector<int>& vec, const int target);
  std::vector<int>* RecvIntVec(int source);
  void RecvIntVec(std::vector<int>& vec, const int source);
  
  void SendDataArrayInt(const ParaMEDMEM::DataArrayInt* da, const int target);
  ParaMEDMEM::DataArrayInt *RecvDataArrayInt(const int source);
  void SendDataArrayDouble(const ParaMEDMEM::DataArrayDouble* da, const int target);
  ParaMEDMEM::DataArrayDouble *RecvDataArrayDouble(const int source);
  
  ParaMEDMEM::MEDCouplingUMesh *CreateEmptyMEDCouplingUMesh();
  
  void TestVectorOfStringMpi();
  void TestMapOfStringIntMpi();
  void TestMapOfStringVectorOfStringMpi();
  void TestDataArrayMpi();
  void TestPersistantMpi0To1(int taille, int nb);
  void TestPersistantMpiRing(int taille, int nb);
  void TestPersistantMpiRingOnCommSplit(int taille, int nb);

  class MyGlobals
  {
  public :
    static int _Verbose;  //0 to 1000 over 200 is debug
    static int _Rank;
    static int _World_Size;
    static int _Randomize;
    static int _Atomize;
    static int _Creates_Boundary_Faces;
    static int _Is0verbose; //trace cout if rank 0 and verbose
    static std::vector<std::string> _File_Names;    //on [iold]
    static std::vector<std::string> _Mesh_Names;    //on [iold]
    static std::vector<std::string> _Field_Descriptions;
    /*! used for descriptions of components of fields for example...*/
    static std::vector<std::string> _General_Informations;
  };
}
#endif
