// Copyright (C) 2007-2016  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
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

#include "MEDPARTITIONER.hxx"

#include "MEDCouplingUMesh.hxx"
#include "BBTree.txx"

#include <string>
#include <vector>
#include <map>

//# define LOCALIZED(message) #message , __FILE__ , __FUNCTION__ , __LINE__

namespace MEDPARTITIONER
{
  MEDPARTITIONER_EXPORT std::string Trim(const std::string& s,const std::string& drop);
  MEDPARTITIONER_EXPORT std::string IntToStr(const int i);
  MEDPARTITIONER_EXPORT std::string DoubleToStr(const double i);
  MEDPARTITIONER_EXPORT int StrToInt(const std::string& s);
  MEDPARTITIONER_EXPORT double StrToDouble(const std::string& s);
  MEDPARTITIONER_EXPORT bool TestArg(const char *arg, const char *argExpected, std::string& argValue);
  MEDPARTITIONER_EXPORT std::vector<int> CreateRandomSize(const int size);
  MEDPARTITIONER_EXPORT void RandomizeAdj(int* xadj, int* adjncy, std::vector<int>& ran, std::vector<int>& vx, std::vector<int>& va);
  MEDPARTITIONER_EXPORT void TestRandomize();
                       
  MEDPARTITIONER_EXPORT std::string ReprVectorOfString(const std::vector<std::string>& vec);
  MEDPARTITIONER_EXPORT std::string ReprVectorOfString(const std::vector<std::string>& vec, const std::string separator);
  MEDPARTITIONER_EXPORT std::string ReprMapOfStringInt(const std::map<std::string,int>& mymap);
  MEDPARTITIONER_EXPORT std::string ReprMapOfStringVectorOfString(const std::map< std::string,std::vector<std::string> >& mymap);
  MEDPARTITIONER_EXPORT std::string ReprFieldDescriptions(const std::vector<std::string>& vec,const  std::string separator);
  
  MEDPARTITIONER_EXPORT std::string SerializeFromString(const std::string& s);
  MEDPARTITIONER_EXPORT std::string SerializeFromVectorOfString(const std::vector<std::string>& vec);
  MEDPARTITIONER_EXPORT std::vector<std::string> DeserializeToVectorOfString(const std::string& str);
  MEDPARTITIONER_EXPORT std::string EraseTagSerialized(const std::string& fromStr, const std::string& tag);
  
  MEDPARTITIONER_EXPORT std::vector<std::string> VectorizeFromMapOfStringInt(const std::map<std::string,int>& mymap);
  MEDPARTITIONER_EXPORT std::map<std::string,int> DevectorizeToMapOfStringInt(const std::vector<std::string>& vec);
  
  MEDPARTITIONER_EXPORT std::vector<std::string> VectorizeFromMapOfStringVectorOfString(const std::map< std::string,std::vector<std::string> >& mymap);
  MEDPARTITIONER_EXPORT std::map< std::string,std::vector<std::string> > DevectorizeToMapOfStringVectorOfString(const std::vector<std::string>& vec);
  
  MEDPARTITIONER_EXPORT std::vector<std::string> SelectTagsInVectorOfString(const std::vector<std::string>& vec, const std::string tag);
  MEDPARTITIONER_EXPORT std::vector<std::string> DeleteDuplicatesInVectorOfString(const std::vector<std::string>& vec);
  MEDPARTITIONER_EXPORT std::map< std::string,std::vector<std::string> > DeleteDuplicatesInMapOfStringVectorOfString(const std::map< std::string,std::vector<std::string> >& mymap);
  
  MEDPARTITIONER_EXPORT std::string Cle1ToStr(const std::string& s, const int inew);
  MEDPARTITIONER_EXPORT void Cle1ToData(const std::string& cle, std::string& s, int& inew);
  
  MEDPARTITIONER_EXPORT std::string Cle2ToStr(const std::string& s,const int inew,const int iold);
  MEDPARTITIONER_EXPORT void Cle2ToData(const std::string& cle, std::string& s, int& inew, int& iold);
  
  MEDPARTITIONER_EXPORT std::string ExtractFromDescription(const std::string& description,const  std::string& tag);
  MEDPARTITIONER_EXPORT void FieldDescriptionToData(const std::string& description,
                              int& idomain, std::string& fileName, std::string& meshName, std::string& fieldName,
                              int& typeField, int& DT, int& IT);
  MEDPARTITIONER_EXPORT void FieldShortDescriptionToData(const std::string& description,
                                   std::string& fieldName, int& typeField, int& entity, int& DT, int& IT);
  
  MEDCoupling::DataArrayInt *CreateDataArrayIntFromVector(const std::vector<int>& v);
  MEDCoupling::DataArrayInt *CreateDataArrayIntFromVector(const std::vector<int>& v, const int nbComponents);
  MEDCoupling::DataArrayDouble *CreateDataArrayDoubleFromVector(const std::vector<double>& v);
  
  MEDCoupling::MEDCouplingUMesh *CreateEmptyMEDCouplingUMesh();

  std::vector<std::string> BrowseFieldDouble(const MEDCoupling::MEDCouplingFieldDouble* fd);
  std::vector<std::string> BrowseAllFields(const std::string& myfile);
  std::vector<std::string> BrowseAllFieldsOnMesh(const std::string& myfile, const std::string& mymesh, const int idomain);
  std::vector<std::string> GetInfosOfField(const char *fileName, const char *meshName, const int idomain );

#ifdef HAVE_MPI
  //not adviced, interblocking, use sendAndReceive
  //void SendVectorOfString(const std::vector<std::string>& vec, const int target);
  //std::vector<std::string> RecvVectorOfString(const int source);
  //TODO void sendRecvVectorOfString(const std::vector<std::string>& vec, const int source, const int target);
  MEDPARTITIONER_EXPORT std::vector<std::string> SendAndReceiveVectorOfString(const std::vector<std::string>& vec, const int source, const int target);
  MEDPARTITIONER_EXPORT std::vector<std::string> AllgathervVectorOfString(const std::vector<std::string>& vec);
  
  void SendDoubleVec(const std::vector<double>& vec, const int target);
  std::vector<double> *RecvDoubleVec(const int source);
  void RecvDoubleVec(std::vector<double>& vec, const int source);
    
  void SendIntVec(const std::vector<int>& vec, const int target);
  std::vector<int>* RecvIntVec(int source);
  void RecvIntVec(std::vector<int>& vec, const int source);
  
  void SendDataArrayInt(const MEDCoupling::DataArrayInt* da, const int target);
  MEDCoupling::DataArrayInt *RecvDataArrayInt(const int source);
  void SendDataArrayDouble(const MEDCoupling::DataArrayDouble* da, const int target);
  MEDCoupling::DataArrayDouble *RecvDataArrayDouble(const int source);

  void TestVectorOfStringMpi();
  void TestMapOfStringIntMpi();
  void TestMapOfStringVectorOfStringMpi();
  void TestDataArrayMpi();
  void TestPersistantMpi0To1(int taille, int nb);
  void TestPersistantMpiRing(int taille, int nb);
  void TestPersistantMpiRingOnCommSplit(int taille, int nb);
#endif

  class MEDPARTITIONER_EXPORT MyGlobals
  {
  public :
    static int _Verbose;  //0 to 1000 over 200 is debug
    static int _Rank;
    static int _World_Size;
    static int _Randomize;
    static int _Atomize;
    static int _Create_Boundary_Faces;
    static int _Create_Joints;
    static int _Is0verbose; //trace cout if rank 0 and verbose
    static std::vector<std::string> _File_Names;    //on [iold]
    static std::vector<std::string> _Mesh_Names;    //on [iold]
    static std::vector<std::string> _Field_Descriptions;
    /*! used for descriptions of components of fields for example...*/
    static std::vector<std::string> _General_Informations;
  };



  /*!
   * \brief Class encapsulating BBTree of dimension given at construction and
   *        providing all features of BBTree
   */
  class BBTreeOfDim
  {
    void * _tree;
    void (BBTreeOfDim::*_PgetElementsAroundPoint)( const double* coordsPtr,
                                                   std::vector<int>& elems ) const;
    void (BBTreeOfDim::*_PgetIntersectingElems)( const double* bb,
                                                 std::vector<int>& elems ) const;

    template< int dim>
    void _getElementsAroundPoint( const double* coordsPtr,
                                  std::vector<int>& elems ) const
    {
      ((BBTree<dim,int>*)_tree)->getElementsAroundPoint( coordsPtr, elems );
    }
    template< int dim>
    void _getIntersectingElems(const double* bb,
                               std::vector<int>& elems) const
    {
      ((BBTree<dim,int>*)_tree)->getIntersectingElems( bb, elems );
    }
  public:

    BBTreeOfDim( int           dim,
                 const double* bbs,
                 int*          elems,
                 int           level,
                 int           nbelems,
                 double        epsilon=1e-12);
    ~BBTreeOfDim();
    void getElementsAroundPoint(const double* coordsPtr, std::vector<int>& elems ) const;
    void getIntersectingElems  (const double* bb,        std::vector<int>& elems)  const;
  };
}
#endif
