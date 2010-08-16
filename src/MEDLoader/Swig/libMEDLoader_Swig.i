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

%module libMEDLoader_Swig

%include std_vector.i
%include std_string.i

%include "MEDCouplingTypemaps.i"
%include "libMEDCoupling_Swig.i"

%{
#include "MEDLoader.hxx"
%}

#if SWIG_VERSION >= 0x010329
%template()  std::vector<std::string>;
#endif

%newobject MEDLoader::ReadUMeshFromFile;
%newobject MEDLoader::ReadFieldDouble;
%newobject MEDLoader::ReadFieldDoubleCell;
%newobject MEDLoader::ReadFieldDoubleNode;
%newobject MEDLoader::ReadFieldDoubleGauss;
%newobject MEDLoader::ReadFieldDoubleGaussNE;

class MEDLoader
{
public:
  static void setEpsilonForNodeComp(double val);
  static void setCompPolicyForCell(int val);
  static void setTooLongStrPolicy(int val);
  static std::vector<std::string> GetMeshNames(const char *fileName);
  static std::vector<std::string> GetMeshGroupsNames(const char *fileName, const char *meshName);
  static std::vector<std::string> GetMeshFamilyNames(const char *fileName, const char *meshName);
  static std::vector<std::string> GetCellFieldNamesOnMesh(const char *fileName, const char *meshName);
  static std::vector<std::string> GetNodeFieldNamesOnMesh(const char *fileName, const char *meshName);
  %extend
     {
       static PyObject *GetFieldIterations(ParaMEDMEM::TypeOfField type, const char *fileName, const char *meshName, const char *fieldName)
       {
         std::vector< std::pair<int,int> > res=MEDLoader::GetFieldIterations(type,fileName,meshName,fieldName);
         PyObject *ret=PyList_New(res.size());
         int rk=0;
         for(std::vector< std::pair<int,int> >::const_iterator iter=res.begin();iter!=res.end();iter++,rk++)
           {
             PyObject *elt=PyTuple_New(2);
             PyTuple_SetItem(elt,0,SWIG_From_int((*iter).first));
             PyTuple_SetItem(elt,1,SWIG_From_int((*iter).second));
             PyList_SetItem(ret,rk,elt);
           }
         return ret;
       }

       static PyObject *GetCellFieldIterations(const char *fileName, const char *meshName, const char *fieldName)
       {
         std::vector< std::pair<int,int> > res=MEDLoader::GetCellFieldIterations(fileName,meshName,fieldName);
         PyObject *ret=PyList_New(res.size());
         int rk=0;
         for(std::vector< std::pair<int,int> >::const_iterator iter=res.begin();iter!=res.end();iter++,rk++)
           {
             PyObject *elt=PyTuple_New(2);
             PyTuple_SetItem(elt,0,SWIG_From_int((*iter).first));
             PyTuple_SetItem(elt,1,SWIG_From_int((*iter).second));
             PyList_SetItem(ret,rk,elt);
           }
         return ret;
       }
       static PyObject *GetNodeFieldIterations(const char *fileName, const char *meshName, const char *fieldName)
       {
         std::vector< std::pair<int,int> > res=MEDLoader::GetNodeFieldIterations(fileName,meshName,fieldName);
         PyObject *ret=PyList_New(res.size());
         int rk=0;
         for(std::vector< std::pair<int,int> >::const_iterator iter=res.begin();iter!=res.end();iter++,rk++)
           {
             PyObject *elt=PyTuple_New(2);
             PyTuple_SetItem(elt,0,SWIG_From_int((*iter).first));
             PyTuple_SetItem(elt,1,SWIG_From_int((*iter).second));
             PyList_SetItem(ret,rk,elt);
           }
         return ret;
       }
     }
  static ParaMEDMEM::MEDCouplingUMesh *ReadUMeshFromFile(const char *fileName, const char *meshName, int meshDimRelToMax=0) throw(INTERP_KERNEL::Exception);
  static ParaMEDMEM::MEDCouplingUMesh *ReadUMeshFromFile(const char *fileName, int meshDimRelToMax=0) throw(INTERP_KERNEL::Exception);
  static ParaMEDMEM::MEDCouplingFieldDouble *ReadFieldDouble(ParaMEDMEM::TypeOfField type, const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order);
  static ParaMEDMEM::MEDCouplingFieldDouble *ReadFieldDoubleCell(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order);
  static ParaMEDMEM::MEDCouplingFieldDouble *ReadFieldDoubleNode(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order);
  static ParaMEDMEM::MEDCouplingFieldDouble *ReadFieldDoubleGauss(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order);
  static ParaMEDMEM::MEDCouplingFieldDouble *ReadFieldDoubleGaussNE(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order);
  static void WriteUMesh(const char *fileName, const ParaMEDMEM::MEDCouplingUMesh *mesh, bool writeFromScratch);
  static void WriteField(const char *fileName, const ParaMEDMEM::MEDCouplingFieldDouble *f, bool writeFromScratch);
  static void WriteFieldUsingAlreadyWrittenMesh(const char *fileName, const ParaMEDMEM::MEDCouplingFieldDouble *f);
};
