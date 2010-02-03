//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
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
%module libMEDCoupling_Swig

#define MEDCOUPLING_EXPORT

%{
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingField.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingTypemaps.i"

using namespace ParaMEDMEM;
using namespace INTERP_KERNEL;
%}

%typemap(out) ParaMEDMEM::MEDCouplingMesh*
{
  $result=convertMesh($1,$owner);
}

%newobject ParaMEDMEM::DataArrayDouble::New;
%newobject ParaMEDMEM::DataArrayInt::New;
%newobject ParaMEDMEM::MEDCouplingUMesh::New;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::New;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::mergeFields;
%newobject ParaMEDMEM::MEDCouplingUMesh::clone;
%newobject ParaMEDMEM::DataArrayDouble::deepCopy;
%newobject ParaMEDMEM::DataArrayDouble::performCpy;
%newobject ParaMEDMEM::DataArrayInt::deepCopy;
%newobject ParaMEDMEM::DataArrayInt::performCpy;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::clone;
%newobject ParaMEDMEM::MEDCouplingMesh::mergeMyselfWith;
%newobject ParaMEDMEM::MEDCouplingUMesh::buildPartOfMySelf;
%newobject ParaMEDMEM::MEDCouplingPointSet::zipCoordsTraducer;
%newobject ParaMEDMEM::MEDCouplingUMesh::getMeasureField;
%newobject ParaMEDMEM::MEDCouplingUMesh::mergeUMeshes;
%feature("unref") DataArrayDouble "$this->decrRef();"
%feature("unref") MEDCouplingUMesh "$this->decrRef();"
%feature("unref") DataArrayInt "$this->decrRef();"
%feature("unref") MEDCouplingFieldDouble "$this->decrRef();"

%ignore ParaMEDMEM::TimeLabel::operator=;
%ignore ParaMEDMEM::MemArray::operator=;
%ignore ParaMEDMEM::MemArray::operator[];

%nodefaultctor;
%include "MEDCouplingTimeLabel.hxx"
%include "MEDCouplingRefCountObject.hxx"
%include "MEDCouplingMesh.hxx"
%include "MEDCouplingPointSet.hxx"
%include "MEDCouplingMemArray.hxx"
%include "MEDCouplingMesh.hxx"
%include "NormalizedUnstructuredMesh.hxx"
%include "MEDCouplingField.hxx"
%include "MEDCouplingNatureOfField.hxx"

namespace ParaMEDMEM
{
  class MEDCouplingUMesh : public ParaMEDMEM::MEDCouplingPointSet
  {
  public:
    static MEDCouplingUMesh *New();
    MEDCouplingUMesh *clone(bool recDeepCpy) const;
    void updateTime();
    void checkCoherency() const throw(INTERP_KERNEL::Exception);
    void setMeshDimension(int meshDim);
    void allocateCells(int nbOfCells);
    void setCoords(DataArrayDouble *coords);
    DataArrayDouble *getCoords() const;
    void finishInsertingCells();
    void setConnectivity(DataArrayInt *conn, DataArrayInt *connIndex, bool isComputingTypes=true);
    DataArrayInt *getNodalConnectivity() const;
    DataArrayInt *getNodalConnectivityIndex() const;
    INTERP_KERNEL::NormalizedCellType getTypeOfCell(int cellId) const;
    int getNumberOfNodesInCell(int cellId) const;
    bool isStructured() const;
    int getNumberOfCells() const;
    int getNumberOfNodes() const;
    int getSpaceDimension() const;
    int getMeshDimension() const;
    int getMeshLength() const;
    //tools
    void zipCoords();
    DataArrayInt *zipCoordsTraducer();
    void getReverseNodalConnectivity(DataArrayInt *revNodal, DataArrayInt *revNodalIndx) const;
    MEDCouplingUMesh *buildPartOfMySelf(const int *start, const int *end, bool keepCoords) const;
    %extend {
      void insertNextCell(INTERP_KERNEL::NormalizedCellType type, int size, PyObject *li)
      {
        int *tmp=convertPyToNewIntArr(li,size);
        self->insertNextCell(type,size,tmp);
        delete [] tmp;
      }
      PyObject *getAllTypes() const
      {
        std::set<INTERP_KERNEL::NormalizedCellType> result=self->getAllTypes();
        std::set<INTERP_KERNEL::NormalizedCellType>::const_iterator iL=result.begin();
        PyObject *res = PyList_New(result.size());
        for (int i=0;iL!=result.end(); i++, iL++)
          PyList_SetItem(res,i,PyInt_FromLong(*iL));
        return res;
      }
      PyObject *mergeNodes(double precision)
      {
        bool ret1;
        DataArrayInt *ret0=self->mergeNodes(precision,ret1);
        PyObject *res = PyList_New(2);
        PyList_SetItem(res,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyList_SetItem(res,1,SWIG_From_bool(ret1));
        return res;
      }
    }
    MEDCouplingFieldDouble *getMeasureField(bool isAbs) const;
    static MEDCouplingUMesh *mergeUMeshes(const MEDCouplingUMesh *mesh1, const MEDCouplingUMesh *mesh2);
  };
}

%extend ParaMEDMEM::DataArrayDouble
 {
   void setValues(PyObject *li, int nbOfTuples, int nbOfElsPerTuple)
   {
     double *tmp=convertPyToNewDblArr2(li);
     self->useArray(tmp,true,CPP_DEALLOC,nbOfTuples,nbOfElsPerTuple);
   }

   PyObject *getValues()
   {
     const double *vals=self->getPointer();
     return convertDblArrToPyList(vals,self->getNbOfElems());
   }
 };

%extend ParaMEDMEM::DataArrayInt
 {
   void setValues(PyObject *li, int nbOfTuples, int nbOfElsPerTuple)
   {
     int *tmp=convertPyToNewIntArr2(li);
     self->useArray(tmp,true,CPP_DEALLOC,nbOfTuples,nbOfElsPerTuple);
   }

   PyObject *getValues()
   {
     const int *vals=self->getPointer();
     return convertIntArrToPyList(vals,self->getNbOfElems());
   }
 };

%include "MEDCouplingField.hxx"

namespace ParaMEDMEM
{
  class MEDCouplingFieldDouble : public ParaMEDMEM::MEDCouplingField
  {
  public:
    static MEDCouplingFieldDouble *New(TypeOfField type);
    MEDCouplingFieldDouble *clone(bool recDeepCpy) const;
    void checkCoherency() const throw(INTERP_KERNEL::Exception);
    double getIJ(int tupleId, int compoId) const;
    void setArray(DataArrayDouble *array);
    DataArrayDouble *getArray() const { return _array; }
    void applyLin(double a, double b, int compoId);
    int getNumberOfComponents() const;
    int getNumberOfTuples() const throw(INTERP_KERNEL::Exception);
    NatureOfField getNature() const { return _nature; }
    void setNature(NatureOfField nat) throw(INTERP_KERNEL::Exception);
    void updateTime();
    static MEDCouplingFieldDouble *mergeFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    %extend {
      void setValues(PyObject *li)
      {
        if(self->getArray()!=0)
          {
            double *tmp=convertPyToNewDblArr2(li);
            int nbTuples=self->getArray()->getNumberOfTuples();
            int nbOfCompo=self->getArray()->getNumberOfComponents();
            self->getArray()->useArray(tmp,true,CPP_DEALLOC,nbTuples,nbOfCompo);
          }
        else
          PyErr_SetString(PyExc_TypeError,"setValuesCpy : field must contain an array behind");
      }
      PyObject *getTime()
      {
        int tmp1,tmp2;
        double tmp0=self->getTime(tmp1,tmp2);
        PyObject *res = PyList_New(3);
        PyList_SetItem(res,0,SWIG_From_double(tmp0));
        PyList_SetItem(res,1,SWIG_From_int(tmp1));
        PyList_SetItem(res,2,SWIG_From_int(tmp2));
        return res;
      }
    }
  };
}
