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

%module MEDCoupling

#define MEDCOUPLING_EXPORT

%include std_vector.i
%include std_string.i

%{
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingExtrudedMesh.hxx"
#include "MEDCouplingCMesh.hxx"
#include "MEDCouplingField.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingGaussLocalization.hxx"
#include "MEDCouplingTypemaps.i"

#include "InterpKernelAutoPtr.hxx"

using namespace ParaMEDMEM;
using namespace INTERP_KERNEL;
%}

%template(ivec) std::vector<int>;
%template(dvec) std::vector<double>;
%template(svec) std::vector<std::string>;

%typemap(out) ParaMEDMEM::MEDCouplingMesh*
{
  $result=convertMesh($1,$owner);
}

%typemap(out) ParaMEDMEM::MEDCouplingPointSet*
{
  $result=convertMesh($1,$owner);
}

#ifdef WITH_NUMPY2
%init %{ import_array(); %}
#endif

%feature("autodoc", "1");
%feature("docstring");

%newobject ParaMEDMEM::MEDCouplingField::buildMeasureField;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::New;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::getArray;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::getEndArray;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::mergeFields;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::meldFields;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::doublyContractedProduct;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::determinant;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::eigenValues;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::eigenVectors;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::inverse;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::trace;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::deviator;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::magnitude;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::maxPerTuple;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::keepSelectedComponents;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::dotFields;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::dot;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::crossProductFields;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::crossProduct;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::maxFields;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::max;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::minFields;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::min;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::getIdsInRange;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::buildSubPart;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::operator+;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::operator-;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::operator*;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::operator/;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::clone;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::cloneWithMesh;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::deepCpy;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::buildNewTimeReprFromThis;
%newobject ParaMEDMEM::DataArrayInt::New;
%newobject ParaMEDMEM::DataArrayInt::convertToDblArr;
%newobject ParaMEDMEM::DataArrayInt::deepCpy;
%newobject ParaMEDMEM::DataArrayInt::performCpy;
%newobject ParaMEDMEM::DataArrayInt::substr;
%newobject ParaMEDMEM::DataArrayInt::changeNbOfComponents;
%newobject ParaMEDMEM::DataArrayInt::keepSelectedComponents;
%newobject ParaMEDMEM::DataArrayInt::selectByTupleId;
%newobject ParaMEDMEM::DataArrayInt::selectByTupleIdSafe;
%newobject ParaMEDMEM::DataArrayInt::renumber;
%newobject ParaMEDMEM::DataArrayInt::renumberR;
%newobject ParaMEDMEM::DataArrayInt::renumberAndReduce;
%newobject ParaMEDMEM::DataArrayInt::invertArrayO2N2N2O;
%newobject ParaMEDMEM::DataArrayInt::invertArrayN2O2O2N;
%newobject ParaMEDMEM::DataArrayInt::getIdsEqual;
%newobject ParaMEDMEM::DataArrayInt::getIdsEqualList;
%newobject ParaMEDMEM::DataArrayInt::aggregate;
%newobject ParaMEDMEM::DataArrayInt::meld;
%newobject ParaMEDMEM::DataArrayInt::fromNoInterlace;
%newobject ParaMEDMEM::DataArrayInt::toNoInterlace;
%newobject ParaMEDMEM::DataArrayInt::buildComplement;
%newobject ParaMEDMEM::DataArrayInt::buildUnion;
%newobject ParaMEDMEM::DataArrayInt::buildSubstraction;
%newobject ParaMEDMEM::DataArrayInt::buildIntersection;
%newobject ParaMEDMEM::DataArrayInt::deltaShiftIndex;
%newobject ParaMEDMEM::DataArrayDouble::New;
%newobject ParaMEDMEM::DataArrayDouble::convertToIntArr;
%newobject ParaMEDMEM::DataArrayDouble::deepCpy;
%newobject ParaMEDMEM::DataArrayDouble::performCpy;
%newobject ParaMEDMEM::DataArrayDouble::aggregate;
%newobject ParaMEDMEM::DataArrayDouble::meld;
%newobject ParaMEDMEM::DataArrayDouble::dot;
%newobject ParaMEDMEM::DataArrayDouble::crossProduct;
%newobject ParaMEDMEM::DataArrayDouble::add;
%newobject ParaMEDMEM::DataArrayDouble::substract;
%newobject ParaMEDMEM::DataArrayDouble::multiply;
%newobject ParaMEDMEM::DataArrayDouble::divide;
%newobject ParaMEDMEM::DataArrayDouble::substr;
%newobject ParaMEDMEM::DataArrayDouble::changeNbOfComponents;
%newobject ParaMEDMEM::DataArrayDouble::keepSelectedComponents;
%newobject ParaMEDMEM::DataArrayDouble::getIdsInRange;
%newobject ParaMEDMEM::DataArrayDouble::selectByTupleId;
%newobject ParaMEDMEM::DataArrayDouble::selectByTupleIdSafe;
%newobject ParaMEDMEM::DataArrayDouble::applyFunc;
%newobject ParaMEDMEM::DataArrayDouble::doublyContractedProduct;
%newobject ParaMEDMEM::DataArrayDouble::determinant;
%newobject ParaMEDMEM::DataArrayDouble::eigenValues;
%newobject ParaMEDMEM::DataArrayDouble::eigenVectors;
%newobject ParaMEDMEM::DataArrayDouble::inverse;
%newobject ParaMEDMEM::DataArrayDouble::trace;
%newobject ParaMEDMEM::DataArrayDouble::deviator;
%newobject ParaMEDMEM::DataArrayDouble::magnitude;
%newobject ParaMEDMEM::DataArrayDouble::maxPerTuple;
%newobject ParaMEDMEM::DataArrayDouble::renumber;
%newobject ParaMEDMEM::DataArrayDouble::renumberR;
%newobject ParaMEDMEM::DataArrayDouble::renumberAndReduce;
%newobject ParaMEDMEM::DataArrayDouble::fromNoInterlace;
%newobject ParaMEDMEM::DataArrayDouble::toNoInterlace;
%newobject ParaMEDMEM::DataArrayDouble::fromPolarToCart;
%newobject ParaMEDMEM::DataArrayDouble::fromCylToCart;
%newobject ParaMEDMEM::DataArrayDouble::fromSpherToCart;
%newobject ParaMEDMEM::MEDCouplingMesh::getCoordinatesAndOwner;
%newobject ParaMEDMEM::MEDCouplingMesh::getBarycenterAndOwner;
%newobject ParaMEDMEM::MEDCouplingMesh::buildOrthogonalField;
%newobject ParaMEDMEM::MEDCouplingMesh::getCellIdsFullyIncludedInNodeIds;
%newobject ParaMEDMEM::MEDCouplingMesh::mergeMyselfWith;
%newobject ParaMEDMEM::MEDCouplingMesh::fillFromAnalytic;
%newobject ParaMEDMEM::MEDCouplingMesh::getMeasureField;
%newobject ParaMEDMEM::MEDCouplingMesh::simplexize;
%newobject ParaMEDMEM::MEDCouplingMesh::buildUnstructured;
%newobject ParaMEDMEM::MEDCouplingPointSet::zipCoordsTraducer;
%newobject ParaMEDMEM::MEDCouplingPointSet::buildBoundaryMesh;
%newobject ParaMEDMEM::MEDCouplingPointSet::mergeNodesArray;
%newobject ParaMEDMEM::MEDCouplingPointSet::buildInstanceFromMeshType;
%newobject ParaMEDMEM::MEDCouplingUMesh::New;
%newobject ParaMEDMEM::MEDCouplingUMesh::clone;
%newobject ParaMEDMEM::MEDCouplingUMesh::zipConnectivityTraducer;
%newobject ParaMEDMEM::MEDCouplingUMesh::buildDescendingConnectivity;
%newobject ParaMEDMEM::MEDCouplingUMesh::buildExtrudedMeshFromThis;
%newobject ParaMEDMEM::MEDCouplingUMesh::mergeUMeshes;
%newobject ParaMEDMEM::MEDCouplingUMesh::mergeUMeshesOnSameCoords;
%newobject ParaMEDMEM::MEDCouplingUMesh::buildNewNumberingFromCommNodesFrmt;
%newobject ParaMEDMEM::MEDCouplingUMesh::rearrange2ConsecutiveCellTypes;
%newobject ParaMEDMEM::MEDCouplingUMesh::convertCellArrayPerGeoType;
%newobject ParaMEDMEM::MEDCouplingUMesh::getRenumArrForConsecutiveCellTypesSpec;
%newobject ParaMEDMEM::MEDCouplingUMesh::buildDirectionVectorField;
%newobject ParaMEDMEM::MEDCouplingUMesh::getEdgeRatioField;
%newobject ParaMEDMEM::MEDCouplingUMesh::getAspectRatioField;
%newobject ParaMEDMEM::MEDCouplingUMesh::getWarpField;
%newobject ParaMEDMEM::MEDCouplingUMesh::getSkewField;
%newobject ParaMEDMEM::MEDCouplingExtrudedMesh::New;
%newobject ParaMEDMEM::MEDCouplingExtrudedMesh::build3DUnstructuredMesh;
%newobject ParaMEDMEM::MEDCouplingCMesh::New;
%feature("unref") DataArrayDouble "$this->decrRef();"
%feature("unref") MEDCouplingPointSet "$this->decrRef();"
%feature("unref") MEDCouplingMesh "$this->decrRef();"
%feature("unref") MEDCouplingUMesh "$this->decrRef();"
%feature("unref") MEDCouplingExtrudedMesh "$this->decrRef();"
%feature("unref") MEDCouplingCMesh "$this->decrRef();"
%feature("unref") DataArrayInt "$this->decrRef();"
%feature("unref") MEDCouplingField "$this->decrRef();"
%feature("unref") MEDCouplingFieldDouble "$this->decrRef();"

%rename(assign) *::operator=;
%ignore ParaMEDMEM::MemArray::operator=;
%ignore ParaMEDMEM::MemArray::operator[];
%ignore ParaMEDMEM::MEDCouplingPointSet::getCoords();
%ignore ParaMEDMEM::MEDCouplingGaussLocalization::pushTinySerializationIntInfo;
%ignore ParaMEDMEM::MEDCouplingGaussLocalization::pushTinySerializationDblInfo;
%ignore ParaMEDMEM::MEDCouplingGaussLocalization::fillWithValues;
%ignore ParaMEDMEM::MEDCouplingGaussLocalization::buildNewInstanceFromTinyInfo;

%rename (Exception) InterpKernelException;
%nodefaultctor;

namespace INTERP_KERNEL
{
  class Exception
  {
  public:
    Exception(const char* what);
    ~Exception() throw ();
    const char *what() const throw ();
  };
}

%include "MEDCouplingTimeLabel.hxx"
%include "MEDCouplingRefCountObject.hxx"

namespace ParaMEDMEM
{
  typedef enum
    {
      UNSTRUCTURED = 5,
      UNSTRUCTURED_DESC = 6,
      CARTESIAN = 7,
      EXTRUDED = 8
    } MEDCouplingMeshType;

  class DataArrayInt;
  class DataArrayDouble;
  class MEDCouplingUMesh;
  class MEDCouplingFieldDouble;

  class MEDCouplingMesh : public RefCountObject, public TimeLabel
  {
  public:
    void setName(const char *name) { _name=name; }
    const char *getName() const { return _name.c_str(); }
    virtual MEDCouplingMeshType getType() const throw(INTERP_KERNEL::Exception) = 0;
    bool isStructured() const throw(INTERP_KERNEL::Exception);
    virtual bool isEqual(const MEDCouplingMesh *other, double prec) const throw(INTERP_KERNEL::Exception);
    virtual bool isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const throw(INTERP_KERNEL::Exception) = 0;
    virtual void copyTinyStringsFrom(const MEDCouplingMesh *other) throw(INTERP_KERNEL::Exception);
    virtual void checkCoherency() const throw(INTERP_KERNEL::Exception) = 0;
    virtual int getNumberOfCells() const throw(INTERP_KERNEL::Exception) = 0;
    virtual int getNumberOfNodes() const throw(INTERP_KERNEL::Exception) = 0;
    virtual int getSpaceDimension() const throw(INTERP_KERNEL::Exception) = 0;
    virtual int getMeshDimension() const throw(INTERP_KERNEL::Exception) = 0;
    virtual DataArrayDouble *getCoordinatesAndOwner() const throw(INTERP_KERNEL::Exception) = 0;
    virtual DataArrayDouble *getBarycenterAndOwner() const throw(INTERP_KERNEL::Exception) = 0;
    virtual int getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType type) const throw(INTERP_KERNEL::Exception) = 0;
    virtual INTERP_KERNEL::NormalizedCellType getTypeOfCell(int cellId) const throw(INTERP_KERNEL::Exception) = 0;
    virtual std::string simpleRepr() const = 0;
    virtual std::string advancedRepr() const = 0;
    // tools
    virtual MEDCouplingFieldDouble *getMeasureField(bool isAbs) const throw(INTERP_KERNEL::Exception) = 0;
    virtual MEDCouplingFieldDouble *getMeasureFieldOnNode(bool isAbs) const throw(INTERP_KERNEL::Exception) = 0;
    virtual MEDCouplingFieldDouble *fillFromAnalytic(TypeOfField t, int nbOfComp, const char *func) const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingFieldDouble *buildOrthogonalField() const throw(INTERP_KERNEL::Exception) = 0;
    virtual MEDCouplingUMesh *buildUnstructured() const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingMesh *mergeMyselfWith(const MEDCouplingMesh *other) const throw(INTERP_KERNEL::Exception) = 0;
    virtual bool areCompatibleForMerge(const MEDCouplingMesh *other) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *simplexize(int policy) throw(INTERP_KERNEL::Exception);
    static MEDCouplingMesh *mergeMeshes(const MEDCouplingMesh *mesh1, const MEDCouplingMesh *mesh2) throw(INTERP_KERNEL::Exception);
    %extend
       {
         std::string __str__() const
         {
           return self->simpleRepr();
         }

         int getCellContainingPoint(PyObject *p, double eps) const throw(INTERP_KERNEL::Exception)
         {
           int sz;
           double *pos=convertPyToNewDblArr2(p,&sz);
           int ret=self->getCellContainingPoint(pos,eps);
           delete [] pos;
           return ret;
         }
         
         void renumberCells(PyObject *li, bool check) throw(INTERP_KERNEL::Exception)
         {
           int size;
           INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
           self->renumberCells(tmp,check);
         }

         PyObject *checkGeoEquivalWith(const MEDCouplingMesh *other, int levOfCheck, double prec) const throw(INTERP_KERNEL::Exception)
         {
           DataArrayInt *cellCor, *nodeCor;
           self->checkGeoEquivalWith(other,levOfCheck,prec,cellCor,nodeCor);
           PyObject *res = PyList_New(2);
           PyList_SetItem(res,0,SWIG_NewPointerObj(SWIG_as_voidptr(cellCor),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, cellCor?SWIG_POINTER_OWN | 0:0 ));
           PyList_SetItem(res,1,SWIG_NewPointerObj(SWIG_as_voidptr(nodeCor),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, nodeCor?SWIG_POINTER_OWN | 0:0 ));
           return res;
         }
         DataArrayInt *getCellIdsFullyIncludedInNodeIds(PyObject *li) const throw(INTERP_KERNEL::Exception)
         {
           int size;
           INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
           DataArrayInt *ret=self->getCellIdsFullyIncludedInNodeIds(tmp,tmp+size);
           return ret;
         }
         PyObject *getNodeIdsOfCell(int cellId) const throw(INTERP_KERNEL::Exception)
         {
           std::vector<int> conn;
           self->getNodeIdsOfCell(cellId,conn);
           return convertIntArrToPyList2(conn);
         }

         PyObject *getCoordinatesOfNode(int nodeId) const throw(INTERP_KERNEL::Exception)
         {
           std::vector<double> coo;
           self->getCoordinatesOfNode(nodeId,coo);
           return convertDblArrToPyList2(coo);
         }

         void scale(PyObject *point, double factor) throw(INTERP_KERNEL::Exception)
         {
           int sz;
           double *p=convertPyToNewDblArr2(point,&sz);
           self->scale(p,factor);
           delete [] p;
         }

         PyObject *getBoundingBox() const throw(INTERP_KERNEL::Exception)
         {
           int spaceDim=self->getSpaceDimension();
           double *tmp=new double[2*spaceDim];
           self->getBoundingBox(tmp);
           PyObject *ret=convertDblArrToPyListOfTuple(tmp,2,spaceDim);
           delete [] tmp;
           return ret;
         }

        PyObject *buildPart(PyObject *li) const throw(INTERP_KERNEL::Exception)
         {
           int size;
           INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
           MEDCouplingMesh *ret=self->buildPart(tmp,tmp+size);
           return convertMesh(ret, SWIG_POINTER_OWN | 0 );
         }

        PyObject *buildPartAndReduceNodes(PyObject *li) const throw(INTERP_KERNEL::Exception)
        {
          int size;
          INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
          DataArrayInt *arr=0;
          MEDCouplingMesh *ret=self->buildPartAndReduceNodes(tmp,tmp+size,arr);
          PyObject *res = PyList_New(2);
          PyObject *obj0=convertMesh(ret, SWIG_POINTER_OWN | 0 );
          PyObject *obj1=SWIG_NewPointerObj(SWIG_as_voidptr(arr),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 );
          PyList_SetItem(res,0,obj0);
          PyList_SetItem(res,1,obj1);
          return res;
        }

         void translate(PyObject *vector) throw(INTERP_KERNEL::Exception)
         {
           int sz;
           INTERP_KERNEL::AutoPtr<double> v=convertPyToNewDblArr2(vector,&sz);
           self->translate(v);
         }

         void rotate(PyObject *center, PyObject *vector, double alpha) throw(INTERP_KERNEL::Exception)
         {
           int sz;
           INTERP_KERNEL::AutoPtr<double> c=convertPyToNewDblArr2(center,&sz);
           if(!c)
             return ;
           INTERP_KERNEL::AutoPtr<double> v=convertPyToNewDblArr2(vector,&sz);
           if(!v)
             { return ; }
           self->rotate(c,v,alpha);
         }
       }
  };
}

%include "MEDCouplingMemArray.hxx"
%include "NormalizedUnstructuredMesh.hxx"
%include "MEDCouplingNatureOfField.hxx"
%include "MEDCouplingTimeDiscretization.hxx"
%include "MEDCouplingGaussLocalization.hxx"

namespace ParaMEDMEM
{
  class MEDCouplingPointSet : public ParaMEDMEM::MEDCouplingMesh
    {
    public:
      void updateTime();
      void setCoords(DataArrayDouble *coords) throw(INTERP_KERNEL::Exception);
      DataArrayDouble *getCoordinatesAndOwner() const throw(INTERP_KERNEL::Exception);
      bool areCoordsEqual(const MEDCouplingPointSet& other, double prec) const throw(INTERP_KERNEL::Exception);
      void zipCoords() throw(INTERP_KERNEL::Exception);
      double getCaracteristicDimension() const throw(INTERP_KERNEL::Exception);
      void changeSpaceDimension(int newSpaceDim, double dftVal=0.) throw(INTERP_KERNEL::Exception);
      void tryToShareSameCoords(const MEDCouplingPointSet& other, double epsilon) throw(INTERP_KERNEL::Exception);
      virtual void tryToShareSameCoordsPermute(const MEDCouplingPointSet& other, double epsilon) throw(INTERP_KERNEL::Exception) = 0;
      static DataArrayDouble *mergeNodesArray(const MEDCouplingPointSet *m1, const MEDCouplingPointSet *m2) throw(INTERP_KERNEL::Exception);
      static MEDCouplingPointSet *buildInstanceFromMeshType(MEDCouplingMeshType type) throw(INTERP_KERNEL::Exception);
      virtual MEDCouplingPointSet *buildBoundaryMesh(bool keepCoords) const throw(INTERP_KERNEL::Exception) = 0;
      virtual bool isEmptyMesh(const std::vector<int>& tinyInfo) const throw(INTERP_KERNEL::Exception) = 0;
      //! size of returned tinyInfo must be always the same.
      void getTinySerializationInformation(std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const throw(INTERP_KERNEL::Exception);
      void resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings) const throw(INTERP_KERNEL::Exception);
      void serialize(DataArrayInt *&a1, DataArrayDouble *&a2) const throw(INTERP_KERNEL::Exception);
      void unserialization(const std::vector<int>& tinyInfo, const DataArrayInt *a1, DataArrayDouble *a2,
                           const std::vector<std::string>& littleStrings) throw(INTERP_KERNEL::Exception);
      virtual void giveElemsInBoundingBox(const INTERP_KERNEL::DirectedBoundingBox& bbox, double eps, std::vector<int>& elems) throw(INTERP_KERNEL::Exception) = 0;
      virtual DataArrayInt *zipCoordsTraducer() throw(INTERP_KERNEL::Exception) = 0;
      %extend 
         {
           std::string __str__() const
           {
             return self->simpleRepr();
           }
           
           PyObject *buildNewNumberingFromCommonNodesFormat(const DataArrayInt *comm, const DataArrayInt *commIndex) const throw(INTERP_KERNEL::Exception)
           {
             int newNbOfNodes;
             DataArrayInt *ret0=self->buildNewNumberingFromCommonNodesFormat(comm,commIndex,newNbOfNodes);
             PyObject *res = PyList_New(2);
             PyList_SetItem(res,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
             PyList_SetItem(res,1,SWIG_From_int(newNbOfNodes));
             return res;
           }
           
           PyObject *findCommonNodes(int limitNodeId, double prec) const throw(INTERP_KERNEL::Exception)
           {
             DataArrayInt *comm, *commIndex;
             self->findCommonNodes(limitNodeId,prec,comm,commIndex);
             PyObject *res = PyList_New(2);
             PyList_SetItem(res,0,SWIG_NewPointerObj(SWIG_as_voidptr(comm),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
             PyList_SetItem(res,1,SWIG_NewPointerObj(SWIG_as_voidptr(commIndex),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
             return res;
           }
           
           PyObject *getCoords() const throw(INTERP_KERNEL::Exception)
           {
             DataArrayDouble *ret1=self->getCoords();
             ret1->incrRef();
             return SWIG_NewPointerObj((void*)ret1,SWIGTYPE_p_ParaMEDMEM__DataArrayDouble,SWIG_POINTER_OWN | 0);
           }
           PyObject *buildPartOfMySelf(PyObject *li, bool keepCoords) const throw(INTERP_KERNEL::Exception)
           {
             int size;
             INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
             MEDCouplingPointSet *ret=self->buildPartOfMySelf(tmp,tmp+size,keepCoords);
             return convertMesh(ret, SWIG_POINTER_OWN | 0 );
           }
           PyObject *buildPartOfMySelfNode(PyObject *li, bool fullyIn) const throw(INTERP_KERNEL::Exception)
           {
             int size;
             INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
             MEDCouplingPointSet *ret=self->buildPartOfMySelfNode(tmp,(int *)tmp+size,fullyIn);
             return convertMesh(ret, SWIG_POINTER_OWN | 0 );
           }
           PyObject *buildFacePartOfMySelfNode(PyObject *li, bool fullyIn) const throw(INTERP_KERNEL::Exception)
           {
             int size;
             INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
             MEDCouplingPointSet *ret=self->buildFacePartOfMySelfNode(tmp,(int *)tmp+size,fullyIn);
             return convertMesh(ret, SWIG_POINTER_OWN | 0 );
           }
           PyObject *findBoundaryNodes() const throw(INTERP_KERNEL::Exception)
           {
             std::vector<int> nodes;
             self->findBoundaryNodes(nodes);
             return convertIntArrToPyList2(nodes);
           }
           void renumberNodes(PyObject *li, int newNbOfNodes) throw(INTERP_KERNEL::Exception)
           {
             int size;
             INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
             self->renumberNodes(tmp,newNbOfNodes);
           }
           void renumberNodes2(PyObject *li, int newNbOfNodes) throw(INTERP_KERNEL::Exception)
           {
             int size;
             INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
             self->renumberNodes2(tmp,newNbOfNodes);
           }
           PyObject *findNodesOnPlane(PyObject *pt, PyObject *vec, double eps) const throw(INTERP_KERNEL::Exception)
             {
               std::vector<int> nodes;
               int sz;
               double *p=convertPyToNewDblArr2(pt,&sz);
               double *v=convertPyToNewDblArr2(vec,&sz);
               self->findNodesOnPlane(p,v,eps,nodes);
               delete [] v;
               delete [] p;
               return convertIntArrToPyList2(nodes);
             }
           PyObject *getNodeIdsNearPoint(PyObject *pt, double eps) const throw(INTERP_KERNEL::Exception)
           {
             int size;
             double *pos=convertPyToNewDblArr2(pt,&size);
             if(size<self->getSpaceDimension())
               {
                 delete [] pos;
                 throw INTERP_KERNEL::Exception("getNodeIdsNearPoint : to tiny array ! must be at least of size SpaceDim !");
               }
             std::vector<int> tmp;
             try
               {
                 tmp=self->getNodeIdsNearPoint(pos,eps);
               }
             catch(INTERP_KERNEL::Exception& e)
               {
                 delete [] pos;
                 throw e;
               }
             delete [] pos;
             return convertIntArrToPyList2(tmp);
           }

           PyObject *getNodeIdsNearPoints(PyObject *pt, int nbOfNodes, double eps) const throw(INTERP_KERNEL::Exception)
           {
             std::vector<int> c,cI;
             int size;
             double *pos=convertPyToNewDblArr2(pt,&size);
             if(size<self->getSpaceDimension()*nbOfNodes)
               {
                 delete [] pos;
                 throw INTERP_KERNEL::Exception("getNodeIdsNearPoints : to tiny array ! must be at least of size SpaceDim*nbOfNodes !");
               }
             try
               {
                 self->getNodeIdsNearPoints(pos,nbOfNodes,eps,c,cI);
               }
             catch(INTERP_KERNEL::Exception& e)
               {
                 delete [] pos;
                 throw e;
               }
             delete [] pos;
             PyObject *ret=PyTuple_New(2);
             PyTuple_SetItem(ret,0,convertIntArrToPyList2(c));
             PyTuple_SetItem(ret,1,convertIntArrToPyList2(cI));
             return ret;
           }

           PyObject *giveElemsInBoundingBox(PyObject *bbox, double eps) throw(INTERP_KERNEL::Exception)
           {
             std::vector<int> elems;
             int size;
             double *tmp=convertPyToNewDblArr2(bbox,&size);
             self->giveElemsInBoundingBox(tmp,eps,elems);
             delete [] tmp;
             return convertIntArrToPyList2(elems);
           }

           static void rotate2DAlg(PyObject *center, double angle, int nbNodes, PyObject *coords) throw(INTERP_KERNEL::Exception)
           {
             int sz;
             double *c=convertPyToNewDblArr2(center,&sz);
             double *coo=convertPyToNewDblArr2(coords,&sz);
             ParaMEDMEM::MEDCouplingPointSet::rotate2DAlg(c,angle,nbNodes,coo);
             for(int i=0;i<sz;i++)
               PyList_SetItem(coords,i,PyFloat_FromDouble(coo[i]));
             delete [] coo;
             delete [] c;
           }
           static void rotate3DAlg(PyObject *center, PyObject *vect, double angle, int nbNodes, PyObject *coords) throw(INTERP_KERNEL::Exception)
           {
             int sz,sz2;
             double *c=convertPyToNewDblArr2(center,&sz);
             double *coo=convertPyToNewDblArr2(coords,&sz);
             double *v=convertPyToNewDblArr2(vect,&sz2);
             ParaMEDMEM::MEDCouplingPointSet::rotate3DAlg(c,v,angle,nbNodes,coo);
             for(int i=0;i<sz;i++)
               PyList_SetItem(coords,i,PyFloat_FromDouble(coo[i]));
             delete [] coo;
             delete [] c;
           }
         }
    };
  
  class MEDCouplingUMesh : public ParaMEDMEM::MEDCouplingPointSet
  {
  public:
    static MEDCouplingUMesh *New();
    static MEDCouplingUMesh *New(const char *meshName, int meshDim);
    MEDCouplingUMesh *clone(bool recDeepCpy) const;
    void updateTime();
    void checkCoherency() const throw(INTERP_KERNEL::Exception);
    void setMeshDimension(int meshDim) throw(INTERP_KERNEL::Exception);
    void allocateCells(int nbOfCells) throw(INTERP_KERNEL::Exception);
    void finishInsertingCells() throw(INTERP_KERNEL::Exception);
    void setConnectivity(DataArrayInt *conn, DataArrayInt *connIndex, bool isComputingTypes=true) throw(INTERP_KERNEL::Exception);
    DataArrayInt *getNodalConnectivity() const throw(INTERP_KERNEL::Exception);
    DataArrayInt *getNodalConnectivityIndex() const throw(INTERP_KERNEL::Exception);
    INTERP_KERNEL::NormalizedCellType getTypeOfCell(int cellId) const throw(INTERP_KERNEL::Exception);
    int getNumberOfNodesInCell(int cellId) const throw(INTERP_KERNEL::Exception);
    int getMeshLength() const throw(INTERP_KERNEL::Exception);
    void computeTypes() throw(INTERP_KERNEL::Exception);
    std::string reprConnectivityOfThis() const throw(INTERP_KERNEL::Exception);
    //tools
    bool checkConsecutiveCellTypes() const throw(INTERP_KERNEL::Exception);
    DataArrayInt *rearrange2ConsecutiveCellTypes() throw(INTERP_KERNEL::Exception);
    DataArrayInt *convertCellArrayPerGeoType(const DataArrayInt *da) const throw(INTERP_KERNEL::Exception);
    DataArrayInt *zipConnectivityTraducer(int compType) throw(INTERP_KERNEL::Exception);
    void getReverseNodalConnectivity(DataArrayInt *revNodal, DataArrayInt *revNodalIndx) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *buildDescendingConnectivity(DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *revDesc, DataArrayInt *revDescIndx) const throw(INTERP_KERNEL::Exception);
    void orientCorrectlyPolyhedrons() throw(INTERP_KERNEL::Exception);
    bool isPresenceOfQuadratic() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *buildDirectionVectorField() const throw(INTERP_KERNEL::Exception);
    bool isContiguous1D() const throw(INTERP_KERNEL::Exception);
    void convertQuadraticCellsToLinear() throw(INTERP_KERNEL::Exception);
    void convertDegeneratedCells() throw(INTERP_KERNEL::Exception);
    bool areOnlySimplexCells() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getEdgeRatioField() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getAspectRatioField() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getWarpField() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getSkewField() const throw(INTERP_KERNEL::Exception);
    static MEDCouplingUMesh *mergeUMeshes(const MEDCouplingUMesh *mesh1, const MEDCouplingUMesh *mesh2) throw(INTERP_KERNEL::Exception);
    static MEDCouplingUMesh *mergeUMeshesOnSameCoords(const MEDCouplingUMesh *mesh1, const MEDCouplingUMesh *mesh2) throw(INTERP_KERNEL::Exception);
    %extend {
      std::string __str__() const
      {
        return self->simpleRepr();
      }
      void insertNextCell(INTERP_KERNEL::NormalizedCellType type, int size, PyObject *li) throw(INTERP_KERNEL::Exception)
      {
        int sz;
        INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&sz);
        self->insertNextCell(type,size,tmp);
      }
      PyObject *getAllTypes() const throw(INTERP_KERNEL::Exception)
      {
        std::set<INTERP_KERNEL::NormalizedCellType> result=self->getAllTypes();
        std::set<INTERP_KERNEL::NormalizedCellType>::const_iterator iL=result.begin();
        PyObject *res = PyList_New(result.size());
        for (int i=0;iL!=result.end(); i++, iL++)
          PyList_SetItem(res,i,PyInt_FromLong(*iL));
        return res;
      }
      PyObject *mergeNodes(double precision) throw(INTERP_KERNEL::Exception)
      {
        bool ret1;
        int ret2;
        DataArrayInt *ret0=self->mergeNodes(precision,ret1,ret2);
        PyObject *res = PyList_New(3);
        PyList_SetItem(res,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyList_SetItem(res,1,SWIG_From_bool(ret1));
        PyList_SetItem(res,2,SWIG_From_int(ret2));
        return res;
      }
      PyObject *mergeNodes2(double precision) throw(INTERP_KERNEL::Exception)
      {
        bool ret1;
        int ret2;
        DataArrayInt *ret0=self->mergeNodes2(precision,ret1,ret2);
        PyObject *res = PyList_New(3);
        PyList_SetItem(res,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyList_SetItem(res,1,SWIG_From_bool(ret1));
        PyList_SetItem(res,2,SWIG_From_int(ret2));
        return res;
      }
      PyObject *checkButterflyCells() throw(INTERP_KERNEL::Exception)
      {
        std::vector<int> cells;
        self->checkButterflyCells(cells);
        return convertIntArrToPyList2(cells);
      }

      PyObject *splitByType() const throw(INTERP_KERNEL::Exception)
      {
        std::vector<MEDCouplingUMesh *> ms=self->splitByType();
        int sz=ms.size();
        PyObject *ret = PyList_New(sz);
        for(int i=0;i<sz;i++)
          PyList_SetItem(ret,i,SWIG_NewPointerObj(SWIG_as_voidptr(ms[i]),SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *keepSpecifiedCells(INTERP_KERNEL::NormalizedCellType type, PyObject *ids) const throw(INTERP_KERNEL::Exception)
      {
        std::vector<int> idsPerGeoType;
        convertPyToNewIntArr3(ids,idsPerGeoType);
        MEDCouplingUMesh *ret=self->keepSpecifiedCells(type,idsPerGeoType);
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 );
      }

      bool checkConsecutiveCellTypesAndOrder(PyObject *li) const throw(INTERP_KERNEL::Exception)
      {
        int sz;
        INTERP_KERNEL::AutoPtr<INTERP_KERNEL::NormalizedCellType> order=(INTERP_KERNEL::NormalizedCellType *)convertPyToNewIntArr2(li,&sz);
        bool ret=self->checkConsecutiveCellTypesAndOrder(order,order+sz);
        return ret;
      }

      DataArrayInt *getRenumArrForConsecutiveCellTypesSpec(PyObject *li) const throw(INTERP_KERNEL::Exception)
      {
        int sz;
        INTERP_KERNEL::AutoPtr<INTERP_KERNEL::NormalizedCellType> order=(INTERP_KERNEL::NormalizedCellType *)convertPyToNewIntArr2(li,&sz);
        DataArrayInt *ret=self->getRenumArrForConsecutiveCellTypesSpec(order,(INTERP_KERNEL::NormalizedCellType *)order+sz);
        return ret;
      }

      PyObject *getCellsContainingPoints(PyObject *p, int nbOfPoints, double eps) const throw(INTERP_KERNEL::Exception)
      {
        int sz;
        double *pos=convertPyToNewDblArr2(p,&sz);
        std::vector<int> elts,eltsIndex;
        self->getCellsContainingPoints(pos,nbOfPoints,eps,elts,eltsIndex);
        delete [] pos;
        PyObject *ret=PyList_New(2);
        PyList_SetItem(ret,0,convertIntArrToPyList2(elts));
        PyList_SetItem(ret,1,convertIntArrToPyList2(eltsIndex));
        return ret;
      }

      PyObject *getCellsContainingPoint(PyObject *p, double eps) const throw(INTERP_KERNEL::Exception)
      {
        int sz;
        double *pos=convertPyToNewDblArr2(p,&sz);
        std::vector<int> elts;
        self->getCellsContainingPoint(pos,eps,elts);
        delete [] pos;
        return convertIntArrToPyList2(elts);
      }

      static PyObject *mergeUMeshesOnSameCoords(PyObject *ms) throw(INTERP_KERNEL::Exception)
      {
        std::vector<ParaMEDMEM::MEDCouplingUMesh *> meshes;
        if(PyList_Check(ms))
          {
            int sz=PyList_Size(ms);
            meshes.resize(sz);
            for(int i=0;i<sz;i++)
              {
                PyObject *o=PyList_GetItem(ms,i);
                void *arg;
                SWIG_ConvertPtr(o,&arg,SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh, 0 |  0 );
                meshes[i]=reinterpret_cast<ParaMEDMEM::MEDCouplingUMesh *>(arg);
              }
          }
        else
          {
            PyErr_SetString(PyExc_TypeError,"mergeUMeshesOnSameCoords : not a list as first parameter");
            PyErr_Print();
            return 0;
          }
        MEDCouplingUMesh *ret=MEDCouplingUMesh::mergeUMeshesOnSameCoords(meshes);
        return convertMesh(ret, SWIG_POINTER_OWN | 0 );
      }

      static PyObject *fuseUMeshesOnSameCoords(PyObject *ms, int compType) throw(INTERP_KERNEL::Exception)
      {
        int sz;
        std::vector<MEDCouplingUMesh *> meshes;
        convertPyObjToVecUMeshes(ms,meshes);
        std::vector<DataArrayInt *> corr;
        MEDCouplingUMesh *um=MEDCouplingUMesh::fuseUMeshesOnSameCoords(meshes,compType,corr);
        sz=corr.size();
        PyObject *ret1=PyList_New(sz);
        for(int i=0;i<sz;i++)
          PyList_SetItem(ret1,i,SWIG_NewPointerObj(SWIG_as_voidptr(corr[i]),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyObject *ret=PyList_New(2);
        PyList_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(um),SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyList_SetItem(ret,1,ret1);
        return ret;
      }

      PyObject *are2DCellsNotCorrectlyOriented(PyObject *vec, bool polyOnly) const throw(INTERP_KERNEL::Exception)
      {
        std::vector<int> cells;
        int sz;
        double *v=convertPyToNewDblArr2(vec,&sz);
        try
          {
            self->are2DCellsNotCorrectlyOriented(v,polyOnly,cells);
          }
        catch(INTERP_KERNEL::Exception& e)
          {
            delete [] v;
            throw e;
          }
        delete [] v;
        return convertIntArrToPyList2(cells);
      }

      void orientCorrectly2DCells(PyObject *vec, bool polyOnly) throw(INTERP_KERNEL::Exception)
      {
        int sz;
        double *v=convertPyToNewDblArr2(vec,&sz);
        try
          {
            self->orientCorrectly2DCells(v,polyOnly);
          }
        catch(INTERP_KERNEL::Exception& e)
          {
            delete [] v;
            throw e;
          }
        delete [] v;
      }
      
      PyObject *arePolyhedronsNotCorrectlyOriented() const throw(INTERP_KERNEL::Exception)
      {
        std::vector<int> cells;
        self->arePolyhedronsNotCorrectlyOriented(cells);
        return convertIntArrToPyList2(cells);
      }

      PyObject *getFastAveragePlaneOfThis() const throw(INTERP_KERNEL::Exception)
      {
        double vec[3];
        double pos[3];
        self->getFastAveragePlaneOfThis(vec,pos);
        double vals[6];
        std::copy(vec,vec+3,vals);
        std::copy(pos,pos+3,vals+3);
        return convertDblArrToPyListOfTuple(vals,3,2);
      }
      
      static MEDCouplingUMesh *mergeUMeshes(PyObject *li) throw(INTERP_KERNEL::Exception)
      {
        std::vector<const ParaMEDMEM::MEDCouplingUMesh *> tmp;
        convertPyObjToVecUMeshesCst(li,tmp);
        return MEDCouplingUMesh::mergeUMeshes(tmp);
      }

      PyObject *areCellsIncludedIn(const MEDCouplingUMesh *other, int compType) const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *ret1;
        bool ret0=self->areCellsIncludedIn(other,compType,ret1);
        PyObject *ret=PyTuple_New(2);
        PyObject *ret0Py=ret0?Py_True:Py_False;
        Py_XINCREF(ret0Py);
        PyTuple_SetItem(ret,0,ret0Py);
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }
    }
    void convertToPolyTypes(const std::vector<int>& cellIdsToConvert) throw(INTERP_KERNEL::Exception);
    void unPolyze() throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *buildExtrudedMeshFromThis(const MEDCouplingUMesh *mesh1D, int policy) throw(INTERP_KERNEL::Exception);
  };

  class MEDCouplingExtrudedMesh : public ParaMEDMEM::MEDCouplingMesh
  {
  public:
    static MEDCouplingExtrudedMesh *New(const MEDCouplingUMesh *mesh3D, MEDCouplingUMesh *mesh2D, int cell2DId) throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *build3DUnstructuredMesh() const throw(INTERP_KERNEL::Exception);
    %extend {
      std::string __str__() const
      {
        return self->simpleRepr();
      }
      PyObject *getMesh2D() const throw(INTERP_KERNEL::Exception)
      {
        MEDCouplingUMesh *ret=self->getMesh2D();
        ret->incrRef();
        return convertMesh(ret, SWIG_POINTER_OWN | 0 );
      }
      PyObject *getMesh1D() const throw(INTERP_KERNEL::Exception)
      {
        MEDCouplingUMesh *ret=self->getMesh1D();
        ret->incrRef();
        return convertMesh(ret, SWIG_POINTER_OWN | 0 );
      }
      PyObject *getMesh3DIds() const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *ret=self->getMesh3DIds();
        ret->incrRef();
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 );
      } 
    }
  };

  class MEDCouplingCMesh : public ParaMEDMEM::MEDCouplingMesh
  {
  public:
    static MEDCouplingCMesh *New();
    void setCoords(DataArrayDouble *coordsX,
                   DataArrayDouble *coordsY=0,
                   DataArrayDouble *coordsZ=0) throw(INTERP_KERNEL::Exception);
    void setCoordsAt(int i, DataArrayDouble *arr) throw(INTERP_KERNEL::Exception);
    %extend {
      std::string __str__() const
      {
        return self->simpleRepr();
      }
    }
  };
}

%extend ParaMEDMEM::DataArray
{
  void copyPartOfStringInfoFrom(const DataArray& other, PyObject *li) throw(INTERP_KERNEL::Exception)
  {
    std::vector<int> tmp;
    convertPyToNewIntArr3(li,tmp);
    self->copyPartOfStringInfoFrom(other,tmp);
  }

  void copyPartOfStringInfoFrom2(PyObject *li, const DataArray& other) throw(INTERP_KERNEL::Exception)
  {
    std::vector<int> tmp;
    convertPyToNewIntArr3(li,tmp);
    self->copyPartOfStringInfoFrom2(tmp,other);
  }
}

%extend ParaMEDMEM::DataArrayDouble
 {
   std::string __str__() const
   {
     return self->repr();
   }

   void setValues(PyObject *li, int nbOfTuples, int nbOfElsPerTuple) throw(INTERP_KERNEL::Exception)
   {
     int sz;
     double *tmp=convertPyToNewDblArr2(li,&sz);
     self->useArray(tmp,true,CPP_DEALLOC,nbOfTuples,nbOfElsPerTuple);
   }

   PyObject *getValues() throw(INTERP_KERNEL::Exception)
   {
     const double *vals=self->getPointer();
     return convertDblArrToPyList(vals,self->getNbOfElems());
   }

   PyObject *getValuesAsTuple() throw(INTERP_KERNEL::Exception)
   {
     const double *vals=self->getPointer();
     int nbOfComp=self->getNumberOfComponents();
     int nbOfTuples=self->getNumberOfTuples();
     return convertDblArrToPyListOfTuple(vals,nbOfComp,nbOfTuples);
   }

   DataArrayDouble *renumber(PyObject *li) throw(INTERP_KERNEL::Exception)
   {
     int size;
     INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
     if(size!=self->getNumberOfTuples())
       {
         throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
       }
     DataArrayDouble *ret=self->renumber(tmp);
     return ret;
   }

   DataArrayDouble *renumberR(PyObject *li) throw(INTERP_KERNEL::Exception)
   {
     int size;
     INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
     if(size!=self->getNumberOfTuples())
       {
         throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
       }
     DataArrayDouble *ret=self->renumberR(tmp);
     return ret;
   }

   DataArrayDouble *renumberAndReduce(PyObject *li, int newNbOfTuple) throw(INTERP_KERNEL::Exception)
   {
     int size;
     INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
     if(size!=self->getNumberOfTuples())
       {
         throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
       }
     DataArrayDouble *ret=self->renumberAndReduce(tmp,newNbOfTuple);
     return ret;
   }

   void renumberInPlace(PyObject *li) throw(INTERP_KERNEL::Exception)
   {
     int size;
     INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
     if(size!=self->getNumberOfTuples())
       {
         throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
       }
     self->renumberInPlace(tmp);
   }

   void renumberInPlaceR(PyObject *li) throw(INTERP_KERNEL::Exception)
   {
     int size;
     INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
     if(size!=self->getNumberOfTuples())
       {
         throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
       }
     self->renumberInPlaceR(tmp);
   }

   DataArrayDouble *selectByTupleId(PyObject *li) const throw(INTERP_KERNEL::Exception)
   {
     int size;
     INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
     DataArrayDouble *ret=self->selectByTupleId(tmp,tmp+size);
     return ret;
   }

   DataArrayDouble *selectByTupleIdSafe(PyObject *li) const throw(INTERP_KERNEL::Exception)
   {
     int size;
     INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
     DataArrayDouble *ret=self->selectByTupleIdSafe(tmp,tmp+size);
     return ret;
   }

   PyObject *getMaxValue() const throw(INTERP_KERNEL::Exception)
   {
     int tmp;
     double r1=self->getMaxValue(tmp);
     PyObject *ret=PyTuple_New(2);
     PyTuple_SetItem(ret,0,PyFloat_FromDouble(r1));
     PyTuple_SetItem(ret,1,PyInt_FromLong(tmp));
     return ret;
   }

   PyObject *getMaxValue2() const throw(INTERP_KERNEL::Exception)
   {
     DataArrayInt *tmp;
     double r1=self->getMaxValue2(tmp);
     PyObject *ret=PyTuple_New(2);
     PyTuple_SetItem(ret,0,PyFloat_FromDouble(r1));
     PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(tmp),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
     return ret;
   }

   PyObject *getMinValue() const throw(INTERP_KERNEL::Exception)
   {
     int tmp;
     double r1=self->getMinValue(tmp);
     PyObject *ret=PyTuple_New(2);
     PyTuple_SetItem(ret,0,PyFloat_FromDouble(r1));
     PyTuple_SetItem(ret,1,PyInt_FromLong(tmp));
     return ret;
   }

   PyObject *getMinValue2() const throw(INTERP_KERNEL::Exception)
   {
     DataArrayInt *tmp;
     double r1=self->getMinValue2(tmp);
     PyObject *ret=PyTuple_New(2);
     PyTuple_SetItem(ret,0,PyFloat_FromDouble(r1));
     PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(tmp),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
     return ret;
   }

   PyObject *accumulate() const throw(INTERP_KERNEL::Exception)
   {
     int sz=self->getNumberOfComponents();
     double *tmp=new double[sz];
     try
       {
         self->accumulate(tmp);
       }
     catch(INTERP_KERNEL::Exception& e)
       {
         delete [] tmp;
         throw e;
       }
     PyObject *ret=convertDblArrToPyList(tmp,sz);
     delete [] tmp;
     return ret;
   }
   
   DataArrayDouble *keepSelectedComponents(PyObject *li) const throw(INTERP_KERNEL::Exception)
   {
     std::vector<int> tmp;
     convertPyToNewIntArr3(li,tmp);
     return self->keepSelectedComponents(tmp);
   }

   void setSelectedComponents(const DataArrayDouble *a, PyObject *li) throw(INTERP_KERNEL::Exception)
   {
     std::vector<int> tmp;
     convertPyToNewIntArr3(li,tmp);
     self->setSelectedComponents(a,tmp);
   }
   
   PyObject *getTuple(int tupleId) throw(INTERP_KERNEL::Exception)
   {
     int sz=self->getNumberOfComponents();
     INTERP_KERNEL::AutoPtr<double> tmp=new double[sz];
     self->getTuple(tupleId,tmp);
     return convertDblArrToPyList(tmp,sz);
   }

   static DataArrayDouble *aggregate(PyObject *li) throw(INTERP_KERNEL::Exception)
   {
     std::vector<const DataArrayDouble *> tmp;
     convertPyObjToVecDataArrayDblCst(li,tmp);
     return DataArrayDouble::aggregate(tmp);
   }

   static DataArrayDouble *meld(PyObject *li) throw(INTERP_KERNEL::Exception)
   {
     std::vector<const DataArrayDouble *> tmp;
     convertPyObjToVecDataArrayDblCst(li,tmp);
     return DataArrayDouble::meld(tmp);
   }
 };

%extend ParaMEDMEM::DataArrayInt
 {
   std::string __str__() const
   {
     return self->repr();
   }

   void setValues(PyObject *li, int nbOfTuples, int nbOfElsPerTuple) throw(INTERP_KERNEL::Exception)
   {
     int size;
     int *tmp=convertPyToNewIntArr2(li,&size);
     self->useArray(tmp,true,CPP_DEALLOC,nbOfTuples,nbOfElsPerTuple);
   }

   PyObject *getValues() throw(INTERP_KERNEL::Exception)
   {
     const int *vals=self->getPointer();
     return convertIntArrToPyList(vals,self->getNbOfElems());
   }

   PyObject *getValuesAsTuple() throw(INTERP_KERNEL::Exception)
   {
     const int *vals=self->getPointer();
     int nbOfComp=self->getNumberOfComponents();
     int nbOfTuples=self->getNumberOfTuples();
     return convertIntArrToPyListOfTuple(vals,nbOfComp,nbOfTuples);
   }

   static PyObject *makePartition(PyObject *gps, int newNb) throw(INTERP_KERNEL::Exception)
   {
     std::vector<DataArrayInt *> groups;
     std::vector< std::vector<int> > fidsOfGroups;
     convertPyObjToVecDataArrayInt(gps,groups);
     ParaMEDMEM::DataArrayInt *ret0=ParaMEDMEM::DataArrayInt::makePartition(groups,newNb,fidsOfGroups);
     PyObject *ret = PyList_New(2);
     PyList_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
     int sz=fidsOfGroups.size();
     PyObject *ret1 = PyList_New(sz);
     for(int i=0;i<sz;i++)
       PyList_SetItem(ret1,i,convertIntArrToPyList2(fidsOfGroups[i]));
     PyList_SetItem(ret,1,ret1);
     return ret;
   }

   void renumberInPlace(PyObject *li) throw(INTERP_KERNEL::Exception)
   {
     int size;
     INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
     if(size!=self->getNumberOfTuples())
       {
         throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
       }
     self->renumberInPlace(tmp);
   }

   void renumberInPlaceR(PyObject *li) throw(INTERP_KERNEL::Exception)
   {
     int size;
     INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
     if(size!=self->getNumberOfTuples())
       {
         throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
       }
     self->renumberInPlaceR(tmp);
   }

   DataArrayInt *renumberAndReduce(PyObject *li, int newNbOfTuple) throw(INTERP_KERNEL::Exception)
   {
     int size;
     INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
     if(size!=self->getNumberOfTuples())
       {
         throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
       }
     DataArrayInt *ret=self->renumberAndReduce(tmp,newNbOfTuple);
     return ret;
   }

   DataArrayInt *renumber(PyObject *li) throw(INTERP_KERNEL::Exception)
   {
     int size;
     INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
     if(size!=self->getNumberOfTuples())
       {
         throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
       }
     DataArrayInt *ret=self->renumber(tmp);
     return ret;
   }

   DataArrayInt *renumberR(PyObject *li) throw(INTERP_KERNEL::Exception)
   {
     int size;
     INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
     if(size!=self->getNumberOfTuples())
       {
         throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
       }
     DataArrayInt *ret=self->renumberR(tmp);
     return ret;
   }

   DataArrayInt *selectByTupleId(PyObject *li) const throw(INTERP_KERNEL::Exception)
   {
     int size;
     INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
     DataArrayInt *ret=self->selectByTupleId(tmp,tmp+size);
     return ret;
   }

   DataArrayInt *selectByTupleIdSafe(PyObject *li) const throw(INTERP_KERNEL::Exception)
   {
     int size;
     INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
     DataArrayInt *ret=self->selectByTupleIdSafe(tmp,tmp+size);
     return ret;
   }

   DataArrayInt *keepSelectedComponents(PyObject *li) const throw(INTERP_KERNEL::Exception)
   {
     std::vector<int> tmp;
     convertPyToNewIntArr3(li,tmp);
     return self->keepSelectedComponents(tmp);
   }

   void setSelectedComponents(const DataArrayInt *a, PyObject *li) throw(INTERP_KERNEL::Exception)
   {
     std::vector<int> tmp;
     convertPyToNewIntArr3(li,tmp);
     self->setSelectedComponents(a,tmp);
   }

   PyObject *getTuple(int tupleId) throw(INTERP_KERNEL::Exception)
   {
     int sz=self->getNumberOfComponents();
     INTERP_KERNEL::AutoPtr<int> tmp=new int[sz];
     self->getTuple(tupleId,tmp);
     return convertIntArrToPyList(tmp,sz);
   }

   static DataArrayInt *meld(PyObject *li) throw(INTERP_KERNEL::Exception)
   {
     std::vector<const DataArrayInt *> tmp;
     convertPyObjToVecDataArrayIntCst(li,tmp);
     return DataArrayInt::meld(tmp);
   }

   PyObject *getMaxValue() const throw(INTERP_KERNEL::Exception)
   {
     int tmp;
     int r1=self->getMaxValue(tmp);
     PyObject *ret=PyTuple_New(2);
     PyTuple_SetItem(ret,0,PyInt_FromLong(r1));
     PyTuple_SetItem(ret,1,PyInt_FromLong(tmp));
     return ret;
   }

   PyObject *getMinValue() const throw(INTERP_KERNEL::Exception)
   {
     int tmp;
     int r1=self->getMinValue(tmp);
     PyObject *ret=PyTuple_New(2);
     PyTuple_SetItem(ret,0,PyInt_FromLong(r1));
     PyTuple_SetItem(ret,1,PyInt_FromLong(tmp));
     return ret;
   }
 };

namespace ParaMEDMEM
{
  class MEDCouplingField : public ParaMEDMEM::RefCountObject, public ParaMEDMEM::TimeLabel
  {
  public:
    virtual void checkCoherency() const throw(INTERP_KERNEL::Exception);
    virtual bool areCompatibleForMerge(const MEDCouplingField *other) const throw(INTERP_KERNEL::Exception);
    virtual bool isEqual(const MEDCouplingField *other, double meshPrec, double valsPrec) const throw(INTERP_KERNEL::Exception);
    virtual bool isEqualWithoutConsideringStr(const MEDCouplingField *other, double meshPrec, double valsPrec) const throw(INTERP_KERNEL::Exception);
    void setMesh(const ParaMEDMEM::MEDCouplingMesh *mesh) throw(INTERP_KERNEL::Exception);
    void setName(const char *name) throw(INTERP_KERNEL::Exception);
    const char *getDescription() const throw(INTERP_KERNEL::Exception);
    void setDescription(const char *desc) throw(INTERP_KERNEL::Exception);
    const char *getName() const throw(INTERP_KERNEL::Exception);
    TypeOfField getTypeOfField() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *buildMeasureField(bool isAbs) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDiscretization *getDiscretization() const throw(INTERP_KERNEL::Exception);
    void setGaussLocalizationOnType(INTERP_KERNEL::NormalizedCellType type, const std::vector<double>& refCoo,
                                    const std::vector<double>& gsCoo, const std::vector<double>& wg) throw(INTERP_KERNEL::Exception);
    void clearGaussLocalizations() throw(INTERP_KERNEL::Exception);
    MEDCouplingGaussLocalization& getGaussLocalization(int locId) throw(INTERP_KERNEL::Exception);
    int getNbOfGaussLocalization() const throw(INTERP_KERNEL::Exception);
    int getGaussLocalizationIdOfOneCell(int cellId) const throw(INTERP_KERNEL::Exception);
    const MEDCouplingGaussLocalization& getGaussLocalization(int locId) const throw(INTERP_KERNEL::Exception);
    %extend {
      PyObject *getMesh() const throw(INTERP_KERNEL::Exception)
      {
        MEDCouplingMesh *ret1=(MEDCouplingMesh *)self->getMesh();
        ret1->incrRef();
        return convertMesh(ret1, SWIG_POINTER_OWN | 0 );
      }

      PyObject *buildSubMeshData(PyObject *li) const throw(INTERP_KERNEL::Exception)
      {
        int size;
        INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
        DataArrayInt *ret1;
        MEDCouplingMesh *ret0=self->buildSubMeshData(tmp,tmp+size,ret1);
        PyObject *res = PyList_New(2);
        PyList_SetItem(res,0,convertMesh(ret0, SWIG_POINTER_OWN | 0 ));
        PyList_SetItem(res,1,SWIG_NewPointerObj((void*)ret1,SWIGTYPE_p_ParaMEDMEM__DataArrayInt,SWIG_POINTER_OWN | 0));
        return res;
      }
      void setGaussLocalizationOnCells(PyObject *li, const std::vector<double>& refCoo,
                                       const std::vector<double>& gsCoo, const std::vector<double>& wg) throw(INTERP_KERNEL::Exception)
      {
        int size;
        INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
        try
          {
            self->setGaussLocalizationOnCells(tmp,tmp+size,refCoo,gsCoo,wg);
          }
        catch(INTERP_KERNEL::Exception& e)
          {
            throw e;
          }
      }
      PyObject *getCellIdsHavingGaussLocalization(int locId) const throw(INTERP_KERNEL::Exception)
      {
        std::vector<int> tmp;
        self->getCellIdsHavingGaussLocalization(locId,tmp);
        return convertIntArrToPyList2(tmp);
      }
    }
  };

  class MEDCouplingFieldDouble : public ParaMEDMEM::MEDCouplingField
  {
  public:
    static MEDCouplingFieldDouble *New(TypeOfField type, TypeOfTimeDiscretization td=NO_TIME);
    void copyTinyStringsFrom(const MEDCouplingFieldDouble *other) throw(INTERP_KERNEL::Exception);
    void copyTinyAttrFrom(const MEDCouplingFieldDouble *other) throw(INTERP_KERNEL::Exception);
    std::string simpleRepr() const;
    std::string advancedRepr() const;
    MEDCouplingFieldDouble *clone(bool recDeepCpy) const;
    MEDCouplingFieldDouble *cloneWithMesh(bool recDeepCpy) const;
    MEDCouplingFieldDouble *deepCpy() const;
    MEDCouplingFieldDouble *buildNewTimeReprFromThis(TypeOfTimeDiscretization td, bool deepCpy) const throw(INTERP_KERNEL::Exception);
    TypeOfTimeDiscretization getTimeDiscretization() const throw(INTERP_KERNEL::Exception);
    void checkCoherency() const throw(INTERP_KERNEL::Exception);
    double getIJ(int tupleId, int compoId) const throw(INTERP_KERNEL::Exception);
    double getIJK(int cellId, int nodeIdInCell, int compoId) const throw(INTERP_KERNEL::Exception);
    void setArray(DataArrayDouble *array) throw(INTERP_KERNEL::Exception);
    void setEndArray(DataArrayDouble *array) throw(INTERP_KERNEL::Exception);
    void setTime(double val, int iteration, int order) throw(INTERP_KERNEL::Exception);
    void setStartTime(double val, int iteration, int order) throw(INTERP_KERNEL::Exception);
    void setEndTime(double val, int iteration, int order) throw(INTERP_KERNEL::Exception);
    void applyLin(double a, double b, int compoId) throw(INTERP_KERNEL::Exception);
    int getNumberOfComponents() const throw(INTERP_KERNEL::Exception);
    int getNumberOfTuples() const throw(INTERP_KERNEL::Exception);
    int getNumberOfValues() const throw(INTERP_KERNEL::Exception);
    NatureOfField getNature() const throw(INTERP_KERNEL::Exception);
    void setNature(NatureOfField nat) throw(INTERP_KERNEL::Exception);
    void setTimeTolerance(double val) throw(INTERP_KERNEL::Exception);
    double getTimeTolerance() const throw(INTERP_KERNEL::Exception);
    void updateTime() throw(INTERP_KERNEL::Exception);
    void changeUnderlyingMesh(const MEDCouplingMesh *other, int levOfCheck, double prec) throw(INTERP_KERNEL::Exception);
    void substractInPlaceDM(const MEDCouplingFieldDouble *f, int levOfCheck, double prec) throw(INTERP_KERNEL::Exception);
    bool mergeNodes(double eps) throw(INTERP_KERNEL::Exception);
    bool mergeNodes2(double eps) throw(INTERP_KERNEL::Exception);
    bool zipCoords() throw(INTERP_KERNEL::Exception);
    bool zipConnectivity(int compType) throw(INTERP_KERNEL::Exception);
    bool simplexize(int policy) throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *doublyContractedProduct() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *determinant() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *eigenValues() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *eigenVectors() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *inverse() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *trace() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *deviator() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *magnitude() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *maxPerTuple() const throw(INTERP_KERNEL::Exception);
    void changeNbOfComponents(int newNbOfComp, double dftValue=0.) throw(INTERP_KERNEL::Exception);
    void sortPerTuple(bool asc) throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble &operator=(double value) throw(INTERP_KERNEL::Exception);
    void fillFromAnalytic(int nbOfComp, const char *func) throw(INTERP_KERNEL::Exception);
    void applyFunc(int nbOfComp, const char *func) throw(INTERP_KERNEL::Exception);
    void applyFunc(int nbOfComp, double val) throw(INTERP_KERNEL::Exception);
    void applyFunc(const char *func) throw(INTERP_KERNEL::Exception);
    void applyFuncFast32(const char *func) throw(INTERP_KERNEL::Exception);
    void applyFuncFast64(const char *func) throw(INTERP_KERNEL::Exception);
    double accumulate(int compId) const throw(INTERP_KERNEL::Exception);
    double getMaxValue() const throw(INTERP_KERNEL::Exception);
    double getMinValue() const throw(INTERP_KERNEL::Exception);
    double getAverageValue() const throw(INTERP_KERNEL::Exception);
    double getWeightedAverageValue() const throw(INTERP_KERNEL::Exception);
    double integral(int compId, bool isWAbs) const throw(INTERP_KERNEL::Exception);
    double normL1(int compId) const throw(INTERP_KERNEL::Exception);
    double normL2(int compId) const throw(INTERP_KERNEL::Exception);
    DataArrayInt *getIdsInRange(double vmin, double vmax) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *buildSubPart(const DataArrayInt *part) const throw(INTERP_KERNEL::Exception);
    static MEDCouplingFieldDouble *mergeFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception);
    static MEDCouplingFieldDouble *meldFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception);
    static MEDCouplingFieldDouble *dotFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *dot(const MEDCouplingFieldDouble& other) const throw(INTERP_KERNEL::Exception);
    static MEDCouplingFieldDouble *crossProductFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *crossProduct(const MEDCouplingFieldDouble& other) const throw(INTERP_KERNEL::Exception);
    static MEDCouplingFieldDouble *maxFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *max(const MEDCouplingFieldDouble& other) const throw(INTERP_KERNEL::Exception);
    static MEDCouplingFieldDouble *minFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *min(const MEDCouplingFieldDouble& other) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *operator+(const MEDCouplingFieldDouble& other) const throw(INTERP_KERNEL::Exception);
    const MEDCouplingFieldDouble &operator+=(const MEDCouplingFieldDouble& other) throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *operator-(const MEDCouplingFieldDouble& other) const throw(INTERP_KERNEL::Exception);
    const MEDCouplingFieldDouble &operator-=(const MEDCouplingFieldDouble& other) throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *operator*(const MEDCouplingFieldDouble& other) const throw(INTERP_KERNEL::Exception);
    const MEDCouplingFieldDouble &operator*=(const MEDCouplingFieldDouble& other) throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *operator/(const MEDCouplingFieldDouble& other) const throw(INTERP_KERNEL::Exception);
    const MEDCouplingFieldDouble &operator/=(const MEDCouplingFieldDouble& other) throw(INTERP_KERNEL::Exception);
    %extend {
      std::string __str__() const
      {
        return self->simpleRepr();
      }

      DataArrayDouble *getArray() const throw(INTERP_KERNEL::Exception)
      {
        DataArrayDouble *ret=self->getArray();
        if(ret)
          ret->incrRef();
        return ret;
      }

      DataArrayDouble *getEndArray() const throw(INTERP_KERNEL::Exception)
      {
        DataArrayDouble *ret=self->getEndArray();
        if(ret)
          ret->incrRef();
        return ret;
      }

      PyObject *getValueOn(PyObject *sl) const throw(INTERP_KERNEL::Exception)
      {
        int sz;
        double *spaceLoc=convertPyToNewDblArr2(sl,&sz);
        sz=self->getNumberOfComponents();
        double *res=new double[sz];
        try
          {
            self->getValueOn(spaceLoc,res);
          }
        catch(INTERP_KERNEL::Exception& e)
          {
            delete [] spaceLoc;
            delete [] res;
            throw e;
          }
        delete [] spaceLoc;
        PyObject *ret=convertDblArrToPyList(res,sz);
        delete [] res;
        return ret;
      }
      PyObject *getValueOn(PyObject *sl, double time) const throw(INTERP_KERNEL::Exception)
      {
        int sz;
        double *spaceLoc=convertPyToNewDblArr2(sl,&sz);
        sz=self->getNumberOfComponents();
        double *res=new double[sz];
        try
          {
            self->getValueOn(spaceLoc,time,res);
          }
        catch(INTERP_KERNEL::Exception& e)
          {
            delete [] spaceLoc;
            delete [] res;
            throw e;
          }
        delete [] spaceLoc;
        PyObject *ret=convertDblArrToPyList(res,sz);
        delete [] res;
        return ret;
      }
      void setValues(PyObject *li) throw(INTERP_KERNEL::Exception)
      {
        if(self->getArray()!=0)
          {
            int sz;
            double *tmp=convertPyToNewDblArr2(li,&sz);
            int nbTuples=self->getArray()->getNumberOfTuples();
            int nbOfCompo=self->getArray()->getNumberOfComponents();
            self->getArray()->useArray(tmp,true,CPP_DEALLOC,nbTuples,nbOfCompo);
          }
        else
          PyErr_SetString(PyExc_TypeError,"setValuesCpy : field must contain an array behind");
      }
      PyObject *getTime() throw(INTERP_KERNEL::Exception)
      {
        int tmp1,tmp2;
        double tmp0=self->getTime(tmp1,tmp2);
        PyObject *res = PyList_New(3);
        PyList_SetItem(res,0,SWIG_From_double(tmp0));
        PyList_SetItem(res,1,SWIG_From_int(tmp1));
        PyList_SetItem(res,2,SWIG_From_int(tmp2));
        return res;
      }

      PyObject *getStartTime() throw(INTERP_KERNEL::Exception)
      {
        int tmp1,tmp2;
        double tmp0=self->getStartTime(tmp1,tmp2);
        PyObject *res = PyList_New(3);
        PyList_SetItem(res,0,SWIG_From_double(tmp0));
        PyList_SetItem(res,1,SWIG_From_int(tmp1));
        PyList_SetItem(res,2,SWIG_From_int(tmp2));
        return res;
      }

      PyObject *getEndTime() throw(INTERP_KERNEL::Exception)
      {
        int tmp1,tmp2;
        double tmp0=self->getEndTime(tmp1,tmp2);
        PyObject *res = PyList_New(3);
        PyList_SetItem(res,0,SWIG_From_double(tmp0));
        PyList_SetItem(res,1,SWIG_From_int(tmp1));
        PyList_SetItem(res,2,SWIG_From_int(tmp2));
        return res;
      }
      PyObject *accumulate() const throw(INTERP_KERNEL::Exception)
      {
        int sz=self->getNumberOfComponents();
        double *tmp=new double[sz];
        self->accumulate(tmp);
        PyObject *ret=convertDblArrToPyList(tmp,sz);
        delete [] tmp;
        return ret;
      }
      PyObject *integral(bool isWAbs) const throw(INTERP_KERNEL::Exception)
      {
        int sz=self->getNumberOfComponents();
        double *tmp=new double[sz];
        self->integral(isWAbs,tmp);
        PyObject *ret=convertDblArrToPyList(tmp,sz);
        delete [] tmp;
        return ret;
      }
      PyObject *normL1() const throw(INTERP_KERNEL::Exception)
      {
        int sz=self->getNumberOfComponents();
        double *tmp=new double[sz];
        self->normL1(tmp);
        PyObject *ret=convertDblArrToPyList(tmp,sz);
        delete [] tmp;
        return ret;
      }
      PyObject *normL2() const throw(INTERP_KERNEL::Exception)
      {
        int sz=self->getNumberOfComponents();
        double *tmp=new double[sz];
        self->normL2(tmp);
        PyObject *ret=convertDblArrToPyList(tmp,sz);
        delete [] tmp;
        return ret;
      }

      void renumberCells(PyObject *li, bool check) throw(INTERP_KERNEL::Exception)
      {
        int size;
        INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
        try
          {
            self->renumberCells(tmp,check);
          }
        catch(INTERP_KERNEL::Exception& e)
          {
            throw e;
          }
      }
      void renumberNodes(PyObject *li) throw(INTERP_KERNEL::Exception)
      {
        int size;
        INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
        try
          {
            self->renumberNodes(tmp);
          }
        catch(INTERP_KERNEL::Exception& e)
          {
            throw e;
          }
      }

      MEDCouplingFieldDouble *buildSubPart(PyObject *li) const throw(INTERP_KERNEL::Exception)
      {
        int size;
        INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
        MEDCouplingFieldDouble *ret=0;
        try
          {
            ret=self->buildSubPart(tmp,tmp+size);
          }
        catch(INTERP_KERNEL::Exception& e)
          {
            throw e;
          }
        return ret;
      }

      PyObject *getMaxValue2() const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *tmp;
        double r1=self->getMaxValue2(tmp);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,PyFloat_FromDouble(r1));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(tmp),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }
      
      PyObject *getMinValue2() const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *tmp;
        double r1=self->getMinValue2(tmp);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,PyFloat_FromDouble(r1));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(tmp),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }
      
      MEDCouplingFieldDouble *keepSelectedComponents(PyObject *li) const throw(INTERP_KERNEL::Exception)
      {
        std::vector<int> tmp;
        convertPyToNewIntArr3(li,tmp);
        return self->keepSelectedComponents(tmp);
      }

      void setSelectedComponents(const MEDCouplingFieldDouble *f, PyObject *li) throw(INTERP_KERNEL::Exception)
      {
        std::vector<int> tmp;
        convertPyToNewIntArr3(li,tmp);
        self->setSelectedComponents(f,tmp);
      }

      static MEDCouplingFieldDouble *mergeFields(PyObject *li) throw(INTERP_KERNEL::Exception)
      {
        std::vector<const MEDCouplingFieldDouble *> tmp;
        convertPyObjToVecFieldDblCst(li,tmp);
        return MEDCouplingFieldDouble::mergeFields(tmp);
      }
    }
  };
}
