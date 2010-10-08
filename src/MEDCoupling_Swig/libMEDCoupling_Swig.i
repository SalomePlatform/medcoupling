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

%module libMEDCoupling_Swig

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

%newobject ParaMEDMEM::DataArrayDouble::New;
%newobject ParaMEDMEM::DataArrayInt::New;
%newobject ParaMEDMEM::DataArrayDouble::convertToIntArr;
%newobject ParaMEDMEM::DataArrayInt::convertToDblArr;
%newobject ParaMEDMEM::MEDCouplingUMesh::New;
%newobject ParaMEDMEM::MEDCouplingField::buildWeightingField;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::New;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::mergeFields;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::doublyContractedProduct;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::determinant;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::eigenValues;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::eigenVectors;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::inverse;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::trace;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::deviator;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::magnitude;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::maxPerTuple;
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
%newobject ParaMEDMEM::MEDCouplingUMesh::clone;
%newobject ParaMEDMEM::DataArrayDouble::deepCopy;
%newobject ParaMEDMEM::DataArrayDouble::performCpy;
%newobject ParaMEDMEM::DataArrayInt::deepCopy;
%newobject ParaMEDMEM::DataArrayInt::performCpy;
%newobject ParaMEDMEM::DataArrayInt::substr;
%newobject ParaMEDMEM::DataArrayInt::changeNbOfComponents;
%newobject ParaMEDMEM::DataArrayInt::selectByTupleId;
%newobject ParaMEDMEM::DataArrayDouble::aggregate;
%newobject ParaMEDMEM::DataArrayDouble::dot;
%newobject ParaMEDMEM::DataArrayDouble::crossProduct;
%newobject ParaMEDMEM::DataArrayDouble::add;
%newobject ParaMEDMEM::DataArrayDouble::substract;
%newobject ParaMEDMEM::DataArrayDouble::multiply;
%newobject ParaMEDMEM::DataArrayDouble::divide;
%newobject ParaMEDMEM::DataArrayDouble::substr;
%newobject ParaMEDMEM::DataArrayDouble::changeNbOfComponents;
%newobject ParaMEDMEM::DataArrayDouble::getIdsInRange;
%newobject ParaMEDMEM::DataArrayDouble::selectByTupleId;
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
%newobject ParaMEDMEM::MEDCouplingFieldDouble::clone;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::cloneWithMesh;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::buildNewTimeReprFromThis;
%newobject ParaMEDMEM::MEDCouplingMesh::getCoordinatesAndOwner;
%newobject ParaMEDMEM::MEDCouplingMesh::getBarycenterAndOwner;
%newobject ParaMEDMEM::MEDCouplingMesh::buildOrthogonalField;
%newobject ParaMEDMEM::MEDCouplingMesh::getCellIdsFullyIncludedInNodeIds;
%newobject ParaMEDMEM::MEDCouplingMesh::buildPart;
%newobject ParaMEDMEM::MEDCouplingMesh::mergeMyselfWith;
%newobject ParaMEDMEM::MEDCouplingMesh::fillFromAnalytic;
%newobject ParaMEDMEM::MEDCouplingPointSet::zipCoordsTraducer;
%newobject ParaMEDMEM::MEDCouplingUMesh::zipConnectivityTraducer;
%newobject ParaMEDMEM::MEDCouplingUMesh::buildDescendingConnectivity;
%newobject ParaMEDMEM::MEDCouplingPointSet::buildBoundaryMesh;
%newobject ParaMEDMEM::MEDCouplingMesh::getMeasureField;
%newobject ParaMEDMEM::MEDCouplingUMesh::buildExtrudedMeshFromThis;
%newobject ParaMEDMEM::MEDCouplingUMesh::mergeUMeshes;
%newobject ParaMEDMEM::MEDCouplingUMesh::buildNewNumberingFromCommNodesFrmt;
%newobject ParaMEDMEM::MEDCouplingUMesh::rearrange2ConsecutiveCellTypes;
%newobject ParaMEDMEM::MEDCouplingUMesh::convertCellArrayPerGeoType;
%newobject ParaMEDMEM::MEDCouplingUMesh::getRenumArrForConsctvCellTypesSpe;
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

%ignore ParaMEDMEM::TimeLabel::operator=;
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
  class MEDCouplingFieldDouble;

  class MEDCouplingMesh : public RefCountObject, public TimeLabel
  {
  public:
    void setName(const char *name) { _name=name; }
    const char *getName() const { return _name.c_str(); }
    virtual MEDCouplingMeshType getType() const = 0;
    bool isStructured() const;
    virtual bool isEqual(const MEDCouplingMesh *other, double prec) const;
    virtual bool isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const = 0;
    virtual void copyTinyStringsFrom(const MEDCouplingMesh *other) throw(INTERP_KERNEL::Exception);
    virtual void checkCoherency() const throw(INTERP_KERNEL::Exception) = 0;
    virtual int getNumberOfCells() const throw(INTERP_KERNEL::Exception) = 0;
    virtual int getNumberOfNodes() const throw(INTERP_KERNEL::Exception) = 0;
    virtual int getSpaceDimension() const throw(INTERP_KERNEL::Exception) = 0;
    virtual int getMeshDimension() const throw(INTERP_KERNEL::Exception) = 0;
    virtual DataArrayDouble *getCoordinatesAndOwner() const = 0;
    virtual DataArrayDouble *getBarycenterAndOwner() const = 0;
    virtual INTERP_KERNEL::NormalizedCellType getTypeOfCell(int cellId) const = 0;
    virtual std::string simpleRepr() const = 0;
    virtual std::string advancedRepr() const = 0;
    // tools
    virtual void getBoundingBox(double *bbox) const = 0;
    virtual MEDCouplingFieldDouble *getMeasureField(bool isAbs) const = 0;
    virtual MEDCouplingFieldDouble *getMeasureFieldOnNode(bool isAbs) const = 0;
    virtual int getCellContainingPoint(const double *pos, double eps) const = 0;
    virtual MEDCouplingFieldDouble *fillFromAnalytic(TypeOfField t, int nbOfComp, const char *func) const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingFieldDouble *buildOrthogonalField() const = 0;
    virtual void rotate(const double *center, const double *vector, double angle) = 0;
    virtual void translate(const double *vector) = 0;
    virtual MEDCouplingMesh *buildPart(const int *start, const int *end) const = 0;
    virtual MEDCouplingMesh *mergeMyselfWith(const MEDCouplingMesh *other) const throw(INTERP_KERNEL::Exception) = 0;
    virtual bool areCompatibleForMerge(const MEDCouplingMesh *other) const;
    static MEDCouplingMesh *mergeMeshes(const MEDCouplingMesh *mesh1, const MEDCouplingMesh *mesh2);
    %extend
       {
         std::string __str__() const
         {
           return self->simpleRepr();
         }

         void renumberCells(PyObject *li, bool check) throw(INTERP_KERNEL::Exception)
         {
           int size;
           int *tmp=convertPyToNewIntArr2(li,&size);
           self->renumberCells(tmp,check);
           delete [] tmp;
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
         DataArrayInt *getCellIdsFullyIncludedInNodeIds(PyObject *li) const
         {
           int size;
           int *tmp=convertPyToNewIntArr2(li,&size);
           DataArrayInt *ret=self->getCellIdsFullyIncludedInNodeIds(tmp,tmp+size);
           delete [] tmp;
           return ret;
         }
         PyObject *getNodeIdsOfCell(int cellId) const
         {
           std::vector<int> conn;
           self->getNodeIdsOfCell(cellId,conn);
           return convertIntArrToPyList2(conn);
         }

         PyObject *getCoordinatesOfNode(int nodeId) const
         {
           std::vector<double> coo;
           self->getCoordinatesOfNode(nodeId,coo);
           return convertDblArrToPyList2(coo);
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
      void setCoords(DataArrayDouble *coords);
      DataArrayDouble *getCoordinatesAndOwner() const;
      bool areCoordsEqual(const MEDCouplingPointSet& other, double prec) const;
      void getBoundingBox(double *bbox) const;
      void zipCoords();
      double getCaracteristicDimension() const;
      void translate(const double *vector);
      void scale(const double *point, double factor);
      void changeSpaceDimension(int newSpaceDim, double dftVal=0.) throw(INTERP_KERNEL::Exception);
      void tryToShareSameCoords(const MEDCouplingPointSet& other, double epsilon) throw(INTERP_KERNEL::Exception);
      virtual void tryToShareSameCoordsPermute(const MEDCouplingPointSet& other, double epsilon) throw(INTERP_KERNEL::Exception) = 0;
      static DataArrayDouble *mergeNodesArray(const MEDCouplingPointSet *m1, const MEDCouplingPointSet *m2);
      static MEDCouplingPointSet *buildInstanceFromMeshType(MEDCouplingMeshType type);
      virtual MEDCouplingPointSet *buildPartOfMySelf(const int *start, const int *end, bool keepCoords) const = 0;
      virtual MEDCouplingPointSet *buildPartOfMySelfNode(const int *start, const int *end, bool fullyIn) const = 0;
      virtual MEDCouplingPointSet *buildFacePartOfMySelfNode(const int *start, const int *end, bool fullyIn) const = 0;
      virtual MEDCouplingPointSet *buildBoundaryMesh(bool keepCoords) const = 0;
      virtual void renumberNodes(const int *newNodeNumbers, int newNbOfNodes);
      virtual bool isEmptyMesh(const std::vector<int>& tinyInfo) const = 0;
      //! size of returned tinyInfo must be always the same.
      void getTinySerializationInformation(std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const;
      void resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings) const;
      void serialize(DataArrayInt *&a1, DataArrayDouble *&a2) const;
      void unserialization(const std::vector<int>& tinyInfo, const DataArrayInt *a1, DataArrayDouble *a2,
                           const std::vector<std::string>& littleStrings);
      virtual void giveElemsInBoundingBox(const double *bbox, double eps, std::vector<int>& elems) = 0;
      virtual void giveElemsInBoundingBox(const INTERP_KERNEL::DirectedBoundingBox& bbox, double eps, std::vector<int>& elems) = 0;
      virtual DataArrayInt *zipCoordsTraducer() = 0;
      %extend 
         {
           std::string __str__() const
           {
             return self->simpleRepr();
           }
           
           PyObject *buildNewNumberingFromCommNodesFrmt(const DataArrayInt *comm, const DataArrayInt *commIndex) const
           {
             int newNbOfNodes;
             DataArrayInt *ret0=self->buildNewNumberingFromCommNodesFrmt(comm,commIndex,newNbOfNodes);
             PyObject *res = PyList_New(2);
             PyList_SetItem(res,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
             PyList_SetItem(res,1,SWIG_From_int(newNbOfNodes));
             return res;
           }
           
           PyObject *findCommonNodes(int limitNodeId, double prec) const
           {
             DataArrayInt *comm, *commIndex;
             self->findCommonNodes(limitNodeId,prec,comm,commIndex);
             PyObject *res = PyList_New(2);
             PyList_SetItem(res,0,SWIG_NewPointerObj(SWIG_as_voidptr(comm),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
             PyList_SetItem(res,1,SWIG_NewPointerObj(SWIG_as_voidptr(commIndex),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
             return res;
           }
           
           PyObject *getCoords() const
           {
             DataArrayDouble *ret1=self->getCoords();
             ret1->incrRef();
             return SWIG_NewPointerObj((void*)ret1,SWIGTYPE_p_ParaMEDMEM__DataArrayDouble,SWIG_POINTER_OWN | 0);
           }
           PyObject *buildPartOfMySelf(PyObject *li, bool keepCoords) const
           {
             int size;
             int *tmp=convertPyToNewIntArr2(li,&size);
             MEDCouplingPointSet *ret=self->buildPartOfMySelf(tmp,tmp+size,keepCoords);
             delete [] tmp;
             return convertMesh(ret, SWIG_POINTER_OWN | 0 );
           }
           PyObject *buildPartOfMySelfNode(PyObject *li, bool fullyIn) const
           {
             int size;
             int *tmp=convertPyToNewIntArr2(li,&size);
             MEDCouplingPointSet *ret=self->buildPartOfMySelfNode(tmp,tmp+size,fullyIn);
             delete [] tmp;
             return convertMesh(ret, SWIG_POINTER_OWN | 0 );
           }
           PyObject *buildFacePartOfMySelfNode(PyObject *li, bool fullyIn) const
           {
             int size;
             int *tmp=convertPyToNewIntArr2(li,&size);
             MEDCouplingPointSet *ret=self->buildFacePartOfMySelfNode(tmp,tmp+size,fullyIn);
             delete [] tmp;
             return convertMesh(ret, SWIG_POINTER_OWN | 0 );
           }
           PyObject *findBoundaryNodes() const
           {
             std::vector<int> nodes;
             self->findBoundaryNodes(nodes);
             return convertIntArrToPyList2(nodes);
           }
           void rotate(PyObject *center, PyObject *vector, double alpha)
           {
             int sz;
             double *c=convertPyToNewDblArr2(center,&sz);
             if(!c)
               return ;
             double *v=convertPyToNewDblArr2(vector,&sz);
             if(!v)
               { delete [] c; return ; }
             self->rotate(c,v,alpha);
             delete [] c;
             delete [] v;
           }
           void translate(PyObject *vector)
           {
             int sz;
             double *v=convertPyToNewDblArr2(vector,&sz);
             self->translate(v);
             delete [] v;
           }
           void scale(PyObject *point, double factor)
           {
             int sz;
             double *p=convertPyToNewDblArr2(point,&sz);
             self->scale(p,factor);
             delete [] p;
           }
           void renumberNodes(PyObject *li, int newNbOfNodes)
           {
             int size;
             int *tmp=convertPyToNewIntArr2(li,&size);
             self->renumberNodes(tmp,newNbOfNodes);
             delete [] tmp;
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
           static void rotate2DAlg(PyObject *center, double angle, int nbNodes, PyObject *coords)
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
           static void rotate3DAlg(PyObject *center, PyObject *vect, double angle, int nbNodes, PyObject *coords)
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
    void allocateCells(int nbOfCells);
    void finishInsertingCells();
    void setConnectivity(DataArrayInt *conn, DataArrayInt *connIndex, bool isComputingTypes=true);
    DataArrayInt *getNodalConnectivity() const;
    DataArrayInt *getNodalConnectivityIndex() const;
    INTERP_KERNEL::NormalizedCellType getTypeOfCell(int cellId) const;
    int getNumberOfNodesInCell(int cellId) const;
    int getMeshLength() const;
    void computeTypes();
    std::string reprConnectivityOfThis() const;
    //tools
    bool checkConsecutiveCellTypes() const;
    DataArrayInt *rearrange2ConsecutiveCellTypes();
    DataArrayInt *convertCellArrayPerGeoType(const DataArrayInt *da) const throw(INTERP_KERNEL::Exception);
    DataArrayInt *zipConnectivityTraducer(int compType);
    void getReverseNodalConnectivity(DataArrayInt *revNodal, DataArrayInt *revNodalIndx) const;
    MEDCouplingUMesh *buildDescendingConnectivity(DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *revDesc, DataArrayInt *revDescIndx) const;
    void orientCorrectlyPolyhedrons() throw(INTERP_KERNEL::Exception);
    bool isPresenceOfQuadratic() const;
    void convertQuadraticCellsToLinear() throw(INTERP_KERNEL::Exception);
    %extend {
      std::string __str__() const
      {
        return self->simpleRepr();
      }

      int getCellContainingPoint(PyObject *p, double eps) const
      {
        int sz;
        double *pos=convertPyToNewDblArr2(p,&sz);
        int ret=self->getCellContainingPoint(pos,eps);
        delete [] pos;
        return ret;
      }
      void insertNextCell(INTERP_KERNEL::NormalizedCellType type, int size, PyObject *li)
      {
        int sz;
        int *tmp=convertPyToNewIntArr2(li,&sz);
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
        int ret2;
        DataArrayInt *ret0=self->mergeNodes(precision,ret1,ret2);
        PyObject *res = PyList_New(3);
        PyList_SetItem(res,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyList_SetItem(res,1,SWIG_From_bool(ret1));
        PyList_SetItem(res,2,SWIG_From_int(ret2));
        return res;
      }
      PyObject *checkButterflyCells()
      {
        std::vector<int> cells;
        self->checkButterflyCells(cells);
        return convertIntArrToPyList2(cells);
      }

      PyObject *splitByType() const
      {
        std::vector<MEDCouplingUMesh *> ms=self->splitByType();
        int sz=ms.size();
        PyObject *ret = PyList_New(sz);
        for(int i=0;i<sz;i++)
          PyList_SetItem(ret,i,SWIG_NewPointerObj(SWIG_as_voidptr(ms[i]),SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *keepSpecifiedCells(INTERP_KERNEL::NormalizedCellType type, PyObject *ids) const
      {
        std::vector<int> idsPerGeoType;
        convertPyToNewIntArr3(ids,idsPerGeoType);
        MEDCouplingUMesh *ret=self->keepSpecifiedCells(type,idsPerGeoType);
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 );
      }

      bool checkConsecutiveCellTypesAndOrder(PyObject *li) const
      {
        int sz;
        INTERP_KERNEL::NormalizedCellType *order=(INTERP_KERNEL::NormalizedCellType *)convertPyToNewIntArr2(li,&sz);
        bool ret=self->checkConsecutiveCellTypesAndOrder(order,order+sz);
        delete [] order;
        return ret;
      }

      DataArrayInt *getRenumArrForConsctvCellTypesSpe(PyObject *li) const
      {
        int sz;
        INTERP_KERNEL::NormalizedCellType *order=(INTERP_KERNEL::NormalizedCellType *)convertPyToNewIntArr2(li,&sz);
        DataArrayInt *ret=self->getRenumArrForConsctvCellTypesSpe(order,order+sz);
        delete [] order;
        return ret;
      }

      PyObject *getCellsContainingPoints(PyObject *p, int nbOfPoints, double eps) const
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

      PyObject *getCellsContainingPoint(PyObject *p, double eps) const
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

      static PyObject *fuseUMeshesOnSameCoords(PyObject *ms, int compType)
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
    }
    void convertToPolyTypes(const std::vector<int>& cellIdsToConvert);
    MEDCouplingUMesh *buildExtrudedMeshFromThis(const MEDCouplingUMesh *mesh1D, int policy);
    static MEDCouplingUMesh *mergeUMeshes(const MEDCouplingUMesh *mesh1, const MEDCouplingUMesh *mesh2) throw(INTERP_KERNEL::Exception);
  };

  class MEDCouplingExtrudedMesh : public ParaMEDMEM::MEDCouplingMesh
  {
  public:
    static MEDCouplingExtrudedMesh *New(const MEDCouplingUMesh *mesh3D, MEDCouplingUMesh *mesh2D, int cell2DId) throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *build3DUnstructuredMesh() const;
    %extend {
      std::string __str__() const
      {
        return self->simpleRepr();
      }
      PyObject *getMesh2D() const
      {
        MEDCouplingUMesh *ret=self->getMesh2D();
        ret->incrRef();
        return convertMesh(ret, SWIG_POINTER_OWN | 0 );
      }
      PyObject *getMesh1D() const
      {
        MEDCouplingUMesh *ret=self->getMesh1D();
        ret->incrRef();
        return convertMesh(ret, SWIG_POINTER_OWN | 0 );
      }
      PyObject *getMesh3DIds() const
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
                   DataArrayDouble *coordsZ=0);
    void setCoordsAt(int i, DataArrayDouble *arr) throw(INTERP_KERNEL::Exception);
    %extend {
      std::string __str__() const
      {
        return self->simpleRepr();
      }
    }
  };
}

%extend ParaMEDMEM::DataArrayDouble
 {
   std::string __str__() const
   {
     return self->repr();
   }

   void setValues(PyObject *li, int nbOfTuples, int nbOfElsPerTuple)
   {
     int sz;
     double *tmp=convertPyToNewDblArr2(li,&sz);
     self->useArray(tmp,true,CPP_DEALLOC,nbOfTuples,nbOfElsPerTuple);
   }

   PyObject *getValues()
   {
     const double *vals=self->getPointer();
     return convertDblArrToPyList(vals,self->getNbOfElems());
   }

   PyObject *getValuesAsTuple()
   {
     const double *vals=self->getPointer();
     int nbOfComp=self->getNumberOfComponents();
     int nbOfTuples=self->getNumberOfTuples();
     return convertDblArrToPyListOfTuple(vals,nbOfComp,nbOfTuples);
   }
 };

%extend ParaMEDMEM::DataArrayInt
 {
   std::string __str__() const
   {
     return self->repr();
   }

   void setValues(PyObject *li, int nbOfTuples, int nbOfElsPerTuple)
   {
     int size;
     int *tmp=convertPyToNewIntArr2(li,&size);
     self->useArray(tmp,true,CPP_DEALLOC,nbOfTuples,nbOfElsPerTuple);
   }

   PyObject *getValues()
   {
     const int *vals=self->getPointer();
     return convertIntArrToPyList(vals,self->getNbOfElems());
   }

   PyObject *getValuesAsTuple()
   {
     const int *vals=self->getPointer();
     int nbOfComp=self->getNumberOfComponents();
     int nbOfTuples=self->getNumberOfTuples();
     return convertIntArrToPyListOfTuple(vals,nbOfComp,nbOfTuples);
   }

   static PyObject *makePartition(PyObject *gps, int newNb)
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
 };

namespace ParaMEDMEM
{
  class MEDCouplingField : public ParaMEDMEM::RefCountObject, public ParaMEDMEM::TimeLabel
  {
  public:
    virtual void checkCoherency() const throw(INTERP_KERNEL::Exception);
    virtual bool areCompatibleForMerge(const MEDCouplingField *other) const;
    virtual bool isEqual(const MEDCouplingField *other, double meshPrec, double valsPrec) const;
    virtual bool isEqualWithoutConsideringStr(const MEDCouplingField *other, double meshPrec, double valsPrec) const;
    void setMesh(const ParaMEDMEM::MEDCouplingMesh *mesh);
    void setName(const char *name);
    const char *getDescription() const;
    void setDescription(const char *desc);
    const char *getName() const;
    TypeOfField getTypeOfField() const;
    MEDCouplingFieldDouble *buildWeightingField(bool isAbs) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDiscretization *getDiscretization() const;
    void setGaussLocalizationOnType(INTERP_KERNEL::NormalizedCellType type, const std::vector<double>& refCoo,
                                    const std::vector<double>& gsCoo, const std::vector<double>& wg) throw(INTERP_KERNEL::Exception);
    void clearGaussLocalizations();
    MEDCouplingGaussLocalization& getGaussLocalization(int locId) throw(INTERP_KERNEL::Exception);
    int getNbOfGaussLocalization() const throw(INTERP_KERNEL::Exception);
    int getGaussLocalizationIdOfOneCell(int cellId) const throw(INTERP_KERNEL::Exception);
    const MEDCouplingGaussLocalization& getGaussLocalization(int locId) const throw(INTERP_KERNEL::Exception);
    %extend {
      PyObject *getMesh() const
      {
        MEDCouplingMesh *ret1=(MEDCouplingMesh *)self->getMesh();
        ret1->incrRef();
        return convertMesh(ret1, SWIG_POINTER_OWN | 0 );
      }

      PyObject *buildSubMeshData(PyObject *li) const
      {
        int size;
        int *tmp=convertPyToNewIntArr2(li,&size);
        DataArrayInt *ret1;
        MEDCouplingMesh *ret0=self->buildSubMeshData(tmp,tmp+size,ret1);
        delete [] tmp;
        PyObject *res = PyList_New(2);
        PyList_SetItem(res,0,convertMesh(ret0, SWIG_POINTER_OWN | 0 ));
        PyList_SetItem(res,1,SWIG_NewPointerObj((void*)ret1,SWIGTYPE_p_ParaMEDMEM__DataArrayInt,SWIG_POINTER_OWN | 0));
        return res;
      }
      void setGaussLocalizationOnCells(PyObject *li, const std::vector<double>& refCoo,
                                       const std::vector<double>& gsCoo, const std::vector<double>& wg) throw(INTERP_KERNEL::Exception)
      {
        int size;
        int *tmp=convertPyToNewIntArr2(li,&size);
        try
          {
            self->setGaussLocalizationOnCells(tmp,tmp+size,refCoo,gsCoo,wg);
          }
        catch(INTERP_KERNEL::Exception& e)
          {
            delete [] tmp;
            throw e;
          }
        delete [] tmp;
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
    std::string simpleRepr() const;
    std::string advancedRepr() const;
    MEDCouplingFieldDouble *clone(bool recDeepCpy) const;
    MEDCouplingFieldDouble *cloneWithMesh(bool recDeepCpy) const;
    MEDCouplingFieldDouble *buildNewTimeReprFromThis(TypeOfTimeDiscretization td, bool deepCpy) const;
    TypeOfTimeDiscretization getTimeDiscretization() const;
    void checkCoherency() const throw(INTERP_KERNEL::Exception);
    double getIJ(int tupleId, int compoId) const;
    double getIJK(int cellId, int nodeIdInCell, int compoId) const;
    void setArray(DataArrayDouble *array) throw(INTERP_KERNEL::Exception);
    void setEndArray(DataArrayDouble *array) throw(INTERP_KERNEL::Exception);
    void setTime(double val, int iteration, int order) throw(INTERP_KERNEL::Exception);
    void setStartTime(double val, int iteration, int order) throw(INTERP_KERNEL::Exception);
    void setEndTime(double val, int iteration, int order) throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getArray() const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getEndArray() const throw(INTERP_KERNEL::Exception);
    void applyLin(double a, double b, int compoId) throw(INTERP_KERNEL::Exception);
    int getNumberOfComponents() const throw(INTERP_KERNEL::Exception);
    int getNumberOfTuples() const throw(INTERP_KERNEL::Exception);
    int getNumberOfValues() const throw(INTERP_KERNEL::Exception);
    NatureOfField getNature() const { return _nature; }
    void setNature(NatureOfField nat) throw(INTERP_KERNEL::Exception);
    void updateTime();
    void changeUnderlyingMesh(const MEDCouplingMesh *other, int levOfCheck, double prec) throw(INTERP_KERNEL::Exception);
    void substractInPlaceDM(const MEDCouplingFieldDouble *f, int levOfCheck, double prec) throw(INTERP_KERNEL::Exception);
    bool mergeNodes(double eps) throw(INTERP_KERNEL::Exception);
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
    void applyFunc(int nbOfComp, const char *func) throw(INTERP_KERNEL::Exception);
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
    static MEDCouplingFieldDouble *dotFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    MEDCouplingFieldDouble *dot(const MEDCouplingFieldDouble& other) const;
    static MEDCouplingFieldDouble *crossProductFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    MEDCouplingFieldDouble *crossProduct(const MEDCouplingFieldDouble& other) const;
    static MEDCouplingFieldDouble *maxFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    MEDCouplingFieldDouble *max(const MEDCouplingFieldDouble& other) const;
    static MEDCouplingFieldDouble *minFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    MEDCouplingFieldDouble *min(const MEDCouplingFieldDouble& other) const;
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
      void setValues(PyObject *li)
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

      PyObject *getStartTime()
      {
        int tmp1,tmp2;
        double tmp0=self->getStartTime(tmp1,tmp2);
        PyObject *res = PyList_New(3);
        PyList_SetItem(res,0,SWIG_From_double(tmp0));
        PyList_SetItem(res,1,SWIG_From_int(tmp1));
        PyList_SetItem(res,2,SWIG_From_int(tmp2));
        return res;
      }

      PyObject *getEndTime()
      {
        int tmp1,tmp2;
        double tmp0=self->getEndTime(tmp1,tmp2);
        PyObject *res = PyList_New(3);
        PyList_SetItem(res,0,SWIG_From_double(tmp0));
        PyList_SetItem(res,1,SWIG_From_int(tmp1));
        PyList_SetItem(res,2,SWIG_From_int(tmp2));
        return res;
      }
      PyObject *accumulate() const
      {
        int sz=self->getNumberOfTuples();
        double *tmp=new double[sz];
        self->accumulate(tmp);
        PyObject *ret=convertDblArrToPyList(tmp,sz);
        delete [] tmp;
        return ret;
      }
      PyObject *integral(bool isWAbs) const
      {
        int sz=self->getNumberOfTuples();
        double *tmp=new double[sz];
        self->integral(isWAbs,tmp);
        PyObject *ret=convertDblArrToPyList(tmp,sz);
        delete [] tmp;
        return ret;
      }
      PyObject *normL1() const throw(INTERP_KERNEL::Exception)
      {
        int sz=self->getNumberOfTuples();
        double *tmp=new double[sz];
        self->normL1(tmp);
        PyObject *ret=convertDblArrToPyList(tmp,sz);
        delete [] tmp;
        return ret;
      }
      PyObject *normL2() const throw(INTERP_KERNEL::Exception)
      {
        int sz=self->getNumberOfTuples();
        double *tmp=new double[sz];
        self->normL2(tmp);
        PyObject *ret=convertDblArrToPyList(tmp,sz);
        delete [] tmp;
        return ret;
      }

      void renumberCells(PyObject *li, bool check) throw(INTERP_KERNEL::Exception)
      {
        int size;
        int *tmp=convertPyToNewIntArr2(li,&size);
        try
          {
            self->renumberCells(tmp,check);
          }
        catch(INTERP_KERNEL::Exception& e)
          {
            delete [] tmp;
            throw e;
          }
        delete [] tmp;
      }
      void renumberNodes(PyObject *li) throw(INTERP_KERNEL::Exception)
      {
        int size;
        int *tmp=convertPyToNewIntArr2(li,&size);
        try
          {
            self->renumberNodes(tmp);
          }
        catch(INTERP_KERNEL::Exception& e)
          {
            delete [] tmp;
            throw e;
          }
        delete [] tmp;
      }

      MEDCouplingFieldDouble *buildSubPart(PyObject *li) const throw(INTERP_KERNEL::Exception)
      {
        int size;
        int *tmp=convertPyToNewIntArr2(li,&size);
        MEDCouplingFieldDouble *ret=0;
        try
          {
            ret=self->buildSubPart(tmp,tmp+size);
          }
        catch(INTERP_KERNEL::Exception& e)
          {
            delete [] tmp;
            throw e;
          }
        delete [] tmp;
        return ret;
      }
    }
  };
}
