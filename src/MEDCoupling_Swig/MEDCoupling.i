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
#include "MEDCouplingFieldTemplate.hxx"
#include "MEDCouplingGaussLocalization.hxx"
#include "MEDCouplingAutoRefCountObjectPtr.hxx"
#include "MEDCouplingMultiFields.hxx"
#include "MEDCouplingFieldOverTime.hxx"
#include "MEDCouplingDefinitionTime.hxx"
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

%typemap(out) ParaMEDMEM::MEDCouplingMultiFields*
{
  $result=convertMultiFields($1,$owner);
}

#ifdef WITH_NUMPY2
%init %{ import_array(); %}
#endif

%feature("autodoc", "1");
%feature("docstring");

%newobject ParaMEDMEM::MEDCouplingFieldDiscretization::getOffsetArr;
%newobject ParaMEDMEM::MEDCouplingField::buildMeasureField;
%newobject ParaMEDMEM::MEDCouplingField::getLocalizationOfDiscr;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::New;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::getArray;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::getEndArray;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::MergeFields;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::MeldFields;
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
%newobject ParaMEDMEM::MEDCouplingFieldDouble::DotFields;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::dot;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::CrossProductFields;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::crossProduct;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::MaxFields;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::max;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::MinFields;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::AddFields;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::SubstractFields;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::MultiplyFields;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::DivideFields;
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
%newobject ParaMEDMEM::MEDCouplingFieldDouble::getValueOnMulti;
%newobject ParaMEDMEM::MEDCouplingFieldTemplate::New;
%newobject ParaMEDMEM::DataArrayInt::New;
%newobject ParaMEDMEM::DataArrayInt::convertToDblArr;
%newobject ParaMEDMEM::DataArrayInt::deepCpy;
%newobject ParaMEDMEM::DataArrayInt::performCpy;
%newobject ParaMEDMEM::DataArrayInt::substr;
%newobject ParaMEDMEM::DataArrayInt::changeNbOfComponents;
%newobject ParaMEDMEM::DataArrayInt::keepSelectedComponents;
%newobject ParaMEDMEM::DataArrayInt::selectByTupleId;
%newobject ParaMEDMEM::DataArrayInt::selectByTupleIdSafe;
%newobject ParaMEDMEM::DataArrayInt::selectByTupleId2;
%newobject ParaMEDMEM::DataArrayInt::checkAndPreparePermutation;
%newobject ParaMEDMEM::DataArrayInt::transformWithIndArrR;
%newobject ParaMEDMEM::DataArrayInt::renumber;
%newobject ParaMEDMEM::DataArrayInt::renumberR;
%newobject ParaMEDMEM::DataArrayInt::renumberAndReduce;
%newobject ParaMEDMEM::DataArrayInt::invertArrayO2N2N2O;
%newobject ParaMEDMEM::DataArrayInt::invertArrayN2O2O2N;
%newobject ParaMEDMEM::DataArrayInt::getIdsEqual;
%newobject ParaMEDMEM::DataArrayInt::getIdsNotEqual;
%newobject ParaMEDMEM::DataArrayInt::getIdsEqualList;
%newobject ParaMEDMEM::DataArrayInt::getIdsNotEqualList;
%newobject ParaMEDMEM::DataArrayInt::Aggregate;
%newobject ParaMEDMEM::DataArrayInt::Meld;
%newobject ParaMEDMEM::DataArrayInt::Add;
%newobject ParaMEDMEM::DataArrayInt::Substract;
%newobject ParaMEDMEM::DataArrayInt::Multiply;
%newobject ParaMEDMEM::DataArrayInt::Divide;
%newobject ParaMEDMEM::DataArrayInt::BuildUnion;
%newobject ParaMEDMEM::DataArrayInt::BuildIntersection;
%newobject ParaMEDMEM::DataArrayInt::fromNoInterlace;
%newobject ParaMEDMEM::DataArrayInt::toNoInterlace;
%newobject ParaMEDMEM::DataArrayInt::buildComplement;
%newobject ParaMEDMEM::DataArrayInt::buildUnion;
%newobject ParaMEDMEM::DataArrayInt::buildSubstraction;
%newobject ParaMEDMEM::DataArrayInt::buildIntersection;
%newobject ParaMEDMEM::DataArrayInt::deltaShiftIndex;
%newobject ParaMEDMEM::DataArrayInt::buildExplicitArrByRanges;
%newobject ParaMEDMEM::DataArrayInt::buildPermutationArr;
%newobject ParaMEDMEM::DataArrayInt::buildPermArrPerLevel;
%newobject ParaMEDMEM::DataArrayInt::__getitem__;
%newobject ParaMEDMEM::DataArrayInt::__add__;
%newobject ParaMEDMEM::DataArrayInt::__radd__;
%newobject ParaMEDMEM::DataArrayInt::__sub__;
%newobject ParaMEDMEM::DataArrayInt::__rsub__;
%newobject ParaMEDMEM::DataArrayInt::__mul__;
%newobject ParaMEDMEM::DataArrayInt::__rmul__;
%newobject ParaMEDMEM::DataArrayInt::__div__;
%newobject ParaMEDMEM::DataArrayInt::__rdiv__;
%newobject ParaMEDMEM::DataArrayInt::__mod__;
%newobject ParaMEDMEM::DataArrayInt::__rmod__;
%newobject ParaMEDMEM::DataArrayDouble::New;
%newobject ParaMEDMEM::DataArrayDouble::convertToIntArr;
%newobject ParaMEDMEM::DataArrayDouble::deepCpy;
%newobject ParaMEDMEM::DataArrayDouble::performCpy;
%newobject ParaMEDMEM::DataArrayDouble::Aggregate;
%newobject ParaMEDMEM::DataArrayDouble::Meld;
%newobject ParaMEDMEM::DataArrayDouble::Dot;
%newobject ParaMEDMEM::DataArrayDouble::CrossProduct;
%newobject ParaMEDMEM::DataArrayDouble::Add;
%newobject ParaMEDMEM::DataArrayDouble::Substract;
%newobject ParaMEDMEM::DataArrayDouble::Multiply;
%newobject ParaMEDMEM::DataArrayDouble::Divide;
%newobject ParaMEDMEM::DataArrayDouble::substr;
%newobject ParaMEDMEM::DataArrayDouble::changeNbOfComponents;
%newobject ParaMEDMEM::DataArrayDouble::keepSelectedComponents;
%newobject ParaMEDMEM::DataArrayDouble::getIdsInRange;
%newobject ParaMEDMEM::DataArrayDouble::selectByTupleId;
%newobject ParaMEDMEM::DataArrayDouble::selectByTupleIdSafe;
%newobject ParaMEDMEM::DataArrayDouble::selectByTupleId2;
%newobject ParaMEDMEM::DataArrayDouble::applyFunc;
%newobject ParaMEDMEM::DataArrayDouble::applyFunc2;
%newobject ParaMEDMEM::DataArrayDouble::applyFunc3;
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
%newobject ParaMEDMEM::DataArrayDouble::__getitem__;
%newobject ParaMEDMEM::DataArrayDouble::__add__;
%newobject ParaMEDMEM::DataArrayDouble::__radd__;
%newobject ParaMEDMEM::DataArrayDouble::__sub__;
%newobject ParaMEDMEM::DataArrayDouble::__rsub__;
%newobject ParaMEDMEM::DataArrayDouble::__mul__;
%newobject ParaMEDMEM::DataArrayDouble::__rmul__;
%newobject ParaMEDMEM::DataArrayDouble::__div__;
%newobject ParaMEDMEM::DataArrayDouble::__rdiv__;
%newobject ParaMEDMEM::MEDCouplingMesh::deepCpy;
%newobject ParaMEDMEM::MEDCouplingMesh::checkTypeConsistencyAndContig;
%newobject ParaMEDMEM::MEDCouplingMesh::getCoordinatesAndOwner;
%newobject ParaMEDMEM::MEDCouplingMesh::getBarycenterAndOwner;
%newobject ParaMEDMEM::MEDCouplingMesh::buildOrthogonalField;
%newobject ParaMEDMEM::MEDCouplingMesh::getCellIdsFullyIncludedInNodeIds;
%newobject ParaMEDMEM::MEDCouplingMesh::mergeMyselfWith;
%newobject ParaMEDMEM::MEDCouplingMesh::fillFromAnalytic;
%newobject ParaMEDMEM::MEDCouplingMesh::fillFromAnalytic2;
%newobject ParaMEDMEM::MEDCouplingMesh::fillFromAnalytic3;
%newobject ParaMEDMEM::MEDCouplingMesh::getMeasureField;
%newobject ParaMEDMEM::MEDCouplingMesh::simplexize;
%newobject ParaMEDMEM::MEDCouplingMesh::buildUnstructured;
%newobject ParaMEDMEM::MEDCouplingMesh::MergeMeshes;
%newobject ParaMEDMEM::MEDCouplingPointSet::zipCoordsTraducer;
%newobject ParaMEDMEM::MEDCouplingPointSet::buildBoundaryMesh;
%newobject ParaMEDMEM::MEDCouplingPointSet::MergeNodesArray;
%newobject ParaMEDMEM::MEDCouplingPointSet::BuildInstanceFromMeshType;
%newobject ParaMEDMEM::MEDCouplingUMesh::New;
%newobject ParaMEDMEM::MEDCouplingUMesh::getNodalConnectivity;
%newobject ParaMEDMEM::MEDCouplingUMesh::getNodalConnectivityIndex;
%newobject ParaMEDMEM::MEDCouplingUMesh::clone;
%newobject ParaMEDMEM::MEDCouplingUMesh::zipConnectivityTraducer;
%newobject ParaMEDMEM::MEDCouplingUMesh::buildDescendingConnectivity;
%newobject ParaMEDMEM::MEDCouplingUMesh::buildExtrudedMesh;
%newobject ParaMEDMEM::MEDCouplingUMesh::MergeUMeshes;
%newobject ParaMEDMEM::MEDCouplingUMesh::MergeUMeshesOnSameCoords;
%newobject ParaMEDMEM::MEDCouplingUMesh::buildNewNumberingFromCommNodesFrmt;
%newobject ParaMEDMEM::MEDCouplingUMesh::rearrange2ConsecutiveCellTypes;
%newobject ParaMEDMEM::MEDCouplingUMesh::sortCellsInMEDFileFrmt;
%newobject ParaMEDMEM::MEDCouplingUMesh::convertCellArrayPerGeoType;
%newobject ParaMEDMEM::MEDCouplingUMesh::getRenumArrForConsecutiveCellTypesSpec;
%newobject ParaMEDMEM::MEDCouplingUMesh::buildDirectionVectorField;
%newobject ParaMEDMEM::MEDCouplingUMesh::getEdgeRatioField;
%newobject ParaMEDMEM::MEDCouplingUMesh::getAspectRatioField;
%newobject ParaMEDMEM::MEDCouplingUMesh::getWarpField;
%newobject ParaMEDMEM::MEDCouplingUMesh::getSkewField;
%newobject ParaMEDMEM::MEDCouplingUMesh::getPartBarycenterAndOwner;
%newobject ParaMEDMEM::MEDCouplingUMesh::getPartMeasureField;
%newobject ParaMEDMEM::MEDCouplingUMesh::buildPartOrthogonalField;
%newobject ParaMEDMEM::MEDCouplingUMesh::keepCellIdsByType;
%newobject ParaMEDMEM::MEDCouplingUMesh::Build0DMeshFromCoords;
%newobject ParaMEDMEM::MEDCouplingUMesh::findCellsIdsOnBoundary;
%newobject ParaMEDMEM::MEDCouplingUMesh::getCellIdsLyingOnNodes;
%newobject ParaMEDMEM::MEDCouplingExtrudedMesh::New;
%newobject ParaMEDMEM::MEDCouplingExtrudedMesh::build3DUnstructuredMesh;
%newobject ParaMEDMEM::MEDCouplingCMesh::New;
%newobject ParaMEDMEM::MEDCouplingCMesh::getCoordsAt;
%newobject ParaMEDMEM::MEDCouplingMultiFields::New;
%newobject ParaMEDMEM::MEDCouplingMultiFields::deepCpy;
%newobject ParaMEDMEM::MEDCouplingFieldOverTime::New;

%feature("unref") DataArrayDouble "$this->decrRef();"
%feature("unref") MEDCouplingPointSet "$this->decrRef();"
%feature("unref") MEDCouplingMesh "$this->decrRef();"
%feature("unref") MEDCouplingUMesh "$this->decrRef();"
%feature("unref") MEDCouplingExtrudedMesh "$this->decrRef();"
%feature("unref") MEDCouplingCMesh "$this->decrRef();"
%feature("unref") DataArrayInt "$this->decrRef();"
%feature("unref") MEDCouplingField "$this->decrRef();"
%feature("unref") MEDCouplingFieldDouble "$this->decrRef();"
%feature("unref") MEDCouplingMultiFields "$this->decrRef();"
%feature("unref") MEDCouplingFieldTemplate "$this->decrRef();"
%feature("unref") MEDCouplingMultiFields "$this->decrRef();"

%rename(assign) *::operator=;
%ignore ParaMEDMEM::RefCountObject::decrRef;
%ignore ParaMEDMEM::MemArray::operator=;
%ignore ParaMEDMEM::MemArray::operator[];
%ignore ParaMEDMEM::MEDCouplingGaussLocalization::pushTinySerializationIntInfo;
%ignore ParaMEDMEM::MEDCouplingGaussLocalization::pushTinySerializationDblInfo;
%ignore ParaMEDMEM::MEDCouplingGaussLocalization::fillWithValues;
%ignore ParaMEDMEM::MEDCouplingGaussLocalization::buildNewInstanceFromTinyInfo;

%nodefaultctor;

%rename (InterpKernelException) INTERP_KERNEL::Exception;

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

  %extend RefCountObject
  {
    std::string getHiddenCppPointer() const
    {
      std::ostringstream oss; oss << "C++ Pointer address is : " << self;
      return oss.str();
    }
  }

  class MEDCouplingMesh : public RefCountObject, public TimeLabel
  {
  public:
    void setName(const char *name);
    const char *getName() const;
    void setDescription(const char *descr);
    const char *getDescription() const;
    void setTime(double val, int iteration, int order);
    void setTimeUnit(const char *unit);
    const char *getTimeUnit() const;
    virtual MEDCouplingMeshType getType() const throw(INTERP_KERNEL::Exception) = 0;
    bool isStructured() const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingMesh *deepCpy() const = 0;
    virtual bool isEqual(const MEDCouplingMesh *other, double prec) const throw(INTERP_KERNEL::Exception);
    virtual bool isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const throw(INTERP_KERNEL::Exception) = 0;
    virtual void copyTinyStringsFrom(const MEDCouplingMesh *other) throw(INTERP_KERNEL::Exception);
    virtual void copyTinyInfoFrom(const MEDCouplingMesh *other) throw(INTERP_KERNEL::Exception);
    virtual void checkCoherency() const throw(INTERP_KERNEL::Exception) = 0;
    virtual void checkCoherency1(double eps=1e-12) const throw(INTERP_KERNEL::Exception) = 0;
    virtual void checkCoherency2(double eps=1e-12) const throw(INTERP_KERNEL::Exception) = 0;
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
    virtual MEDCouplingFieldDouble *fillFromAnalytic2(TypeOfField t, int nbOfComp, const char *func) const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingFieldDouble *fillFromAnalytic3(TypeOfField t, int nbOfComp, const std::vector<std::string>& varsOrder, const char *func) const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingFieldDouble *buildOrthogonalField() const throw(INTERP_KERNEL::Exception) = 0;
    virtual MEDCouplingUMesh *buildUnstructured() const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingMesh *mergeMyselfWith(const MEDCouplingMesh *other) const throw(INTERP_KERNEL::Exception) = 0;
    virtual bool areCompatibleForMerge(const MEDCouplingMesh *other) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *simplexize(int policy) throw(INTERP_KERNEL::Exception);
    static MEDCouplingMesh *MergeMeshes(const MEDCouplingMesh *mesh1, const MEDCouplingMesh *mesh2) throw(INTERP_KERNEL::Exception);
    %extend
       {
         std::string __str__() const
         {
           return self->simpleRepr();
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

         int getCellContainingPoint(PyObject *p, double eps) const throw(INTERP_KERNEL::Exception)
         {
           int sz;
           double *pos=convertPyToNewDblArr2(p,&sz);
           int ret=self->getCellContainingPoint(pos,eps);
           delete [] pos;
           return ret;
         }

         PyObject *getCellsContainingPoints(PyObject *p, int nbOfPoints, double eps) const throw(INTERP_KERNEL::Exception)
         {
           int sz;
           INTERP_KERNEL::AutoPtr<double> pos=convertPyToNewDblArr2(p,&sz);
           std::vector<int> elts,eltsIndex;
           self->getCellsContainingPoints(pos,nbOfPoints,eps,elts,eltsIndex);
           MEDCouplingAutoRefCountObjectPtr<DataArrayInt> d0=DataArrayInt::New();
           MEDCouplingAutoRefCountObjectPtr<DataArrayInt> d1=DataArrayInt::New();
           d0->alloc(elts.size(),1);
           d1->alloc(eltsIndex.size(),1);
           std::copy(elts.begin(),elts.end(),d0->getPointer());
           std::copy(eltsIndex.begin(),eltsIndex.end(),d1->getPointer());
           PyObject *ret=PyTuple_New(2);
           PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(d0),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
           PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(d1),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
           d0->incrRef();
           d1->incrRef();
           return ret;
         }

         PyObject *getCellsContainingPoint(PyObject *p, double eps) const throw(INTERP_KERNEL::Exception)
         {
           int sz;
           INTERP_KERNEL::AutoPtr<double> pos=convertPyToNewDblArr2(p,&sz);
           std::vector<int> elts;
           self->getCellsContainingPoint(pos,eps,elts);
           return convertIntArrToPyList2(elts);
         }
         
         void renumberCells(PyObject *li, bool check) throw(INTERP_KERNEL::Exception)
         {
           void *da=0;
           int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 |  0 );
           if (!SWIG_IsOK(res1))
             {
               int size;
               INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
               self->renumberCells(tmp,check);
             }
           else
             {
               DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
               if(!da2)
                   throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
               da2->checkAllocated();
               self->renumberCells(da2->getConstPointer(),check);
             }
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
           void *da=0;
           int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 |  0 );
           if (!SWIG_IsOK(res1))
             {
               int size;
               INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
               return self->getCellIdsFullyIncludedInNodeIds(tmp,((const int *)tmp)+size);
             }
           else
             {
               DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
               if(!da2)
                 throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
               da2->checkAllocated();
               return self->getCellIdsFullyIncludedInNodeIds(da2->getConstPointer(),da2->getConstPointer()+da2->getNbOfElems());
             }
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
           void *da=0;
           int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 |  0 );
           if (!SWIG_IsOK(res1))
             {
               int size;
               INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
               MEDCouplingMesh *ret=self->buildPart(tmp,((const int *)tmp)+size);
               return convertMesh(ret, SWIG_POINTER_OWN | 0 );
             }
           else
             {
               DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
               if(!da2)
                 throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
               da2->checkAllocated();
               MEDCouplingMesh *ret=self->buildPart(da2->getConstPointer(),da2->getConstPointer()+da2->getNbOfElems());
               ret->setName(da2->getName().c_str());
               return convertMesh(ret, SWIG_POINTER_OWN | 0 );
             }
         }
        
        PyObject *buildPartAndReduceNodes(PyObject *li) const throw(INTERP_KERNEL::Exception)
        {
          void *da=0;
          DataArrayInt *arr=0;
          MEDCouplingMesh *ret=0;
           int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 |  0 );
           if (!SWIG_IsOK(res1))
             {
               int size;
               INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
               ret=self->buildPartAndReduceNodes(tmp,((const int *)tmp)+size,arr);
             }
           else
             {
               DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
               if(!da2)
                 throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
               da2->checkAllocated();
               ret=self->buildPartAndReduceNodes(da2->getConstPointer(),da2->getConstPointer()+da2->getNbOfElems(),arr);
               ret->setName(da2->getName().c_str());
             }
          PyObject *res = PyList_New(2);
          PyObject *obj0=convertMesh(ret, SWIG_POINTER_OWN | 0 );
          PyObject *obj1=SWIG_NewPointerObj(SWIG_as_voidptr(arr),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 );
          PyList_SetItem(res,0,obj0);
          PyList_SetItem(res,1,obj1);
          return res;
        }

        DataArrayInt *checkTypeConsistencyAndContig(PyObject *li, PyObject *li2) const throw(INTERP_KERNEL::Exception)
        {
          std::vector<int> code;
          std::vector<const DataArrayInt *> idsPerType;
          convertPyObjToVecDataArrayIntCst(li2,idsPerType);
          convertPyToNewIntArr3(li,code);
          return self->checkTypeConsistencyAndContig(code,idsPerType);
        }

        PyObject *splitProfilePerType(const DataArrayInt *profile) const throw(INTERP_KERNEL::Exception)
        {
          std::vector<int> code;
          std::vector<DataArrayInt *> idsInPflPerType;
          std::vector<DataArrayInt *> idsPerType;
          self->splitProfilePerType(profile,code,idsInPflPerType,idsPerType);
          PyObject *ret=PyTuple_New(3);
          PyTuple_SetItem(ret,0,convertIntArrToPyList2(code));
          PyObject *ret1=PyList_New(idsInPflPerType.size());
          for(std::size_t j=0;j<idsInPflPerType.size();j++)
            PyList_SetItem(ret1,j,SWIG_NewPointerObj(SWIG_as_voidptr(idsInPflPerType[j]),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
          PyTuple_SetItem(ret,1,ret1);
          int n=idsPerType.size();
          PyObject *ret2=PyList_New(n);
          for(int i=0;i<n;i++)
            PyList_SetItem(ret2,i,SWIG_NewPointerObj(SWIG_as_voidptr(idsPerType[i]),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
          PyTuple_SetItem(ret,2,ret2);
          return ret;
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

         PyObject *getAllGeoTypes() const throw(INTERP_KERNEL::Exception)
         {
           std::set<INTERP_KERNEL::NormalizedCellType> result=self->getAllGeoTypes();
           std::set<INTERP_KERNEL::NormalizedCellType>::const_iterator iL=result.begin();
           PyObject *res=PyList_New(result.size());
           for(int i=0;iL!=result.end(); i++, iL++)
             PyList_SetItem(res,i,PyInt_FromLong(*iL));
           return res;
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
      void updateTime() const;
      void setCoords(const DataArrayDouble *coords) throw(INTERP_KERNEL::Exception);
      DataArrayDouble *getCoordinatesAndOwner() const throw(INTERP_KERNEL::Exception);
      bool areCoordsEqual(const MEDCouplingPointSet& other, double prec) const throw(INTERP_KERNEL::Exception);
      void zipCoords() throw(INTERP_KERNEL::Exception);
      double getCaracteristicDimension() const throw(INTERP_KERNEL::Exception);
      void changeSpaceDimension(int newSpaceDim, double dftVal=0.) throw(INTERP_KERNEL::Exception);
      void tryToShareSameCoords(const MEDCouplingPointSet& other, double epsilon) throw(INTERP_KERNEL::Exception);
      virtual void tryToShareSameCoordsPermute(const MEDCouplingPointSet& other, double epsilon) throw(INTERP_KERNEL::Exception) = 0;
      static DataArrayDouble *MergeNodesArray(const MEDCouplingPointSet *m1, const MEDCouplingPointSet *m2) throw(INTERP_KERNEL::Exception);
      static MEDCouplingPointSet *BuildInstanceFromMeshType(MEDCouplingMeshType type) throw(INTERP_KERNEL::Exception);
      virtual MEDCouplingPointSet *buildBoundaryMesh(bool keepCoords) const throw(INTERP_KERNEL::Exception) = 0;
      virtual bool isEmptyMesh(const std::vector<int>& tinyInfo) const throw(INTERP_KERNEL::Exception) = 0;
      //! size of returned tinyInfo must be always the same.
      void getTinySerializationInformation(std::vector<double>& tinyInfoD, std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const throw(INTERP_KERNEL::Exception);
      void resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings) const throw(INTERP_KERNEL::Exception);
      void serialize(DataArrayInt *&a1, DataArrayDouble *&a2) const throw(INTERP_KERNEL::Exception);
      void unserialization(const std::vector<double>& tinyInfoD, const std::vector<int>& tinyInfo, const DataArrayInt *a1, DataArrayDouble *a2,
                           const std::vector<std::string>& littleStrings) throw(INTERP_KERNEL::Exception);
      virtual void getCellsInBoundingBox(const INTERP_KERNEL::DirectedBoundingBox& bbox, double eps, std::vector<int>& elems) throw(INTERP_KERNEL::Exception) = 0;
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
           
           PyObject *getCoords() throw(INTERP_KERNEL::Exception)
           {
             DataArrayDouble *ret1=self->getCoords();
             ret1->incrRef();
             return SWIG_NewPointerObj((void*)ret1,SWIGTYPE_p_ParaMEDMEM__DataArrayDouble,SWIG_POINTER_OWN | 0);
           }
           PyObject *buildPartOfMySelf(PyObject *li, bool keepCoords) const throw(INTERP_KERNEL::Exception)
           {
             void *da=0;
             int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 |  0 );
             if (!SWIG_IsOK(res1))
               {
                 int size;
                 INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
                 MEDCouplingPointSet *ret=self->buildPartOfMySelf(tmp,((const int *)tmp)+size,keepCoords);
                 return convertMesh(ret, SWIG_POINTER_OWN | 0 );
               }
             else
               {
                 DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
                 if(!da2)
                   throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
                 da2->checkAllocated();
                 MEDCouplingPointSet *ret=self->buildPartOfMySelf(da2->getConstPointer(),da2->getConstPointer()+da2->getNbOfElems(),keepCoords);
                 ret->setName(da2->getName().c_str());
                 return convertMesh(ret, SWIG_POINTER_OWN | 0 );
               }
           }
           PyObject *buildPartOfMySelfNode(PyObject *li, bool fullyIn) const throw(INTERP_KERNEL::Exception)
           {
             void *da=0;
             int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 |  0 );
             if (!SWIG_IsOK(res1))
               {
                 int size;
                 INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
                 MEDCouplingPointSet *ret=self->buildPartOfMySelfNode(tmp,((const int *)tmp)+size,fullyIn);
                 return convertMesh(ret, SWIG_POINTER_OWN | 0 );
               }
             else
               {
                 DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
                 if(!da2)
                   throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
                 da2->checkAllocated();
                 MEDCouplingPointSet *ret=self->buildPartOfMySelfNode(da2->getConstPointer(),da2->getConstPointer()+da2->getNbOfElems(),fullyIn);
                 ret->setName(da2->getName().c_str());
                 return convertMesh(ret, SWIG_POINTER_OWN | 0 );
               }
           }
           PyObject *buildFacePartOfMySelfNode(PyObject *li, bool fullyIn) const throw(INTERP_KERNEL::Exception)
           {
             void *da=0;
             int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 |  0 );
             if (!SWIG_IsOK(res1))
               {
                 int size;
                 INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
                 MEDCouplingPointSet *ret=self->buildFacePartOfMySelfNode(tmp,((const int *)tmp)+size,fullyIn);
                 return convertMesh(ret, SWIG_POINTER_OWN | 0 );
               }
             else
               {
                 DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
                 if(!da2)
                   throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
                 da2->checkAllocated();
                 MEDCouplingPointSet *ret=self->buildFacePartOfMySelfNode(da2->getConstPointer(),da2->getConstPointer()+da2->getNbOfElems(),fullyIn);
                 ret->setName(da2->getName().c_str());
                 return convertMesh(ret, SWIG_POINTER_OWN | 0 );
               }
           }
           PyObject *findBoundaryNodes() const throw(INTERP_KERNEL::Exception)
           {
             std::vector<int> nodes;
             self->findBoundaryNodes(nodes);
             return convertIntArrToPyList2(nodes);
           }
           void renumberNodes(PyObject *li, int newNbOfNodes) throw(INTERP_KERNEL::Exception)
           {
             void *da=0;
             int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 |  0 );
             if (!SWIG_IsOK(res1))
               {
                 int size;
                 INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
                 self->renumberNodes(tmp,newNbOfNodes);
               }
             else
               {
                 DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
                 if(!da2)
                   throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
                 da2->checkAllocated();
                 self->renumberNodes(da2->getConstPointer(),newNbOfNodes);
               }
           }
           void renumberNodes2(PyObject *li, int newNbOfNodes) throw(INTERP_KERNEL::Exception)
           {
             void *da=0;
             int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 |  0 );
             if (!SWIG_IsOK(res1))
               {
                 int size;
                 INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
                 self->renumberNodes2(tmp,newNbOfNodes);
               }
             else
               {
                 DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
                 if(!da2)
                   throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
                 da2->checkAllocated();
                 self->renumberNodes2(da2->getConstPointer(),newNbOfNodes);
               }
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
             INTERP_KERNEL::AutoPtr<double> pos=convertPyToNewDblArr2(pt,&size);
             if(size<self->getSpaceDimension())
               throw INTERP_KERNEL::Exception("getNodeIdsNearPoint : to tiny array ! must be at least of size SpaceDim !");
             std::vector<int> tmp=self->getNodeIdsNearPoint(pos,eps);
             return convertIntArrToPyList2(tmp);
           }

           PyObject *getNodeIdsNearPoints(PyObject *pt, int nbOfNodes, double eps) const throw(INTERP_KERNEL::Exception)
           {
             std::vector<int> c,cI;
             int size;
             INTERP_KERNEL::AutoPtr<double> pos=convertPyToNewDblArr2(pt,&size);
             if(size<self->getSpaceDimension()*nbOfNodes)
               throw INTERP_KERNEL::Exception("getNodeIdsNearPoints : to tiny array ! must be at least of size SpaceDim*nbOfNodes !");
             self->getNodeIdsNearPoints(pos,nbOfNodes,eps,c,cI);
             PyObject *ret=PyTuple_New(2);
             MEDCouplingAutoRefCountObjectPtr<DataArrayInt> d0=DataArrayInt::New();
             MEDCouplingAutoRefCountObjectPtr<DataArrayInt> d1=DataArrayInt::New();
             d0->alloc(c.size(),1);
             d1->alloc(cI.size(),1);
             std::copy(c.begin(),c.end(),d0->getPointer());
             std::copy(cI.begin(),cI.end(),d1->getPointer());
             PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(d0),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
             PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(d1),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
             d0->incrRef();
             d1->incrRef();
             return ret;
           }

           PyObject *getCellsInBoundingBox(PyObject *bbox, double eps) throw(INTERP_KERNEL::Exception)
           {
             std::vector<int> elems;
             int size;
             double *tmp=convertPyToNewDblArr2(bbox,&size);
             self->getCellsInBoundingBox(tmp,eps,elems);
             delete [] tmp;
             return convertIntArrToPyList2(elems);
           }

           static void Rotate2DAlg(PyObject *center, double angle, int nbNodes, PyObject *coords) throw(INTERP_KERNEL::Exception)
           {
             int sz;
             double *c=convertPyToNewDblArr2(center,&sz);
             double *coo=convertPyToNewDblArr2(coords,&sz);
             ParaMEDMEM::MEDCouplingPointSet::Rotate2DAlg(c,angle,nbNodes,coo);
             for(int i=0;i<sz;i++)
               PyList_SetItem(coords,i,PyFloat_FromDouble(coo[i]));
             delete [] coo;
             delete [] c;
           }
           static void Rotate3DAlg(PyObject *center, PyObject *vect, double angle, int nbNodes, PyObject *coords) throw(INTERP_KERNEL::Exception)
           {
             int sz,sz2;
             double *c=convertPyToNewDblArr2(center,&sz);
             double *coo=convertPyToNewDblArr2(coords,&sz);
             double *v=convertPyToNewDblArr2(vect,&sz2);
             ParaMEDMEM::MEDCouplingPointSet::Rotate3DAlg(c,v,angle,nbNodes,coo);
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
    void updateTime() const;
    void checkCoherency() const throw(INTERP_KERNEL::Exception);
    void setMeshDimension(int meshDim) throw(INTERP_KERNEL::Exception);
    void allocateCells(int nbOfCells) throw(INTERP_KERNEL::Exception);
    void finishInsertingCells() throw(INTERP_KERNEL::Exception);
    void setConnectivity(DataArrayInt *conn, DataArrayInt *connIndex, bool isComputingTypes=true) throw(INTERP_KERNEL::Exception);
    INTERP_KERNEL::NormalizedCellType getTypeOfCell(int cellId) const throw(INTERP_KERNEL::Exception);
    int getNumberOfNodesInCell(int cellId) const throw(INTERP_KERNEL::Exception);
    int getMeshLength() const throw(INTERP_KERNEL::Exception);
    void computeTypes() throw(INTERP_KERNEL::Exception);
    std::string reprConnectivityOfThis() const throw(INTERP_KERNEL::Exception);
    //tools
    DataArrayInt *findCellsIdsOnBoundary() const throw(INTERP_KERNEL::Exception);
    bool checkConsecutiveCellTypes() const throw(INTERP_KERNEL::Exception);
    DataArrayInt *rearrange2ConsecutiveCellTypes() throw(INTERP_KERNEL::Exception);
    DataArrayInt *sortCellsInMEDFileFrmt() throw(INTERP_KERNEL::Exception);
    DataArrayInt *convertCellArrayPerGeoType(const DataArrayInt *da) const throw(INTERP_KERNEL::Exception);
    DataArrayInt *zipConnectivityTraducer(int compType) throw(INTERP_KERNEL::Exception);
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
    static MEDCouplingUMesh *Build0DMeshFromCoords(DataArrayDouble *da) throw(INTERP_KERNEL::Exception);
    static MEDCouplingUMesh *MergeUMeshes(const MEDCouplingUMesh *mesh1, const MEDCouplingUMesh *mesh2) throw(INTERP_KERNEL::Exception);
    static MEDCouplingUMesh *MergeUMeshesOnSameCoords(const MEDCouplingUMesh *mesh1, const MEDCouplingUMesh *mesh2) throw(INTERP_KERNEL::Exception);
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
      DataArrayInt *getNodalConnectivity() throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *ret=self->getNodalConnectivity();
        if(ret)
          ret->incrRef();
        return ret;
      }
      DataArrayInt *getNodalConnectivityIndex() throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *ret=self->getNodalConnectivityIndex();
        if(ret)
          ret->incrRef();
        return ret;
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
        int size;
        INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(ids,&size);
        MEDCouplingUMesh *ret=self->keepSpecifiedCells(type,tmp,tmp+size);
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

      PyObject *getLevArrPerCellTypes(PyObject *li) const throw(INTERP_KERNEL::Exception)
      {
        int sz;
        INTERP_KERNEL::AutoPtr<INTERP_KERNEL::NormalizedCellType> order=(INTERP_KERNEL::NormalizedCellType *)convertPyToNewIntArr2(li,&sz);
        DataArrayInt *tmp0,*tmp1=0;
        tmp0=self->getLevArrPerCellTypes(order,(INTERP_KERNEL::NormalizedCellType *)order+sz,tmp1);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(tmp0),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(tmp1),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      static PyObject *MergeUMeshesOnSameCoords(PyObject *ms) throw(INTERP_KERNEL::Exception)
      {
        std::vector<const ParaMEDMEM::MEDCouplingUMesh *> meshes;
        convertPyObjToVecUMeshesCst(ms,meshes);
        MEDCouplingUMesh *ret=MEDCouplingUMesh::MergeUMeshesOnSameCoords(meshes);
        return convertMesh(ret, SWIG_POINTER_OWN | 0 );
      }

      static PyObject *FuseUMeshesOnSameCoords(PyObject *ms, int compType) throw(INTERP_KERNEL::Exception)
      {
        int sz;
        std::vector<const MEDCouplingUMesh *> meshes;
        convertPyObjToVecUMeshesCst(ms,meshes);
        std::vector<DataArrayInt *> corr;
        MEDCouplingUMesh *um=MEDCouplingUMesh::FuseUMeshesOnSameCoords(meshes,compType,corr);
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
      
      static MEDCouplingUMesh *MergeUMeshes(PyObject *li) throw(INTERP_KERNEL::Exception)
      {
        std::vector<const ParaMEDMEM::MEDCouplingUMesh *> tmp;
        convertPyObjToVecUMeshesCst(li,tmp);
        return MEDCouplingUMesh::MergeUMeshes(tmp);
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

      PyObject *buildDescendingConnectivity() const throw(INTERP_KERNEL::Exception)
      {
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> d0=DataArrayInt::New();
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> d1=DataArrayInt::New();
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> d2=DataArrayInt::New();
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> d3=DataArrayInt::New();
        MEDCouplingUMesh *m=self->buildDescendingConnectivity(d0,d1,d2,d3);
        PyObject *ret=PyTuple_New(5);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(m),SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(d0),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(d1),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,3,SWIG_NewPointerObj(SWIG_as_voidptr(d2),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,4,SWIG_NewPointerObj(SWIG_as_voidptr(d3),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        d0->incrRef();
        d1->incrRef();
        d2->incrRef();
        d3->incrRef();
        return ret;
      }
      
      PyObject *emulateMEDMEMBDC(const MEDCouplingUMesh *nM1LevMesh)
      {
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> d0=DataArrayInt::New();
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> d1=DataArrayInt::New();
        DataArrayInt *d2,*d3,*d4,*dd5;
        MEDCouplingUMesh *mOut=self->emulateMEDMEMBDC(nM1LevMesh,d0,d1,d2,d3,d4,dd5);
        PyObject *ret=PyTuple_New(7);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(mOut),SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr((DataArrayInt *)d0),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr((DataArrayInt *)d1),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,3,SWIG_NewPointerObj(SWIG_as_voidptr(d2),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,4,SWIG_NewPointerObj(SWIG_as_voidptr(d3),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,5,SWIG_NewPointerObj(SWIG_as_voidptr(d4),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,6,SWIG_NewPointerObj(SWIG_as_voidptr(dd5),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        d0->incrRef();
        d1->incrRef();
        return ret;
      }

      PyObject *getReverseNodalConnectivity() const throw(INTERP_KERNEL::Exception)
      {
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> d0=DataArrayInt::New();
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> d1=DataArrayInt::New();
        self->getReverseNodalConnectivity(d0,d1);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(d0),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(d1),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        d0->incrRef();
        d1->incrRef();
        return ret;
      }

      DataArrayDouble *getPartBarycenterAndOwner(DataArrayInt *da) const throw(INTERP_KERNEL::Exception)
      {
        if(!da)
          throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
        da->checkAllocated();
        return self->getPartBarycenterAndOwner(da->getConstPointer(),da->getConstPointer()+da->getNbOfElems());
      }

      DataArrayDouble *getPartMeasureField(bool isAbs, DataArrayInt *da) const throw(INTERP_KERNEL::Exception)
      {
        if(!da)
          throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
        da->checkAllocated();
        return self->getPartMeasureField(isAbs,da->getConstPointer(),da->getConstPointer()+da->getNbOfElems());
      }

      MEDCouplingFieldDouble *buildPartOrthogonalField(DataArrayInt *da) const throw(INTERP_KERNEL::Exception)
      {
        if(!da)
          throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
        da->checkAllocated();
        return self->buildPartOrthogonalField(da->getConstPointer(),da->getConstPointer()+da->getNbOfElems());
      }

      PyObject *getTypesOfPart(DataArrayInt *da) const throw(INTERP_KERNEL::Exception)
      {
        if(!da)
          throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
        da->checkAllocated();
        std::set<INTERP_KERNEL::NormalizedCellType> result=self->getTypesOfPart(da->getConstPointer(),da->getConstPointer()+da->getNbOfElems());
        std::set<INTERP_KERNEL::NormalizedCellType>::const_iterator iL=result.begin();
        PyObject *res = PyList_New(result.size());
        for (int i=0;iL!=result.end(); i++, iL++)
          PyList_SetItem(res,i,PyInt_FromLong(*iL));
        return res;
      }

      DataArrayInt *keepCellIdsByType(INTERP_KERNEL::NormalizedCellType type, DataArrayInt *da) const throw(INTERP_KERNEL::Exception)
      {
        if(!da)
          throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
        da->checkAllocated();
        DataArrayInt *ret=self->keepCellIdsByType(type,da->getConstPointer(),da->getConstPointer()+da->getNbOfElems());
        ret->setName(da->getName().c_str());
        return ret;
      }

      DataArrayInt *getCellIdsLyingOnNodes(PyObject *li, bool fullyIn) const throw(INTERP_KERNEL::Exception)
      {
        void *da=0;
        int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 |  0 );
        if (!SWIG_IsOK(res1))
          {
            int size;
            INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
            return self->getCellIdsLyingOnNodes(tmp,((const int *)tmp)+size,fullyIn);
          }
        else
          {
            DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
            if(!da2)
              throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
            da2->checkAllocated();
            return self->getCellIdsLyingOnNodes(da2->getConstPointer(),da2->getConstPointer()+da2->getNbOfElems(),fullyIn);
          }
      }
    }
    void convertToPolyTypes(const std::vector<int>& cellIdsToConvert) throw(INTERP_KERNEL::Exception);
    void convertAllToPoly();
    void unPolyze() throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *buildExtrudedMesh(const MEDCouplingUMesh *mesh1D, int policy) throw(INTERP_KERNEL::Exception);
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
    void setCoords(const DataArrayDouble *coordsX,
                   const DataArrayDouble *coordsY=0,
                   const DataArrayDouble *coordsZ=0) throw(INTERP_KERNEL::Exception);
    void setCoordsAt(int i, const DataArrayDouble *arr) throw(INTERP_KERNEL::Exception);
    %extend {
      std::string __str__() const
      {
        return self->simpleRepr();
      }
      DataArrayDouble *getCoordsAt(int i) throw(INTERP_KERNEL::Exception)
      {
        DataArrayDouble *ret=self->getCoordsAt(i);
        if(ret)
          ret->incrRef();
        return ret;
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
     double *tmp=new double[nbOfTuples*nbOfElsPerTuple];
     try
       {
         fillArrayWithPyListDbl(li,tmp,nbOfTuples*nbOfElsPerTuple,0.);
       }
     catch(INTERP_KERNEL::Exception& e)
       {
         delete [] tmp;
         throw e;
       }
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
     void *da=0;
     int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 |  0 );
     if (!SWIG_IsOK(res1))
       {
         int size;
         INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
         if(size!=self->getNumberOfTuples())
           {
             throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
           }
         return self->renumber(tmp);
       }
     else
       {
         DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
         if(!da2)
           throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
         da2->checkAllocated();
         int size=self->getNumberOfTuples();
         if(size!=self->getNumberOfTuples())
           {
             throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
           }
         return self->renumber(da2->getConstPointer());
       }
   }

   DataArrayDouble *renumberR(PyObject *li) throw(INTERP_KERNEL::Exception)
   {
     void *da=0;
     int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 |  0 );
     if (!SWIG_IsOK(res1))
       {
         int size;
         INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
         if(size!=self->getNumberOfTuples())
           {
             throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
           }
         return self->renumberR(tmp);
       }
     else
       {
         DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
         if(!da2)
           throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
         da2->checkAllocated();
         int size=self->getNumberOfTuples();
         if(size!=self->getNumberOfTuples())
           {
             throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
           }
         return self->renumberR(da2->getConstPointer());
       }
   }

   DataArrayDouble *renumberAndReduce(PyObject *li, int newNbOfTuple) throw(INTERP_KERNEL::Exception)
   {
     void *da=0;
     int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 |  0 );
     if (!SWIG_IsOK(res1))
       {
         int size;
         INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
         if(size!=self->getNumberOfTuples())
           {
             throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
           }
         return self->renumberAndReduce(tmp,newNbOfTuple);
       }
     else
       {
         DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
         if(!da2)
           throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
         da2->checkAllocated();
         int size=self->getNumberOfTuples();
         if(size!=self->getNumberOfTuples())
           {
             throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
           }
         return self->renumberAndReduce(da2->getConstPointer(),newNbOfTuple);
       }
   }

   void renumberInPlace(PyObject *li) throw(INTERP_KERNEL::Exception)
   {
     void *da=0;
     int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 |  0 );
     if (!SWIG_IsOK(res1))
       {
         int size;
         INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
         if(size!=self->getNumberOfTuples())
           {
             throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
           }
         self->renumberInPlace(tmp);
       }
     else
       {
         DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
         if(!da2)
           throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
         da2->checkAllocated();
         int size=self->getNumberOfTuples();
         if(size!=self->getNumberOfTuples())
           {
             throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
           }
         self->renumberInPlace(da2->getConstPointer());
       }
   }

   void renumberInPlaceR(PyObject *li) throw(INTERP_KERNEL::Exception)
   {
     void *da=0;
     int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 |  0 );
     if (!SWIG_IsOK(res1))
       {
         int size;
         INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
         if(size!=self->getNumberOfTuples())
           {
             throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
           }
         self->renumberInPlaceR(tmp);
       }
     else
       {
         DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
         if(!da2)
           throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
         da2->checkAllocated();
         int size=self->getNumberOfTuples();
         if(size!=self->getNumberOfTuples())
           {
             throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
           }
         self->renumberInPlaceR(da2->getConstPointer());
       }
   }

   DataArrayDouble *selectByTupleId(PyObject *li) const throw(INTERP_KERNEL::Exception)
   {
     void *da=0;
     int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 |  0 );
     if (!SWIG_IsOK(res1))
       {
         int size;
         INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
         return self->selectByTupleId(tmp,tmp+size);
       }
     else
       {
         DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
         if(!da2)
           throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
         da2->checkAllocated();
         return self->selectByTupleId(da2->getConstPointer(),da2->getConstPointer()+da2->getNbOfElems());
       }
   }

   DataArrayDouble *selectByTupleIdSafe(PyObject *li) const throw(INTERP_KERNEL::Exception)
   {
     void *da=0;
     int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 |  0 );
     if (!SWIG_IsOK(res1))
       {
         int size;
         INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
         return self->selectByTupleIdSafe(tmp,tmp+size);
       }
     else
       {
         DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
         if(!da2)
           throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
         da2->checkAllocated();
         return self->selectByTupleIdSafe(da2->getConstPointer(),da2->getConstPointer()+da2->getNbOfElems());
       }
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

   static DataArrayDouble *Aggregate(PyObject *li) throw(INTERP_KERNEL::Exception)
   {
     std::vector<const DataArrayDouble *> tmp;
     convertPyObjToVecDataArrayDblCst(li,tmp);
     return DataArrayDouble::Aggregate(tmp);
   }

   static DataArrayDouble *Meld(PyObject *li) throw(INTERP_KERNEL::Exception)
   {
     std::vector<const DataArrayDouble *> tmp;
     convertPyObjToVecDataArrayDblCst(li,tmp);
     return DataArrayDouble::Meld(tmp);
   }

   DataArrayDouble *__getitem__(PyObject *obj) throw(INTERP_KERNEL::Exception)
   {
     const char msg[]="Unexpected situation in __getitem__ !";
     self->checkAllocated();
     int nbOfTuples=self->getNumberOfTuples();
     int nbOfComponents=self->getNumberOfComponents();
     int it1,ic1;
     std::vector<int> vt1,vc1;
     std::pair<int, std::pair<int,int> > pt1,pc1;
     DataArrayInt *dt1=0,*dc1=0;
     int sw;
     convertObjToPossibleCpp3(obj,nbOfTuples,nbOfComponents,sw,it1,ic1,vt1,vc1,pt1,pc1,dt1,dc1);
     MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret;
     switch(sw)
       {
       case 1:
         return self->selectByTupleIdSafe(&it1,&it1+1);
       case 2:
         return self->selectByTupleIdSafe(&vt1[0],&vt1[0]+vt1.size());
       case 3:
         return self->selectByTupleId2(pt1.first,pt1.second.first,pt1.second.second);
       case 4:
         return self->selectByTupleIdSafe(dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems());
       case 5:
         {
           ret=self->selectByTupleIdSafe(&it1,&it1+1);
           std::vector<int> v2(1,ic1);
           return ret->keepSelectedComponents(v2);
         }
       case 6:
         {
           ret=self->selectByTupleIdSafe(&vt1[0],&vt1[0]+vt1.size());
           std::vector<int> v2(1,ic1);
           return ret->keepSelectedComponents(v2);
         }
       case 7:
         {
           ret=self->selectByTupleId2(pt1.first,pt1.second.first,pt1.second.second);
           std::vector<int> v2(1,ic1);
           return ret->keepSelectedComponents(v2);
         }
       case 8:
         {
           ret=self->selectByTupleIdSafe(dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems());
           std::vector<int> v2(1,ic1);
           return ret->keepSelectedComponents(v2);
         }
       case 9:
         {
           ret=self->selectByTupleIdSafe(&it1,&it1+1);
           return ret->keepSelectedComponents(vc1);
         }
       case 10:
         {
           ret=self->selectByTupleIdSafe(&vt1[0],&vt1[0]+vt1.size());
           return ret->keepSelectedComponents(vc1);
         }
       case 11:
         {
           ret=self->selectByTupleId2(pt1.first,pt1.second.first,pt1.second.second);
           return ret->keepSelectedComponents(vc1);
         }
       case 12:
         {
           ret=self->selectByTupleIdSafe(dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems());
           return ret->keepSelectedComponents(vc1);
         }
       case 13:
         {
           ret=self->selectByTupleIdSafe(&it1,&it1+1);
           int nbOfComp=(pc1.second.first-1-pc1.first)/pc1.second.second+1;
           std::vector<int> v2(nbOfComp);
           for(int i=0;i<nbOfComp;i++)
             v2[i]=pc1.first+i*pc1.second.second;
           return ret->keepSelectedComponents(v2);
         }
       case 14:
         {
           ret=self->selectByTupleIdSafe(&vt1[0],&vt1[0]+vt1.size());
           int nbOfComp=(pc1.second.first-1-pc1.first)/pc1.second.second+1;
           std::vector<int> v2(nbOfComp);
           for(int i=0;i<nbOfComp;i++)
             v2[i]=pc1.first+i*pc1.second.second;
           return ret->keepSelectedComponents(v2);
         }
       case 15:
         {
           ret=self->selectByTupleId2(pt1.first,pt1.second.first,pt1.second.second);
           int nbOfComp=(pc1.second.first-1-pc1.first)/pc1.second.second+1;
           std::vector<int> v2(nbOfComp);
           for(int i=0;i<nbOfComp;i++)
             v2[i]=pc1.first+i*pc1.second.second;
           return ret->keepSelectedComponents(v2);
         }
       case 16:
         {
           ret=self->selectByTupleIdSafe(dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems());
           int nbOfComp=(pc1.second.first-1-pc1.first)/pc1.second.second+1;
           std::vector<int> v2(nbOfComp);
           for(int i=0;i<nbOfComp;i++)
             v2[i]=pc1.first+i*pc1.second.second;
           return ret->keepSelectedComponents(v2);
         }
       default:
         throw INTERP_KERNEL::Exception(msg);
       }
   }

   DataArrayDouble *__setitem__(PyObject *obj, PyObject *value) throw(INTERP_KERNEL::Exception)
   {
     self->checkAllocated();
     const char msg[]="Unexpected situation in __setitem__ !";
     int nbOfTuples=self->getNumberOfTuples();
     int nbOfComponents=self->getNumberOfComponents();
     int sw1,sw2;
     double i1;
     std::vector<double> v1;
     DataArrayDouble *d1=0;
     convertObjToPossibleCpp4(value,sw1,i1,v1,d1);
     int it1,ic1;
     std::vector<int> vt1,vc1;
     std::pair<int, std::pair<int,int> > pt1,pc1;
     DataArrayInt *dt1=0,*dc1=0;
     convertObjToPossibleCpp3(obj,nbOfTuples,nbOfComponents,sw2,it1,ic1,vt1,vc1,pt1,pc1,dt1,dc1);
     MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> tmp;
     switch(sw2)
       {
       case 1:
         {
           switch(sw1)
             {
             case 1:
               self->setPartOfValuesSimple1(i1,it1,it1+1,1,0,nbOfComponents,1);
               return self;
             case 2:
               tmp=DataArrayDouble::New();
               tmp->useArray(&v1[0],false,CPP_DEALLOC,v1.size(),1);
               self->setPartOfValues1(tmp,it1,it1+1,1,0,nbOfComponents,1,false);
               return self;
             case 3:
               self->setPartOfValues1(d1,it1,it1+1,1,0,nbOfComponents,1);
               return self;
             default:
               throw INTERP_KERNEL::Exception(msg);
             }
           break;
         }
       case 2:
         {
           switch(sw1)
             {
             case 1:
               self->setPartOfValuesSimple3(i1,&vt1[0],&vt1[0]+vt1.size(),0,nbOfComponents,1);
               return self;
             case 2:
               tmp=DataArrayDouble::New();
               tmp->useArray(&v1[0],false,CPP_DEALLOC,v1.size(),1);
               self->setPartOfValues3(tmp,&vt1[0],&vt1[0]+vt1.size(),0,nbOfComponents,1,false);
               return self;
             case 3:
               self->setPartOfValues3(d1,&vt1[0],&vt1[0]+vt1.size(),0,nbOfComponents,1);
               return self;
             default:
               throw INTERP_KERNEL::Exception(msg);
             }
           break;
         }
       case 3:
         {
           switch(sw1)
             {
             case 1:
               self->setPartOfValuesSimple1(i1,pt1.first,pt1.second.first,pt1.second.second,0,nbOfComponents,1);
               return self;
             case 2:
               tmp=DataArrayDouble::New();
               tmp->useArray(&v1[0],false,CPP_DEALLOC,v1.size(),1);
               self->setPartOfValues1(tmp,pt1.first,pt1.second.first,pt1.second.second,0,nbOfComponents,1,false);
               return self;
             case 3:
               self->setPartOfValues1(d1,pt1.first,pt1.second.first,pt1.second.second,0,nbOfComponents,1);
               return self;
             default:
               throw INTERP_KERNEL::Exception(msg);
             }
           break;
         }
       case 4:
         {
           switch(sw1)
             {
             case 1:
               self->setPartOfValuesSimple3(i1,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),0,nbOfComponents,1);
               return self;
             case 2:
               tmp=DataArrayDouble::New();
               tmp->useArray(&v1[0],false,CPP_DEALLOC,v1.size(),1);
               self->setPartOfValues3(tmp,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),0,nbOfComponents,1,false);
               return self;
             case 3:
               self->setPartOfValues3(d1,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),0,nbOfComponents,1);
               return self;
             default:
               throw INTERP_KERNEL::Exception(msg);
             }
           break;
         }
       case 5:
         {
           switch(sw1)
             {
             case 1:
               self->setPartOfValuesSimple1(i1,it1,it1+1,1,ic1,ic1+1,1);
               return self;
             case 2:
               tmp=DataArrayDouble::New();
               tmp->useArray(&v1[0],false,CPP_DEALLOC,v1.size(),1);
               self->setPartOfValues1(tmp,it1,it1+1,1,ic1,ic1+1,1,false);
               return self;
             case 3:
               self->setPartOfValues1(d1,it1,it1+1,1,ic1,ic1+1,1);
               return self;
             default:
               throw INTERP_KERNEL::Exception(msg);
             }
           break;
         }
       case 6:
         {
           switch(sw1)
             {
             case 1:
               self->setPartOfValuesSimple3(i1,&vt1[0],&vt1[0]+vt1.size(),ic1,ic1+1,1);
               return self;
             case 2:
               tmp=DataArrayDouble::New();
               tmp->useArray(&v1[0],false,CPP_DEALLOC,v1.size(),1);
               self->setPartOfValues3(tmp,&vt1[0],&vt1[0]+vt1.size(),ic1,ic1+1,1,false);
               return self;
             case 3:
               self->setPartOfValues3(d1,&vt1[0],&vt1[0]+vt1.size(),ic1,ic1+1,1);
               return self;
             default:
               throw INTERP_KERNEL::Exception(msg);
             }
           break;
         }
       case 7:
         {
           switch(sw1)
             {
             case 1:
               self->setPartOfValuesSimple1(i1,pt1.first,pt1.second.first,pt1.second.second,ic1,ic1+1,1);
               return self;
             case 2:
               tmp=DataArrayDouble::New();
               tmp->useArray(&v1[0],false,CPP_DEALLOC,v1.size(),1);
               self->setPartOfValues1(tmp,pt1.first,pt1.second.first,pt1.second.second,ic1,ic1+1,1,false);
               return self;
             case 3:
               self->setPartOfValues1(d1,pt1.first,pt1.second.first,pt1.second.second,ic1,ic1+1,1);
               return self;
             default:
               throw INTERP_KERNEL::Exception(msg);
             }
           break;
         }
       case 8:
         {
           switch(sw1)
             {
             case 1:
               self->setPartOfValuesSimple3(i1,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),ic1,ic1+1,1);
               return self;
             case 2:
               tmp=DataArrayDouble::New();
               tmp->useArray(&v1[0],false,CPP_DEALLOC,v1.size(),1);
               self->setPartOfValues3(tmp,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),ic1,ic1+1,1,false);
               return self;
             case 3:
               self->setPartOfValues3(d1,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),ic1,ic1+1,1);
               return self;
             default:
               throw INTERP_KERNEL::Exception(msg);
             }
           break;
         }
       case 9:
         {
           switch(sw1)
             {
             case 1:
               self->setPartOfValuesSimple2(i1,&it1,&it1+1,&vc1[0],&vc1[0]+vc1.size());
               return self;
             case 2:
               tmp=DataArrayDouble::New();
               tmp->useArray(&v1[0],false,CPP_DEALLOC,v1.size(),1);
               self->setPartOfValues2(tmp,&it1,&it1+1,&vc1[0],&vc1[0]+vc1.size(),false);
               return self;
             case 3:
               self->setPartOfValues2(d1,&it1,&it1+1,&vc1[0],&vc1[0]+vc1.size());
               return self;
             default:
               throw INTERP_KERNEL::Exception(msg);
             }
           break;
         }
       case 10:
         {
           switch(sw1)
             {
             case 1:
               self->setPartOfValuesSimple2(i1,&vt1[0],&vt1[0]+vt1.size(),&vc1[0],&vc1[0]+vc1.size());
               return self;
             case 2:
               tmp=DataArrayDouble::New();
               tmp->useArray(&v1[0],false,CPP_DEALLOC,v1.size(),1);
               self->setPartOfValues2(tmp,&vt1[0],&vt1[0]+vt1.size(),&vc1[0],&vc1[0]+vc1.size(),false);
               return self;
             case 3:
               self->setPartOfValues2(d1,&vt1[0],&vt1[0]+vt1.size(),&vc1[0],&vc1[0]+vc1.size());
               return self;
             default:
               throw INTERP_KERNEL::Exception(msg);
             }
           break;
         }
       case 11:
         {
           int bb=pt1.first;
           int ee=pt1.second.first;
           int ss=pt1.second.second;
           if(ee<bb || ss<=0)
             throw INTERP_KERNEL::Exception("Invalid slice in tuple selection");
           int nbOfE=(ee-bb)/ss;
           std::vector<int> nv(nbOfE);
           for(int jj=0;jj<nbOfE;jj++)
             nv[jj]=bb+jj*ss;
           switch(sw1)
             {
             case 1:
               self->setPartOfValuesSimple2(i1,&nv[0],&nv[0]+nv.size(),&vc1[0],&vc1[0]+vc1.size());
               return self;
             case 2:
               tmp=DataArrayDouble::New();
               tmp->useArray(&v1[0],false,CPP_DEALLOC,v1.size(),1);
               self->setPartOfValues2(tmp,&nv[0],&nv[0]+nv.size(),&vc1[0],&vc1[0]+vc1.size(),false);
               return self;
             case 3:
               self->setPartOfValues2(d1,&nv[0],&nv[0]+nv.size(),&vc1[0],&vc1[0]+vc1.size());
               return self;
             default:
               throw INTERP_KERNEL::Exception(msg);
             }
           break;
         }
       case 12:
         {
           switch(sw1)
             {
             case 1:
               self->setPartOfValuesSimple2(i1,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),&vc1[0],&vc1[0]+vc1.size());
               return self;
             case 2:
               tmp=DataArrayDouble::New();
               tmp->useArray(&v1[0],false,CPP_DEALLOC,v1.size(),1);
               self->setPartOfValues2(tmp,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),&vc1[0],&vc1[0]+vc1.size(),false);
               return self;
             case 3:
               self->setPartOfValues2(d1,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),&vc1[0],&vc1[0]+vc1.size());
               return self;
             default:
               throw INTERP_KERNEL::Exception(msg);
             }
           break;
         }
       case 13:
         {
           switch(sw1)
             {
             case 1:
               self->setPartOfValuesSimple1(i1,it1,it1+1,1,pc1.first,pc1.second.first,pc1.second.second);
               return self;
             case 2:
               tmp=DataArrayDouble::New();
               tmp->useArray(&v1[0],false,CPP_DEALLOC,v1.size(),1);
               self->setPartOfValues1(tmp,it1,it1+1,1,pc1.first,pc1.second.first,pc1.second.second,false);
               return self;
             case 3:
               self->setPartOfValues1(d1,it1,it1+1,1,pc1.first,pc1.second.first,pc1.second.second);
               return self;
             default:
               throw INTERP_KERNEL::Exception(msg);
             }
           break;
         }
       case 14:
         {
           switch(sw1)
             {
             case 1:
               self->setPartOfValuesSimple3(i1,&vt1[0],&vt1[0]+vt1.size(),pc1.first,pc1.second.first,pc1.second.second);
               return self;
             case 2:
               tmp=DataArrayDouble::New();
               tmp->useArray(&v1[0],false,CPP_DEALLOC,v1.size(),1);
               self->setPartOfValues3(tmp,&vt1[0],&vt1[0]+vt1.size(),pc1.first,pc1.second.first,pc1.second.second,false);
               return self;
             case 3:
               self->setPartOfValues3(d1,&vt1[0],&vt1[0]+vt1.size(),pc1.first,pc1.second.first,pc1.second.second);
               return self;
             default:
               throw INTERP_KERNEL::Exception(msg);
             }
           break;
         }
       case 15:
         {
           switch(sw1)
             {
             case 1:
               self->setPartOfValuesSimple1(i1,pt1.first,pt1.second.first,pt1.second.second,pc1.first,pc1.second.first,pc1.second.second);
               return self;
             case 2:
                 tmp=DataArrayDouble::New();
                 tmp->useArray(&v1[0],false,CPP_DEALLOC,v1.size(),1);
                 self->setPartOfValues1(tmp,pt1.first,pt1.second.first,pt1.second.second,pc1.first,pc1.second.first,pc1.second.second,false);
                 return self;
             case 3:
               self->setPartOfValues1(d1,pt1.first,pt1.second.first,pt1.second.second,pc1.first,pc1.second.first,pc1.second.second);
               return self;
             default:
               throw INTERP_KERNEL::Exception(msg);
             }
           break;
         }
       case 16:
         {
           switch(sw1)
             {
             case 1:
               self->setPartOfValuesSimple3(i1,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),pc1.first,pc1.second.first,pc1.second.second);
               return self;
             case 2:
               tmp=DataArrayDouble::New();
               tmp->useArray(&v1[0],false,CPP_DEALLOC,v1.size(),1);
               self->setPartOfValues3(tmp,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),pc1.first,pc1.second.first,pc1.second.second,false);
               return self;
             case 3:
               self->setPartOfValues3(d1,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),pc1.first,pc1.second.first,pc1.second.second);
               return self;
             default:
               throw INTERP_KERNEL::Exception(msg);
             }
           break;
         }
       default:
         throw INTERP_KERNEL::Exception(msg);
       }
     return self;
   }

   DataArrayDouble *__add__(PyObject *obj) throw(INTERP_KERNEL::Exception)
   {
     const char msg[]="Unexpected situation in __add__ !";
     double val;
     DataArrayDouble *a;
     int sw;
     convertObjToPossibleCpp5(obj,sw,val,a);
     switch(sw)
       {
       case 1:
         {
           MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=self->deepCpy();
           ret->applyLin(1.,val);
           ret->incrRef();
           return ret;
         }
       case 2:
         {
           return DataArrayDouble::Add(self,a);
         }
       default:
         throw INTERP_KERNEL::Exception(msg);
       }
   }

   DataArrayDouble *__radd__(PyObject *obj) throw(INTERP_KERNEL::Exception)
   {
     const char msg[]="Unexpected situation in __radd__ !";
     double val;
     DataArrayDouble *a;
     int sw;
     convertObjToPossibleCpp5(obj,sw,val,a);
     switch(sw)
       {
       case 1:
         {
           MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=self->deepCpy();
           ret->applyLin(1.,val);
           ret->incrRef();
           return ret;
         }
       default:
         throw INTERP_KERNEL::Exception(msg);
       }
   }

   DataArrayDouble *operator+=(PyObject *obj) throw(INTERP_KERNEL::Exception)
   {
     const char msg[]="Unexpected situation in __iadd__ !";
     double val;
     DataArrayDouble *a;
     int sw;
     convertObjToPossibleCpp5(obj,sw,val,a);
     switch(sw)
       {
       case 1:
         {
           self->applyLin(1.,val);
           return self;
         }
       case 2:
         {
           self->addEqual(a);
           return self;
         }
       default:
         throw INTERP_KERNEL::Exception(msg);
       }
   }

   DataArrayDouble *__sub__(PyObject *obj) throw(INTERP_KERNEL::Exception)
   {
     const char msg[]="Unexpected situation in __sub__ !";
     double val;
     DataArrayDouble *a;
     int sw;
     convertObjToPossibleCpp5(obj,sw,val,a);
     switch(sw)
       {
       case 1:
         {
           MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=self->deepCpy();
           ret->applyLin(1.,-val);
           ret->incrRef();
           return ret;
         }
       case 2:
         {
           return DataArrayDouble::Substract(self,a);
         }
       default:
         throw INTERP_KERNEL::Exception(msg);
       }
   }

   DataArrayDouble *__rsub__(PyObject *obj) throw(INTERP_KERNEL::Exception)
   {
     const char msg[]="Unexpected situation in __rsub__ !";
     double val;
     DataArrayDouble *a;
     int sw;
     convertObjToPossibleCpp5(obj,sw,val,a);
     switch(sw)
       {
       case 1:
         {
           MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=self->deepCpy();
           ret->applyLin(-1.,val);
           ret->incrRef();
           return ret;
         }
       default:
         throw INTERP_KERNEL::Exception(msg);
       }
   }

   DataArrayDouble *operator-=(PyObject *obj) throw(INTERP_KERNEL::Exception)
   {
     const char msg[]="Unexpected situation in __isub__ !";
     double val;
     DataArrayDouble *a;
     int sw;
     convertObjToPossibleCpp5(obj,sw,val,a);
     switch(sw)
       {
       case 1:
         {
           self->applyLin(1.,-val);
           return self;
         }
       case 2:
         {
           self->substractEqual(a);
           return self;
         }
       default:
         throw INTERP_KERNEL::Exception(msg);
       }
   }

   DataArrayDouble *__mul__(PyObject *obj) throw(INTERP_KERNEL::Exception)
   {
     const char msg[]="Unexpected situation in __mul__ !";
     double val;
     DataArrayDouble *a;
     int sw;
     convertObjToPossibleCpp5(obj,sw,val,a);
     switch(sw)
       {
       case 1:
         {
           MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=self->deepCpy();
           ret->applyLin(val,0.);
           ret->incrRef();
           return ret;
         }
       case 2:
         {
           return DataArrayDouble::Multiply(self,a);
         }
       default:
         throw INTERP_KERNEL::Exception(msg);
       }
   }

   DataArrayDouble *__rmul__(PyObject *obj) throw(INTERP_KERNEL::Exception)
   {
     const char msg[]="Unexpected situation in __rmul__ !";
     double val;
     DataArrayDouble *a;
     int sw;
     convertObjToPossibleCpp5(obj,sw,val,a);
     switch(sw)
       {
       case 1:
         {
           MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=self->deepCpy();
           ret->applyLin(val,0.);
           ret->incrRef();
           return ret;
         }
       default:
         throw INTERP_KERNEL::Exception(msg);
       }
   }

   DataArrayDouble *operator*=(PyObject *obj) throw(INTERP_KERNEL::Exception)
   {
     const char msg[]="Unexpected situation in __imul__ !";
     double val;
     DataArrayDouble *a;
     int sw;
     convertObjToPossibleCpp5(obj,sw,val,a);
     switch(sw)
       {
       case 1:
         {
           self->applyLin(val,0.);
           return self;
         }
       case 2:
         {
           self->multiplyEqual(a);
           return self;
         }
       default:
         throw INTERP_KERNEL::Exception(msg);
       }
   }

   DataArrayDouble *__div__(PyObject *obj) throw(INTERP_KERNEL::Exception)
   {
     const char msg[]="Unexpected situation in __div__ !";
     double val;
     DataArrayDouble *a;
     int sw;
     convertObjToPossibleCpp5(obj,sw,val,a);
     switch(sw)
       {
       case 1:
         {
           if(val==0.)
             throw INTERP_KERNEL::Exception("DataArrayDouble::__div__ : trying to divide by zero !");
           MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=self->deepCpy();
           ret->applyLin(1/val,0.);
           ret->incrRef();
           return ret;
         }
       case 2:
         {
           return DataArrayDouble::Divide(self,a);
         }
       default:
         throw INTERP_KERNEL::Exception(msg);
       }
   }

   DataArrayDouble *__rdiv__(PyObject *obj) throw(INTERP_KERNEL::Exception)
   {
     const char msg[]="Unexpected situation in __rdiv__ !";
     double val;
     DataArrayDouble *a;
     int sw;
     convertObjToPossibleCpp5(obj,sw,val,a);
     switch(sw)
       {
       case 1:
         {
           MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=self->deepCpy();
           ret->applyInv(val);
           ret->incrRef();
           return ret;
         }
       default:
         throw INTERP_KERNEL::Exception(msg);
       }
   }

   DataArrayDouble *operator/=(PyObject *obj) throw(INTERP_KERNEL::Exception)
   {
     const char msg[]="Unexpected situation in __imul__ !";
     double val;
     DataArrayDouble *a;
     int sw;
     convertObjToPossibleCpp5(obj,sw,val,a);
     switch(sw)
       {
       case 1:
         {
           if(val==0.)
             throw INTERP_KERNEL::Exception("DataArrayDouble::__div__ : trying to divide by zero !");
           self->applyLin(1./val,0.);
           return self;
         }
       case 2:
         {
           self->divideEqual(a);
           return self;
         }
       default:
         throw INTERP_KERNEL::Exception(msg);
       }
   }
 };

%extend ParaMEDMEM::DataArrayInt
 {
   std::string __str__() const
   {
     return self->repr();
   }

   PyObject *getDifferentValues(bool val) const throw(INTERP_KERNEL::Exception)
   {
     std::set<int> ret=self->getDifferentValues();
     return convertIntArrToPyList3(ret);
   }

   void setValues(PyObject *li, int nbOfTuples, int nbOfElsPerTuple) throw(INTERP_KERNEL::Exception)
   {
     int *tmp=new int[nbOfTuples*nbOfElsPerTuple];
     try
       {
         fillArrayWithPyListInt(li,tmp,nbOfTuples*nbOfElsPerTuple,0.);
       }
     catch(INTERP_KERNEL::Exception& e)
       {
         delete [] tmp;
         throw e;
       }
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

   static PyObject *MakePartition(PyObject *gps, int newNb) throw(INTERP_KERNEL::Exception)
   {
     std::vector<const DataArrayInt *> groups;
     std::vector< std::vector<int> > fidsOfGroups;
     convertPyObjToVecDataArrayIntCst(gps,groups);
     ParaMEDMEM::DataArrayInt *ret0=ParaMEDMEM::DataArrayInt::MakePartition(groups,newNb,fidsOfGroups);
     PyObject *ret = PyList_New(2);
     PyList_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
     int sz=fidsOfGroups.size();
     PyObject *ret1 = PyList_New(sz);
     for(int i=0;i<sz;i++)
       PyList_SetItem(ret1,i,convertIntArrToPyList2(fidsOfGroups[i]));
     PyList_SetItem(ret,1,ret1);
     return ret;
   }

   void transformWithIndArr(PyObject *li)
   {
     void *da=0;
     int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 |  0 );
     if (!SWIG_IsOK(res1))
       {
         int size;
         INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
         self->transformWithIndArr(tmp,tmp+size);
       }
     else
       {
         DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
         self->transformWithIndArr(da2->getConstPointer(),da2->getConstPointer()+da2->getNbOfElems());
       }
   }

   PyObject *splitByValueRange(PyObject *li) const throw(INTERP_KERNEL::Exception)
   {
     DataArrayInt *ret0=0,*ret1=0,*ret2=0;
     void *da=0;
     int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 |  0 );
     if (!SWIG_IsOK(res1))
       {
         int size;
         INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
         self->splitByValueRange(tmp,(int *)tmp+size,ret0,ret1,ret2);
       }
     else
       {
         DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
         if(!da2)
           throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
         da2->checkAllocated();
         int size=self->getNumberOfTuples();
         self->splitByValueRange(da2->getConstPointer(),da2->getConstPointer()+size,ret0,ret1,ret2);
       }
     PyObject *ret = PyList_New(3);
     PyList_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
     PyList_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
     PyList_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(ret2),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
     return ret;
   }

   DataArrayInt *transformWithIndArrR(PyObject *li) const
   {
     void *da=0;
     int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 |  0 );
     if (!SWIG_IsOK(res1))
       {
         int size;
         INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
         return self->transformWithIndArrR(tmp,tmp+size);
       }
     else
       {
         DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
         return self->transformWithIndArrR(da2->getConstPointer(),da2->getConstPointer()+da2->getNbOfElems());
       }
   }

   void renumberInPlace(PyObject *li) throw(INTERP_KERNEL::Exception)
   {
     void *da=0;
     int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 |  0 );
     if (!SWIG_IsOK(res1))
       {
         int size;
         INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
         if(size!=self->getNumberOfTuples())
           {
             throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
           }
         self->renumberInPlace(tmp);
       }
     else
       {
         DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
         if(!da2)
           throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
         da2->checkAllocated();
         int size=self->getNumberOfTuples();
         if(size!=self->getNumberOfTuples())
           {
             throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
           }
         self->renumberInPlace(da2->getConstPointer());
       }
   }

   void renumberInPlaceR(PyObject *li) throw(INTERP_KERNEL::Exception)
   {
     void *da=0;
     int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 |  0 );
     if (!SWIG_IsOK(res1))
       {
         int size;
         INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
         if(size!=self->getNumberOfTuples())
           {
             throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
           }
         self->renumberInPlaceR(tmp);
       }
     else
       {
         DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
         if(!da2)
           throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
         da2->checkAllocated();
         int size=self->getNumberOfTuples();
         if(size!=self->getNumberOfTuples())
           {
             throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
           }
         self->renumberInPlaceR(da2->getConstPointer());
       }
   }

   DataArrayInt *renumberAndReduce(PyObject *li, int newNbOfTuple) throw(INTERP_KERNEL::Exception)
   {
     void *da=0;
     int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 |  0 );
     if (!SWIG_IsOK(res1))
       {
         int size;
         INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
         if(size!=self->getNumberOfTuples())
           {
             throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
           }
         return self->renumberAndReduce(tmp,newNbOfTuple);
       }
     else
       {
         DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
         if(!da2)
           throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
         da2->checkAllocated();
         int size=self->getNumberOfTuples();
         if(size!=self->getNumberOfTuples())
           {
             throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
           }
         return self->renumberAndReduce(da2->getConstPointer(),newNbOfTuple);
       }
   }

   DataArrayInt *renumber(PyObject *li) throw(INTERP_KERNEL::Exception)
   {
     void *da=0;
     int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 |  0 );
     if (!SWIG_IsOK(res1))
       {
         int size;
         INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
         if(size!=self->getNumberOfTuples())
           {
             throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
           }
         return self->renumber(tmp);
       }
     else
       {
         DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
         if(!da2)
           throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
         da2->checkAllocated();
         int size=self->getNumberOfTuples();
         if(size!=self->getNumberOfTuples())
           {
             throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
           }
         return self->renumber(da2->getConstPointer());
       }
   }

   DataArrayInt *renumberR(PyObject *li) throw(INTERP_KERNEL::Exception)
   {
     void *da=0;
     int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 |  0 );
     if (!SWIG_IsOK(res1))
       {
         int size;
         INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
         if(size!=self->getNumberOfTuples())
           {
             throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
           }
         return self->renumberR(tmp);
       }
     else
       {
         DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
         if(!da2)
           throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
         da2->checkAllocated();
         int size=self->getNumberOfTuples();
         if(size!=self->getNumberOfTuples())
           {
             throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
           }
         return self->renumberR(da2->getConstPointer());
       }
   }

   DataArrayInt *selectByTupleId(PyObject *li) const throw(INTERP_KERNEL::Exception)
   {
     void *da=0;
     int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 |  0 );
     if (!SWIG_IsOK(res1))
       {
         int size;
         INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
         return self->selectByTupleId(tmp,tmp+size);
       }
     else
       {
         DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
         if(!da2)
          throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
         da2->checkAllocated();
         return self->selectByTupleId(da2->getConstPointer(),da2->getConstPointer()+da2->getNbOfElems());
       }
   }

   DataArrayInt *selectByTupleIdSafe(PyObject *li) const throw(INTERP_KERNEL::Exception)
   {
     void *da=0;
     int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 |  0 );
     if (!SWIG_IsOK(res1))
       {
         int size;
         INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
         return self->selectByTupleIdSafe(tmp,tmp+size);
       }
     else
       {
         DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
         if(!da2)
          throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
         da2->checkAllocated();
         return self->selectByTupleIdSafe(da2->getConstPointer(),da2->getConstPointer()+da2->getNbOfElems());
       }
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

   PyObject *changeSurjectiveFormat(int targetNb) const throw(INTERP_KERNEL::Exception)
   {
     DataArrayInt *arr=0;
     DataArrayInt *arrI=0;
     self->changeSurjectiveFormat(targetNb,arr,arrI);
     PyObject *res = PyList_New(2);
     PyList_SetItem(res,0,SWIG_NewPointerObj((void*)arr,SWIGTYPE_p_ParaMEDMEM__DataArrayInt,SWIG_POINTER_OWN | 0));
     PyList_SetItem(res,1,SWIG_NewPointerObj((void*)arrI,SWIGTYPE_p_ParaMEDMEM__DataArrayInt,SWIG_POINTER_OWN | 0));
     return res;
   }

   static DataArrayInt *Meld(PyObject *li) throw(INTERP_KERNEL::Exception)
   {
     std::vector<const DataArrayInt *> tmp;
     convertPyObjToVecDataArrayIntCst(li,tmp);
     return DataArrayInt::Meld(tmp);
   }

   static DataArrayInt *Aggregate(PyObject *li) throw(INTERP_KERNEL::Exception)
   {
     std::vector<const DataArrayInt *> tmp;
     convertPyObjToVecDataArrayIntCst(li,tmp);
     return DataArrayInt::Aggregate(tmp);
   }

   static DataArrayInt *BuildUnion(PyObject *li) throw(INTERP_KERNEL::Exception)
   {
     std::vector<const DataArrayInt *> tmp;
     convertPyObjToVecDataArrayIntCst(li,tmp);
     return DataArrayInt::BuildUnion(tmp);
   }

   static DataArrayInt *BuildIntersection(PyObject *li) throw(INTERP_KERNEL::Exception)
   {
     std::vector<const DataArrayInt *> tmp;
     convertPyObjToVecDataArrayIntCst(li,tmp);
     return DataArrayInt::BuildIntersection(tmp);
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

   DataArrayInt *__getitem__(PyObject *obj) throw(INTERP_KERNEL::Exception)
   {
     const char msg[]="Unexpected situation in __getitem__ !";
     self->checkAllocated();
     int nbOfTuples=self->getNumberOfTuples();
     int nbOfComponents=self->getNumberOfComponents();
     int it1,ic1;
     std::vector<int> vt1,vc1;
     std::pair<int, std::pair<int,int> > pt1,pc1;
     DataArrayInt *dt1=0,*dc1=0;
     int sw;
     convertObjToPossibleCpp3(obj,nbOfTuples,nbOfComponents,sw,it1,ic1,vt1,vc1,pt1,pc1,dt1,dc1);
     MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret;
     switch(sw)
       {
       case 1:
         return self->selectByTupleIdSafe(&it1,&it1+1);
       case 2:
         return self->selectByTupleIdSafe(&vt1[0],&vt1[0]+vt1.size());
       case 3:
         return self->selectByTupleId2(pt1.first,pt1.second.first,pt1.second.second);
       case 4:
         return self->selectByTupleIdSafe(dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems());
       case 5:
         {
           ret=self->selectByTupleIdSafe(&it1,&it1+1);
           std::vector<int> v2(1,ic1);
           return ret->keepSelectedComponents(v2);
         }
       case 6:
         {
           ret=self->selectByTupleIdSafe(&vt1[0],&vt1[0]+vt1.size());
           std::vector<int> v2(1,ic1);
           return ret->keepSelectedComponents(v2);
         }
       case 7:
         {
           ret=self->selectByTupleId2(pt1.first,pt1.second.first,pt1.second.second);
           std::vector<int> v2(1,ic1);
           return ret->keepSelectedComponents(v2);
         }
       case 8:
         {
           ret=self->selectByTupleIdSafe(dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems());
           std::vector<int> v2(1,ic1);
           return ret->keepSelectedComponents(v2);
         }
       case 9:
         {
           ret=self->selectByTupleIdSafe(&it1,&it1+1);
           return ret->keepSelectedComponents(vc1);
         }
       case 10:
         {
           ret=self->selectByTupleIdSafe(&vt1[0],&vt1[0]+vt1.size());
           return ret->keepSelectedComponents(vc1);
         }
       case 11:
         {
           ret=self->selectByTupleId2(pt1.first,pt1.second.first,pt1.second.second);
           return ret->keepSelectedComponents(vc1);
         }
       case 12:
         {
           ret=self->selectByTupleIdSafe(dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems());
           return ret->keepSelectedComponents(vc1);
         }
       case 13:
         {
           ret=self->selectByTupleIdSafe(&it1,&it1+1);
           int nbOfComp=(pc1.second.first-1-pc1.first)/pc1.second.second+1;
           std::vector<int> v2(nbOfComp);
           for(int i=0;i<nbOfComp;i++)
             v2[i]=pc1.first+i*pc1.second.second;
           return ret->keepSelectedComponents(v2);
         }
       case 14:
         {
           ret=self->selectByTupleIdSafe(&vt1[0],&vt1[0]+vt1.size());
           int nbOfComp=(pc1.second.first-1-pc1.first)/pc1.second.second+1;
           std::vector<int> v2(nbOfComp);
           for(int i=0;i<nbOfComp;i++)
             v2[i]=pc1.first+i*pc1.second.second;
           return ret->keepSelectedComponents(v2);
         }
       case 15:
         {
           ret=self->selectByTupleId2(pt1.first,pt1.second.first,pt1.second.second);
           int nbOfComp=(pc1.second.first-1-pc1.first)/pc1.second.second+1;
           std::vector<int> v2(nbOfComp);
           for(int i=0;i<nbOfComp;i++)
             v2[i]=pc1.first+i*pc1.second.second;
           return ret->keepSelectedComponents(v2);
         }
       case 16:
         {
           ret=self->selectByTupleIdSafe(dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems());
           int nbOfComp=(pc1.second.first-1-pc1.first)/pc1.second.second+1;
           std::vector<int> v2(nbOfComp);
           for(int i=0;i<nbOfComp;i++)
             v2[i]=pc1.first+i*pc1.second.second;
           return ret->keepSelectedComponents(v2);
         }
       default:
         throw INTERP_KERNEL::Exception(msg);
       }
   }

   DataArrayInt *__setitem__(PyObject *obj, PyObject *value) throw(INTERP_KERNEL::Exception)
   {
     self->checkAllocated();
     const char msg[]="Unexpected situation in __setitem__ !";
     int nbOfTuples=self->getNumberOfTuples();
     int nbOfComponents=self->getNumberOfComponents();
     int sw1,sw2;
     int i1;
     std::vector<int> v1;
     DataArrayInt *d1=0;
     convertObjToPossibleCpp1(value,sw1,i1,v1,d1);
     int it1,ic1;
     std::vector<int> vt1,vc1;
     std::pair<int, std::pair<int,int> > pt1,pc1;
     DataArrayInt *dt1=0,*dc1=0;
     convertObjToPossibleCpp3(obj,nbOfTuples,nbOfComponents,sw2,it1,ic1,vt1,vc1,pt1,pc1,dt1,dc1);
     MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmp;
     switch(sw2)
       {
       case 1:
         {
           switch(sw1)
             {
             case 1:
               self->setPartOfValuesSimple1(i1,it1,it1+1,1,0,nbOfComponents,1);
               return self;
             case 2:
               tmp=DataArrayInt::New();
               tmp->useArray(&v1[0],false,CPP_DEALLOC,v1.size(),1);
               self->setPartOfValues1(tmp,it1,it1+1,1,0,nbOfComponents,1,false);
               return self;
             case 3:
               self->setPartOfValues1(d1,it1,it1+1,1,0,nbOfComponents,1);
               return self;
             default:
               throw INTERP_KERNEL::Exception(msg);
             }
           break;
         }
       case 2:
         {
           switch(sw1)
             {
             case 1:
               self->setPartOfValuesSimple3(i1,&vt1[0],&vt1[0]+vt1.size(),0,nbOfComponents,1);
               return self;
             case 2:
               tmp=DataArrayInt::New();
               tmp->useArray(&v1[0],false,CPP_DEALLOC,v1.size(),1);
               self->setPartOfValues3(tmp,&vt1[0],&vt1[0]+vt1.size(),0,nbOfComponents,1,false);
               return self;
             case 3:
               self->setPartOfValues3(d1,&vt1[0],&vt1[0]+vt1.size(),0,nbOfComponents,1);
               return self;
             default:
               throw INTERP_KERNEL::Exception(msg);
             }
           break;
         }
       case 3:
         {
           switch(sw1)
             {
             case 1:
               self->setPartOfValuesSimple1(i1,pt1.first,pt1.second.first,pt1.second.second,0,nbOfComponents,1);
               return self;
             case 2:
               tmp=DataArrayInt::New();
               tmp->useArray(&v1[0],false,CPP_DEALLOC,v1.size(),1);
               self->setPartOfValues1(tmp,pt1.first,pt1.second.first,pt1.second.second,0,nbOfComponents,1,false);
               return self;
             case 3:
               self->setPartOfValues1(d1,pt1.first,pt1.second.first,pt1.second.second,0,nbOfComponents,1);
               return self;
             default:
               throw INTERP_KERNEL::Exception(msg);
             }
           break;
         }
       case 4:
         {
           switch(sw1)
             {
             case 1:
               self->setPartOfValuesSimple3(i1,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),0,nbOfComponents,1);
               return self;
             case 2:
               tmp=DataArrayInt::New();
               tmp->useArray(&v1[0],false,CPP_DEALLOC,v1.size(),1);
               self->setPartOfValues3(tmp,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),0,nbOfComponents,1,false);
               return self;
             case 3:
               self->setPartOfValues3(d1,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),0,nbOfComponents,1);
               return self;
             default:
               throw INTERP_KERNEL::Exception(msg);
             }
           break;
         }
       case 5:
         {
           switch(sw1)
             {
             case 1:
               self->setPartOfValuesSimple1(i1,it1,it1+1,1,ic1,ic1+1,1);
               return self;
             case 2:
               tmp=DataArrayInt::New();
               tmp->useArray(&v1[0],false,CPP_DEALLOC,v1.size(),1);
               self->setPartOfValues1(tmp,it1,it1+1,1,ic1,ic1+1,1,false);
               return self;
             case 3:
               self->setPartOfValues1(d1,it1,it1+1,1,ic1,ic1+1,1);
               return self;
             default:
               throw INTERP_KERNEL::Exception(msg);
             }
           break;
         }
       case 6:
         {
           switch(sw1)
             {
             case 1:
               self->setPartOfValuesSimple3(i1,&vt1[0],&vt1[0]+vt1.size(),ic1,ic1+1,1);
               return self;
             case 2:
               tmp=DataArrayInt::New();
               tmp->useArray(&v1[0],false,CPP_DEALLOC,v1.size(),1);
               self->setPartOfValues3(tmp,&vt1[0],&vt1[0]+vt1.size(),ic1,ic1+1,1,false);
               return self;
             case 3:
               self->setPartOfValues3(d1,&vt1[0],&vt1[0]+vt1.size(),ic1,ic1+1,1);
               return self;
             default:
               throw INTERP_KERNEL::Exception(msg);
             }
           break;
         }
       case 7:
         {
           switch(sw1)
             {
             case 1:
               self->setPartOfValuesSimple1(i1,pt1.first,pt1.second.first,pt1.second.second,ic1,ic1+1,1);
               return self;
             case 2:
               tmp=DataArrayInt::New();
               tmp->useArray(&v1[0],false,CPP_DEALLOC,v1.size(),1);
               self->setPartOfValues1(tmp,pt1.first,pt1.second.first,pt1.second.second,ic1,ic1+1,1,false);
               return self;
             case 3:
               self->setPartOfValues1(d1,pt1.first,pt1.second.first,pt1.second.second,ic1,ic1+1,1);
               return self;
             default:
               throw INTERP_KERNEL::Exception(msg);
             }
           break;
         }
       case 8:
         {
           switch(sw1)
             {
             case 1:
               self->setPartOfValuesSimple3(i1,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),ic1,ic1+1,1);
               return self;
             case 2:
               tmp=DataArrayInt::New();
               tmp->useArray(&v1[0],false,CPP_DEALLOC,v1.size(),1);
               self->setPartOfValues3(tmp,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),ic1,ic1+1,1,false);
               return self;
             case 3:
               self->setPartOfValues3(d1,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),ic1,ic1+1,1);
               return self;
             default:
               throw INTERP_KERNEL::Exception(msg);
             }
           break;
         }
       case 9:
         {
           switch(sw1)
             {
             case 1:
               self->setPartOfValuesSimple2(i1,&it1,&it1+1,&vc1[0],&vc1[0]+vc1.size());
               return self;
             case 2:
               tmp=DataArrayInt::New();
               tmp->useArray(&v1[0],false,CPP_DEALLOC,v1.size(),1);
               self->setPartOfValues2(tmp,&it1,&it1+1,&vc1[0],&vc1[0]+vc1.size(),false);
               return self;
             case 3:
               self->setPartOfValues2(d1,&it1,&it1+1,&vc1[0],&vc1[0]+vc1.size());
               return self;
             default:
               throw INTERP_KERNEL::Exception(msg);
             }
           break;
         }
       case 10:
         {
           switch(sw1)
             {
             case 1:
               self->setPartOfValuesSimple2(i1,&vt1[0],&vt1[0]+vt1.size(),&vc1[0],&vc1[0]+vc1.size());
               return self;
             case 2:
               tmp=DataArrayInt::New();
               tmp->useArray(&v1[0],false,CPP_DEALLOC,v1.size(),1);
               self->setPartOfValues2(tmp,&vt1[0],&vt1[0]+vt1.size(),&vc1[0],&vc1[0]+vc1.size(),false);
               return self;
             case 3:
               self->setPartOfValues2(d1,&vt1[0],&vt1[0]+vt1.size(),&vc1[0],&vc1[0]+vc1.size());
               return self;
             default:
               throw INTERP_KERNEL::Exception(msg);
             }
           break;
         }
       case 11:
         {
           int bb=pt1.first;
           int ee=pt1.second.first;
           int ss=pt1.second.second;
           if(ee<bb || ss<=0)
             throw INTERP_KERNEL::Exception("Invalid slice in tuple selection");
           int nbOfE=(ee-bb)/ss;
           std::vector<int> nv(nbOfE);
           for(int jj=0;jj<nbOfE;jj++)
             nv[jj]=bb+jj*ss;
           switch(sw1)
             {
             case 1:
               self->setPartOfValuesSimple2(i1,&nv[0],&nv[0]+nv.size(),&vc1[0],&vc1[0]+vc1.size());
               return self;
             case 2:
               tmp=DataArrayInt::New();
               tmp->useArray(&v1[0],false,CPP_DEALLOC,v1.size(),1);
               self->setPartOfValues2(tmp,&nv[0],&nv[0]+nv.size(),&vc1[0],&vc1[0]+vc1.size(),false);
               return self;
             case 3:
               self->setPartOfValues2(d1,&nv[0],&nv[0]+nv.size(),&vc1[0],&vc1[0]+vc1.size());
               return self;
             default:
               throw INTERP_KERNEL::Exception(msg);
             }
           break;
         }
       case 12:
         {
           switch(sw1)
             {
             case 1:
               self->setPartOfValuesSimple2(i1,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),&vc1[0],&vc1[0]+vc1.size());
               return self;
             case 2:
               tmp=DataArrayInt::New();
               tmp->useArray(&v1[0],false,CPP_DEALLOC,v1.size(),1);
               self->setPartOfValues2(tmp,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),&vc1[0],&vc1[0]+vc1.size(),false);
               return self;
             case 3:
               self->setPartOfValues2(d1,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),&vc1[0],&vc1[0]+vc1.size());
               return self;
             default:
               throw INTERP_KERNEL::Exception(msg);
             }
           break;
         }
       case 13:
         {
           switch(sw1)
             {
             case 1:
               self->setPartOfValuesSimple1(i1,it1,it1+1,1,pc1.first,pc1.second.first,pc1.second.second);
               return self;
             case 2:
               tmp=DataArrayInt::New();
               tmp->useArray(&v1[0],false,CPP_DEALLOC,v1.size(),1);
               self->setPartOfValues1(tmp,it1,it1+1,1,pc1.first,pc1.second.first,pc1.second.second,false);
               return self;
             case 3:
               self->setPartOfValues1(d1,it1,it1+1,1,pc1.first,pc1.second.first,pc1.second.second);
               return self;
             default:
               throw INTERP_KERNEL::Exception(msg);
             }
           break;
         }
       case 14:
         {
           switch(sw1)
             {
             case 1:
               self->setPartOfValuesSimple3(i1,&vt1[0],&vt1[0]+vt1.size(),pc1.first,pc1.second.first,pc1.second.second);
               return self;
             case 2:
               tmp=DataArrayInt::New();
               tmp->useArray(&v1[0],false,CPP_DEALLOC,v1.size(),1);
               self->setPartOfValues3(tmp,&vt1[0],&vt1[0]+vt1.size(),pc1.first,pc1.second.first,pc1.second.second,false);
               return self;
             case 3:
               self->setPartOfValues3(d1,&vt1[0],&vt1[0]+vt1.size(),pc1.first,pc1.second.first,pc1.second.second);
               return self;
             default:
               throw INTERP_KERNEL::Exception(msg);
             }
           break;
         }
       case 15:
         {
           switch(sw1)
             {
             case 1:
               self->setPartOfValuesSimple1(i1,pt1.first,pt1.second.first,pt1.second.second,pc1.first,pc1.second.first,pc1.second.second);
               return self;
             case 2:
                 tmp=DataArrayInt::New();
                 tmp->useArray(&v1[0],false,CPP_DEALLOC,v1.size(),1);
                 self->setPartOfValues1(tmp,pt1.first,pt1.second.first,pt1.second.second,pc1.first,pc1.second.first,pc1.second.second,false);
                 return self;
             case 3:
               self->setPartOfValues1(d1,pt1.first,pt1.second.first,pt1.second.second,pc1.first,pc1.second.first,pc1.second.second);
               return self;
             default:
               throw INTERP_KERNEL::Exception(msg);
             }
           break;
         }
       case 16:
         {
           switch(sw1)
             {
             case 1:
               self->setPartOfValuesSimple3(i1,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),pc1.first,pc1.second.first,pc1.second.second);
               return self;
             case 2:
               tmp=DataArrayInt::New();
               tmp->useArray(&v1[0],false,CPP_DEALLOC,v1.size(),1);
               self->setPartOfValues3(tmp,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),pc1.first,pc1.second.first,pc1.second.second,false);
               return self;
             case 3:
               self->setPartOfValues3(d1,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),pc1.first,pc1.second.first,pc1.second.second);
               return self;
             default:
               throw INTERP_KERNEL::Exception(msg);
             }
           break;
         }
       default:
         throw INTERP_KERNEL::Exception(msg);
       }
     return self;
   }

   DataArrayInt *__add__(PyObject *obj) throw(INTERP_KERNEL::Exception)
   {
     const char msg[]="Unexpected situation in __add__ !";
     int val;
     DataArrayInt *a;
     std::vector<int> aa;
     int sw;
     convertObjToPossibleCpp1(obj,sw,val,aa,a);
     switch(sw)
       {
       case 1:
         {
           MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=self->deepCpy();
           ret->applyLin(1,val);
           ret->incrRef();
           return ret;
         }
       case 3:
         {
           return DataArrayInt::Add(self,a);
         }
       default:
         throw INTERP_KERNEL::Exception(msg);
       }
   }

   DataArrayInt *__radd__(PyObject *obj) throw(INTERP_KERNEL::Exception)
   {
     const char msg[]="Unexpected situation in __radd__ !";
     int val;
     DataArrayInt *a;
     std::vector<int> aa;
     int sw;
     convertObjToPossibleCpp1(obj,sw,val,aa,a);
     switch(sw)
       {
       case 1:
         {
           MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=self->deepCpy();
           ret->applyLin(1,val);
           ret->incrRef();
           return ret;
         }
       default:
         throw INTERP_KERNEL::Exception(msg);
       }
   }

   DataArrayInt *operator+=(PyObject *obj) throw(INTERP_KERNEL::Exception)
   {
     const char msg[]="Unexpected situation in __iadd__ !";
     int val;
     DataArrayInt *a;
     std::vector<int> aa;
     int sw;
     convertObjToPossibleCpp1(obj,sw,val,aa,a);
     switch(sw)
       {
       case 1:
         {
           self->applyLin(1,val);
           return self;
         }
       case 3:
         {
           self->addEqual(a);
           return self;
         }
       default:
         throw INTERP_KERNEL::Exception(msg);
       }
   }

   DataArrayInt *__sub__(PyObject *obj) throw(INTERP_KERNEL::Exception)
   {
     const char msg[]="Unexpected situation in __sub__ !";
     int val;
     DataArrayInt *a;
     std::vector<int> aa;
     int sw;
     convertObjToPossibleCpp1(obj,sw,val,aa,a);
     switch(sw)
       {
       case 1:
         {
           MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=self->deepCpy();
           ret->applyLin(1,-val);
           ret->incrRef();
           return ret;
         }
       case 3:
         {
           return DataArrayInt::Substract(self,a);
         }
       default:
         throw INTERP_KERNEL::Exception(msg);
       }
   }

   DataArrayInt *__rsub__(PyObject *obj) throw(INTERP_KERNEL::Exception)
   {
     const char msg[]="Unexpected situation in __rsub__ !";
     int val;
     DataArrayInt *a;
     std::vector<int> aa;
     int sw;
     convertObjToPossibleCpp1(obj,sw,val,aa,a);
     switch(sw)
       {
       case 1:
         {
           MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=self->deepCpy();
           ret->applyLin(-1,val);
           ret->incrRef();
           return ret;
         }
       default:
         throw INTERP_KERNEL::Exception(msg);
       }
   }

   DataArrayInt *operator-=(PyObject *obj) throw(INTERP_KERNEL::Exception)
   {
     const char msg[]="Unexpected situation in __isub__ !";
     int val;
     DataArrayInt *a;
     std::vector<int> aa;
     int sw;
     convertObjToPossibleCpp1(obj,sw,val,aa,a);
     switch(sw)
       {
       case 1:
         {
           self->applyLin(1.,-val);
           return self;
         }
       case 3:
         {
           self->substractEqual(a);
           return self;
         }
       default:
         throw INTERP_KERNEL::Exception(msg);
       }
   }

   DataArrayInt *__mul__(PyObject *obj) throw(INTERP_KERNEL::Exception)
   {
     const char msg[]="Unexpected situation in __mul__ !";
     int val;
     DataArrayInt *a;
     std::vector<int> aa;
     int sw;
     convertObjToPossibleCpp1(obj,sw,val,aa,a);
     switch(sw)
       {
       case 1:
         {
           MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=self->deepCpy();
           ret->applyLin(val,0);
           ret->incrRef();
           return ret;
         }
       case 3:
         {
           return DataArrayInt::Multiply(self,a);
         }
       default:
         throw INTERP_KERNEL::Exception(msg);
       }
   }

   DataArrayInt *__rmul__(PyObject *obj) throw(INTERP_KERNEL::Exception)
   {
     const char msg[]="Unexpected situation in __rmul__ !";
     int val;
     DataArrayInt *a;
     std::vector<int> aa;
     int sw;
     convertObjToPossibleCpp1(obj,sw,val,aa,a);
     switch(sw)
       {
       case 1:
         {
           MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=self->deepCpy();
           ret->applyLin(val,0);
           ret->incrRef();
           return ret;
         }
       default:
         throw INTERP_KERNEL::Exception(msg);
       }
   }

   DataArrayInt *operator*=(PyObject *obj) throw(INTERP_KERNEL::Exception)
   {
     const char msg[]="Unexpected situation in __imul__ !";
     int val;
     DataArrayInt *a;
     std::vector<int> aa;
     int sw;
     convertObjToPossibleCpp1(obj,sw,val,aa,a);
     switch(sw)
       {
       case 1:
         {
           self->applyLin(val,0);
           return self;
         }
       case 3:
         {
           self->multiplyEqual(a);
           return self;
         }
       default:
         throw INTERP_KERNEL::Exception(msg);
       }
   }

   DataArrayInt *__div__(PyObject *obj) throw(INTERP_KERNEL::Exception)
   {
     const char msg[]="Unexpected situation in __div__ !";
     int val;
     DataArrayInt *a;
     std::vector<int> aa;
     int sw;
     convertObjToPossibleCpp1(obj,sw,val,aa,a);
     switch(sw)
       {
       case 1:
         {
           MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=self->deepCpy();
           ret->applyDivideBy(val);
           ret->incrRef();
           return ret;
         }
       case 3:
         {
           return DataArrayInt::Divide(self,a);
         }
       default:
         throw INTERP_KERNEL::Exception(msg);
       }
   }

   DataArrayInt *__rdiv__(PyObject *obj) throw(INTERP_KERNEL::Exception)
   {
     const char msg[]="Unexpected situation in __rdiv__ !";
     int val;
     DataArrayInt *a;
     std::vector<int> aa;
     int sw;
     convertObjToPossibleCpp1(obj,sw,val,aa,a);
     switch(sw)
       {
       case 1:
         {
           MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=self->deepCpy();
           ret->applyInv(val);
           ret->incrRef();
           return ret;
         }
       default:
         throw INTERP_KERNEL::Exception(msg);
       }
   }

   DataArrayInt *operator/=(PyObject *obj) throw(INTERP_KERNEL::Exception)
   {
     const char msg[]="Unexpected situation in __imul__ !";
     int val;
     DataArrayInt *a;
     std::vector<int> aa;
     int sw;
     convertObjToPossibleCpp1(obj,sw,val,aa,a);
     switch(sw)
       {
       case 1:
         {
           self->applyDivideBy(val);
           return self;
         }
       case 3:
         {
           self->divideEqual(a);
           return self;
         }
       default:
         throw INTERP_KERNEL::Exception(msg);
       }
   }

   DataArrayInt *__mod__(PyObject *obj) throw(INTERP_KERNEL::Exception)
   {
     const char msg[]="Unexpected situation in __mod__ !";
     int val;
     DataArrayInt *a;
     std::vector<int> aa;
     int sw;
     convertObjToPossibleCpp1(obj,sw,val,aa,a);
     switch(sw)
       {
       case 1:
         {
           MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=self->deepCpy();
           ret->applyModulus(val);
           ret->incrRef();
           return ret;
         }
       case 3:
         {
           return DataArrayInt::Modulus(self,a);
         }
       default:
         throw INTERP_KERNEL::Exception(msg);
       }
   }

   DataArrayInt *__rmod__(PyObject *obj) throw(INTERP_KERNEL::Exception)
   {
     const char msg[]="Unexpected situation in __rmod__ !";
     int val;
     DataArrayInt *a;
     std::vector<int> aa;
     int sw;
     convertObjToPossibleCpp1(obj,sw,val,aa,a);
     switch(sw)
       {
       case 1:
         {
           MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=self->deepCpy();
           ret->applyRModulus(val);
           ret->incrRef();
           return ret;
         }
       case 3:
         {
           return DataArrayInt::Modulus(self,a);
         }
       default:
         throw INTERP_KERNEL::Exception(msg);
       }
   }

   DataArrayInt *operator%=(PyObject *obj) throw(INTERP_KERNEL::Exception)
   {
     const char msg[]="Unexpected situation in __imod__ !";
     int val;
     DataArrayInt *a;
     std::vector<int> aa;
     int sw;
     convertObjToPossibleCpp1(obj,sw,val,aa,a);
     switch(sw)
       {
       case 1:
         {
           self->applyModulus(val);
           return self;
         }
       case 3:
         {
           self->modulusEqual(a);
           return self;
         }
       default:
         throw INTERP_KERNEL::Exception(msg);
       }
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
    NatureOfField getNature() const throw(INTERP_KERNEL::Exception);
    virtual void setNature(NatureOfField nat) throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getLocalizationOfDiscr() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *buildMeasureField(bool isAbs) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDiscretization *getDiscretization() const throw(INTERP_KERNEL::Exception);
    int getNumberOfTuplesExpected() const throw(INTERP_KERNEL::Exception);
    int getNumberOfMeshPlacesExpected() const throw(INTERP_KERNEL::Exception);
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
        if(ret1)
          ret1->incrRef();
        return convertMesh(ret1, SWIG_POINTER_OWN | 0 );
      }

      PyObject *buildSubMeshData(PyObject *li) const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *ret1=0;
        MEDCouplingMesh *ret0=0;
        void *da=0;
        int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 |  0 );
        if (!SWIG_IsOK(res1))
          {
            int size;
            INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
            ret0=self->buildSubMeshData(tmp,tmp+size,ret1);
          }
        else
          {
            DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
            if(!da2)
              throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
            da2->checkAllocated();
            ret0=self->buildSubMeshData(da2->getConstPointer(),da2->getConstPointer()+da2->getNbOfElems(),ret1);
          }
        PyObject *res = PyList_New(2);
        PyList_SetItem(res,0,convertMesh(ret0, SWIG_POINTER_OWN | 0 ));
        PyList_SetItem(res,1,SWIG_NewPointerObj((void*)ret1,SWIGTYPE_p_ParaMEDMEM__DataArrayInt,SWIG_POINTER_OWN | 0));
        return res;
      }

      void setGaussLocalizationOnCells(PyObject *li, const std::vector<double>& refCoo,
                                       const std::vector<double>& gsCoo, const std::vector<double>& wg) throw(INTERP_KERNEL::Exception)
      {
        void *da=0;
        int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 |  0 );
        if (!SWIG_IsOK(res1))
          {
            int size;
            INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
            self->setGaussLocalizationOnCells(tmp,((int *)tmp)+size,refCoo,gsCoo,wg);
          }
        else
          {
            DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
            if(!da2)
              throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
            da2->checkAllocated();
            self->setGaussLocalizationOnCells(da2->getConstPointer(),da2->getConstPointer()+da2->getNbOfElems(),refCoo,gsCoo,wg);
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
    static MEDCouplingFieldDouble *New(const MEDCouplingFieldTemplate *ft, TypeOfTimeDiscretization td=NO_TIME);
    void setTimeUnit(const char *unit);
    const char *getTimeUnit() const;
    void copyTinyStringsFrom(const MEDCouplingFieldDouble *other) throw(INTERP_KERNEL::Exception);
    void copyTinyAttrFrom(const MEDCouplingFieldDouble *other) throw(INTERP_KERNEL::Exception);
    std::string simpleRepr() const;
    std::string advancedRepr() const;
    MEDCouplingFieldDouble *clone(bool recDeepCpy) const;
    MEDCouplingFieldDouble *cloneWithMesh(bool recDeepCpy) const;
    MEDCouplingFieldDouble *deepCpy() const;
    MEDCouplingFieldDouble *buildNewTimeReprFromThis(TypeOfTimeDiscretization td, bool deepCpy) const throw(INTERP_KERNEL::Exception);
    TypeOfTimeDiscretization getTimeDiscretization() const throw(INTERP_KERNEL::Exception);
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
    void setTimeTolerance(double val) throw(INTERP_KERNEL::Exception);
    double getTimeTolerance() const throw(INTERP_KERNEL::Exception);
    void setIteration(int it) throw(INTERP_KERNEL::Exception);
    void setEndIteration(int it) throw(INTERP_KERNEL::Exception);
    void setOrder(int order) throw(INTERP_KERNEL::Exception);
    void setEndOrder(int order) throw(INTERP_KERNEL::Exception);
    void setTimeValue(double val) throw(INTERP_KERNEL::Exception);
    void setEndTimeValue(double val) throw(INTERP_KERNEL::Exception);
    void updateTime() const throw(INTERP_KERNEL::Exception);
    void changeUnderlyingMesh(const MEDCouplingMesh *other, int levOfCheck, double prec) throw(INTERP_KERNEL::Exception);
    void substractInPlaceDM(const MEDCouplingFieldDouble *f, int levOfCheck, double prec) throw(INTERP_KERNEL::Exception);
    bool mergeNodes(double eps, double epsOnVals=1e-15) throw(INTERP_KERNEL::Exception);
    bool mergeNodes2(double eps, double epsOnVals=1e-15) throw(INTERP_KERNEL::Exception);
    bool zipCoords(double epsOnVals=1e-15) throw(INTERP_KERNEL::Exception);
    bool zipConnectivity(int compType,double epsOnVals=1e-15) throw(INTERP_KERNEL::Exception);
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
    void applyFunc2(int nbOfComp, const char *func) throw(INTERP_KERNEL::Exception);
    void applyFunc3(int nbOfComp, const std::vector<std::string>& varsOrder, const char *func) throw(INTERP_KERNEL::Exception);
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
    static MEDCouplingFieldDouble *MergeFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception);
    static MEDCouplingFieldDouble *MeldFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception);
    static MEDCouplingFieldDouble *DotFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *dot(const MEDCouplingFieldDouble& other) const throw(INTERP_KERNEL::Exception);
    static MEDCouplingFieldDouble *CrossProductFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *crossProduct(const MEDCouplingFieldDouble& other) const throw(INTERP_KERNEL::Exception);
    static MEDCouplingFieldDouble *MaxFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *max(const MEDCouplingFieldDouble& other) const throw(INTERP_KERNEL::Exception);
    static MEDCouplingFieldDouble *MinFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception);
    static MEDCouplingFieldDouble *AddFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception);
    static MEDCouplingFieldDouble *SubstractFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception);
    static MEDCouplingFieldDouble *MultiplyFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception);
    static MEDCouplingFieldDouble *DivideFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2) throw(INTERP_KERNEL::Exception);
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

      DataArrayDouble *getArray() throw(INTERP_KERNEL::Exception)
      {
        DataArrayDouble *ret=self->getArray();
        if(ret)
          ret->incrRef();
        return ret;
      }

      PyObject *getArrays() const throw(INTERP_KERNEL::Exception)
      {
        std::vector<DataArrayDouble *> arrs=self->getArrays();
        for(std::vector<DataArrayDouble *>::iterator it=arrs.begin();it!=arrs.end();it++)
          if(*it)
            (*it)->incrRef();
        int sz=arrs.size();
        PyObject *ret=PyTuple_New(sz);
        for(int i=0;i<sz;i++)
          {
            if(arrs[i])
              PyTuple_SetItem(ret,i,SWIG_NewPointerObj(SWIG_as_voidptr(arrs[i]),SWIGTYPE_p_ParaMEDMEM__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
            else
              PyTuple_SetItem(ret,i,SWIG_NewPointerObj(SWIG_as_voidptr(0),SWIGTYPE_p_ParaMEDMEM__DataArrayDouble, 0 | 0 ));
          }
        return ret;
      }

      void setArrays(PyObject *ls) throw(INTERP_KERNEL::Exception)
      {
        std::vector<const DataArrayDouble *> tmp;
        convertPyObjToVecDataArrayDblCst(ls,tmp);
        int sz=tmp.size();
        std::vector<DataArrayDouble *> arrs(sz);
        for(int i=0;i<sz;i++)
          arrs[i]=const_cast<DataArrayDouble *>(tmp[i]);
        self->setArrays(arrs);
      }

      DataArrayDouble *getEndArray() throw(INTERP_KERNEL::Exception)
      {
        DataArrayDouble *ret=self->getEndArray();
        if(ret)
          ret->incrRef();
        return ret;
      }

      PyObject *getValueOn(PyObject *sl) const throw(INTERP_KERNEL::Exception)
      {
        int sz;
        INTERP_KERNEL::AutoPtr<double> spaceLoc=convertPyToNewDblArr2(sl,&sz);
        sz=self->getNumberOfComponents();
        INTERP_KERNEL::AutoPtr<double> res=new double[sz];
        self->getValueOn(spaceLoc,res);
        return convertDblArrToPyList(res,sz);
      }

      DataArrayDouble *getValueOnMulti(PyObject *li) const throw(INTERP_KERNEL::Exception)
      {
        void *da=0;
        int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayDouble, 0 |  0 );
        if (!SWIG_IsOK(res1))
          {
            int size;
            INTERP_KERNEL::AutoPtr<double> tmp=convertPyToNewDblArr2(li,&size);
            int spaceDim=self->getMesh()->getSpaceDimension();
            int nbOfPoints=size/spaceDim;
            if(size%spaceDim!=0)
              {
                throw INTERP_KERNEL::Exception("Invalid list length ! Must be a multiple of self.getMesh().getSpaceDimension() !");
              }
            return self->getValueOnMulti(tmp,nbOfPoints);
          }
        else
          {
            DataArrayDouble *da2=reinterpret_cast< DataArrayDouble * >(da);
            if(!da2)
              throw INTERP_KERNEL::Exception("Not null DataArrayDouble instance expected !");
            da2->checkAllocated();
            int size=da2->getNumberOfTuples();
            int nbOfCompo=da2->getNumberOfComponents();
            if(nbOfCompo!=self->getMesh()->getSpaceDimension())
              {
                throw INTERP_KERNEL::Exception("Invalid DataArrayDouble nb of components ! Expected same as self.getMesh().getSpaceDimension() !");
              }
            return self->getValueOnMulti(da2->getConstPointer(),size);
          }
      }

      PyObject *getValueOn(PyObject *sl, double time) const throw(INTERP_KERNEL::Exception)
      {
        int sz;
        INTERP_KERNEL::AutoPtr<double> spaceLoc=convertPyToNewDblArr2(sl,&sz);
        sz=self->getNumberOfComponents();
        INTERP_KERNEL::AutoPtr<double> res=new double[sz];
        self->getValueOn(spaceLoc,time,res);
        return convertDblArrToPyList(res,sz);
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
        INTERP_KERNEL::AutoPtr<double> tmp=new double[sz];
        self->accumulate(tmp);
        PyObject *ret=convertDblArrToPyList(tmp,sz);
        return ret;
      }
      PyObject *integral(bool isWAbs) const throw(INTERP_KERNEL::Exception)
      {
        int sz=self->getNumberOfComponents();
        INTERP_KERNEL::AutoPtr<double> tmp=new double[sz];
        self->integral(isWAbs,tmp);
        PyObject *ret=convertDblArrToPyList(tmp,sz);
        return ret;
      }
      PyObject *normL1() const throw(INTERP_KERNEL::Exception)
      {
        int sz=self->getNumberOfComponents();
        INTERP_KERNEL::AutoPtr<double> tmp=new double[sz];
        self->normL1(tmp);
        PyObject *ret=convertDblArrToPyList(tmp,sz);
        return ret;
      }
      PyObject *normL2() const throw(INTERP_KERNEL::Exception)
      {
        int sz=self->getNumberOfComponents();
        INTERP_KERNEL::AutoPtr<double> tmp=new double[sz];
        self->normL2(tmp);
        PyObject *ret=convertDblArrToPyList(tmp,sz);
        return ret;
      }
      void renumberCells(PyObject *li, bool check) throw(INTERP_KERNEL::Exception)
      {
        void *da=0;
        int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 |  0 );
        if (!SWIG_IsOK(res1))
          {
            int size;
            INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
            self->renumberCells(tmp,check);
          }
        else
          {
            DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
            if(!da2)
              throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
            da2->checkAllocated();
            self->renumberCells(da2->getConstPointer(),check);
          }
      }
      void renumberNodes(PyObject *li) throw(INTERP_KERNEL::Exception)
      {
        void *da=0;
        int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 |  0 );
        if (!SWIG_IsOK(res1))
          {
            int size;
            INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
            self->renumberNodes(tmp);
          }
        else
          {
            DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
            if(!da2)
              throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
            da2->checkAllocated();
            self->renumberNodes(da2->getConstPointer());
          }
      }

      MEDCouplingFieldDouble *buildSubPart(PyObject *li) const throw(INTERP_KERNEL::Exception)
      {
        void *da=0;
        int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 |  0 );
        if (!SWIG_IsOK(res1))
          {
            int size;
            INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
            return self->buildSubPart(tmp,((const int *)tmp)+size);
          }
        else
          {
            DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
            if(!da2)
              throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
            da2->checkAllocated();
            return self->buildSubPart(da2->getConstPointer(),da2->getConstPointer()+da2->getNbOfElems());
          }
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

      static MEDCouplingFieldDouble *MergeFields(PyObject *li) throw(INTERP_KERNEL::Exception)
      {
        std::vector<const MEDCouplingFieldDouble *> tmp;
        convertPyObjToVecFieldDblCst(li,tmp);
        return MEDCouplingFieldDouble::MergeFields(tmp);
      }
    }
  };

  class MEDCouplingFieldTemplate : public ParaMEDMEM::MEDCouplingField
  {
  public:
    static MEDCouplingFieldTemplate *New(const MEDCouplingFieldDouble *f) throw(INTERP_KERNEL::Exception);
    static MEDCouplingFieldTemplate *New(TypeOfField type);
    std::string simpleRepr() const;
    std::string advancedRepr() const;
    %extend
       {
         std::string __str__() const
           {
             return self->simpleRepr();
           }
       }
  };

  class MEDCouplingMultiFields : public RefCountObject, public TimeLabel
  {
  public:
    int getNumberOfFields() const;
    MEDCouplingMultiFields *deepCpy() const;
    virtual std::string simpleRepr() const;
    virtual std::string advancedRepr() const;
    virtual bool isEqual(const MEDCouplingMultiFields *other, double meshPrec, double valsPrec) const;
    virtual bool isEqualWithoutConsideringStr(const MEDCouplingMultiFields *other, double meshPrec, double valsPrec) const;
    virtual void checkCoherency() const throw(INTERP_KERNEL::Exception);
    %extend
       {
         std::string __str__() const
         {
           return self->simpleRepr();
         }
         static MEDCouplingMultiFields *New(PyObject *li) throw(INTERP_KERNEL::Exception)
         {
           std::vector<const ParaMEDMEM::MEDCouplingFieldDouble *> tmp;
           convertPyObjToVecFieldDblCst(li,tmp);
           int sz=tmp.size();
           std::vector<MEDCouplingFieldDouble *> fs(sz);
           for(int i=0;i<sz;i++)
             fs[i]=const_cast<MEDCouplingFieldDouble *>(tmp[i]);
           return MEDCouplingMultiFields::New(fs);
         }
         PyObject *getFields() const
         {
           std::vector<const MEDCouplingFieldDouble *> fields=self->getFields();
           int sz=fields.size();
           PyObject *res = PyList_New(sz);
           for(int i=0;i<sz;i++)
             {
               if(fields[i])
                 {
                   fields[i]->incrRef();
                   PyList_SetItem(res,i,SWIG_NewPointerObj(SWIG_as_voidptr(fields[i]),SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble, SWIG_POINTER_OWN | 0 ));
                 }
               else
                 {
                   PyList_SetItem(res,i,SWIG_NewPointerObj(SWIG_as_voidptr(0),SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble, 0 ));
                 }
             }
           return res;
         }
         PyObject *getFieldAtPos(int id) const throw(INTERP_KERNEL::Exception)
         {
           const MEDCouplingFieldDouble *ret=self->getFieldAtPos(id);
           if(ret)
             {
               ret->incrRef();
               return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble, SWIG_POINTER_OWN | 0 );
             }
           else
             return SWIG_NewPointerObj(SWIG_as_voidptr(0),SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble, 0 );
         }
         PyObject *getMeshes() const throw(INTERP_KERNEL::Exception)
         {
           std::vector<MEDCouplingMesh *> ms=self->getMeshes();
           int sz=ms.size();
           PyObject *res = PyList_New(sz);
           for(int i=0;i<sz;i++)
             {
               if(ms[i])
                 {
                   ms[i]->incrRef();
                   PyList_SetItem(res,i,convertMesh(ms[i], SWIG_POINTER_OWN | 0 ));
                 }
               else
                 {
                   PyList_SetItem(res,i,SWIG_NewPointerObj(SWIG_as_voidptr(0),SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh, 0 ));
                 }
             }
           return res;
         }
         PyObject *getDifferentMeshes() const throw(INTERP_KERNEL::Exception)
         {
           std::vector<int> refs;
           std::vector<MEDCouplingMesh *> ms=self->getDifferentMeshes(refs);
           int sz=ms.size();
           PyObject *res = PyList_New(sz);
           for(int i=0;i<sz;i++)
             {
               if(ms[i])
                 {
                   ms[i]->incrRef();
                   PyList_SetItem(res,i,convertMesh(ms[i], SWIG_POINTER_OWN | 0 ));
                 }
               else
                 {
                   PyList_SetItem(res,i,SWIG_NewPointerObj(SWIG_as_voidptr(0),SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh, 0 ));
                 }
             }
           //
           PyObject *ret=PyTuple_New(2);
           PyTuple_SetItem(ret,0,res);
           PyTuple_SetItem(ret,1,convertIntArrToPyList2(refs));
           return ret;
         }
         PyObject *getArrays() const throw(INTERP_KERNEL::Exception)
         {
           std::vector<DataArrayDouble *> ms=self->getArrays();
           int sz=ms.size();
           PyObject *res = PyList_New(sz);
           for(int i=0;i<sz;i++)
             {
               if(ms[i])
                 {
                   ms[i]->incrRef();
                   PyList_SetItem(res,i,SWIG_NewPointerObj(SWIG_as_voidptr(ms[i]),SWIGTYPE_p_ParaMEDMEM__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
                 }
               else
                 {
                   PyList_SetItem(res,i,SWIG_NewPointerObj(SWIG_as_voidptr(0),SWIGTYPE_p_ParaMEDMEM__DataArrayDouble, 0 ));
                 }
             }
           return res;
         }
         PyObject *getDifferentArrays() const throw(INTERP_KERNEL::Exception)
         {
           std::vector< std::vector<int> > refs;
           std::vector<DataArrayDouble *> ms=self->getDifferentArrays(refs);
           int sz=ms.size();
           PyObject *res = PyList_New(sz);
           PyObject *res2 = PyList_New(sz);
           for(int i=0;i<sz;i++)
             {
               if(ms[i])
                 {
                   ms[i]->incrRef();
                   PyList_SetItem(res,i,SWIG_NewPointerObj(SWIG_as_voidptr(ms[i]),SWIGTYPE_p_ParaMEDMEM__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
                 }
               else
                 {
                   PyList_SetItem(res,i,SWIG_NewPointerObj(SWIG_as_voidptr(0),SWIGTYPE_p_ParaMEDMEM__DataArrayDouble, 0 ));
                 }
               PyList_SetItem(res2,i,convertIntArrToPyList2(refs[i]));
             }
           //
           PyObject *ret=PyTuple_New(2);
           PyTuple_SetItem(ret,0,res);
           PyTuple_SetItem(ret,1,res2);
           return ret;
         }
       }
  };
  
  class MEDCouplingDefinitionTime
  {
  public:
    MEDCouplingDefinitionTime();
    void assign(const MEDCouplingDefinitionTime& other);
    bool isEqual(const MEDCouplingDefinitionTime& other) const;
    double getTimeResolution() const;
    std::vector<double> getHotSpotsTime() const;
    %extend
      {
        std::string __str__() const
          {
            std::ostringstream oss;
            self->appendRepr(oss);
            return oss.str();
          }

        PyObject *getIdsOnTimeRight(double tm) const throw(INTERP_KERNEL::Exception)
        {
          int meshId,arrId,arrIdInField,fieldId;
          self->getIdsOnTimeRight(tm,meshId,arrId,arrIdInField,fieldId);
          PyObject *res=PyList_New(4);
          PyList_SetItem(res,0,PyInt_FromLong(meshId));
          PyList_SetItem(res,1,PyInt_FromLong(arrId));
          PyList_SetItem(res,2,PyInt_FromLong(arrIdInField));
          PyList_SetItem(res,3,PyInt_FromLong(fieldId));
          return res;
        }

        PyObject *getIdsOnTimeLeft(double tm) const throw(INTERP_KERNEL::Exception)
        {
          int meshId,arrId,arrIdInField,fieldId;
          self->getIdsOnTimeLeft(tm,meshId,arrId,arrIdInField,fieldId);
          PyObject *res=PyList_New(4);
          PyList_SetItem(res,0,PyInt_FromLong(meshId));
          PyList_SetItem(res,1,PyInt_FromLong(arrId));
          PyList_SetItem(res,2,PyInt_FromLong(arrIdInField));
          PyList_SetItem(res,3,PyInt_FromLong(fieldId));
          return res;
        }
      }
  };

  class MEDCouplingFieldOverTime : public MEDCouplingMultiFields
  {
  public:
    double getTimeTolerance() const throw(INTERP_KERNEL::Exception);
    MEDCouplingDefinitionTime getDefinitionTimeZone() const;
    
    %extend
      {
        std::string __str__() const
          {
            return self->simpleRepr();
          }
        static MEDCouplingFieldOverTime *New(PyObject *li) throw(INTERP_KERNEL::Exception)
        {
          std::vector<const ParaMEDMEM::MEDCouplingFieldDouble *> tmp;
          convertPyObjToVecFieldDblCst(li,tmp);
           int sz=tmp.size();
           std::vector<MEDCouplingFieldDouble *> fs(sz);
           for(int i=0;i<sz;i++)
             fs[i]=const_cast<MEDCouplingFieldDouble *>(tmp[i]);
           return MEDCouplingFieldOverTime::New(fs);
         }
      }
  };
}
