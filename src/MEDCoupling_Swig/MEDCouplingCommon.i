// Copyright (C) 2007-2013  CEA/DEN, EDF R&D
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
#include "MEDCouplingCurveLinearMesh.hxx"
#include "MEDCouplingField.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingFieldTemplate.hxx"
#include "MEDCouplingGaussLocalization.hxx"
#include "MEDCouplingAutoRefCountObjectPtr.hxx"
#include "MEDCouplingMultiFields.hxx"
#include "MEDCouplingFieldOverTime.hxx"
#include "MEDCouplingDefinitionTime.hxx"
#include "MEDCouplingFieldDiscretization.hxx"
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

%typemap(out) ParaMEDMEM::MEDCouplingStructuredMesh*
{
  $result=convertMesh($1,$owner);
}

%typemap(out) ParaMEDMEM::MEDCouplingFieldDiscretization*
{
  $result=convertFieldDiscretization($1,$owner);
}

%typemap(out) ParaMEDMEM::MEDCouplingMultiFields*
{
  $result=convertMultiFields($1,$owner);
}

%typemap(out) ParaMEDMEM::DataArray*
{
  $result=convertDataArray($1,$owner);
}

%typemap(out) ParaMEDMEM::DataArrayChar*
{
  $result=convertDataArrayChar($1,$owner);
}

#ifdef WITH_NUMPY
%init %{ import_array(); %}
#endif

%feature("autodoc", "1");
%feature("docstring");

%newobject ParaMEDMEM::MEDCouplingFieldDiscretization::New;
%newobject ParaMEDMEM::MEDCouplingFieldDiscretization::getOffsetArr;
%newobject ParaMEDMEM::MEDCouplingFieldDiscretization::deepCpy;
%newobject ParaMEDMEM::MEDCouplingFieldDiscretization::clone;
%newobject ParaMEDMEM::MEDCouplingFieldDiscretization::clonePart;
%newobject ParaMEDMEM::MEDCouplingFieldDiscretization::clonePartRange;
%newobject ParaMEDMEM::MEDCouplingFieldDiscretization::getMeasureField;
%newobject ParaMEDMEM::MEDCouplingFieldDiscretization::getOffsetArr;
%newobject ParaMEDMEM::MEDCouplingFieldDiscretization::getLocalizationOfDiscValues;
%newobject ParaMEDMEM::MEDCouplingFieldDiscretization::getValueOnMulti;
%newobject ParaMEDMEM::MEDCouplingFieldDiscretization::computeTupleIdsToSelectFromCellIds;
%newobject ParaMEDMEM::MEDCouplingFieldDiscretization::buildSubMeshData;
%newobject ParaMEDMEM::MEDCouplingField::buildMeasureField;
%newobject ParaMEDMEM::MEDCouplingField::getLocalizationOfDiscr;
%newobject ParaMEDMEM::MEDCouplingField::computeTupleIdsToSelectFromCellIds;
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
%newobject ParaMEDMEM::MEDCouplingFieldDouble::extractSlice3D;
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
%newobject ParaMEDMEM::MEDCouplingFieldDouble::negate;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::getIdsInRange;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::buildSubPart;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::buildSubPartRange;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::__getitem__;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::__neg__;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::__add__;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::__sub__;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::__mul__;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::__div__;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::__pow__;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::__radd__;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::__rsub__;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::__rmul__;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::__rdiv__;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::clone;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::cloneWithMesh;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::deepCpy;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::buildNewTimeReprFromThis;
%newobject ParaMEDMEM::MEDCouplingFieldDouble::getValueOnMulti;
%newobject ParaMEDMEM::MEDCouplingFieldTemplate::New;
%newobject ParaMEDMEM::DataArray::selectByTupleRanges;
%newobject ParaMEDMEM::DataArrayInt::New;
%newobject ParaMEDMEM::DataArrayInt::__iter__;
%newobject ParaMEDMEM::DataArrayInt::convertToDblArr;
%newobject ParaMEDMEM::DataArrayInt::deepCpy;
%newobject ParaMEDMEM::DataArrayInt::performCpy;
%newobject ParaMEDMEM::DataArrayInt::substr;
%newobject ParaMEDMEM::DataArrayInt::changeNbOfComponents;
%newobject ParaMEDMEM::DataArrayInt::accumulatePerChunck;
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
%newobject ParaMEDMEM::DataArrayInt::invertArrayO2N2N2OBis;
%newobject ParaMEDMEM::DataArrayInt::getIdsEqual;
%newobject ParaMEDMEM::DataArrayInt::getIdsNotEqual;
%newobject ParaMEDMEM::DataArrayInt::getIdsEqualList;
%newobject ParaMEDMEM::DataArrayInt::getIdsNotEqualList;
%newobject ParaMEDMEM::DataArrayInt::negate;
%newobject ParaMEDMEM::DataArrayInt::getIdsInRange;
%newobject ParaMEDMEM::DataArrayInt::Aggregate;
%newobject ParaMEDMEM::DataArrayInt::Meld;
%newobject ParaMEDMEM::DataArrayInt::Add;
%newobject ParaMEDMEM::DataArrayInt::Substract;
%newobject ParaMEDMEM::DataArrayInt::Multiply;
%newobject ParaMEDMEM::DataArrayInt::Divide;
%newobject ParaMEDMEM::DataArrayInt::Pow;
%newobject ParaMEDMEM::DataArrayInt::BuildUnion;
%newobject ParaMEDMEM::DataArrayInt::BuildIntersection;
%newobject ParaMEDMEM::DataArrayInt::Range;
%newobject ParaMEDMEM::DataArrayInt::fromNoInterlace;
%newobject ParaMEDMEM::DataArrayInt::toNoInterlace;
%newobject ParaMEDMEM::DataArrayInt::buildComplement;
%newobject ParaMEDMEM::DataArrayInt::buildUnion;
%newobject ParaMEDMEM::DataArrayInt::buildSubstraction;
%newobject ParaMEDMEM::DataArrayInt::buildSubstractionOptimized;
%newobject ParaMEDMEM::DataArrayInt::buildIntersection;
%newobject ParaMEDMEM::DataArrayInt::buildUnique;
%newobject ParaMEDMEM::DataArrayInt::deltaShiftIndex;
%newobject ParaMEDMEM::DataArrayInt::buildExplicitArrByRanges;
%newobject ParaMEDMEM::DataArrayInt::findRangeIdForEachTuple;
%newobject ParaMEDMEM::DataArrayInt::findIdInRangeForEachTuple;
%newobject ParaMEDMEM::DataArrayInt::duplicateEachTupleNTimes;
%newobject ParaMEDMEM::DataArrayInt::buildPermutationArr;
%newobject ParaMEDMEM::DataArrayInt::buildPermArrPerLevel;
%newobject ParaMEDMEM::DataArrayInt::getDifferentValues;
%newobject ParaMEDMEM::DataArrayInt::__neg__;
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
%newobject ParaMEDMEM::DataArrayInt::__pow__;
%newobject ParaMEDMEM::DataArrayInt::__rpow__;
%newobject ParaMEDMEM::DataArrayIntTuple::buildDAInt;
%newobject ParaMEDMEM::DataArrayChar::deepCpy;
%newobject ParaMEDMEM::DataArrayChar::convertToIntArr;
%newobject ParaMEDMEM::DataArrayChar::renumber;
%newobject ParaMEDMEM::DataArrayChar::renumberR;
%newobject ParaMEDMEM::DataArrayChar::renumberAndReduce;
%newobject ParaMEDMEM::DataArrayChar::selectByTupleIdSafe;
%newobject ParaMEDMEM::DataArrayChar::selectByTupleId2;
%newobject ParaMEDMEM::DataArrayChar::changeNbOfComponents;
%newobject ParaMEDMEM::DataArrayChar::getIdsEqual;
%newobject ParaMEDMEM::DataArrayChar::getIdsNotEqual;
%newobject ParaMEDMEM::DataArrayChar::Aggregate;
%newobject ParaMEDMEM::DataArrayChar::Meld;
%newobject ParaMEDMEM::DataArrayByte::New;
%newobject ParaMEDMEM::DataArrayByte::__iter__;
%newobject ParaMEDMEM::DataArrayByte::performCpy;
%newobject ParaMEDMEM::DataArrayByteTuple::buildDAByte;
%newobject ParaMEDMEM::DataArrayChar::substr;
%newobject ParaMEDMEM::DataArrayAsciiChar::New;
%newobject ParaMEDMEM::DataArrayAsciiChar::__iter__;
%newobject ParaMEDMEM::DataArrayAsciiChar::performCpy;
%newobject ParaMEDMEM::DataArrayAsciiCharTuple::buildDAAsciiChar;
%newobject ParaMEDMEM::DataArrayDouble::New;
%newobject ParaMEDMEM::DataArrayDouble::__iter__;
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
%newobject ParaMEDMEM::DataArrayDouble::Pow;
%newobject ParaMEDMEM::DataArrayDouble::substr;
%newobject ParaMEDMEM::DataArrayDouble::changeNbOfComponents;
%newobject ParaMEDMEM::DataArrayDouble::accumulatePerChunck;
%newobject ParaMEDMEM::DataArrayDouble::getIdsInRange;
%newobject ParaMEDMEM::DataArrayDouble::selectByTupleId;
%newobject ParaMEDMEM::DataArrayDouble::selectByTupleIdSafe;
%newobject ParaMEDMEM::DataArrayDouble::selectByTupleId2;
%newobject ParaMEDMEM::DataArrayDouble::negate;
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
%newobject ParaMEDMEM::DataArrayDouble::computeBBoxPerTuple;
%newobject ParaMEDMEM::DataArrayDouble::buildEuclidianDistanceDenseMatrix;
%newobject ParaMEDMEM::DataArrayDouble::buildEuclidianDistanceDenseMatrixWith;
%newobject ParaMEDMEM::DataArrayDouble::renumber;
%newobject ParaMEDMEM::DataArrayDouble::renumberR;
%newobject ParaMEDMEM::DataArrayDouble::renumberAndReduce;
%newobject ParaMEDMEM::DataArrayDouble::fromNoInterlace;
%newobject ParaMEDMEM::DataArrayDouble::toNoInterlace;
%newobject ParaMEDMEM::DataArrayDouble::fromPolarToCart;
%newobject ParaMEDMEM::DataArrayDouble::fromCylToCart;
%newobject ParaMEDMEM::DataArrayDouble::fromSpherToCart;
%newobject ParaMEDMEM::DataArrayDouble::getDifferentValues;
%newobject ParaMEDMEM::DataArrayDouble::findClosestTupleId;
%newobject ParaMEDMEM::DataArrayDouble::duplicateEachTupleNTimes;
%newobject ParaMEDMEM::DataArrayDouble::__neg__;
%newobject ParaMEDMEM::DataArrayDouble::__radd__;
%newobject ParaMEDMEM::DataArrayDouble::__rsub__;
%newobject ParaMEDMEM::DataArrayDouble::__rmul__;
%newobject ParaMEDMEM::DataArrayDouble::__rdiv__;
%newobject ParaMEDMEM::DataArrayDouble::__pow__;
%newobject ParaMEDMEM::DataArrayDouble::__rpow__;
%newobject ParaMEDMEM::DataArrayDoubleTuple::buildDADouble;
%newobject ParaMEDMEM::MEDCouplingMesh::deepCpy;
%newobject ParaMEDMEM::MEDCouplingMesh::checkDeepEquivalOnSameNodesWith;
%newobject ParaMEDMEM::MEDCouplingMesh::checkTypeConsistencyAndContig;
%newobject ParaMEDMEM::MEDCouplingMesh::computeNbOfNodesPerCell;
%newobject ParaMEDMEM::MEDCouplingMesh::computeNbOfFacesPerCell;
%newobject ParaMEDMEM::MEDCouplingMesh::buildPartRange;
%newobject ParaMEDMEM::MEDCouplingMesh::giveCellsWithType;
%newobject ParaMEDMEM::MEDCouplingMesh::getCoordinatesAndOwner;
%newobject ParaMEDMEM::MEDCouplingMesh::getBarycenterAndOwner;
%newobject ParaMEDMEM::MEDCouplingMesh::computeIsoBarycenterOfNodesPerCell;
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
%newobject ParaMEDMEM::MEDCouplingPointSet::getCellsInBoundingBox;
%newobject ParaMEDMEM::MEDCouplingPointSet::findBoundaryNodes;
%newobject ParaMEDMEM::MEDCouplingPointSet::buildBoundaryMesh;
%newobject ParaMEDMEM::MEDCouplingPointSet::MergeNodesArray;
%newobject ParaMEDMEM::MEDCouplingPointSet::buildPartOfMySelf2;
%newobject ParaMEDMEM::MEDCouplingPointSet::BuildInstanceFromMeshType;
%newobject ParaMEDMEM::MEDCouplingUMesh::New;
%newobject ParaMEDMEM::MEDCouplingUMesh::getNodalConnectivity;
%newobject ParaMEDMEM::MEDCouplingUMesh::getNodalConnectivityIndex;
%newobject ParaMEDMEM::MEDCouplingUMesh::clone;
%newobject ParaMEDMEM::MEDCouplingUMesh::__iter__;
%newobject ParaMEDMEM::MEDCouplingUMesh::__getitem__;
%newobject ParaMEDMEM::MEDCouplingUMesh::cellsByType;
%newobject ParaMEDMEM::MEDCouplingUMesh::zipConnectivityTraducer;
%newobject ParaMEDMEM::MEDCouplingUMesh::buildDescendingConnectivity;
%newobject ParaMEDMEM::MEDCouplingUMesh::buildDescendingConnectivity2;
%newobject ParaMEDMEM::MEDCouplingUMesh::explode3DMeshTo1D;
%newobject ParaMEDMEM::MEDCouplingUMesh::buildExtrudedMesh;
%newobject ParaMEDMEM::MEDCouplingUMesh::buildSpreadZonesWithPoly;
%newobject ParaMEDMEM::MEDCouplingUMesh::MergeUMeshes;
%newobject ParaMEDMEM::MEDCouplingUMesh::MergeUMeshesOnSameCoords;
%newobject ParaMEDMEM::MEDCouplingUMesh::ComputeSpreadZoneGradually;
%newobject ParaMEDMEM::MEDCouplingUMesh::ComputeSpreadZoneGraduallyFromSeed;
%newobject ParaMEDMEM::MEDCouplingUMesh::buildNewNumberingFromCommNodesFrmt;
%newobject ParaMEDMEM::MEDCouplingUMesh::rearrange2ConsecutiveCellTypes;
%newobject ParaMEDMEM::MEDCouplingUMesh::sortCellsInMEDFileFrmt;
%newobject ParaMEDMEM::MEDCouplingUMesh::getRenumArrForMEDFileFrmt;
%newobject ParaMEDMEM::MEDCouplingUMesh::convertCellArrayPerGeoType;
%newobject ParaMEDMEM::MEDCouplingUMesh::computeFetchedNodeIds;
%newobject ParaMEDMEM::MEDCouplingUMesh::getRenumArrForConsecutiveCellTypesSpec;
%newobject ParaMEDMEM::MEDCouplingUMesh::buildDirectionVectorField;
%newobject ParaMEDMEM::MEDCouplingUMesh::convertLinearCellsToQuadratic;
%newobject ParaMEDMEM::MEDCouplingUMesh::getEdgeRatioField;
%newobject ParaMEDMEM::MEDCouplingUMesh::getAspectRatioField;
%newobject ParaMEDMEM::MEDCouplingUMesh::getWarpField;
%newobject ParaMEDMEM::MEDCouplingUMesh::getSkewField;
%newobject ParaMEDMEM::MEDCouplingUMesh::getPartBarycenterAndOwner;
%newobject ParaMEDMEM::MEDCouplingUMesh::getPartMeasureField;
%newobject ParaMEDMEM::MEDCouplingUMesh::buildPartOrthogonalField;
%newobject ParaMEDMEM::MEDCouplingUMesh::keepCellIdsByType;
%newobject ParaMEDMEM::MEDCouplingUMesh::Build0DMeshFromCoords;
%newobject ParaMEDMEM::MEDCouplingUMesh::findAndCorrectBadOriented3DExtrudedCells;
%newobject ParaMEDMEM::MEDCouplingUMesh::findAndCorrectBadOriented3DCells;
%newobject ParaMEDMEM::MEDCouplingUMesh::findCellIdsOnBoundary;
%newobject ParaMEDMEM::MEDCouplingUMesh::computeSkin;
%newobject ParaMEDMEM::MEDCouplingUMesh::getCellIdsLyingOnNodes;
%newobject ParaMEDMEM::MEDCouplingUMesh::buildSetInstanceFromThis;
%newobject ParaMEDMEM::MEDCouplingUMesh::getCellIdsCrossingPlane;
%newobject ParaMEDMEM::MEDCouplingUMesh::convexEnvelop2D;
%newobject ParaMEDMEM::MEDCouplingUMesh::ComputeRangesFromTypeDistribution;
%newobject ParaMEDMEM::MEDCouplingUMeshCellByTypeEntry::__iter__;
%newobject ParaMEDMEM::MEDCouplingUMeshCellEntry::__iter__;
%newobject ParaMEDMEM::MEDCouplingExtrudedMesh::New;
%newobject ParaMEDMEM::MEDCouplingExtrudedMesh::build3DUnstructuredMesh;
%newobject ParaMEDMEM::MEDCouplingCMesh::New;
%newobject ParaMEDMEM::MEDCouplingCMesh::clone;
%newobject ParaMEDMEM::MEDCouplingCMesh::getCoordsAt;
%newobject ParaMEDMEM::MEDCouplingCurveLinearMesh::New;
%newobject ParaMEDMEM::MEDCouplingCurveLinearMesh::clone;
%newobject ParaMEDMEM::MEDCouplingCurveLinearMesh::getCoords;
%newobject ParaMEDMEM::MEDCouplingMultiFields::New;
%newobject ParaMEDMEM::MEDCouplingMultiFields::deepCpy;
%newobject ParaMEDMEM::MEDCouplingFieldOverTime::New;

%feature("unref") DataArray "$this->decrRef();"
%feature("unref") DataArrayDouble "$this->decrRef();"
%feature("unref") MEDCouplingPointSet "$this->decrRef();"
%feature("unref") MEDCouplingMesh "$this->decrRef();"
%feature("unref") MEDCouplingUMesh "$this->decrRef();"
%feature("unref") MEDCouplingExtrudedMesh "$this->decrRef();"
%feature("unref") MEDCouplingCMesh "$this->decrRef();"
%feature("unref") DataArrayInt "$this->decrRef();"
%feature("unref") DataArrayChar "$this->decrRef();"
%feature("unref") DataArrayAsciiChar "$this->decrRef();"
%feature("unref") DataArrayByte "$this->decrRef();"
%feature("unref") MEDCouplingField "$this->decrRef();"
%feature("unref") MEDCouplingFieldDiscretizationP0 "$this->decrRef();"
%feature("unref") MEDCouplingFieldDiscretizationP1 "$this->decrRef();"
%feature("unref") MEDCouplingFieldDiscretizationGauss "$this->decrRef();"
%feature("unref") MEDCouplingFieldDiscretizationGaussNE "$this->decrRef();"
%feature("unref") MEDCouplingFieldDiscretizationKriging "$this->decrRef();"
%feature("unref") MEDCouplingFieldDouble "$this->decrRef();"
%feature("unref") MEDCouplingMultiFields "$this->decrRef();"
%feature("unref") MEDCouplingFieldTemplate "$this->decrRef();"
%feature("unref") MEDCouplingMultiFields "$this->decrRef();"

%rename(assign) *::operator=;
%ignore ParaMEDMEM::MEDCouplingVersionMajMinRel;
%ignore ParaMEDMEM::RefCountObject::decrRef;
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
    %extend
    {
      std::string __str__() const
        {
          return std::string(self->what());
        }
    }
  };
}

%include "MEDCouplingTimeLabel.hxx"

namespace ParaMEDMEM
{
  typedef enum
    {
      C_DEALLOC = 2,
      CPP_DEALLOC = 3
    } DeallocType;

  typedef enum
    {
      ON_CELLS = 0,
      ON_NODES = 1,
      ON_GAUSS_PT = 2,
      ON_GAUSS_NE = 3,
      ON_NODES_KR = 4
    } TypeOfField;

  typedef enum
    {
      NO_TIME = 4,
      ONE_TIME = 5,
      LINEAR_TIME = 6,
      CONST_ON_TIME_INTERVAL = 7
    } TypeOfTimeDiscretization;

  const char *MEDCouplingVersionStr();
  int MEDCouplingVersion();
  int MEDCouplingSizeOfVoidStar();
  PyObject *MEDCouplingVersionMajMinRel()
  {
    int tmp0=0,tmp1=0,tmp2=0;
    MEDCouplingVersionMajMinRel(tmp0,tmp1,tmp2);
    PyObject *res = PyList_New(3);
    PyList_SetItem(res,0,SWIG_From_int(tmp0));
    PyList_SetItem(res,1,SWIG_From_int(tmp1));
    PyList_SetItem(res,2,SWIG_From_int(tmp2));
    return res;
  }

  class RefCountObject
  {
  protected:
    RefCountObject();
    RefCountObject(const RefCountObject& other);
    ~RefCountObject();
  public:
    bool decrRef() const;
    void incrRef() const;
    virtual std::size_t getHeapMemorySize() const;
  };
}

%inline
{
  bool MEDCouplingHasNumPyBindings()
  {
#ifdef WITH_NUMPY
    return true;
#else
    return false;
#endif
  }

  std::string MEDCouplingCompletionScript() throw(INTERP_KERNEL::Exception)
  {
    static const char script[]="import rlcompleter,readline\nreadline.parse_and_bind('tab:complete')";
    std::ostringstream oss; oss << "MEDCouplingCompletionScript : error when trying to activate completion ! readline not present ?\nScript is :\n" << script;
    if(PyRun_SimpleString(script)!=0)
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    return std::string(script);
  }
}

%include "MEDCouplingMemArray.i"

namespace ParaMEDMEM
{
  typedef enum
    {
      UNSTRUCTURED = 5,
      UNSTRUCTURED_DESC = 6,
      CARTESIAN = 7,
      EXTRUDED = 8,
      CURVE_LINEAR = 9
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

  %extend MEDCouplingGaussLocalization
  {
    std::string __str__() const throw(INTERP_KERNEL::Exception)
    {
      return self->getStringRepr();
    }

    std::string __repr__() const throw(INTERP_KERNEL::Exception)
    {
      std::ostringstream oss; oss << "MEDCouplingGaussLocalization C++ instance at " << self << "." << std::endl;
      oss << self->getStringRepr();
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
    virtual MEDCouplingMeshType getType() const throw(INTERP_KERNEL::Exception);
    bool isStructured() const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingMesh *deepCpy() const;
    virtual bool isEqual(const MEDCouplingMesh *other, double prec) const throw(INTERP_KERNEL::Exception);
    virtual bool isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const throw(INTERP_KERNEL::Exception);
    virtual void copyTinyStringsFrom(const MEDCouplingMesh *other) throw(INTERP_KERNEL::Exception);
    virtual void copyTinyInfoFrom(const MEDCouplingMesh *other) throw(INTERP_KERNEL::Exception);
    virtual void checkCoherency() const throw(INTERP_KERNEL::Exception);
    virtual void checkCoherency1(double eps=1e-12) const throw(INTERP_KERNEL::Exception);
    virtual void checkCoherency2(double eps=1e-12) const throw(INTERP_KERNEL::Exception);
    virtual int getNumberOfCells() const throw(INTERP_KERNEL::Exception);
    virtual int getNumberOfNodes() const throw(INTERP_KERNEL::Exception);
    virtual int getSpaceDimension() const throw(INTERP_KERNEL::Exception);
    virtual int getMeshDimension() const throw(INTERP_KERNEL::Exception);
    virtual DataArrayDouble *getCoordinatesAndOwner() const throw(INTERP_KERNEL::Exception);
    virtual DataArrayDouble *getBarycenterAndOwner() const throw(INTERP_KERNEL::Exception);
    virtual DataArrayDouble *computeIsoBarycenterOfNodesPerCell() const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *giveCellsWithType(INTERP_KERNEL::NormalizedCellType type) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *computeNbOfNodesPerCell() const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *computeNbOfFacesPerCell() const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingMesh *buildPartRange(int beginCellIds, int endCellIds, int stepCellIds) const throw(INTERP_KERNEL::Exception);
    virtual int getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType type) const throw(INTERP_KERNEL::Exception);
    virtual INTERP_KERNEL::NormalizedCellType getTypeOfCell(int cellId) const throw(INTERP_KERNEL::Exception);
    virtual std::string simpleRepr() const throw(INTERP_KERNEL::Exception);
    virtual std::string advancedRepr() const throw(INTERP_KERNEL::Exception);
    void writeVTK(const char *fileName) const throw(INTERP_KERNEL::Exception);
    // tools
    virtual MEDCouplingFieldDouble *getMeasureField(bool isAbs) const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingFieldDouble *getMeasureFieldOnNode(bool isAbs) const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingFieldDouble *fillFromAnalytic(TypeOfField t, int nbOfComp, const char *func) const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingFieldDouble *fillFromAnalytic2(TypeOfField t, int nbOfComp, const char *func) const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingFieldDouble *fillFromAnalytic3(TypeOfField t, int nbOfComp, const std::vector<std::string>& varsOrder, const char *func) const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingFieldDouble *buildOrthogonalField() const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingUMesh *buildUnstructured() const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingMesh *mergeMyselfWith(const MEDCouplingMesh *other) const throw(INTERP_KERNEL::Exception);
    virtual bool areCompatibleForMerge(const MEDCouplingMesh *other) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *simplexize(int policy) throw(INTERP_KERNEL::Exception);
    static MEDCouplingMesh *MergeMeshes(const MEDCouplingMesh *mesh1, const MEDCouplingMesh *mesh2) throw(INTERP_KERNEL::Exception);
    static bool IsStaticGeometricType(INTERP_KERNEL::NormalizedCellType type) throw(INTERP_KERNEL::Exception);
    static bool IsLinearGeometricType(INTERP_KERNEL::NormalizedCellType type) throw(INTERP_KERNEL::Exception);
    static INTERP_KERNEL::NormalizedCellType GetCorrespondingPolyType(INTERP_KERNEL::NormalizedCellType type) throw(INTERP_KERNEL::Exception);
    static int GetNumberOfNodesOfGeometricType(INTERP_KERNEL::NormalizedCellType type) throw(INTERP_KERNEL::Exception);
    static int GetDimensionOfGeometricType(INTERP_KERNEL::NormalizedCellType type) throw(INTERP_KERNEL::Exception);
    static const char *GetReprOfGeometricType(INTERP_KERNEL::NormalizedCellType type) throw(INTERP_KERNEL::Exception);
    %extend
       {
         std::string __str__() const throw(INTERP_KERNEL::Exception)
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
           double val;
           DataArrayDouble *a;
           DataArrayDoubleTuple *aa;
           std::vector<double> bb;
           int sw;
           int spaceDim=self->getSpaceDimension();
           const char msg[]="Python wrap of MEDCouplingMesh::getCellContainingPoint : ";
           const double *pos=convertObjToPossibleCpp5_Safe(p,sw,val,a,aa,bb,msg,1,spaceDim,true);
           return self->getCellContainingPoint(pos,eps);
         }

         PyObject *getCellsContainingPoints(PyObject *p, int nbOfPoints, double eps) const throw(INTERP_KERNEL::Exception)
         {
           double val;
           DataArrayDouble *a;
           DataArrayDoubleTuple *aa;
           std::vector<double> bb;
           int sw;
           int spaceDim=self->getSpaceDimension();
           const char msg[]="Python wrap of MEDCouplingMesh::getCellsContainingPoint : ";
           const double *pos=convertObjToPossibleCpp5_Safe(p,sw,val,a,aa,bb,msg,nbOfPoints,spaceDim,true);
           std::vector<int> elts,eltsIndex;
           self->getCellsContainingPoints(pos,nbOfPoints,eps,elts,eltsIndex);
           MEDCouplingAutoRefCountObjectPtr<DataArrayInt> d0=DataArrayInt::New();
           MEDCouplingAutoRefCountObjectPtr<DataArrayInt> d1=DataArrayInt::New();
           d0->alloc(elts.size(),1);
           d1->alloc(eltsIndex.size(),1);
           std::copy(elts.begin(),elts.end(),d0->getPointer());
           std::copy(eltsIndex.begin(),eltsIndex.end(),d1->getPointer());
           PyObject *ret=PyTuple_New(2);
           PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(d0.retn()),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
           PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(d1.retn()),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
           return ret;
         }

         PyObject *getCellsContainingPoints(PyObject *p, double eps) const throw(INTERP_KERNEL::Exception)
         {
           std::vector<int> elts,eltsIndex;
           int spaceDim=self->getSpaceDimension();
           void *da=0;
           int res1=SWIG_ConvertPtr(p,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayDouble, 0 |  0 );
           if (!SWIG_IsOK(res1))
             {
               int size;
               INTERP_KERNEL::AutoCPtr<double> tmp=convertPyToNewDblArr2(p,&size);
               int nbOfPoints=size/spaceDim;
               if(size%spaceDim!=0)
                 {
                   throw INTERP_KERNEL::Exception("MEDCouplingMesh::getCellsContainingPoints : Invalid list length ! Must be a multiple of self.getSpaceDimension() !");
                 }
               self->getCellsContainingPoints(tmp,nbOfPoints,eps,elts,eltsIndex);
             }
           else
             {
               DataArrayDouble *da2=reinterpret_cast< DataArrayDouble * >(da);
               if(!da2)
                 throw INTERP_KERNEL::Exception("MEDCouplingMesh::getCellsContainingPoints : Not null DataArrayDouble instance expected !");
               da2->checkAllocated();
               int size=da2->getNumberOfTuples();
               int nbOfCompo=da2->getNumberOfComponents();
               if(nbOfCompo!=spaceDim)
                 {
                   throw INTERP_KERNEL::Exception("MEDCouplingMesh::getCellsContainingPoints : Invalid DataArrayDouble nb of components ! Expected same as self.getSpaceDimension() !");
                 }
               self->getCellsContainingPoints(da2->getConstPointer(),size,eps,elts,eltsIndex);
             }
           MEDCouplingAutoRefCountObjectPtr<DataArrayInt> d0=DataArrayInt::New();
           MEDCouplingAutoRefCountObjectPtr<DataArrayInt> d1=DataArrayInt::New();
           d0->alloc(elts.size(),1);
           d1->alloc(eltsIndex.size(),1);
           std::copy(elts.begin(),elts.end(),d0->getPointer());
           std::copy(eltsIndex.begin(),eltsIndex.end(),d1->getPointer());
           PyObject *ret=PyTuple_New(2);
           PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(d0.retn()),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
           PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(d1.retn()),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
           return ret;
         }

         PyObject *getCellsContainingPoint(PyObject *p, double eps) const throw(INTERP_KERNEL::Exception)
         {
           double val;
           DataArrayDouble *a;
           DataArrayDoubleTuple *aa;
           std::vector<double> bb;
           int sw;
           int spaceDim=self->getSpaceDimension();
           const char msg[]="Python wrap of MEDCouplingUMesh::getCellsContainingPoint : ";
           const double *pos=convertObjToPossibleCpp5_Safe(p,sw,val,a,aa,bb,msg,1,spaceDim,true);
           std::vector<int> elts;
           self->getCellsContainingPoint(pos,eps,elts);
           DataArrayInt *ret=DataArrayInt::New();
           ret->alloc((int)elts.size(),1);
           std::copy(elts.begin(),elts.end(),ret->getPointer());
           return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 );
         }
         
         void renumberCells(PyObject *li, bool check=true) throw(INTERP_KERNEL::Exception)
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

         PyObject *checkDeepEquivalWith(const MEDCouplingMesh *other, int cellCompPol, double prec) const throw(INTERP_KERNEL::Exception)
         {
           DataArrayInt *cellCor=0,*nodeCor=0;
           self->checkDeepEquivalWith(other,cellCompPol,prec,cellCor,nodeCor);
           PyObject *res = PyList_New(2);
           PyList_SetItem(res,0,SWIG_NewPointerObj(SWIG_as_voidptr(cellCor),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, cellCor?SWIG_POINTER_OWN | 0:0 ));
           PyList_SetItem(res,1,SWIG_NewPointerObj(SWIG_as_voidptr(nodeCor),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, nodeCor?SWIG_POINTER_OWN | 0:0 ));
           return res;
         }
         
         DataArrayInt *checkDeepEquivalOnSameNodesWith(const MEDCouplingMesh *other, int cellCompPol, double prec) const throw(INTERP_KERNEL::Exception)
         {
           DataArrayInt *cellCor=0;
           self->checkDeepEquivalOnSameNodesWith(other,cellCompPol,prec,cellCor);
           return cellCor;
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
           double val;
           DataArrayDouble *a;
           DataArrayDoubleTuple *aa;
           std::vector<double> bb;
           int sw;
           int spaceDim=self->getSpaceDimension();
           const char msg[]="Python wrap of MEDCouplingPointSet::scale : ";
           const double *pointPtr=convertObjToPossibleCpp5_Safe(point,sw,val,a,aa,bb,msg,1,spaceDim,true);
           self->scale(pointPtr,factor);
         }

         PyObject *getBoundingBox() const throw(INTERP_KERNEL::Exception)
         {
           int spaceDim=self->getSpaceDimension();
           INTERP_KERNEL::AutoPtr<double> tmp=new double[2*spaceDim];
           self->getBoundingBox(tmp);
           PyObject *ret=convertDblArrToPyListOfTuple(tmp,2,spaceDim);
           return ret;
         }

         PyObject *isEqualIfNotWhy(const MEDCouplingMesh *other, double prec) const throw(INTERP_KERNEL::Exception)
         {
           std::string ret1;
           bool ret0=self->isEqualIfNotWhy(other,prec,ret1);
           PyObject *ret=PyTuple_New(2);
           PyObject *ret0Py=ret0?Py_True:Py_False;
           Py_XINCREF(ret0Py);
           PyTuple_SetItem(ret,0,ret0Py);
           PyTuple_SetItem(ret,1,PyString_FromString(ret1.c_str()));
           return ret;
         }

         PyObject *buildPart(PyObject *li) const throw(INTERP_KERNEL::Exception)
         {
           int szArr,sw,iTypppArr;
           std::vector<int> stdvecTyyppArr;
           const int *tmp=convertObjToPossibleCpp1_Safe(li,sw,szArr,iTypppArr,stdvecTyyppArr);
           MEDCouplingMesh *ret=self->buildPart(tmp,tmp+szArr);
           if(sw==3)//DataArrayInt
             { 
               void *argp; SWIG_ConvertPtr(li,&argp,SWIGTYPE_p_ParaMEDMEM__DataArrayInt,0|0);
               DataArrayInt *argpt=reinterpret_cast< ParaMEDMEM::DataArrayInt * >(argp);
               std::string name=argpt->getName();
               if(!name.empty())
                 ret->setName(name.c_str());
             }
           return convertMesh(ret, SWIG_POINTER_OWN | 0 );
         }
        
         PyObject *buildPartAndReduceNodes(PyObject *li) const throw(INTERP_KERNEL::Exception)
         {
           int szArr,sw,iTypppArr;
           std::vector<int> stdvecTyyppArr;
           DataArrayInt *arr=0;
           const int *tmp=convertObjToPossibleCpp1_Safe(li,sw,szArr,iTypppArr,stdvecTyyppArr);
           MEDCouplingMesh *ret=self->buildPartAndReduceNodes(tmp,tmp+szArr,arr);
           if(sw==3)//DataArrayInt
             { 
               void *argp; SWIG_ConvertPtr(li,&argp,SWIGTYPE_p_ParaMEDMEM__DataArrayInt,0|0);
               DataArrayInt *argpt=reinterpret_cast< ParaMEDMEM::DataArrayInt * >(argp);
               std::string name=argpt->getName();
               if(!name.empty())
                 ret->setName(name.c_str());
             }
           //
           PyObject *res = PyList_New(2);
           PyObject *obj0=convertMesh(ret, SWIG_POINTER_OWN | 0 );
           PyObject *obj1=SWIG_NewPointerObj(SWIG_as_voidptr(arr),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 );
           PyList_SetItem(res,0,obj0);
           PyList_SetItem(res,1,obj1);
           return res;
         }

         PyObject *buildPartRangeAndReduceNodes(int beginCellIds, int endCellIds, int stepCellIds) const throw(INTERP_KERNEL::Exception)
         {
           int a,b,c;
           DataArrayInt *arr=0;
           MEDCouplingMesh *ret=self->buildPartRangeAndReduceNodes(beginCellIds,endCellIds,stepCellIds,a,b,c,arr);
           PyObject *res = PyTuple_New(2);
           PyObject *obj0=convertMesh(ret, SWIG_POINTER_OWN | 0 );
           PyObject *obj1=0;
           if(arr)
             obj1=SWIG_NewPointerObj(SWIG_as_voidptr(arr),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 );
           else
             obj1=PySlice_New(PyInt_FromLong(a),PyInt_FromLong(b),PyInt_FromLong(b));
           PyTuple_SetItem(res,0,obj0);
           PyTuple_SetItem(res,1,obj1);
           return res;
         }

        PyObject *getDistributionOfTypes() const throw(INTERP_KERNEL::Exception)
        {
          std::vector<int> vals=self->getDistributionOfTypes();
          PyObject *ret=PyList_New((int)vals.size()/3);
          for(int j=0;j<(int)vals.size()/3;j++)
             {
               PyObject *ret1=PyList_New(3);
               PyList_SetItem(ret1,0,SWIG_From_int(vals[3*j]));
               PyList_SetItem(ret1,1,SWIG_From_int(vals[3*j+1]));
               PyList_SetItem(ret1,2,SWIG_From_int(vals[3*j+2]));
               PyList_SetItem(ret,j,ret1);
             }
          return ret;
        }

        DataArrayInt *checkTypeConsistencyAndContig(PyObject *li, PyObject *li2) const throw(INTERP_KERNEL::Exception)
        {
          std::vector<int> code;
          std::vector<const DataArrayInt *> idsPerType;
          convertFromPyObjVectorOfObj<const ParaMEDMEM::DataArrayInt *>(li2,SWIGTYPE_p_ParaMEDMEM__DataArrayInt,"DataArrayInt",idsPerType);
          convertPyToNewIntArr4(li,1,3,code);
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
          double val;
          DataArrayDouble *a;
          DataArrayDoubleTuple *aa;
          std::vector<double> bb;
          int sw;
          int spaceDim=self->getSpaceDimension();
          const char msg[]="Python wrap of MEDCouplingPointSet::translate : ";
          const double *vectorPtr=convertObjToPossibleCpp5_Safe(vector,sw,val,a,aa,bb,msg,1,spaceDim,true);
          self->translate(vectorPtr);
        }

         void rotate(PyObject *center, double alpha) throw(INTERP_KERNEL::Exception)
         {
           const char msg[]="Python wrap of MEDCouplingPointSet::rotate : ";
           double val;
           DataArrayDouble *a;
           DataArrayDoubleTuple *aa;
           std::vector<double> bb;
           int sw;
           int spaceDim=self->getSpaceDimension();
           const double *centerPtr=convertObjToPossibleCpp5_Safe(center,sw,val,a,aa,bb,msg,1,spaceDim,true);
           self->rotate(centerPtr,0,alpha);
         }

         void rotate(PyObject *center, PyObject *vector, double alpha) throw(INTERP_KERNEL::Exception)
         {
           const char msg[]="Python wrap of MEDCouplingPointSet::rotate : ";
           double val,val2;
           DataArrayDouble *a,*a2;
           DataArrayDoubleTuple *aa,*aa2;
           std::vector<double> bb,bb2;
           int sw;
           int spaceDim=self->getSpaceDimension();
           const double *centerPtr=convertObjToPossibleCpp5_Safe(center,sw,val,a,aa,bb,msg,1,spaceDim,true);
           const double *vectorPtr=convertObjToPossibleCpp5_Safe(vector,sw,val2,a2,aa2,bb2,msg,1,spaceDim,false);//vectorPtr can be null in case of space dim 2
           self->rotate(centerPtr,vectorPtr,alpha);
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
         
         static MEDCouplingMesh *MergeMeshes(PyObject *li) throw(INTERP_KERNEL::Exception)
         {
            std::vector<const ParaMEDMEM::MEDCouplingMesh *> tmp;
            convertFromPyObjVectorOfObj<const ParaMEDMEM::MEDCouplingMesh *>(li,SWIGTYPE_p_ParaMEDMEM__MEDCouplingMesh,"MEDCouplingMesh",tmp);
            return MEDCouplingMesh::MergeMeshes(tmp);
         }
       }
  };
}

%include "NormalizedUnstructuredMesh.hxx"
%include "MEDCouplingNatureOfField.hxx"
%include "MEDCouplingTimeDiscretization.hxx"
%include "MEDCouplingGaussLocalization.hxx"
%include "MEDCouplingFieldDiscretization.hxx"

%ignore ParaMEDMEM::MEDCouplingFieldDiscretization::clonePart;
%ignore ParaMEDMEM::MEDCouplingFieldDiscretization::buildSubMeshDataRange;
%ignore ParaMEDMEM::MEDCouplingFieldDiscretizationPerCell::getArrayOfDiscIds;

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
      void recenterForMaxPrecision(double eps) throw(INTERP_KERNEL::Exception);
      void changeSpaceDimension(int newSpaceDim, double dftVal=0.) throw(INTERP_KERNEL::Exception);
      void tryToShareSameCoords(const MEDCouplingPointSet& other, double epsilon) throw(INTERP_KERNEL::Exception);
      virtual MEDCouplingPointSet *buildPartOfMySelf2(int start, int end, int step) const throw(INTERP_KERNEL::Exception);
      virtual void tryToShareSameCoordsPermute(const MEDCouplingPointSet& other, double epsilon) throw(INTERP_KERNEL::Exception);
      static DataArrayDouble *MergeNodesArray(const MEDCouplingPointSet *m1, const MEDCouplingPointSet *m2) throw(INTERP_KERNEL::Exception);
      static MEDCouplingPointSet *BuildInstanceFromMeshType(MEDCouplingMeshType type) throw(INTERP_KERNEL::Exception);
      virtual MEDCouplingPointSet *buildBoundaryMesh(bool keepCoords) const throw(INTERP_KERNEL::Exception);
      virtual bool isEmptyMesh(const std::vector<int>& tinyInfo) const throw(INTERP_KERNEL::Exception);
      virtual DataArrayInt *getCellsInBoundingBox(const INTERP_KERNEL::DirectedBoundingBox& bbox, double eps) throw(INTERP_KERNEL::Exception);
      virtual DataArrayInt *zipCoordsTraducer() throw(INTERP_KERNEL::Exception);
      virtual DataArrayInt *findBoundaryNodes() const;
      %extend 
         {
           std::string __str__() const throw(INTERP_KERNEL::Exception)
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
           
           PyObject *findCommonNodes(double prec, int limitTupleId=-1) const throw(INTERP_KERNEL::Exception)
           {
             DataArrayInt *comm, *commIndex;
             self->findCommonNodes(prec,limitTupleId,comm,commIndex);
             PyObject *res = PyList_New(2);
             PyList_SetItem(res,0,SWIG_NewPointerObj(SWIG_as_voidptr(comm),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
             PyList_SetItem(res,1,SWIG_NewPointerObj(SWIG_as_voidptr(commIndex),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
             return res;
           }
           
           PyObject *getCoords() throw(INTERP_KERNEL::Exception)
           {
             DataArrayDouble *ret1=self->getCoords();
             if (ret1)
                ret1->incrRef();
             return SWIG_NewPointerObj((void*)ret1,SWIGTYPE_p_ParaMEDMEM__DataArrayDouble,SWIG_POINTER_OWN | 0);
           }
           
           PyObject *buildPartOfMySelf(PyObject *li, bool keepCoords=true) const throw(INTERP_KERNEL::Exception)
           {
             int szArr,sw,iTypppArr;
             std::vector<int> stdvecTyyppArr;
             const int *tmp=convertObjToPossibleCpp1_Safe(li,sw,szArr,iTypppArr,stdvecTyyppArr);
             MEDCouplingPointSet *ret=self->buildPartOfMySelf(tmp,tmp+szArr,keepCoords);
             if(sw==3)//DataArrayInt
               { 
                 void *argp; SWIG_ConvertPtr(li,&argp,SWIGTYPE_p_ParaMEDMEM__DataArrayInt,0|0);
                 DataArrayInt *argpt=reinterpret_cast< ParaMEDMEM::DataArrayInt * >(argp);
                 std::string name=argpt->getName();
                 if(!name.empty())
                   ret->setName(name.c_str());
               }
             return convertMesh(ret, SWIG_POINTER_OWN | 0 );
           }
           
           PyObject *buildPartOfMySelfNode(PyObject *li, bool fullyIn) const throw(INTERP_KERNEL::Exception)
           {
             int szArr,sw,iTypppArr;
             std::vector<int> stdvecTyyppArr;
             const int *tmp=convertObjToPossibleCpp1_Safe(li,sw,szArr,iTypppArr,stdvecTyyppArr);
             MEDCouplingPointSet *ret=self->buildPartOfMySelfNode(tmp,tmp+szArr,fullyIn);
             if(sw==3)//DataArrayInt
               { 
                 void *argp; SWIG_ConvertPtr(li,&argp,SWIGTYPE_p_ParaMEDMEM__DataArrayInt,0|0);
                 DataArrayInt *argpt=reinterpret_cast< ParaMEDMEM::DataArrayInt * >(argp);
                 std::string name=argpt->getName();
                 if(!name.empty())
                   ret->setName(name.c_str());
               }
             return convertMesh(ret, SWIG_POINTER_OWN | 0 );
           }

           PyObject *buildFacePartOfMySelfNode(PyObject *li, bool fullyIn) const throw(INTERP_KERNEL::Exception)
           {
             int szArr,sw,iTypppArr;
             std::vector<int> stdvecTyyppArr;
             const int *tmp=convertObjToPossibleCpp1_Safe(li,sw,szArr,iTypppArr,stdvecTyyppArr);
             MEDCouplingPointSet *ret=self->buildFacePartOfMySelfNode(tmp,tmp+szArr,fullyIn);
             if(sw==3)//DataArrayInt
               { 
                 void *argp; SWIG_ConvertPtr(li,&argp,SWIGTYPE_p_ParaMEDMEM__DataArrayInt,0|0);
                 DataArrayInt *argpt=reinterpret_cast< ParaMEDMEM::DataArrayInt * >(argp);
                 std::string name=argpt->getName();
                 if(!name.empty())
                   ret->setName(name.c_str());
               }
             return convertMesh(ret, SWIG_POINTER_OWN | 0 );
           }

           void renumberNodes(PyObject *li, int newNbOfNodes) throw(INTERP_KERNEL::Exception)
           {
             int szArr,sw,iTypppArr;
             std::vector<int> stdvecTyyppArr;
             const int *tmp=convertObjToPossibleCpp1_Safe(li,sw,szArr,iTypppArr,stdvecTyyppArr);
             self->renumberNodes(tmp,newNbOfNodes);
           }

           void renumberNodes2(PyObject *li, int newNbOfNodes) throw(INTERP_KERNEL::Exception)
           {
             int szArr,sw,iTypppArr;
             std::vector<int> stdvecTyyppArr;
             const int *tmp=convertObjToPossibleCpp1_Safe(li,sw,szArr,iTypppArr,stdvecTyyppArr);
             self->renumberNodes2(tmp,newNbOfNodes);
           }

           PyObject *findNodesOnLine(PyObject *pt, PyObject *vec, double eps) const throw(INTERP_KERNEL::Exception)
             {
               int spaceDim=self->getSpaceDimension();
               double val,val2;
               DataArrayDouble *a,*a2;
               DataArrayDoubleTuple *aa,*aa2;
               std::vector<double> bb,bb2;
               int sw;
               const char msg[]="Python wrap of MEDCouplingPointSet::findNodesOnLine : 1st paramater for point.";
               const char msg2[]="Python wrap of MEDCouplingPointSet::findNodesOnLine : 2nd paramater for vector.";
               const double *p=convertObjToPossibleCpp5_Safe(pt,sw,val,a,aa,bb,msg,1,spaceDim,true);
               const double *v=convertObjToPossibleCpp5_Safe(vec,sw,val2,a2,aa2,bb2,msg2,1,spaceDim,true);
               std::vector<int> nodes;
               self->findNodesOnLine(p,v,eps,nodes);
               DataArrayInt *ret=DataArrayInt::New();
               ret->alloc((int)nodes.size(),1);
               std::copy(nodes.begin(),nodes.end(),ret->getPointer());
               return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 );
             }
           PyObject *findNodesOnPlane(PyObject *pt, PyObject *vec, double eps) const throw(INTERP_KERNEL::Exception)
             {
               int spaceDim=self->getSpaceDimension();
               double val,val2;
               DataArrayDouble *a,*a2;
               DataArrayDoubleTuple *aa,*aa2;
               std::vector<double> bb,bb2;
               int sw;
               const char msg[]="Python wrap of MEDCouplingPointSet::findNodesOnPlane : 1st paramater for point.";
               const char msg2[]="Python wrap of MEDCouplingPointSet::findNodesOnPlane : 2nd paramater for vector.";
               const double *p=convertObjToPossibleCpp5_Safe(pt,sw,val,a,aa,bb,msg,1,spaceDim,true);
               const double *v=convertObjToPossibleCpp5_Safe(vec,sw,val2,a2,aa2,bb2,msg2,1,spaceDim,true);
               std::vector<int> nodes;
               self->findNodesOnPlane(p,v,eps,nodes);
               DataArrayInt *ret=DataArrayInt::New();
               ret->alloc((int)nodes.size(),1);
               std::copy(nodes.begin(),nodes.end(),ret->getPointer());
               return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 );
             }
           
           PyObject *getNodeIdsNearPoint(PyObject *pt, double eps) const throw(INTERP_KERNEL::Exception)
           {
             double val;
             DataArrayDouble *a;
             DataArrayDoubleTuple *aa;
             std::vector<double> bb;
             int sw;
             int spaceDim=self->getSpaceDimension();
             const char msg[]="Python wrap of MEDCouplingPointSet::getNodeIdsNearPoint : ";
             const double *pos=convertObjToPossibleCpp5_Safe(pt,sw,val,a,aa,bb,msg,1,spaceDim,true);
             DataArrayInt *ret=self->getNodeIdsNearPoint(pos,eps);
             return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 );
           }

           PyObject *getNodeIdsNearPoints(PyObject *pt, int nbOfPoints, double eps) const throw(INTERP_KERNEL::Exception)
           {
             DataArrayInt *c=0,*cI=0;
             //
             double val;
             DataArrayDouble *a;
             DataArrayDoubleTuple *aa;
             std::vector<double> bb;
             int sw;
             int spaceDim=self->getSpaceDimension();
             const char msg[]="Python wrap of MEDCouplingPointSet::getNodeIdsNearPoints : ";
             const double *pos=convertObjToPossibleCpp5_Safe(pt,sw,val,a,aa,bb,msg,nbOfPoints,spaceDim,true);
             self->getNodeIdsNearPoints(pos,nbOfPoints,eps,c,cI);
             PyObject *ret=PyTuple_New(2);
             PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(c),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
             PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(cI),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
             return ret;
           }

           PyObject *getNodeIdsNearPoints(PyObject *pt, double eps) const throw(INTERP_KERNEL::Exception)
           {
             DataArrayInt *c=0,*cI=0;
             int spaceDim=self->getSpaceDimension();
             double val;
             DataArrayDouble *a;
             DataArrayDoubleTuple *aa;
             std::vector<double> bb;
             int sw;
             int nbOfTuples=-1;
             const double *ptPtr=convertObjToPossibleCpp5_Safe2(pt,sw,val,a,aa,bb,"Python wrap of MEDCouplingUMesh::getNodeIdsNearPoints",spaceDim,true,nbOfTuples);
             self->getNodeIdsNearPoints(ptPtr,nbOfTuples,eps,c,cI);
             //
             PyObject *ret=PyTuple_New(2);
             PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(c),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
             PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(cI),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
             return ret;
           }

           PyObject *getCellsInBoundingBox(PyObject *bbox, double eps) const throw(INTERP_KERNEL::Exception)
           {
             double val;
             DataArrayDouble *a;
             DataArrayDoubleTuple *aa;
             std::vector<double> bb;
             int sw;
             int spaceDim=self->getSpaceDimension();
             const char msg[]="Python wrap of MEDCouplingPointSet::getCellsInBoundingBox : ";
             const double *tmp=convertObjToPossibleCpp5_Safe(bbox,sw,val,a,aa,bb,msg,spaceDim,2,true);
             //
             DataArrayInt *elems=self->getCellsInBoundingBox(tmp,eps);
             return SWIG_NewPointerObj(SWIG_as_voidptr(elems),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 );
           }

           void duplicateNodesInCoords(PyObject *li) throw(INTERP_KERNEL::Exception)
           {
             int sw;
             int singleVal;
             std::vector<int> multiVal;
             std::pair<int, std::pair<int,int> > slic;
             ParaMEDMEM::DataArrayInt *daIntTyypp=0;
             convertObjToPossibleCpp2(li,self->getNumberOfNodes(),sw,singleVal,multiVal,slic,daIntTyypp);
             switch(sw)
               {
               case 1:
                 return self->duplicateNodesInCoords(&singleVal,&singleVal+1);
               case 2:
                 return self->duplicateNodesInCoords(&multiVal[0],&multiVal[0]+multiVal.size());
               case 4:
                 return self->duplicateNodesInCoords(daIntTyypp->begin(),daIntTyypp->end());
               default:
                 throw INTERP_KERNEL::Exception("MEDCouplingPointSet::duplicateNodesInCoords : unrecognized type entered, expected list of int, tuple of int or DataArrayInt !");
               }
           }

           static void Rotate2DAlg(PyObject *center, double angle, int nbNodes, PyObject *coords) throw(INTERP_KERNEL::Exception)
           {
             int sz;
             INTERP_KERNEL::AutoCPtr<double> c=convertPyToNewDblArr2(center,&sz);
             INTERP_KERNEL::AutoCPtr<double> coo=convertPyToNewDblArr2(coords,&sz);
             ParaMEDMEM::MEDCouplingPointSet::Rotate2DAlg(c,angle,nbNodes,coo);
             for(int i=0;i<sz;i++)
               PyList_SetItem(coords,i,PyFloat_FromDouble(coo[i]));
           }
           
           static void Rotate2DAlg(PyObject *center, double angle, PyObject *coords) throw(INTERP_KERNEL::Exception)
           {
             int sz;
             INTERP_KERNEL::AutoCPtr<double> c=convertPyToNewDblArr2(center,&sz);
             int sw,nbNodes=0;
             double val0;  ParaMEDMEM::DataArrayDouble *val1=0; ParaMEDMEM::DataArrayDoubleTuple *val2=0;
             std::vector<double> val3;
             const double *coo=convertObjToPossibleCpp5_Safe2(coords,sw,val0,val1,val2,val3,
                                                            "Rotate2DAlg",2,true,nbNodes);
             if(sw!=2 && sw!=3)
               throw INTERP_KERNEL::Exception("Invalid call to MEDCouplingPointSet::Rotate2DAlg : try another overload method !");
             ParaMEDMEM::MEDCouplingPointSet::Rotate2DAlg(c,angle,nbNodes,const_cast<double *>(coo));
           }
           
           static void Rotate3DAlg(PyObject *center, PyObject *vect, double angle, int nbNodes, PyObject *coords) throw(INTERP_KERNEL::Exception)
           {
             int sz,sz2;
             INTERP_KERNEL::AutoCPtr<double> c=convertPyToNewDblArr2(center,&sz);
             INTERP_KERNEL::AutoCPtr<double> coo=convertPyToNewDblArr2(coords,&sz);
             INTERP_KERNEL::AutoCPtr<double> v=convertPyToNewDblArr2(vect,&sz2);
             ParaMEDMEM::MEDCouplingPointSet::Rotate3DAlg(c,v,angle,nbNodes,coo);
             for(int i=0;i<sz;i++)
               PyList_SetItem(coords,i,PyFloat_FromDouble(coo[i]));
           }
           
           static void Rotate3DAlg(PyObject *center, PyObject *vect, double angle, PyObject *coords) throw(INTERP_KERNEL::Exception)
           {
             int sz,sz2;
             INTERP_KERNEL::AutoCPtr<double> c=convertPyToNewDblArr2(center,&sz);
             int sw,nbNodes=0;
             double val0;  ParaMEDMEM::DataArrayDouble *val1=0; ParaMEDMEM::DataArrayDoubleTuple *val2=0;
             std::vector<double> val3;
             const double *coo=convertObjToPossibleCpp5_Safe2(coords,sw,val0,val1,val2,val3,
                                                            "Rotate3DAlg",3,true,nbNodes);
             if(sw!=2 && sw!=3)
               throw INTERP_KERNEL::Exception("Invalid call to MEDCouplingPointSet::Rotate3DAlg : try another overload method !");
             INTERP_KERNEL::AutoCPtr<double> v=convertPyToNewDblArr2(vect,&sz2);
             ParaMEDMEM::MEDCouplingPointSet::Rotate3DAlg(c,v,angle,nbNodes,const_cast<double *>(coo));
           }
         }
    };

  class MEDCouplingUMeshCell
  {
  public:
    INTERP_KERNEL::NormalizedCellType getType() const;
    %extend
      {
        std::string __str__() const throw(INTERP_KERNEL::Exception)
        {
          return self->repr();
        }

        PyObject *getAllConn() const throw(INTERP_KERNEL::Exception)
        {
          int ret2;
          const int *r=self->getAllConn(ret2);
          PyObject *ret=PyTuple_New(ret2);
          for(int i=0;i<ret2;i++)
            PyTuple_SetItem(ret,i,PyInt_FromLong(r[i]));
          return ret;
        }
      }
  };

  class MEDCouplingUMeshCellIterator
  {
  public:
    %extend
      {
        PyObject *next()
        {
          MEDCouplingUMeshCell *ret=self->nextt();
          if(ret)
            return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMeshCell,0|0);
          else
            {
              PyErr_SetString(PyExc_StopIteration,"No more data.");
              return 0;
            }
        }
      }
  };

  class MEDCouplingUMeshCellByTypeIterator
  {
  public:
    ~MEDCouplingUMeshCellByTypeIterator();
    %extend
      {
        PyObject *next()
        {
          MEDCouplingUMeshCellEntry *ret=self->nextt();
          if(ret)
            return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMeshCellEntry,SWIG_POINTER_OWN | 0);
          else
            {
              PyErr_SetString(PyExc_StopIteration,"No more data.");
              return 0;
            }
        }
      }
  };

  class MEDCouplingUMeshCellByTypeEntry
  {
  public:
    ~MEDCouplingUMeshCellByTypeEntry();
    %extend
      {
        MEDCouplingUMeshCellByTypeIterator *__iter__()
        {
          return self->iterator();
        }
      }
  };

  class MEDCouplingUMeshCellEntry
  {
  public:
    INTERP_KERNEL::NormalizedCellType getType() const;
    int getNumberOfElems() const;
    %extend
      {
        MEDCouplingUMeshCellIterator *__iter__()
        {
          return self->iterator();
        }
      }
  };
  
  class MEDCouplingUMesh : public ParaMEDMEM::MEDCouplingPointSet
  {
  public:
    static MEDCouplingUMesh *New() throw(INTERP_KERNEL::Exception);
    static MEDCouplingUMesh *New(const char *meshName, int meshDim) throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *clone(bool recDeepCpy) const;
    void updateTime() const;
    void checkCoherency() const throw(INTERP_KERNEL::Exception);
    void setMeshDimension(int meshDim) throw(INTERP_KERNEL::Exception);
    void allocateCells(int nbOfCells=0) throw(INTERP_KERNEL::Exception);
    void finishInsertingCells() throw(INTERP_KERNEL::Exception);
    MEDCouplingUMeshCellByTypeEntry *cellsByType() throw(INTERP_KERNEL::Exception);
    void setConnectivity(DataArrayInt *conn, DataArrayInt *connIndex, bool isComputingTypes=true) throw(INTERP_KERNEL::Exception);
    INTERP_KERNEL::NormalizedCellType getTypeOfCell(int cellId) const throw(INTERP_KERNEL::Exception);
    void setPartOfMySelf2(int start, int end, int step, const MEDCouplingUMesh& otherOnSameCoordsThanThis) throw(INTERP_KERNEL::Exception);
    int getNumberOfNodesInCell(int cellId) const throw(INTERP_KERNEL::Exception);
    int getMeshLength() const throw(INTERP_KERNEL::Exception);
    void computeTypes() throw(INTERP_KERNEL::Exception);
    std::string reprConnectivityOfThis() const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *buildSetInstanceFromThis(int spaceDim) const throw(INTERP_KERNEL::Exception);
    //tools
    void shiftNodeNumbersInConn(int delta) throw(INTERP_KERNEL::Exception);
    std::vector<bool> getQuadraticStatus() const throw(INTERP_KERNEL::Exception);
    DataArrayInt *findCellIdsOnBoundary() const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *computeSkin() const throw(INTERP_KERNEL::Exception);
    bool checkConsecutiveCellTypes() const throw(INTERP_KERNEL::Exception);
    bool checkConsecutiveCellTypesForMEDFileFrmt() const throw(INTERP_KERNEL::Exception);
    DataArrayInt *rearrange2ConsecutiveCellTypes() throw(INTERP_KERNEL::Exception);
    DataArrayInt *sortCellsInMEDFileFrmt() throw(INTERP_KERNEL::Exception);
    DataArrayInt *getRenumArrForMEDFileFrmt() const throw(INTERP_KERNEL::Exception);
    DataArrayInt *convertCellArrayPerGeoType(const DataArrayInt *da) const throw(INTERP_KERNEL::Exception);
    DataArrayInt *computeFetchedNodeIds() const throw(INTERP_KERNEL::Exception);
    DataArrayInt *zipConnectivityTraducer(int compType, int startCellId=0) throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *buildDescendingConnectivity(DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *revDesc, DataArrayInt *revDescIndx) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *buildDescendingConnectivity2(DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *revDesc, DataArrayInt *revDescIndx) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *explode3DMeshTo1D(DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *revDesc, DataArrayInt *revDescIndx) const throw(INTERP_KERNEL::Exception);
    void orientCorrectlyPolyhedrons() throw(INTERP_KERNEL::Exception);
    bool isPresenceOfQuadratic() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *buildDirectionVectorField() const throw(INTERP_KERNEL::Exception);
    bool isContiguous1D() const throw(INTERP_KERNEL::Exception);
    void tessellate2D(double eps) throw(INTERP_KERNEL::Exception);
    void tessellate2DCurve(double eps) throw(INTERP_KERNEL::Exception);
    void convertQuadraticCellsToLinear() throw(INTERP_KERNEL::Exception);
    DataArrayInt *convertLinearCellsToQuadratic(int conversionType=0) throw(INTERP_KERNEL::Exception);
    void convertDegeneratedCells() throw(INTERP_KERNEL::Exception);
    bool areOnlySimplexCells() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getEdgeRatioField() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getAspectRatioField() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getWarpField() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getSkewField() const throw(INTERP_KERNEL::Exception);
    DataArrayInt *convexEnvelop2D() throw(INTERP_KERNEL::Exception);
    std::string cppRepr() const throw(INTERP_KERNEL::Exception);
    DataArrayInt *findAndCorrectBadOriented3DExtrudedCells() throw(INTERP_KERNEL::Exception);
    DataArrayInt *findAndCorrectBadOriented3DCells() throw(INTERP_KERNEL::Exception);
    static MEDCouplingUMesh *Build0DMeshFromCoords(DataArrayDouble *da) throw(INTERP_KERNEL::Exception);
    static MEDCouplingUMesh *MergeUMeshes(const MEDCouplingUMesh *mesh1, const MEDCouplingUMesh *mesh2) throw(INTERP_KERNEL::Exception);
    static MEDCouplingUMesh *MergeUMeshesOnSameCoords(const MEDCouplingUMesh *mesh1, const MEDCouplingUMesh *mesh2) throw(INTERP_KERNEL::Exception);
    static DataArrayInt *ComputeSpreadZoneGradually(const DataArrayInt *arrIn, const DataArrayInt *arrIndxIn) throw(INTERP_KERNEL::Exception);
    static DataArrayInt *ComputeRangesFromTypeDistribution(const std::vector<int>& code) throw(INTERP_KERNEL::Exception);
    %extend {
      MEDCouplingUMesh() throw(INTERP_KERNEL::Exception)
      {
        return MEDCouplingUMesh::New();
      }
      
      MEDCouplingUMesh(const char *meshName, int meshDim) throw(INTERP_KERNEL::Exception)
      {
        return MEDCouplingUMesh::New(meshName,meshDim);
      }
      
      std::string __str__() const throw(INTERP_KERNEL::Exception)
      {
        return self->simpleRepr();
      }
      
      std::string __repr__() const throw(INTERP_KERNEL::Exception)
      {
        std::ostringstream oss;
        self->reprQuickOverview(oss);
        return oss.str();
      }
      
      MEDCouplingUMeshCellIterator *__iter__() throw(INTERP_KERNEL::Exception)
      {
        return self->cellIterator();
      }

      MEDCouplingPointSet *__getitem__(PyObject *listOrDataArrI) throw(INTERP_KERNEL::Exception)
      {
        int sw;
        int singleVal;
        std::vector<int> multiVal;
        std::pair<int, std::pair<int,int> > slic;
        ParaMEDMEM::DataArrayInt *daIntTyypp=0;
        int nbc=self->getNumberOfCells();
        convertObjToPossibleCpp2(listOrDataArrI,nbc,sw,singleVal,multiVal,slic,daIntTyypp);
        switch(sw)
          {
          case 1:
            {
              if(singleVal>=nbc)
                {
                  std::ostringstream oss;
                  oss << "Requesting for cell id " << singleVal << " having only " << nbc << " cells !";
                  throw INTERP_KERNEL::Exception(oss.str().c_str());
                }
              if(singleVal>=0)
                return self->buildPartOfMySelf(&singleVal,&singleVal+1,true);
              else
                {
                  if(nbc+singleVal>0)
                    {
                      int tmp=nbc+singleVal;
                      return self->buildPartOfMySelf(&tmp,&tmp+1,true);
                    }
                  else
                    {
                      std::ostringstream oss;
                      oss << "Requesting for cell id " << singleVal << " having only " << nbc << " cells !";
                      throw INTERP_KERNEL::Exception(oss.str().c_str());
                    }
                }
            }
          case 2:
            {
              return static_cast<MEDCouplingUMesh *>(self->buildPartOfMySelf(&multiVal[0],&multiVal[0]+multiVal.size(),true));
            }
          case 3:
            {
              return self->buildPartOfMySelf2(slic.first,slic.second.first,slic.second.second,true);
            }
          case 4:
            {
              if(!daIntTyypp)
                throw INTERP_KERNEL::Exception("MEDCouplingUMesh::__getitem__ : null instance has been given in input !");
              daIntTyypp->checkAllocated();
              return self->buildPartOfMySelf(daIntTyypp->begin(),daIntTyypp->end(),true);
            }
          default:
            throw INTERP_KERNEL::Exception("MEDCouplingUMesh::__getitem__ : unrecognized type in input ! Possibilities are : int, list or tuple of int DataArrayInt instance !");
          }
      }
      
      void setPartOfMySelf(PyObject *li, const MEDCouplingUMesh& otherOnSameCoordsThanThis) throw(INTERP_KERNEL::Exception)
      {
        int sw;
        int singleVal;
        std::vector<int> multiVal;
        std::pair<int, std::pair<int,int> > slic;
        ParaMEDMEM::DataArrayInt *daIntTyypp=0;
        int nbc=self->getNumberOfCells();
        convertObjToPossibleCpp2(li,nbc,sw,singleVal,multiVal,slic,daIntTyypp);
        switch(sw)
          {
          case 1:
            {
              if(singleVal>=nbc)
                {
                  std::ostringstream oss;
                  oss << "Requesting for cell id " << singleVal << " having only " << nbc << " cells !";
                  throw INTERP_KERNEL::Exception(oss.str().c_str());
                }
              if(singleVal>=0)
                {
                  self->setPartOfMySelf(&singleVal,&singleVal+1,otherOnSameCoordsThanThis);
                  break;
                }
              else
                {
                  if(nbc+singleVal>0)
                    {
                      int tmp=nbc+singleVal;
                      self->setPartOfMySelf(&tmp,&tmp+1,otherOnSameCoordsThanThis);
                      break;
                    }
                  else
                    {
                      std::ostringstream oss;
                      oss << "Requesting for cell id " << singleVal << " having only " << nbc << " cells !";
                      throw INTERP_KERNEL::Exception(oss.str().c_str());
                    }
                }
            }
          case 2:
            {
              self->setPartOfMySelf(&multiVal[0],&multiVal[0]+multiVal.size(),otherOnSameCoordsThanThis);
              break;
            }
          case 4:
            {
              if(!daIntTyypp)
                throw INTERP_KERNEL::Exception("MEDCouplingUMesh::setPartOfMySelf : null instance has been given in input !");
              daIntTyypp->checkAllocated();
              self->setPartOfMySelf(daIntTyypp->begin(),daIntTyypp->end(),otherOnSameCoordsThanThis);
              break;
            }
          default:
            throw INTERP_KERNEL::Exception("MEDCouplingUMesh::setPartOfMySelf : unrecognized type in input ! Possibilities are : int, list or tuple of int DataArrayInt instance !");
          }
      }

      void __setitem__(PyObject *li, const MEDCouplingUMesh& otherOnSameCoordsThanThis) throw(INTERP_KERNEL::Exception)
      {
        int sw;
        int singleVal;
        std::vector<int> multiVal;
        std::pair<int, std::pair<int,int> > slic;
        ParaMEDMEM::DataArrayInt *daIntTyypp=0;
        int nbc=self->getNumberOfCells();
        convertObjToPossibleCpp2(li,nbc,sw,singleVal,multiVal,slic,daIntTyypp);
        switch(sw)
          {
          case 1:
            {
              if(singleVal>=nbc)
                {
                  std::ostringstream oss;
                  oss << "Requesting for cell id " << singleVal << " having only " << nbc << " cells !";
                  throw INTERP_KERNEL::Exception(oss.str().c_str());
                }
              if(singleVal>=0)
                {
                  self->setPartOfMySelf(&singleVal,&singleVal+1,otherOnSameCoordsThanThis);
                  break;
                }
              else
                {
                  if(nbc+singleVal>0)
                    {
                      int tmp=nbc+singleVal;
                      self->setPartOfMySelf(&tmp,&tmp+1,otherOnSameCoordsThanThis);
                      break;
                    }
                  else
                    {
                      std::ostringstream oss;
                      oss << "Requesting for cell id " << singleVal << " having only " << nbc << " cells !";
                      throw INTERP_KERNEL::Exception(oss.str().c_str());
                    }
                }
            }
          case 2:
            {
              self->setPartOfMySelf(&multiVal[0],&multiVal[0]+multiVal.size(),otherOnSameCoordsThanThis);
              break;
            }
          case 3:
            {
              self->setPartOfMySelf2(slic.first,slic.second.first,slic.second.second,otherOnSameCoordsThanThis);
              break;
            }
          case 4:
            {
              if(!daIntTyypp)
                throw INTERP_KERNEL::Exception("MEDCouplingUMesh::__setitem__ : null instance has been given in input !");
              daIntTyypp->checkAllocated();
              self->setPartOfMySelf(daIntTyypp->begin(),daIntTyypp->end(),otherOnSameCoordsThanThis);
              break;
            }
          default:
            throw INTERP_KERNEL::Exception("MEDCouplingUMesh::__setitem__ : unrecognized type in input ! Possibilities are : int, list or tuple of int, slice, DataArrayInt instance !");
          }
      }

      void insertNextCell(INTERP_KERNEL::NormalizedCellType type, int size, PyObject *li) throw(INTERP_KERNEL::Exception)
      {
        int szArr,sw,iTypppArr;
        std::vector<int> stdvecTyyppArr;
        const int *tmp=convertObjToPossibleCpp1_Safe(li,sw,szArr,iTypppArr,stdvecTyyppArr);
        if(size>szArr)
          {
            std::ostringstream oss; oss << "Wrap of MEDCouplingUMesh::insertNextCell : request of connectivity with length " << size << " whereas the length of input is " << szArr << " !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
        self->insertNextCell(type,size,tmp);
      }

      void insertNextCell(INTERP_KERNEL::NormalizedCellType type, PyObject *li) throw(INTERP_KERNEL::Exception)
      {
        int szArr,sw,iTypppArr;
        std::vector<int> stdvecTyyppArr;
        const int *tmp=convertObjToPossibleCpp1_Safe(li,sw,szArr,iTypppArr,stdvecTyyppArr);
        self->insertNextCell(type,szArr,tmp);
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
      
      static PyObject *ComputeSpreadZoneGraduallyFromSeed(PyObject *seed, const DataArrayInt *arrIn, const DataArrayInt *arrIndxIn, int nbOfDepthPeeling=-1) throw(INTERP_KERNEL::Exception)
      {
        int szArr,sw,iTypppArr;
        std::vector<int> stdvecTyyppArr;
        const int *seedPtr=convertObjToPossibleCpp1_Safe(seed,sw,szArr,iTypppArr,stdvecTyyppArr);
        int nbOfDepthPeelingPerformed=0;
        DataArrayInt *ret0=MEDCouplingUMesh::ComputeSpreadZoneGraduallyFromSeed(seedPtr,seedPtr+szArr,arrIn,arrIndxIn,nbOfDepthPeeling,nbOfDepthPeelingPerformed);
        PyObject *res=PyTuple_New(2);
        PyTuple_SetItem(res,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(res,1,PyInt_FromLong(nbOfDepthPeelingPerformed));
        return res;
      }

      PyObject *findCommonCells(int compType, int startCellId=0) const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *v0=0,*v1=0;
        self->findCommonCells(compType,startCellId,v0,v1);
        PyObject *res = PyList_New(2);
        PyList_SetItem(res,0,SWIG_NewPointerObj(SWIG_as_voidptr(v0),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyList_SetItem(res,1,SWIG_NewPointerObj(SWIG_as_voidptr(v1),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return res;
      }

      static PyObject *FindCommonCellsAlg(int compType, int startCellId, const DataArrayInt *nodal, const DataArrayInt *nodalI, const DataArrayInt *revNodal, const DataArrayInt *revNodalI) throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *v0=0,*v1=0;
        MEDCouplingUMesh::FindCommonCellsAlg(compType,startCellId,nodal,nodalI,revNodal,revNodalI,v0,v1);
        PyObject *res = PyList_New(2);
        PyList_SetItem(res,0,SWIG_NewPointerObj(SWIG_as_voidptr(v0),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyList_SetItem(res,1,SWIG_NewPointerObj(SWIG_as_voidptr(v1),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return res;
      }
      
      PyObject *distanceToPoint(PyObject *point) const throw(INTERP_KERNEL::Exception)
      {
        double val;
        DataArrayDouble *a;
        DataArrayDoubleTuple *aa;
        std::vector<double> bb;
        int sw;
        int nbOfCompo=self->getSpaceDimension();
        const double *pt=convertObjToPossibleCpp5_Safe(point,sw,val,a,aa,bb,"Python wrap of MEDCouplingUMesh::distanceToPoint",1,nbOfCompo,true);
        //
        int cellId=-1;
        double ret0=self->distanceToPoint(pt,pt+nbOfCompo,cellId);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,PyFloat_FromDouble(ret0));
        PyTuple_SetItem(ret,1,PyInt_FromLong(cellId));
        return ret;
      }

      PyObject *distanceToPoints(const DataArrayDouble *pts) const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *ret1=0;
        DataArrayDouble *ret0=self->distanceToPoints(pts,ret1);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
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
      PyObject *checkButterflyCells(double eps=1e-12) throw(INTERP_KERNEL::Exception)
      {
        std::vector<int> cells;
        self->checkButterflyCells(cells,eps);
        DataArrayInt *ret=DataArrayInt::New();
        ret->alloc((int)cells.size(),1);
        std::copy(cells.begin(),cells.end(),ret->getPointer());
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 );
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

      PyObject *partitionBySpreadZone() const throw(INTERP_KERNEL::Exception)
      {
        std::vector<DataArrayInt *> retCpp=self->partitionBySpreadZone();
        int sz=retCpp.size();
        PyObject *ret=PyList_New(sz);
        for(int i=0;i<sz;i++)
          PyList_SetItem(ret,i,SWIG_NewPointerObj(SWIG_as_voidptr(retCpp[i]),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
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

      PyObject *findNodesToDuplicate(const MEDCouplingUMesh& otherDimM1OnSameCoords) const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *tmp0=0,*tmp1=0,*tmp2=0;
        self->findNodesToDuplicate(otherDimM1OnSameCoords,tmp0,tmp1,tmp2);
        PyObject *ret=PyTuple_New(3);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(tmp0),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(tmp1),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(tmp2),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *findCellIdsLyingOn(const MEDCouplingUMesh& otherDimM1OnSameCoords) const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *tmp0=0,*tmp1=0;
        self->findCellIdsLyingOn(otherDimM1OnSameCoords,tmp0,tmp1);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(tmp0),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(tmp1),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      void duplicateNodes(PyObject *li) throw(INTERP_KERNEL::Exception)
      {
        int sw;
        int singleVal;
        std::vector<int> multiVal;
        std::pair<int, std::pair<int,int> > slic;
        ParaMEDMEM::DataArrayInt *daIntTyypp=0;
        convertObjToPossibleCpp2(li,self->getNumberOfNodes(),sw,singleVal,multiVal,slic,daIntTyypp);
        switch(sw)
          {
          case 1:
            return self->duplicateNodes(&singleVal,&singleVal+1);
          case 2:
            return self->duplicateNodes(&multiVal[0],&multiVal[0]+multiVal.size());
          case 4:
            return self->duplicateNodes(daIntTyypp->begin(),daIntTyypp->end());
          default:
            throw INTERP_KERNEL::Exception("MEDCouplingUMesh::duplicateNodes : unrecognized type entered, expected list of int, tuple of int or DataArrayInt !");
          }
      }

      void duplicateNodesInConn(PyObject *li, int offset) throw(INTERP_KERNEL::Exception)
      {
        int sw;
        int singleVal;
        std::vector<int> multiVal;
        std::pair<int, std::pair<int,int> > slic;
        ParaMEDMEM::DataArrayInt *daIntTyypp=0;
        convertObjToPossibleCpp2(li,self->getNumberOfNodes(),sw,singleVal,multiVal,slic,daIntTyypp);
        switch(sw)
          {
          case 1:
            return self->duplicateNodesInConn(&singleVal,&singleVal+1,offset);
          case 2:
            return self->duplicateNodesInConn(&multiVal[0],&multiVal[0]+multiVal.size(),offset);
          case 4:
            return self->duplicateNodesInConn(daIntTyypp->begin(),daIntTyypp->end(),offset);
          default:
            throw INTERP_KERNEL::Exception("MEDCouplingUMesh::duplicateNodesInConn : unrecognized type entered, expected list of int, tuple of int or DataArrayInt !");
          }
      }

      void renumberNodesInConn(PyObject *li) throw(INTERP_KERNEL::Exception)
      {
        void *da=0;
        int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 | 0 );
        if (!SWIG_IsOK(res1))
          {
            int size;
            INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
            self->renumberNodesInConn(tmp);
          }
        else
          {
            DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
            if(!da2)
              throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
            da2->checkAllocated();
            self->renumberNodesInConn(da2->getConstPointer());
          }
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

      static PyObject *AggregateSortedByTypeMeshesOnSameCoords(PyObject *ms) throw(INTERP_KERNEL::Exception)
      {
        std::vector<const ParaMEDMEM::MEDCouplingUMesh *> meshes;
        convertFromPyObjVectorOfObj<const ParaMEDMEM::MEDCouplingUMesh *>(ms,SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh,"MEDCouplingUMesh",meshes);
        DataArrayInt *ret1=0,*ret2=0;
        MEDCouplingUMesh *ret0=MEDCouplingUMesh::AggregateSortedByTypeMeshesOnSameCoords(meshes,ret1,ret2);
        PyObject *ret=PyTuple_New(3);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(ret2),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      static PyObject *MergeUMeshesOnSameCoords(PyObject *ms) throw(INTERP_KERNEL::Exception)
      {
        std::vector<const ParaMEDMEM::MEDCouplingUMesh *> meshes;
        convertFromPyObjVectorOfObj<const ParaMEDMEM::MEDCouplingUMesh *>(ms,SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh,"MEDCouplingUMesh",meshes);
        MEDCouplingUMesh *ret=MEDCouplingUMesh::MergeUMeshesOnSameCoords(meshes);
        return convertMesh(ret, SWIG_POINTER_OWN | 0 );
      }

      static PyObject *FuseUMeshesOnSameCoords(PyObject *ms, int compType) throw(INTERP_KERNEL::Exception)
      {
        int sz;
        std::vector<const MEDCouplingUMesh *> meshes;
        convertFromPyObjVectorOfObj<const ParaMEDMEM::MEDCouplingUMesh *>(ms,SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh,"MEDCouplingUMesh",meshes);
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

      static void PutUMeshesOnSameAggregatedCoords(PyObject *ms) throw(INTERP_KERNEL::Exception)
      {
        std::vector<MEDCouplingUMesh *> meshes;
        convertFromPyObjVectorOfObj<ParaMEDMEM::MEDCouplingUMesh *>(ms,SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh,"MEDCouplingUMesh",meshes);
        MEDCouplingUMesh::PutUMeshesOnSameAggregatedCoords(meshes);
      }

      static void MergeNodesOnUMeshesSharingSameCoords(PyObject *ms, double eps) throw(INTERP_KERNEL::Exception)
      {
        std::vector<MEDCouplingUMesh *> meshes;
        convertFromPyObjVectorOfObj<ParaMEDMEM::MEDCouplingUMesh *>(ms,SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh,"MEDCouplingUMesh",meshes);
        MEDCouplingUMesh::MergeNodesOnUMeshesSharingSameCoords(meshes,eps);
      }

      static bool RemoveIdsFromIndexedArrays(PyObject *li, DataArrayInt *arr, DataArrayInt *arrIndx, int offsetForRemoval=0) throw(INTERP_KERNEL::Exception)
      {
        int sw;
        int singleVal;
        std::vector<int> multiVal;
        std::pair<int, std::pair<int,int> > slic;
        ParaMEDMEM::DataArrayInt *daIntTyypp=0;
        if(!arrIndx)
          throw INTERP_KERNEL::Exception("MEDCouplingUMesh::RemoveIdsFromIndexedArrays : null pointer as arrIndex !");
        convertObjToPossibleCpp2(li,arrIndx->getNumberOfTuples()-1,sw,singleVal,multiVal,slic,daIntTyypp);
        switch(sw)
          {
          case 1:
            return MEDCouplingUMesh::RemoveIdsFromIndexedArrays(&singleVal,&singleVal+1,arr,arrIndx,offsetForRemoval);
          case 2:
            return MEDCouplingUMesh::RemoveIdsFromIndexedArrays(&multiVal[0],&multiVal[0]+multiVal.size(),arr,arrIndx,offsetForRemoval);
          case 4:
            return MEDCouplingUMesh::RemoveIdsFromIndexedArrays(daIntTyypp->begin(),daIntTyypp->end(),arr,arrIndx,offsetForRemoval);
          default:
            throw INTERP_KERNEL::Exception("MEDCouplingUMesh::RemoveIdsFromIndexedArrays : unrecognized type entered, expected list of int, tuple of int or DataArrayInt !");
          }
      }
      
      static PyObject *ExtractFromIndexedArrays(PyObject *li, const DataArrayInt *arrIn, const DataArrayInt *arrIndxIn) throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *arrOut=0,*arrIndexOut=0;
        int sw;
        int singleVal;
        std::vector<int> multiVal;
        std::pair<int, std::pair<int,int> > slic;
        ParaMEDMEM::DataArrayInt *daIntTyypp=0;
        if(!arrIndxIn)
          throw INTERP_KERNEL::Exception("MEDCouplingUMesh::ExtractFromIndexedArrays : null pointer as arrIndxIn !");
        convertObjToPossibleCpp2(li,arrIndxIn->getNumberOfTuples()-1,sw,singleVal,multiVal,slic,daIntTyypp);
        switch(sw)
          {
          case 1:
            {
              MEDCouplingUMesh::ExtractFromIndexedArrays(&singleVal,&singleVal+1,arrIn,arrIndxIn,arrOut,arrIndexOut);
              break;
            }
          case 2:
            {
              MEDCouplingUMesh::ExtractFromIndexedArrays(&multiVal[0],&multiVal[0]+multiVal.size(),arrIn,arrIndxIn,arrOut,arrIndexOut);
              break;
            }
          case 4:
            {
              MEDCouplingUMesh::ExtractFromIndexedArrays(daIntTyypp->begin(),daIntTyypp->end(),arrIn,arrIndxIn,arrOut,arrIndexOut);
              break;
            }
          default:
            throw INTERP_KERNEL::Exception("MEDCouplingUMesh::ExtractFromIndexedArrays : unrecognized type entered, expected list of int, tuple of int or DataArrayInt !");
          }
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(arrOut),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(arrIndexOut),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      static PyObject *SetPartOfIndexedArrays(PyObject *li,
                                              const DataArrayInt *arrIn, const DataArrayInt *arrIndxIn,
                                              const DataArrayInt *srcArr, const DataArrayInt *srcArrIndex) throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *arrOut=0,*arrIndexOut=0;
        int sw;
        int singleVal;
        std::vector<int> multiVal;
        std::pair<int, std::pair<int,int> > slic;
        ParaMEDMEM::DataArrayInt *daIntTyypp=0;
        if(!arrIndxIn)
          throw INTERP_KERNEL::Exception("MEDCouplingUMesh::SetPartOfIndexedArrays : null pointer as arrIndex !");
        convertObjToPossibleCpp2(li,arrIndxIn->getNumberOfTuples()-1,sw,singleVal,multiVal,slic,daIntTyypp);
        switch(sw)
          {
          case 1:
            {
              MEDCouplingUMesh::SetPartOfIndexedArrays(&singleVal,&singleVal+1,arrIn,arrIndxIn,srcArr,srcArrIndex,arrOut,arrIndexOut);
              break;
            }
          case 2:
            {
              MEDCouplingUMesh::SetPartOfIndexedArrays(&multiVal[0],&multiVal[0]+multiVal.size(),arrIn,arrIndxIn,srcArr,srcArrIndex,arrOut,arrIndexOut);
              break;
            }
          case 4:
            {
              MEDCouplingUMesh::SetPartOfIndexedArrays(daIntTyypp->begin(),daIntTyypp->end(),arrIn,arrIndxIn,srcArr,srcArrIndex,arrOut,arrIndexOut);
              break;
            }
          default:
            throw INTERP_KERNEL::Exception("MEDCouplingUMesh::SetPartOfIndexedArrays : unrecognized type entered, expected list of int, tuple of int or DataArrayInt !");
          }
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(arrOut),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(arrIndexOut),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      static void SetPartOfIndexedArraysSameIdx(PyObject *li, DataArrayInt *arrIn, const DataArrayInt *arrIndxIn,
                                                const DataArrayInt *srcArr, const DataArrayInt *srcArrIndex) throw(INTERP_KERNEL::Exception)
      {
        int sw;
        int singleVal;
        std::vector<int> multiVal;
        std::pair<int, std::pair<int,int> > slic;
        ParaMEDMEM::DataArrayInt *daIntTyypp=0;
        if(!arrIndxIn)
          throw INTERP_KERNEL::Exception("MEDCouplingUMesh::SetPartOfIndexedArraysSameIdx : null pointer as arrIndex !");
        convertObjToPossibleCpp2(li,arrIndxIn->getNumberOfTuples()-1,sw,singleVal,multiVal,slic,daIntTyypp);
        switch(sw)
          {
          case 1:
            {
              MEDCouplingUMesh::SetPartOfIndexedArraysSameIdx(&singleVal,&singleVal+1,arrIn,arrIndxIn,srcArr,srcArrIndex);
              break;
            }
          case 2:
            {
              MEDCouplingUMesh::SetPartOfIndexedArraysSameIdx(&multiVal[0],&multiVal[0]+multiVal.size(),arrIn,arrIndxIn,srcArr,srcArrIndex);
              break;
            }
          case 4:
            {
              MEDCouplingUMesh::SetPartOfIndexedArraysSameIdx(daIntTyypp->begin(),daIntTyypp->end(),arrIn,arrIndxIn,srcArr,srcArrIndex);
              break;
            }
          default:
            throw INTERP_KERNEL::Exception("MEDCouplingUMesh::SetPartOfIndexedArraysSameIdx : unrecognized type entered, expected list of int, tuple of int or DataArrayInt !");
          }
      }

      PyObject *are2DCellsNotCorrectlyOriented(PyObject *vec, bool polyOnly) const throw(INTERP_KERNEL::Exception)
      {
        double val;
        DataArrayDouble *a;
        DataArrayDoubleTuple *aa;
        std::vector<double> bb;
        int sw;
        int spaceDim=self->getSpaceDimension();
        const char msg[]="Python wrap of MEDCouplingUMesh::are2DCellsNotCorrectlyOriented : ";
        const double *v=convertObjToPossibleCpp5_Safe(vec,sw,val,a,aa,bb,msg,1,spaceDim,true);
        //
        std::vector<int> cells;
        self->are2DCellsNotCorrectlyOriented(v,polyOnly,cells);
        DataArrayInt *ret=DataArrayInt::New();
        ret->alloc((int)cells.size(),1);
        std::copy(cells.begin(),cells.end(),ret->getPointer());
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 );
      }

      void orientCorrectly2DCells(PyObject *vec, bool polyOnly) throw(INTERP_KERNEL::Exception)
      {
        double val;
        DataArrayDouble *a;
        DataArrayDoubleTuple *aa;
        std::vector<double> bb;
        int sw;
        int spaceDim=self->getSpaceDimension();
        const char msg[]="Python wrap of MEDCouplingUMesh::orientCorrectly2DCells : ";
        const double *v=convertObjToPossibleCpp5_Safe(vec,sw,val,a,aa,bb,msg,1,spaceDim,true);
        self->orientCorrectly2DCells(v,polyOnly);
      }
      
      PyObject *arePolyhedronsNotCorrectlyOriented() const throw(INTERP_KERNEL::Exception)
      {
        std::vector<int> cells;
        self->arePolyhedronsNotCorrectlyOriented(cells);
        DataArrayInt *ret=DataArrayInt::New();
        ret->alloc((int)cells.size(),1);
        std::copy(cells.begin(),cells.end(),ret->getPointer());
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 );
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
        convertFromPyObjVectorOfObj<const ParaMEDMEM::MEDCouplingUMesh *>(li,SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh,"MEDCouplingUMesh",tmp);
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

      PyObject *areCellsIncludedIn2(const MEDCouplingUMesh *other) const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *ret1;
        bool ret0=self->areCellsIncludedIn2(other,ret1);
        PyObject *ret=PyTuple_New(2);
        PyObject *ret0Py=ret0?Py_True:Py_False;
        Py_XINCREF(ret0Py);
        PyTuple_SetItem(ret,0,ret0Py);
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *explode3DMeshTo1D() const throw(INTERP_KERNEL::Exception)
      {
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> d0=DataArrayInt::New();
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> d1=DataArrayInt::New();
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> d2=DataArrayInt::New();
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> d3=DataArrayInt::New();
        MEDCouplingUMesh *m=self->explode3DMeshTo1D(d0,d1,d2,d3);
        PyObject *ret=PyTuple_New(5);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(m),SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(d0.retn()),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(d1.retn()),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,3,SWIG_NewPointerObj(SWIG_as_voidptr(d2.retn()),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,4,SWIG_NewPointerObj(SWIG_as_voidptr(d3.retn()),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
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
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(d0.retn()),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(d1.retn()),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,3,SWIG_NewPointerObj(SWIG_as_voidptr(d2.retn()),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,4,SWIG_NewPointerObj(SWIG_as_voidptr(d3.retn()),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *buildDescendingConnectivity2() const throw(INTERP_KERNEL::Exception)
      {
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> d0=DataArrayInt::New();
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> d1=DataArrayInt::New();
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> d2=DataArrayInt::New();
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> d3=DataArrayInt::New();
        MEDCouplingUMesh *m=self->buildDescendingConnectivity2(d0,d1,d2,d3);
        PyObject *ret=PyTuple_New(5);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(m),SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(d0.retn()),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(d1.retn()),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,3,SWIG_NewPointerObj(SWIG_as_voidptr(d2.retn()),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,4,SWIG_NewPointerObj(SWIG_as_voidptr(d3.retn()),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }
      
      PyObject *computeNeighborsOfCells() const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *neighbors=0,*neighborsIdx=0;
        self->computeNeighborsOfCells(neighbors,neighborsIdx);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(neighbors),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(neighborsIdx),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      static PyObject *ComputeNeighborsOfCellsAdv(const DataArrayInt *desc, const DataArrayInt *descI, const DataArrayInt *revDesc, const DataArrayInt *revDescI) throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *neighbors=0,*neighborsIdx=0;
        MEDCouplingUMesh::ComputeNeighborsOfCellsAdv(desc,descI,revDesc,revDescI,neighbors,neighborsIdx);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(neighbors),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(neighborsIdx),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
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
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(d0.retn()),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(d1.retn()),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,3,SWIG_NewPointerObj(SWIG_as_voidptr(d2),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,4,SWIG_NewPointerObj(SWIG_as_voidptr(d3),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,5,SWIG_NewPointerObj(SWIG_as_voidptr(d4),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,6,SWIG_NewPointerObj(SWIG_as_voidptr(dd5),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *getReverseNodalConnectivity() const throw(INTERP_KERNEL::Exception)
      {
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> d0=DataArrayInt::New();
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> d1=DataArrayInt::New();
        self->getReverseNodalConnectivity(d0,d1);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(d0.retn()),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(d1.retn()),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
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

      PyObject *getNodeIdsInUse() const throw(INTERP_KERNEL::Exception)
      {
        int ret1=-1;
        DataArrayInt *ret0=self->getNodeIdsInUse(ret1);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,PyInt_FromLong(ret1));
        return ret;
      }

      static PyObject *Intersect2DMeshes(const MEDCouplingUMesh *m1, const MEDCouplingUMesh *m2, double eps) throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *cellNb1=0,*cellNb2=0;
        MEDCouplingUMesh *mret=MEDCouplingUMesh::Intersect2DMeshes(m1,m2,eps,cellNb1,cellNb2);
        PyObject *ret=PyTuple_New(3);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(mret),SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(cellNb1),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(cellNb2),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *buildSlice3D(PyObject *origin, PyObject *vec, double eps) const throw(INTERP_KERNEL::Exception)
      {
        int spaceDim=self->getSpaceDimension();
        if(spaceDim!=3)
          throw INTERP_KERNEL::Exception("Python wrap of MEDCouplingUMesh::buildSlice3D : works only for spaceDim 3 !");
        double val,val2;
        DataArrayDouble *a,*a2;
        DataArrayDoubleTuple *aa,*aa2;
        std::vector<double> bb,bb2;
        int sw;
        const char msg[]="Python wrap of MEDCouplingUMesh::buildSlice3D : 1st paramater for origin.";
        const char msg2[]="Python wrap of MEDCouplingUMesh::buildSlice3D : 2nd paramater for vector.";
        const double *orig=convertObjToPossibleCpp5_Safe(origin,sw,val,a,aa,bb,msg,1,spaceDim,true);
        const double *vect=convertObjToPossibleCpp5_Safe(vec,sw,val2,a2,aa2,bb2,msg2,1,spaceDim,true);
        //
        DataArrayInt *cellIds=0;
        MEDCouplingUMesh *ret0=self->buildSlice3D(orig,vect,eps,cellIds);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(cellIds),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *buildSlice3DSurf(PyObject *origin, PyObject *vec, double eps) const throw(INTERP_KERNEL::Exception)
      {
        int spaceDim=self->getSpaceDimension();
        if(spaceDim!=3)
          throw INTERP_KERNEL::Exception("Python wrap of MEDCouplingUMesh::buildSlice3DSurf : works only for spaceDim 3 !");
        double val,val2;
        DataArrayDouble *a,*a2;
        DataArrayDoubleTuple *aa,*aa2;
        std::vector<double> bb,bb2;
        int sw;
        const char msg[]="Python wrap of MEDCouplingUMesh::buildSlice3DSurf : 1st paramater for origin.";
        const char msg2[]="Python wrap of MEDCouplingUMesh::buildSlice3DSurf : 2nd paramater for vector.";
        const double *orig=convertObjToPossibleCpp5_Safe(origin,sw,val,a,aa,bb,msg,1,spaceDim,true);
        const double *vect=convertObjToPossibleCpp5_Safe(vec,sw,val2,a2,aa2,bb2,msg2,1,spaceDim,true);
        //
        DataArrayInt *cellIds=0;
        MEDCouplingUMesh *ret0=self->buildSlice3DSurf(orig,vect,eps,cellIds);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(cellIds),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      DataArrayInt *getCellIdsCrossingPlane(PyObject *origin, PyObject *vec, double eps) const throw(INTERP_KERNEL::Exception)
      {
        int spaceDim=self->getSpaceDimension();
        if(spaceDim!=3)
          throw INTERP_KERNEL::Exception("Python wrap of MEDCouplingUMesh::getCellIdsCrossingPlane : works only for spaceDim 3 !");
        double val,val2;
        DataArrayDouble *a,*a2;
        DataArrayDoubleTuple *aa,*aa2;
        std::vector<double> bb,bb2;
        int sw;
        const char msg[]="Python wrap of MEDCouplingUMesh::getCellIdsCrossingPlane : 1st paramater for origin.";
        const char msg2[]="Python wrap of MEDCouplingUMesh::getCellIdsCrossingPlane : 2nd paramater for vector.";
        const double *orig=convertObjToPossibleCpp5_Safe(origin,sw,val,a,aa,bb,msg,1,spaceDim,true);
        const double *vect=convertObjToPossibleCpp5_Safe(vec,sw,val2,a2,aa2,bb2,msg2,1,spaceDim,true);
        return self->getCellIdsCrossingPlane(orig,vect,eps);
      }

      void convertToPolyTypes(PyObject *li) throw(INTERP_KERNEL::Exception)
      {
        int sw;
        int pos1;
        std::vector<int> pos2;
        DataArrayInt *pos3=0;
        DataArrayIntTuple *pos4=0;
        convertObjToPossibleCpp1(li,sw,pos1,pos2,pos3,pos4);
        switch(sw)
          {
          case 1:
            {
              self->convertToPolyTypes(&pos1,&pos1+1);
              return;
            }
          case 2:
            {
              if(pos2.empty())
                return;
              self->convertToPolyTypes(&pos2[0],&pos2[0]+pos2.size());
              return ;
            }
          case 3:
            {
              self->convertToPolyTypes(pos3->begin(),pos3->end());
              return ;
            }
          default:
            throw INTERP_KERNEL::Exception("MEDCouplingUMesh::convertToPolyTypes : unexpected input array type recognized !");
          }
      }
    }
    void convertAllToPoly();
    void convertExtrudedPolyhedra() throw(INTERP_KERNEL::Exception);
    bool unPolyze() throw(INTERP_KERNEL::Exception);
    void simplifyPolyhedra(double eps) throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *buildSpreadZonesWithPoly() const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *buildExtrudedMesh(const MEDCouplingUMesh *mesh1D, int policy) throw(INTERP_KERNEL::Exception);
  };

  class MEDCouplingExtrudedMesh : public ParaMEDMEM::MEDCouplingMesh
  {
  public:
    static MEDCouplingExtrudedMesh *New(const MEDCouplingUMesh *mesh3D, const MEDCouplingUMesh *mesh2D, int cell2DId) throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *build3DUnstructuredMesh() const throw(INTERP_KERNEL::Exception);
    void updateTime() const throw(INTERP_KERNEL::Exception);
    %extend {
      MEDCouplingExtrudedMesh(const MEDCouplingUMesh *mesh3D, const MEDCouplingUMesh *mesh2D, int cell2DId) throw(INTERP_KERNEL::Exception)
      {
        return MEDCouplingExtrudedMesh::New(mesh3D,mesh2D,cell2DId);
      }
      
      std::string __str__() const throw(INTERP_KERNEL::Exception)
      {
        return self->simpleRepr();
      }

      std::string __repr__() const throw(INTERP_KERNEL::Exception)
      {
        std::ostringstream oss;
        self->reprQuickOverview(oss);
        return oss.str();
      }
      
      PyObject *getMesh2D() const throw(INTERP_KERNEL::Exception)
      {
        MEDCouplingUMesh *ret=self->getMesh2D();
        if(ret)
          ret->incrRef();
        return convertMesh(ret, SWIG_POINTER_OWN | 0 );
      }
      PyObject *getMesh1D() const throw(INTERP_KERNEL::Exception)
      {
        MEDCouplingUMesh *ret=self->getMesh1D();
        if(ret)
          ret->incrRef();
        return convertMesh(ret, SWIG_POINTER_OWN | 0 );
      }
      PyObject *getMesh3DIds() const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *ret=self->getMesh3DIds();
        if(ret)
          ret->incrRef();
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 );
      } 
    }
  };

  class MEDCouplingStructuredMesh : public ParaMEDMEM::MEDCouplingMesh
  {
  public:
    void updateTime() const throw(INTERP_KERNEL::Exception);
    int getCellIdFromPos(int i, int j, int k) const throw(INTERP_KERNEL::Exception);
    int getNodeIdFromPos(int i, int j, int k) const throw(INTERP_KERNEL::Exception);
  };

  class MEDCouplingCMesh : public ParaMEDMEM::MEDCouplingStructuredMesh
  {
  public:
    static MEDCouplingCMesh *New();
    static MEDCouplingCMesh *New(const char *meshName);
    MEDCouplingCMesh *clone(bool recDeepCpy) const;
    void setCoords(const DataArrayDouble *coordsX,
                   const DataArrayDouble *coordsY=0,
                   const DataArrayDouble *coordsZ=0) throw(INTERP_KERNEL::Exception);
    void setCoordsAt(int i, const DataArrayDouble *arr) throw(INTERP_KERNEL::Exception);
    %extend {
      MEDCouplingCMesh()
      {
        return MEDCouplingCMesh::New();
      }
      MEDCouplingCMesh(const char *meshName)
      {
        return MEDCouplingCMesh::New(meshName);
      }
      std::string __str__() const throw(INTERP_KERNEL::Exception)
      {
        return self->simpleRepr();
      }
      std::string __repr__() const throw(INTERP_KERNEL::Exception)
      {
        std::ostringstream oss;
        self->reprQuickOverview(oss);
        return oss.str();
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

  class MEDCouplingCurveLinearMesh : public ParaMEDMEM::MEDCouplingStructuredMesh
  {
  public:
    static MEDCouplingCurveLinearMesh *New();
    static MEDCouplingCurveLinearMesh *New(const char *meshName);
    MEDCouplingCurveLinearMesh *clone(bool recDeepCpy) const;
    void setCoords(const DataArrayDouble *coords) throw(INTERP_KERNEL::Exception);
    std::vector<int> getNodeGridStructure() const throw(INTERP_KERNEL::Exception);
    %extend {
      MEDCouplingCurveLinearMesh()
      {
        return MEDCouplingCurveLinearMesh::New();
      }
      MEDCouplingCurveLinearMesh(const char *meshName)
      {
        return MEDCouplingCurveLinearMesh::New(meshName);
      }
      std::string __str__() const throw(INTERP_KERNEL::Exception) 
      {
        return self->simpleRepr();
      }
      std::string __repr__() const throw(INTERP_KERNEL::Exception)
      {
        std::ostringstream oss;
        self->reprQuickOverview(oss);
        return oss.str();
      }
      DataArrayDouble *getCoords() throw(INTERP_KERNEL::Exception)
      {
        DataArrayDouble *ret=self->getCoords();
        if(ret)
          ret->incrRef();
        return ret;
      }
      void setNodeGridStructure(PyObject *gridStruct) throw(INTERP_KERNEL::Exception)
      {
        int szArr,sw,iTypppArr;
        std::vector<int> stdvecTyyppArr;
        const int *tmp=convertObjToPossibleCpp1_Safe(gridStruct,sw,szArr,iTypppArr,stdvecTyyppArr);
        self->setNodeGridStructure(tmp,tmp+szArr);
      }
    }
  };
  
}

%extend ParaMEDMEM::MEDCouplingFieldDiscretization
{
  MEDCouplingFieldDiscretization *clonePart(PyObject *li)
  {
    int sz=0,sw=-1,val1=-1;
    std::vector<int> val2;
    const int *inp=convertObjToPossibleCpp1_Safe(li,sw,sz,val1,val2);
    return self->clonePart(inp,inp+sz);
  }

  PyObject *buildSubMeshDataRange(const MEDCouplingMesh *mesh, int beginCellIds, int endCellIds, int stepCellIds, int& beginOut, int& endOut, int& stepOut, DataArrayInt *&di) const throw(INTERP_KERNEL::Exception)
  {
    DataArrayInt *ret1=0;
    int bb,ee,ss;
    MEDCouplingMesh *ret0=self->buildSubMeshDataRange(mesh,begin,end,step,bb,ee,ss,ret1);
    PyObject *res=PyTuple_New(2);
    PyTuple_SetItem(res,0,convertMesh(ret0, SWIG_POINTER_OWN | 0 ));
    if(ret1)
      PyTuple_SetItem(res,1,SWIG_NewPointerObj((void*)ret1,SWIGTYPE_p_ParaMEDMEM__DataArrayInt,SWIG_POINTER_OWN | 0));
    else
      {
        PyObject *res1=PySlice_New(PyInt_FromLong(bb),PyInt_FromLong(ee),PyInt_FromLong(ss));
        PyTuple_SetItem(res,1,res1);
      }
    return res;
  }
  
  PyObject *computeMeshRestrictionFromTupleIds(const MEDCouplingMesh *mesh, PyObject *tupleIds) const throw(INTERP_KERNEL::Exception)
  {
    std::vector<int> vVal; int iVal=-1;
    int sz=-1,sw=0;
    const int *tupleIdsBg=convertObjToPossibleCpp1_Safe(tupleIds,sw,sz,iVal,vVal);
    if(sw==0)
      throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretization::computeMeshRestrictionFromTupleIds : none parameter in input !");
    DataArrayInt *ret0=0,*ret1=0;
    self->computeMeshRestrictionFromTupleIds(mesh,tupleIdsBg,tupleIdsBg+sz,ret0,ret1);
    PyObject *pyRet=PyTuple_New(2);
    PyTuple_SetItem(pyRet,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
    PyTuple_SetItem(pyRet,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
    return pyRet;
  }
}

%extend ParaMEDMEM::MEDCouplingFieldDiscretizationP0
{
  PyObject *computeMeshRestrictionFromTupleIds(const MEDCouplingMesh *mesh, PyObject *tupleIds) const throw(INTERP_KERNEL::Exception)
  { return ParaMEDMEM_MEDCouplingFieldDiscretization_computeMeshRestrictionFromTupleIds__SWIG_1(self,mesh,tupleIds); }
}

%extend ParaMEDMEM::MEDCouplingFieldDiscretizationOnNodes
{
  PyObject *computeMeshRestrictionFromTupleIds(const MEDCouplingMesh *mesh, PyObject *tupleIds) const throw(INTERP_KERNEL::Exception)
  { return ParaMEDMEM_MEDCouplingFieldDiscretization_computeMeshRestrictionFromTupleIds__SWIG_1(self,mesh,tupleIds); }
}

%extend ParaMEDMEM::MEDCouplingFieldDiscretizationGauss
{
  PyObject *computeMeshRestrictionFromTupleIds(const MEDCouplingMesh *mesh, PyObject *tupleIds) const throw(INTERP_KERNEL::Exception)
  { return ParaMEDMEM_MEDCouplingFieldDiscretization_computeMeshRestrictionFromTupleIds__SWIG_1(self,mesh,tupleIds); }
}

%extend ParaMEDMEM::MEDCouplingFieldDiscretizationGaussNE
{
  PyObject *computeMeshRestrictionFromTupleIds(const MEDCouplingMesh *mesh, PyObject *tupleIds) const throw(INTERP_KERNEL::Exception)
  { return ParaMEDMEM_MEDCouplingFieldDiscretization_computeMeshRestrictionFromTupleIds__SWIG_1(self,mesh,tupleIds); }
}

%extend ParaMEDMEM::MEDCouplingFieldDiscretizationPerCell
{
  PyObject *getArrayOfDiscIds() const
  {
    DataArrayInt *ret=const_cast<DataArrayInt *>(self->getArrayOfDiscIds());
    if(ret)
      ret->incrRef();
    return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 );
  }

  PyObject *splitIntoSingleGaussDicrPerCellType() const throw(INTERP_KERNEL::Exception)
  {
    std::vector<int> ret1;
    std::vector<DataArrayInt *> ret0=self->splitIntoSingleGaussDicrPerCellType(ret1);
    std::size_t sz=ret0.size();
    PyObject *pyRet=PyTuple_New(2);
    PyObject *pyRet0=PyList_New((int)sz);
    PyObject *pyRet1=PyList_New((int)sz);
    for(std::size_t i=0;i<sz;i++)
      {
        PyList_SetItem(pyRet0,i,SWIG_NewPointerObj(SWIG_as_voidptr(ret0[i]),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyList_SetItem(pyRet1,i,PyInt_FromLong(ret1[i]));
      }
    PyTuple_SetItem(pyRet,0,pyRet0);
    PyTuple_SetItem(pyRet,1,pyRet1);
    return pyRet;
  }
}

%extend ParaMEDMEM::MEDCouplingFieldDiscretizationKriging
{
  PyObject *computeVectorOfCoefficients(const MEDCouplingMesh *mesh, const DataArrayDouble *arr) const
  {
    int ret1;
    DataArrayDouble *ret0=self->computeVectorOfCoefficients(mesh,arr,ret1);
    PyObject *ret=PyTuple_New(2);
    PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
    PyTuple_SetItem(ret,1,PyInt_FromLong(ret1));
    return ret;
  }
}

namespace ParaMEDMEM
{
  class MEDCouplingField : public ParaMEDMEM::RefCountObject, public ParaMEDMEM::TimeLabel
  {
  public:
    virtual void checkCoherency() const throw(INTERP_KERNEL::Exception);
    virtual bool areCompatibleForMerge(const MEDCouplingField *other) const throw(INTERP_KERNEL::Exception);
    virtual bool isEqual(const MEDCouplingField *other, double meshPrec, double valsPrec) const throw(INTERP_KERNEL::Exception);
    virtual bool isEqualWithoutConsideringStr(const MEDCouplingField *other, double meshPrec, double valsPrec) const throw(INTERP_KERNEL::Exception);
    virtual void copyTinyStringsFrom(const MEDCouplingField *other) throw(INTERP_KERNEL::Exception);
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
    int getNumberOfTuplesExpected() const throw(INTERP_KERNEL::Exception);
    int getNumberOfMeshPlacesExpected() const throw(INTERP_KERNEL::Exception);
    void setGaussLocalizationOnType(INTERP_KERNEL::NormalizedCellType type, const std::vector<double>& refCoo,
                                    const std::vector<double>& gsCoo, const std::vector<double>& wg) throw(INTERP_KERNEL::Exception);
    void clearGaussLocalizations() throw(INTERP_KERNEL::Exception);
    MEDCouplingGaussLocalization& getGaussLocalization(int locId) throw(INTERP_KERNEL::Exception);
    int getNbOfGaussLocalization() const throw(INTERP_KERNEL::Exception);
    int getGaussLocalizationIdOfOneCell(int cellId) const throw(INTERP_KERNEL::Exception);
    const MEDCouplingGaussLocalization& getGaussLocalization(int locId) const throw(INTERP_KERNEL::Exception);
    int getGaussLocalizationIdOfOneType(INTERP_KERNEL::NormalizedCellType type) const throw(INTERP_KERNEL::Exception);
    void setDiscretization(MEDCouplingFieldDiscretization *newDisc);
    %extend {
      PyObject *getMesh() const throw(INTERP_KERNEL::Exception)
      {
        MEDCouplingMesh *ret1=const_cast<MEDCouplingMesh *>(self->getMesh());
        if(ret1)
          ret1->incrRef();
        return convertMesh(ret1,SWIG_POINTER_OWN | 0 );
      }

      PyObject *getDiscretization() throw(INTERP_KERNEL::Exception)
      {
        MEDCouplingFieldDiscretization *ret=self->getDiscretization();
        if(ret)
          ret->incrRef();
        return convertFieldDiscretization(ret,SWIG_POINTER_OWN | 0 );
      }

      PyObject *getGaussLocalizationIdsOfOneType(INTERP_KERNEL::NormalizedCellType type) const throw(INTERP_KERNEL::Exception)
      {
        std::set<int> ret=self->getGaussLocalizationIdsOfOneType(type);
        return convertIntArrToPyList3(ret);
      }

      PyObject *isEqualIfNotWhy(const MEDCouplingField *other, double meshPrec, double valsPrec) const throw(INTERP_KERNEL::Exception)
      {
        std::string ret1;
        bool ret0=self->isEqualIfNotWhy(other,meshPrec,valsPrec,ret1);
        PyObject *ret=PyTuple_New(2);
        PyObject *ret0Py=ret0?Py_True:Py_False;
        Py_XINCREF(ret0Py);
        PyTuple_SetItem(ret,0,ret0Py);
        PyTuple_SetItem(ret,1,PyString_FromString(ret1.c_str()));
        return ret;
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

      PyObject *buildSubMeshDataRange(int begin, int end, int step) const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *ret1=0;
        int bb,ee,ss;
        MEDCouplingMesh *ret0=self->buildSubMeshDataRange(begin,end,step,bb,ee,ss,ret1);
        PyObject *res=PyTuple_New(2);
        PyTuple_SetItem(res,0,convertMesh(ret0, SWIG_POINTER_OWN | 0 ));
        if(ret1)
          PyTuple_SetItem(res,1,SWIG_NewPointerObj((void*)ret1,SWIGTYPE_p_ParaMEDMEM__DataArrayInt,SWIG_POINTER_OWN | 0));
        else
          {
            PyObject *res1=PySlice_New(PyInt_FromLong(bb),PyInt_FromLong(ee),PyInt_FromLong(ss));
            PyTuple_SetItem(res,1,res1);
          }
        return res;
      }

      DataArrayInt *computeTupleIdsToSelectFromCellIds(PyObject *li) const
      {
        int sw;
        int pos1;
        std::vector<int> pos2;
        DataArrayInt *pos3=0;
        DataArrayIntTuple *pos4=0;
        convertObjToPossibleCpp1(li,sw,pos1,pos2,pos3,pos4);
        switch(sw)
          {
          case 1:
            {
              return self->computeTupleIdsToSelectFromCellIds(&pos1,&pos1+1);
            }
          case 2:
            {
              return self->computeTupleIdsToSelectFromCellIds(&pos2[0],&pos2[0]+pos2.size());
            }
          case 3:
            {
              return self->computeTupleIdsToSelectFromCellIds(pos3->begin(),pos3->end());
            }
          default:
            throw INTERP_KERNEL::Exception("MEDCouplingField::computeTupleIdsToSelectFromCellIds : unexpected input array type recognized !");
          }
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
        DataArrayInt *ret=DataArrayInt::New();
        ret->alloc((int)tmp.size(),1);
        std::copy(tmp.begin(),tmp.end(),ret->getPointer());
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 );
      }
    }
  };
  
  class MEDCouplingFieldTemplate : public ParaMEDMEM::MEDCouplingField
  {
  public:
    static MEDCouplingFieldTemplate *New(const MEDCouplingFieldDouble& f) throw(INTERP_KERNEL::Exception);
    static MEDCouplingFieldTemplate *New(TypeOfField type);
    std::string simpleRepr() const throw(INTERP_KERNEL::Exception);
    std::string advancedRepr() const throw(INTERP_KERNEL::Exception);
    void updateTime() const;
    %extend
       {
         MEDCouplingFieldTemplate(const MEDCouplingFieldDouble& f) throw(INTERP_KERNEL::Exception)
         {
           return MEDCouplingFieldTemplate::New(f);
         }
         
         MEDCouplingFieldTemplate(TypeOfField type) throw(INTERP_KERNEL::Exception)
         {
           return MEDCouplingFieldTemplate::New(type);
         }
         
         std::string __str__() const throw(INTERP_KERNEL::Exception)
         {
           return self->simpleRepr();
         }
         
         std::string __repr__() const throw(INTERP_KERNEL::Exception)
         {
           std::ostringstream oss;
           self->reprQuickOverview(oss);
           return oss.str();
         }
       }
  };
  
  class MEDCouplingFieldDouble : public ParaMEDMEM::MEDCouplingField
  {
  public:
    static MEDCouplingFieldDouble *New(TypeOfField type, TypeOfTimeDiscretization td=ONE_TIME);
    static MEDCouplingFieldDouble *New(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td=ONE_TIME);
    void setTimeUnit(const char *unit);
    const char *getTimeUnit() const;
    void synchronizeTimeWithSupport() throw(INTERP_KERNEL::Exception);
    void copyTinyAttrFrom(const MEDCouplingFieldDouble *other) throw(INTERP_KERNEL::Exception);
    void copyAllTinyAttrFrom(const MEDCouplingFieldDouble *other) throw(INTERP_KERNEL::Exception);
    std::string simpleRepr() const throw(INTERP_KERNEL::Exception);
    std::string advancedRepr() const throw(INTERP_KERNEL::Exception);
    void writeVTK(const char *fileName) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *clone(bool recDeepCpy) const;
    MEDCouplingFieldDouble *cloneWithMesh(bool recDeepCpy) const;
    MEDCouplingFieldDouble *deepCpy() const;
    MEDCouplingFieldDouble *buildNewTimeReprFromThis(TypeOfTimeDiscretization td, bool deepCpy) const throw(INTERP_KERNEL::Exception);
    TypeOfTimeDiscretization getTimeDiscretization() const throw(INTERP_KERNEL::Exception);
    double getIJ(int tupleId, int compoId) const throw(INTERP_KERNEL::Exception);
    double getIJK(int cellId, int nodeIdInCell, int compoId) const throw(INTERP_KERNEL::Exception);
    void synchronizeTimeWithMesh() throw(INTERP_KERNEL::Exception);
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
    void changeUnderlyingMesh(const MEDCouplingMesh *other, int levOfCheck, double precOnMesh, double eps=1e-15) throw(INTERP_KERNEL::Exception);
    void substractInPlaceDM(const MEDCouplingFieldDouble *f, int levOfCheck, double precOnMesh, double eps=1e-15) throw(INTERP_KERNEL::Exception);
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
    void fillFromAnalytic2(int nbOfComp, const char *func) throw(INTERP_KERNEL::Exception);
    void fillFromAnalytic3(int nbOfComp, const std::vector<std::string>& varsOrder, const char *func) throw(INTERP_KERNEL::Exception);
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
    double norm2() const throw(INTERP_KERNEL::Exception);
    double normMax() const throw(INTERP_KERNEL::Exception);
    //do not put a default value to isWAbs because confusion in python with overloaded getWeightedAverageValue method
    double getWeightedAverageValue(int compId, bool isWAbs) const throw(INTERP_KERNEL::Exception);
    double integral(int compId, bool isWAbs) const throw(INTERP_KERNEL::Exception);
    double normL1(int compId) const throw(INTERP_KERNEL::Exception);
    double normL2(int compId) const throw(INTERP_KERNEL::Exception);
    DataArrayInt *getIdsInRange(double vmin, double vmax) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *buildSubPartRange(int begin, int end, int step) const throw(INTERP_KERNEL::Exception);
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
    MEDCouplingFieldDouble *negate() const throw(INTERP_KERNEL::Exception);
    %extend {
      MEDCouplingFieldDouble(TypeOfField type, TypeOfTimeDiscretization td=ONE_TIME)
      {
        return MEDCouplingFieldDouble::New(type,td);
      }

      MEDCouplingFieldDouble(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td=ONE_TIME)
      {
        return MEDCouplingFieldDouble::New(ft,td);
      }

      std::string __str__() const throw(INTERP_KERNEL::Exception)
      {
        return self->simpleRepr();
      }

      std::string __repr__() const throw(INTERP_KERNEL::Exception)
      {
        std::ostringstream oss;
        self->reprQuickOverview(oss);
        return oss.str();
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
        convertFromPyObjVectorOfObj<const DataArrayDouble *>(ls,SWIGTYPE_p_ParaMEDMEM__DataArrayDouble,"DataArrayDouble",tmp);
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
        double val;
        DataArrayDouble *a;
        DataArrayDoubleTuple *aa;
        std::vector<double> bb;
        int sw;
        const MEDCouplingMesh *mesh=self->getMesh();
        if(!mesh)
          throw INTERP_KERNEL::Exception("Python wrap of MEDCouplingFieldDouble::getValueOn : no underlying mesh !");
        int spaceDim=mesh->getSpaceDimension();
        const char msg[]="Python wrap of MEDCouplingFieldDouble::getValueOn : ";
        const double *spaceLoc=convertObjToPossibleCpp5_Safe(sl,sw,val,a,aa,bb,msg,1,spaceDim,true);
        //
        int sz=self->getNumberOfComponents();
        INTERP_KERNEL::AutoPtr<double> res=new double[sz];
        self->getValueOn(spaceLoc,res);
        return convertDblArrToPyList(res,sz);
      }

       PyObject *getValueOnPos(int i, int j, int k) const throw(INTERP_KERNEL::Exception)
       {
         int sz=self->getNumberOfComponents();
         INTERP_KERNEL::AutoPtr<double> res=new double[sz];
         self->getValueOnPos(i,j,k,res);
         return convertDblArrToPyList(res,sz);
       }

      DataArrayDouble *getValueOnMulti(PyObject *li) const throw(INTERP_KERNEL::Exception)
      {
        void *da=0;
        int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayDouble, 0 |  0 );
        if (!SWIG_IsOK(res1))
          {
            int size;
            INTERP_KERNEL::AutoCPtr<double> tmp=convertPyToNewDblArr2(li,&size);
            const MEDCouplingMesh *mesh=self->getMesh();
            if(!mesh)
              throw INTERP_KERNEL::Exception("Python wrap MEDCouplingFieldDouble::getValueOnMulti : lying on a null mesh !");
            int spaceDim=mesh->getSpaceDimension();
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
            const MEDCouplingMesh *mesh=self->getMesh();
            if(!mesh)
              throw INTERP_KERNEL::Exception("Python wrap MEDCouplingFieldDouble::getValueOnMulti : lying on a null mesh !");
            if(nbOfCompo!=mesh->getSpaceDimension())
              {
                throw INTERP_KERNEL::Exception("Invalid DataArrayDouble nb of components ! Expected same as self.getMesh().getSpaceDimension() !");
              }
            return self->getValueOnMulti(da2->getConstPointer(),size);
          }
      }

      PyObject *getValueOn(PyObject *sl, double time) const throw(INTERP_KERNEL::Exception)
      {
        double val;
        DataArrayDouble *a;
        DataArrayDoubleTuple *aa;
        std::vector<double> bb;
        int sw;
        const MEDCouplingMesh *mesh=self->getMesh();
        if(!mesh)
          throw INTERP_KERNEL::Exception("Python wrap of MEDCouplingFieldDouble::getValueOn : no underlying mesh !");
        int spaceDim=mesh->getSpaceDimension();
        const char msg[]="Python wrap of MEDCouplingFieldDouble::getValueOn : ";
        const double *spaceLoc=convertObjToPossibleCpp5_Safe(sl,sw,val,a,aa,bb,msg,1,spaceDim,true);
        //
        //
        int sz=self->getNumberOfComponents();
        INTERP_KERNEL::AutoPtr<double> res=new double[sz];
        self->getValueOn(spaceLoc,time,res);
        return convertDblArrToPyList(res,sz);
      }

      void setValues(PyObject *li, PyObject *nbOfTuples=0, PyObject *nbOfComp=0) throw(INTERP_KERNEL::Exception)
      {
        if(self->getArray()!=0)
          ParaMEDMEM_DataArrayDouble_setValues__SWIG_0(self->getArray(),li,nbOfTuples,nbOfComp);
        else
          {
            MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr=DataArrayDouble::New();
            ParaMEDMEM_DataArrayDouble_setValues__SWIG_0(arr,li,nbOfTuples,nbOfComp);
            self->setArray(arr);
          }
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
        return convertDblArrToPyList(tmp,sz);
      }
      PyObject *integral(bool isWAbs) const throw(INTERP_KERNEL::Exception)
      {
        int sz=self->getNumberOfComponents();
        INTERP_KERNEL::AutoPtr<double> tmp=new double[sz];
        self->integral(isWAbs,tmp);
        return convertDblArrToPyList(tmp,sz);
      }
      PyObject *getWeightedAverageValue(bool isWAbs=true) const throw(INTERP_KERNEL::Exception)
      {
        int sz=self->getNumberOfComponents();
        INTERP_KERNEL::AutoPtr<double> tmp=new double[sz];
        self->getWeightedAverageValue(tmp,isWAbs);
        return convertDblArrToPyList(tmp,sz);
      }
      PyObject *normL1() const throw(INTERP_KERNEL::Exception)
      {
        int sz=self->getNumberOfComponents();
        INTERP_KERNEL::AutoPtr<double> tmp=new double[sz];
        self->normL1(tmp);
        return convertDblArrToPyList(tmp,sz);
      }
      PyObject *normL2() const throw(INTERP_KERNEL::Exception)
      {
        int sz=self->getNumberOfComponents();
        INTERP_KERNEL::AutoPtr<double> tmp=new double[sz];
        self->normL2(tmp);
        return convertDblArrToPyList(tmp,sz);
      }
      void renumberCells(PyObject *li, bool check=true) throw(INTERP_KERNEL::Exception)
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
      void renumberNodes(PyObject *li, double eps=1e-15) throw(INTERP_KERNEL::Exception)
      {
        void *da=0;
        int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_ParaMEDMEM__DataArrayInt, 0 |  0 );
        if (!SWIG_IsOK(res1))
          {
            int size;
            INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
            self->renumberNodes(tmp,eps);
          }
        else
          {
            DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
            if(!da2)
              throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
            da2->checkAllocated();
            self->renumberNodes(da2->getConstPointer(),eps);
          }
      }

      MEDCouplingFieldDouble *buildSubPart(PyObject *li) const throw(INTERP_KERNEL::Exception)
      {
        int sw;
        int singleVal;
        std::vector<int> multiVal;
        std::pair<int, std::pair<int,int> > slic;
        ParaMEDMEM::DataArrayInt *daIntTyypp=0;
        const MEDCouplingMesh *mesh=self->getMesh();
        if(!mesh)
          throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::buildSubPart : field lies on a null mesh !");
        int nbc=mesh->getNumberOfCells();
        convertObjToPossibleCpp2(li,nbc,sw,singleVal,multiVal,slic,daIntTyypp);
        switch(sw)
          {
          case 1:
            {
              if(singleVal>=nbc)
                {
                  std::ostringstream oss;
                  oss << "Requesting for cell id " << singleVal << " having only " << nbc << " cells !";
                  throw INTERP_KERNEL::Exception(oss.str().c_str());
                }
              if(singleVal>=0)
                return self->buildSubPart(&singleVal,&singleVal+1);
              else
                {
                  if(nbc+singleVal>0)
                    {
                      int tmp=nbc+singleVal;
                      return self->buildSubPart(&tmp,&tmp+1);
                    }
                  else
                    {
                      std::ostringstream oss;
                      oss << "Requesting for cell id " << singleVal << " having only " << nbc << " cells !";
                      throw INTERP_KERNEL::Exception(oss.str().c_str());
                    }
                }
            }
          case 2:
            {
              return self->buildSubPart(&multiVal[0],&multiVal[0]+multiVal.size());
            }
          case 3:
            {
              return self->buildSubPartRange(slic.first,slic.second.first,slic.second.second);
            }
          case 4:
            {
              if(!daIntTyypp)
                throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::buildSubPart : null instance has been given in input !");
              daIntTyypp->checkAllocated();
              return self->buildSubPart(daIntTyypp->begin(),daIntTyypp->end());
            }
          default:
            throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::buildSubPart : unrecognized type in input ! Possibilities are : int, list or tuple of int DataArrayInt instance !");
          }
      }

      MEDCouplingFieldDouble *__getitem__(PyObject *li) const throw(INTERP_KERNEL::Exception)
      {
        const char msg[]="MEDCouplingFieldDouble::__getitem__ : invalid call  Available API are : \n-myField[dataArrayInt]\n-myField[slice]\n-myField[pythonListOfCellIds]\n-myField[integer]\n-myField[dataArrayInt,1]\n-myField[slice,1]\n-myField[pythonListOfCellIds,1]\n-myField[integer,1]\n";
        if(PyTuple_Check(li))
          {
            Py_ssize_t sz=PyTuple_Size(li);
            if(sz!=2)
              throw INTERP_KERNEL::Exception(msg);
            PyObject *elt0=PyTuple_GetItem(li,0),*elt1=PyTuple_GetItem(li,1);
            int sw;
            int singleVal;
            std::vector<int> multiVal;
            std::pair<int, std::pair<int,int> > slic;
            ParaMEDMEM::DataArrayInt *daIntTyypp=0;
            if(!self->getArray())
              throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::__getitem__ : no array set on field to deduce number of components !");
            try
              { convertObjToPossibleCpp2(elt1,self->getArray()->getNumberOfComponents(),sw,singleVal,multiVal,slic,daIntTyypp); }
            catch(INTERP_KERNEL::Exception& e)
              { std::ostringstream oss; oss << "MEDCouplingFieldDouble::__getitem__ : invalid type in 2nd parameter (compo) !" << e.what(); throw INTERP_KERNEL::Exception(oss.str().c_str()); }
            MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret0=ParaMEDMEM_MEDCouplingFieldDouble_buildSubPart(self,elt0);
            DataArrayDouble *ret0Arr=ret0->getArray();
            if(!ret0Arr)
              throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::__getitem__ : no array exists to apply restriction on component on it !");
            switch(sw)
              {
              case 1:
                {
                  std::vector<int> v2(1,singleVal);
                  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> aarr=static_cast<DataArrayDouble *>(ret0Arr->keepSelectedComponents(v2));
                  ret0->setArray(aarr);
                  return ret0.retn();
                }
              case 2:
                {
                  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> aarr=static_cast<DataArrayDouble *>(ret0Arr->keepSelectedComponents(multiVal));
                  ret0->setArray(aarr);
                  return ret0.retn();
                }
              case 3:
                {
                  int nbOfComp=DataArray::GetNumberOfItemGivenBESRelative(slic.first,slic.second.first,slic.second.second,"MEDCouplingFieldDouble::__getitem__ : invalid range in 2nd parameter (components) !");
                  std::vector<int> v2(nbOfComp);
                  for(int i=0;i<nbOfComp;i++)
                    v2[i]=slic.first+i*slic.second.second;
                  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> aarr=static_cast<DataArrayDouble *>(ret0Arr->keepSelectedComponents(v2));
                  ret0->setArray(aarr);
                  return ret0.retn();
                }
              default:
                throw INTERP_KERNEL::Exception(msg);
              }
            
          }
        else
          return ParaMEDMEM_MEDCouplingFieldDouble_buildSubPart(self,li);
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

      MEDCouplingFieldDouble *extractSlice3D(PyObject *origin, PyObject *vec, double eps) const throw(INTERP_KERNEL::Exception)
      {
        double val,val2;
        DataArrayDouble *a,*a2;
        DataArrayDoubleTuple *aa,*aa2;
        std::vector<double> bb,bb2;
        int sw;
        int spaceDim=3;
        const char msg[]="Python wrap of MEDCouplingFieldDouble::extractSlice3D : 1st paramater for origin.";
        const char msg2[]="Python wrap of MEDCouplingFieldDouble::extractSlice3D : 2nd paramater for vector.";
        const double *orig=convertObjToPossibleCpp5_Safe(origin,sw,val,a,aa,bb,msg,1,spaceDim,true);
        const double *vect=convertObjToPossibleCpp5_Safe(vec,sw,val2,a2,aa2,bb2,msg2,1,spaceDim,true);
        //
        return self->extractSlice3D(orig,vect,eps);
      }

      MEDCouplingFieldDouble *__add__(PyObject *obj) throw(INTERP_KERNEL::Exception)
      {
        return ParaMEDMEM_MEDCouplingFieldDouble___add__Impl(self,obj);
      }

      MEDCouplingFieldDouble *__radd__(PyObject *obj) throw(INTERP_KERNEL::Exception)
      {
        return ParaMEDMEM_MEDCouplingFieldDouble___radd__Impl(self,obj);
      }

      MEDCouplingFieldDouble *__sub__(PyObject *obj) throw(INTERP_KERNEL::Exception)
      {
        const char msg[]="Unexpected situation in MEDCouplingFieldDouble.__sub__ ! Expecting a not null MEDCouplingFieldDouble or DataArrayDouble or DataArrayDoubleTuple instance, or a list of double, or a double.";
        const char msg2[]="in MEDCouplingFieldDouble.__sub__ : self field has no Array of values set !";
        void *argp;
        //
        if(SWIG_IsOK(SWIG_ConvertPtr(obj,&argp,SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble,0|0)))
          {
            MEDCouplingFieldDouble *other=reinterpret_cast< ParaMEDMEM::MEDCouplingFieldDouble * >(argp);
            if(other)
              return (*self)-(*other);
            else
              throw INTERP_KERNEL::Exception(msg);
          }
        //
        double val;
        DataArrayDouble *a;
        DataArrayDoubleTuple *aa;
        std::vector<double> bb;
        int sw;
        convertObjToPossibleCpp5(obj,sw,val,a,aa,bb);
        switch(sw)
          {
          case 1:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=self->getArray()->deepCpy();
              ret->applyLin(1.,-val);
              MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(ret);
              return ret2.retn();
            }
          case 2:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=DataArrayDouble::Substract(self->getArray(),a);
              MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(ret);
              return ret2.retn();
            }
          case 3:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> aaa=aa->buildDADouble(1,self->getNumberOfComponents());
              MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=DataArrayDouble::Substract(self->getArray(),aaa);
              MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(ret);
              return ret2.retn();
            }
          case 4:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> aaa=DataArrayDouble::New(); aaa->useArray(&bb[0],false,CPP_DEALLOC,1,(int)bb.size());
              MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=DataArrayDouble::Substract(self->getArray(),aaa);
              MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(ret);
              return ret2.retn();
            }
          default:
            { throw INTERP_KERNEL::Exception(msg); }
          }
      }

      MEDCouplingFieldDouble *__rsub__(PyObject *obj) throw(INTERP_KERNEL::Exception)
      {
        return ParaMEDMEM_MEDCouplingFieldDouble___rsub__Impl(self,obj);
      }

      MEDCouplingFieldDouble *__mul__(PyObject *obj) throw(INTERP_KERNEL::Exception)
      {
        return ParaMEDMEM_MEDCouplingFieldDouble___mul__Impl(self,obj);
      }

      MEDCouplingFieldDouble *__rmul__(PyObject *obj) throw(INTERP_KERNEL::Exception)
      {
        return ParaMEDMEM_MEDCouplingFieldDouble___rmul__Impl(self,obj);
      }

      MEDCouplingFieldDouble *__div__(PyObject *obj) throw(INTERP_KERNEL::Exception)
      {
        const char msg[]="Unexpected situation in MEDCouplingFieldDouble.__div__ ! Expecting a not null MEDCouplingFieldDouble or DataArrayDouble or DataArrayDoubleTuple instance, or a list of double, or a double.";
        const char msg2[]="in MEDCouplingFieldDouble.__div__ : self field has no Array of values set !";
        void *argp;
        //
        if(SWIG_IsOK(SWIG_ConvertPtr(obj,&argp,SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble,0|0)))
          {
            MEDCouplingFieldDouble *other=reinterpret_cast< ParaMEDMEM::MEDCouplingFieldDouble * >(argp);
            if(other)
              return (*self)/(*other);
            else
              throw INTERP_KERNEL::Exception(msg);
          }
        //
        double val;
        DataArrayDouble *a;
        DataArrayDoubleTuple *aa;
        std::vector<double> bb;
        int sw;
        convertObjToPossibleCpp5(obj,sw,val,a,aa,bb);
        switch(sw)
          {
          case 1:
            {
              if(val==0.)
                throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble.__div__ : trying to divide by zero !");
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=self->getArray()->deepCpy();
              ret->applyLin(1./val,0);
              MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(ret);
              return ret2.retn();
            }
          case 2:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=DataArrayDouble::Divide(self->getArray(),a);
              MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(ret);
              return ret2.retn();
            }
          case 3:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> aaa=aa->buildDADouble(1,self->getNumberOfComponents());
              MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=DataArrayDouble::Divide(self->getArray(),aaa);
              MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(ret);
              return ret2.retn();
            }
          case 4:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> aaa=DataArrayDouble::New(); aaa->useArray(&bb[0],false,CPP_DEALLOC,1,(int)bb.size());
              MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=DataArrayDouble::Divide(self->getArray(),aaa);
              MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(ret);
              return ret2.retn();
            }
          default:
            { throw INTERP_KERNEL::Exception(msg); }
          }
      }

      MEDCouplingFieldDouble *__rdiv__(PyObject *obj) throw(INTERP_KERNEL::Exception)
      {
        return ParaMEDMEM_MEDCouplingFieldDouble___rdiv__Impl(self,obj);
      }

      MEDCouplingFieldDouble *__pow__(PyObject *obj) throw(INTERP_KERNEL::Exception)
      {
        const char msg[]="Unexpected situation in MEDCouplingFieldDouble.__pow__ ! Expecting a not null MEDCouplingFieldDouble or DataArrayDouble or DataArrayDoubleTuple instance, or a list of double, or a double.";
        const char msg2[]="in MEDCouplingFieldDouble.__pow__ : self field has no Array of values set !";
        void *argp;
        //
        if(SWIG_IsOK(SWIG_ConvertPtr(obj,&argp,SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble,0|0)))
          {
            MEDCouplingFieldDouble *other=reinterpret_cast< ParaMEDMEM::MEDCouplingFieldDouble * >(argp);
            if(other)
              return (*self)^(*other);
            else
              throw INTERP_KERNEL::Exception(msg);
          }
        //
        double val;
        DataArrayDouble *a;
        DataArrayDoubleTuple *aa;
        std::vector<double> bb;
        int sw;
        convertObjToPossibleCpp5(obj,sw,val,a,aa,bb);
        switch(sw)
          {
          case 1:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=self->getArray()->deepCpy();
              ret->applyPow(val);
              MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(ret);
              return ret2.retn();
            }
          case 2:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=DataArrayDouble::Pow(self->getArray(),a);
              MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(ret);
              return ret2.retn();
            }
          case 3:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> aaa=aa->buildDADouble(1,self->getNumberOfComponents());
              MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=DataArrayDouble::Pow(self->getArray(),aaa);
              MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(ret);
              return ret2.retn();
            }
          case 4:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> aaa=DataArrayDouble::New(); aaa->useArray(&bb[0],false,CPP_DEALLOC,1,(int)bb.size());
              MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=DataArrayDouble::Pow(self->getArray(),aaa);
              MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(ret);
              return ret2.retn();
            }
          default:
            { throw INTERP_KERNEL::Exception(msg); }
          }
      }

      MEDCouplingFieldDouble *__neg__() const throw(INTERP_KERNEL::Exception)
      {
        return self->negate();
      }

      PyObject *___iadd___(PyObject *trueSelf, PyObject *obj) throw(INTERP_KERNEL::Exception)
      {
        const char msg[]="Unexpected situation in MEDCouplingFieldDouble.__iadd__ ! Expecting a not null MEDCouplingFieldDouble or DataArrayDouble or DataArrayDoubleTuple instance, or a list of double, or a double.";
        const char msg2[]="in MEDCouplingFieldDouble.__iadd__ : self field has no Array of values set !";
        void *argp;
        //
        if(SWIG_IsOK(SWIG_ConvertPtr(obj,&argp,SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble,0|0)))
          {
            MEDCouplingFieldDouble *other=reinterpret_cast< ParaMEDMEM::MEDCouplingFieldDouble * >(argp);
            if(other)
              {
                *self+=*other;
                Py_XINCREF(trueSelf);
                return trueSelf;
              }
            else
              throw INTERP_KERNEL::Exception(msg);
          }
        //
        double val;
        DataArrayDouble *a;
        DataArrayDoubleTuple *aa;
        std::vector<double> bb;
        int sw;
        convertObjToPossibleCpp5(obj,sw,val,a,aa,bb);
        switch(sw)
          {
          case 1:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              self->getArray()->applyLin(1.,val);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 2:
            {
              MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(a);
              *self+=*ret2;
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 3:
            {
              MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> aaa=aa->buildDADouble(1,self->getNumberOfComponents());
              MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(aaa);
              *self+=*ret2;
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 4:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> aaa=DataArrayDouble::New(); aaa->useArray(&bb[0],false,CPP_DEALLOC,1,(int)bb.size());
              self->getArray()->addEqual(aaa);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          default:
            { throw INTERP_KERNEL::Exception(msg); }
          }
      }

      PyObject *___isub___(PyObject *trueSelf, PyObject *obj) throw(INTERP_KERNEL::Exception)
      {
        const char msg[]="Unexpected situation in MEDCouplingFieldDouble.__isub__ ! Expecting a not null MEDCouplingFieldDouble or DataArrayDouble or DataArrayDoubleTuple instance, or a list of double, or a double.";
        const char msg2[]="in MEDCouplingFieldDouble.__isub__ : self field has no Array of values set !";
        void *argp;
        //
        if(SWIG_IsOK(SWIG_ConvertPtr(obj,&argp,SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble,0|0)))
          {
            MEDCouplingFieldDouble *other=reinterpret_cast< ParaMEDMEM::MEDCouplingFieldDouble * >(argp);
            if(other)
              {
                *self-=*other;
                Py_XINCREF(trueSelf);
                return trueSelf;
              }
            else
              throw INTERP_KERNEL::Exception(msg);
          }
        //
        double val;
        DataArrayDouble *a;
        DataArrayDoubleTuple *aa;
        std::vector<double> bb;
        int sw;
        convertObjToPossibleCpp5(obj,sw,val,a,aa,bb);
        switch(sw)
          {
          case 1:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              self->getArray()->applyLin(1.,-val);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 2:
            {
              MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(a);
              *self-=*ret2;
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 3:
            {
              MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> aaa=aa->buildDADouble(1,self->getNumberOfComponents());
              MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(aaa);
              *self-=*ret2;
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 4:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> aaa=DataArrayDouble::New(); aaa->useArray(&bb[0],false,CPP_DEALLOC,1,(int)bb.size());
              self->getArray()->substractEqual(aaa);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          default:
            { throw INTERP_KERNEL::Exception(msg); }
          }
      }

      PyObject *___imul___(PyObject *trueSelf, PyObject *obj) throw(INTERP_KERNEL::Exception)
      {
        const char msg[]="Unexpected situation in MEDCouplingFieldDouble.__imul__ ! Expecting a not null MEDCouplingFieldDouble or DataArrayDouble or DataArrayDoubleTuple instance, or a list of double, or a double.";
        const char msg2[]="in MEDCouplingFieldDouble.__imul__ : self field has no Array of values set !";
        void *argp;
        //
        if(SWIG_IsOK(SWIG_ConvertPtr(obj,&argp,SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble,0|0)))
          {
            MEDCouplingFieldDouble *other=reinterpret_cast< ParaMEDMEM::MEDCouplingFieldDouble * >(argp);
            if(other)
              {
                *self*=*other;
                Py_XINCREF(trueSelf);
                return trueSelf;
              }
            else
              throw INTERP_KERNEL::Exception(msg);
          }
        //
        double val;
        DataArrayDouble *a;
        DataArrayDoubleTuple *aa;
        std::vector<double> bb;
        int sw;
        convertObjToPossibleCpp5(obj,sw,val,a,aa,bb);
        switch(sw)
          {
          case 1:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              self->getArray()->applyLin(val,0);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 2:
            {
              MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(a);
              *self*=*ret2;
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 3:
            {
              MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> aaa=aa->buildDADouble(1,self->getNumberOfComponents());
              MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(aaa);
              *self*=*ret2;
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 4:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> aaa=DataArrayDouble::New(); aaa->useArray(&bb[0],false,CPP_DEALLOC,1,(int)bb.size());
              self->getArray()->multiplyEqual(aaa);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          default:
            { throw INTERP_KERNEL::Exception(msg); }
          }
      }

      PyObject *___idiv___(PyObject *trueSelf, PyObject *obj) throw(INTERP_KERNEL::Exception)
      {
        const char msg[]="Unexpected situation in MEDCouplingFieldDouble.__idiv__ ! Expecting a not null MEDCouplingFieldDouble or DataArrayDouble or DataArrayDoubleTuple instance, or a list of double, or a double.";
        const char msg2[]="in MEDCouplingFieldDouble.__idiv__ : self field has no Array of values set !";
        void *argp;
        //
        if(SWIG_IsOK(SWIG_ConvertPtr(obj,&argp,SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble,0|0)))
          {
            MEDCouplingFieldDouble *other=reinterpret_cast< ParaMEDMEM::MEDCouplingFieldDouble * >(argp);
            if(other)
              {
                *self/=*other;
                Py_XINCREF(trueSelf);
                return trueSelf;
              }
            else
              throw INTERP_KERNEL::Exception(msg);
          }
        //
        double val;
        DataArrayDouble *a;
        DataArrayDoubleTuple *aa;
        std::vector<double> bb;
        int sw;
        convertObjToPossibleCpp5(obj,sw,val,a,aa,bb);
        switch(sw)
          {
          case 1:
            {
              if(val==0.)
                throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble.__idiv__ : trying to divide by zero !");
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              self->getArray()->applyLin(1./val,0);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 2:
            {
              MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(a);
              *self/=*ret2;
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 3:
            {
              MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> aaa=aa->buildDADouble(1,self->getNumberOfComponents());
              MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(aaa);
              *self/=*ret2;
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 4:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> aaa=DataArrayDouble::New(); aaa->useArray(&bb[0],false,CPP_DEALLOC,1,(int)bb.size());
              self->getArray()->divideEqual(aaa);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          default:
            { throw INTERP_KERNEL::Exception(msg); }
          }
      }

      PyObject *___ipow___(PyObject *trueSelf, PyObject *obj) throw(INTERP_KERNEL::Exception)
      {
        const char msg[]="Unexpected situation in MEDCouplingFieldDouble.__ipow__ ! Expecting a not null MEDCouplingFieldDouble or DataArrayDouble or DataArrayDoubleTuple instance, or a list of double, or a double.";
        const char msg2[]="in MEDCouplingFieldDouble.__ipow__ : self field has no Array of values set !";
        void *argp;
        //
        if(SWIG_IsOK(SWIG_ConvertPtr(obj,&argp,SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble,0|0)))
          {
            MEDCouplingFieldDouble *other=reinterpret_cast< ParaMEDMEM::MEDCouplingFieldDouble * >(argp);
            if(other)
              {
                *self^=*other;
                Py_XINCREF(trueSelf);
                return trueSelf;
              }
            else
              throw INTERP_KERNEL::Exception(msg);
          }
        //
        double val;
        DataArrayDouble *a;
        DataArrayDoubleTuple *aa;
        std::vector<double> bb;
        int sw;
        convertObjToPossibleCpp5(obj,sw,val,a,aa,bb);
        switch(sw)
          {
          case 1:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              self->getArray()->applyPow(val);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 2:
            {
              MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(a);
              *self^=*ret2;
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 3:
            {
              MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> aaa=aa->buildDADouble(1,self->getNumberOfComponents());
              MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(aaa);
              *self^=*ret2;
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 4:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> aaa=DataArrayDouble::New(); aaa->useArray(&bb[0],false,CPP_DEALLOC,1,(int)bb.size());
              self->getArray()->powEqual(aaa);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          default:
            { throw INTERP_KERNEL::Exception(msg); }
          }
      }

      static MEDCouplingFieldDouble *MergeFields(PyObject *li) throw(INTERP_KERNEL::Exception)
      {
        std::vector<const MEDCouplingFieldDouble *> tmp;
        convertFromPyObjVectorOfObj<const ParaMEDMEM::MEDCouplingFieldDouble *>(li,SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble,"MEDCouplingFieldDouble",tmp);
        return MEDCouplingFieldDouble::MergeFields(tmp);
      }

      static void WriteVTK(const char *fileName, PyObject *li) throw(INTERP_KERNEL::Exception)
      {
        std::vector<const MEDCouplingFieldDouble *> tmp;
        convertFromPyObjVectorOfObj<const ParaMEDMEM::MEDCouplingFieldDouble *>(li,SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble,"MEDCouplingFieldDouble",tmp);
        MEDCouplingFieldDouble::WriteVTK(fileName,tmp);
      }
    }
  };

  class MEDCouplingMultiFields : public RefCountObject, public TimeLabel
  {
  public:
    int getNumberOfFields() const;
    MEDCouplingMultiFields *deepCpy() const;
    virtual std::string simpleRepr() const throw(INTERP_KERNEL::Exception);
    virtual std::string advancedRepr() const throw(INTERP_KERNEL::Exception);
    virtual bool isEqual(const MEDCouplingMultiFields *other, double meshPrec, double valsPrec) const;
    virtual bool isEqualWithoutConsideringStr(const MEDCouplingMultiFields *other, double meshPrec, double valsPrec) const;
    virtual void checkCoherency() const throw(INTERP_KERNEL::Exception);
    void updateTime() const throw(INTERP_KERNEL::Exception);
    %extend
       {
         std::string __str__() const throw(INTERP_KERNEL::Exception)
         {
           return self->simpleRepr();
         }
         static MEDCouplingMultiFields *New(PyObject *li) throw(INTERP_KERNEL::Exception)
         {
           std::vector<const ParaMEDMEM::MEDCouplingFieldDouble *> tmp;
           convertFromPyObjVectorOfObj<const ParaMEDMEM::MEDCouplingFieldDouble *>(li,SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble,"MEDCouplingFieldDouble",tmp);
           int sz=tmp.size();
           std::vector<MEDCouplingFieldDouble *> fs(sz);
           for(int i=0;i<sz;i++)
             fs[i]=const_cast<MEDCouplingFieldDouble *>(tmp[i]);
           return MEDCouplingMultiFields::New(fs);
         }
         MEDCouplingMultiFields(PyObject *li) throw(INTERP_KERNEL::Exception)
         {
           std::vector<const ParaMEDMEM::MEDCouplingFieldDouble *> tmp;
           convertFromPyObjVectorOfObj<const ParaMEDMEM::MEDCouplingFieldDouble *>(li,SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble,"MEDCouplingFieldDouble",tmp);
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
        std::string __str__() const throw(INTERP_KERNEL::Exception)
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
        MEDCouplingFieldOverTime(PyObject *li) throw(INTERP_KERNEL::Exception)
          {
            std::vector<const ParaMEDMEM::MEDCouplingFieldDouble *> tmp;
            convertFromPyObjVectorOfObj<const ParaMEDMEM::MEDCouplingFieldDouble *>(li,SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble,"MEDCouplingFieldDouble",tmp);
            int sz=tmp.size();
            std::vector<MEDCouplingFieldDouble *> fs(sz);
            for(int i=0;i<sz;i++)
              fs[i]=const_cast<MEDCouplingFieldDouble *>(tmp[i]);
            return MEDCouplingFieldOverTime::New(fs);
          }
        std::string __str__() const throw(INTERP_KERNEL::Exception)
          {
            return self->simpleRepr();
          }
        static MEDCouplingFieldOverTime *New(PyObject *li) throw(INTERP_KERNEL::Exception)
        {
          std::vector<const ParaMEDMEM::MEDCouplingFieldDouble *> tmp;
          convertFromPyObjVectorOfObj<const ParaMEDMEM::MEDCouplingFieldDouble *>(li,SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble,"MEDCouplingFieldDouble",tmp);
           int sz=tmp.size();
           std::vector<MEDCouplingFieldDouble *> fs(sz);
           for(int i=0;i<sz;i++)
             fs[i]=const_cast<MEDCouplingFieldDouble *>(tmp[i]);
           return MEDCouplingFieldOverTime::New(fs);
         }
      }
  };
}

%pythoncode %{
import os
__filename=os.environ.get('PYTHONSTARTUP')
if __filename and os.path.isfile(__filename):
  execfile(__filename)
  pass
%}
