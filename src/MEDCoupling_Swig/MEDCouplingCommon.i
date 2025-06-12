// Copyright (C) 2017-2025  CEA, EDF
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
// Author : Anthony Geay (EDF R&D)

#ifdef WITH_DOCSTRINGS
%include MEDCoupling_doc.i
#endif

%include std_vector.i
%include std_string.i

%{
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingMemArray.txx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingMappedExtrudedMesh.hxx"
#include "MEDCouplingCMesh.hxx"
#include "MEDCouplingIMesh.hxx"
#include "MEDCouplingMap.txx"
#include "MEDCouplingCurveLinearMesh.hxx"
#include "MEDCoupling1GTUMesh.hxx"
#include "MEDCouplingField.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingFieldInt32.hxx"
#include "MEDCouplingFieldInt64.hxx"
#include "MEDCouplingFieldFloat.hxx"
#include "MEDCouplingFieldTemplate.hxx"
#include "MEDCouplingGaussLocalization.hxx"
#include "MCAuto.hxx"
#include "MEDCouplingMultiFields.hxx"
#include "MEDCouplingFieldOverTime.hxx"
#include "MEDCouplingDefinitionTime.hxx"
#include "MEDCouplingFieldDiscretization.hxx"
#include "MEDCouplingFieldDiscretizationOnNodesFE.hxx"
#include "MEDCouplingCartesianAMRMesh.hxx"
#include "MEDCouplingAMRAttribute.hxx"
#include "MEDCouplingMatrix.hxx"
#include "MEDCouplingPartDefinition.hxx"
#include "MEDCouplingSkyLineArray.hxx"
#include "MEDCouplingTypemaps.i"

#include "InterpKernelAutoPtr.hxx"
#include "BoxSplittingOptions.hxx"

using namespace MEDCoupling;
using namespace INTERP_KERNEL;

%}

%template(dvec) std::vector<double>;
%template(svec) std::vector<std::string>;

//%include stdint.i

#ifndef MEDCOUPLING_USE_64BIT_IDS
//typedef std::int32_t mcIdType;
typedef int mcIdType;
typedef DataArrayInt32 DataArrayIdType;
%template(ivec) std::vector<int>;
%template(i64vec) std::vector<int64_t>;
#else
//typedef std::int64_t mcIdType;
typedef DataArrayInt64 DataArrayIdType;
#ifdef WIN32
%template(ivec) std::vector<long long>;
typedef long long mcIdType;
#else
%template(ivec) std::vector<long>;
typedef long int mcIdType;
#endif
%template(i32vec) std::vector<int>;
#endif
#ifdef WIN32
typedef long long mcPyPtrType;
#else
typedef long mcPyPtrType;
#endif

////////////////////
%typemap(out) MEDCoupling::MEDCouplingMesh*
{
  $result=convertMesh($1,$owner);
}

%typemap(out) MEDCouplingMesh*
{
  $result=convertMesh($1,$owner);
}
//$$$$$$$$$$$$$$$$$$

////////////////////
%typemap(out) MEDCoupling::MEDCouplingPointSet*
{
  $result=convertMesh($1,$owner);
}

%typemap(out) MEDCouplingPointSet*
{
  $result=convertMesh($1,$owner);
}
//$$$$$$$$$$$$$$$$$$

////////////////////
%typemap(out) MEDCouplingCartesianAMRPatchGen*
{
  $result=convertCartesianAMRPatch($1,$owner);
}
//$$$$$$$$$$$$$$$$$$

////////////////////
%typemap(out) MEDCouplingCartesianAMRMeshGen*
{
  $result=convertCartesianAMRMesh($1,$owner);
}
//$$$$$$$$$$$$$$$$$$

////////////////////
%typemap(out) MEDCouplingDataForGodFather*
{
  $result=convertDataForGodFather($1,$owner);
}
//$$$$$$$$$$$$$$$$$$

////////////////////
%typemap(out) MEDCoupling::MEDCoupling1GTUMesh*
{
  $result=convertMesh($1,$owner);
}

%typemap(out) MEDCoupling1GTUMesh*
{
  $result=convertMesh($1,$owner);
}
//$$$$$$$$$$$$$$$$$$

////////////////////
%typemap(out) MEDCoupling::MEDCouplingStructuredMesh*
{
  $result=convertMesh($1,$owner);
}

%typemap(out) MEDCouplingStructuredMesh*
{
  $result=convertMesh($1,$owner);
}
//$$$$$$$$$$$$$$$$$$

////////////////////
%typemap(out) MEDCoupling::MEDCouplingFieldDiscretization*
{
  $result=convertFieldDiscretization($1,$owner);
}

%typemap(out) MEDCouplingFieldDiscretization*
{
  $result=convertFieldDiscretization($1,$owner);
}
//$$$$$$$$$$$$$$$$$$

////////////////////
%typemap(out) MEDCoupling::MEDCouplingField*
{
  $result=convertField($1,$owner);
}

%typemap(out) MEDCouplingField*
{
  $result=convertField($1,$owner);
}
//$$$$$$$$$$$$$$$$$$

////////////////////
%typemap(out) MEDCoupling::MEDCouplingMultiFields*
{
  $result=convertMultiFields($1,$owner);
}

%typemap(out) MEDCouplingMultiFields*
{
  $result=convertMultiFields($1,$owner);
}
//$$$$$$$$$$$$$$$$$$

////////////////////
%typemap(out) MEDCoupling::PartDefinition*
{
  $result=convertPartDefinition($1,$owner);
}

%typemap(out) PartDefinition*
{
  $result=convertPartDefinition($1,$owner);
}
//$$$$$$$$$$$$$$$$$$

#ifdef WITH_NUMPY
%init %{ import_array(); %}
#endif

%init %{ initializeMe(); %}

%feature("autodoc", "1");
%feature("docstring");

%newobject MEDCoupling::MEDCouplingFieldDiscretizationOnNodesFE::getCooInRefElement;
%newobject MEDCoupling::MEDCouplingFieldDiscretizationOnNodesFE::getClosestCooInRefElement;
%newobject MEDCoupling::MEDCouplingField::buildMeasureField;
%newobject MEDCoupling::MEDCouplingField::getLocalizationOfDiscr;
%newobject MEDCoupling::MEDCouplingField::computeTupleIdsToSelectFromCellIds;
%newobject MEDCoupling::MEDCouplingFieldDouble::New;
%newobject MEDCoupling::MEDCouplingFieldDouble::getArray;
%newobject MEDCoupling::MEDCouplingFieldDouble::getEndArray;
%newobject MEDCoupling::MEDCouplingFieldDouble::MergeFields;
%newobject MEDCoupling::MEDCouplingFieldDouble::MeldFields;
%newobject MEDCoupling::MEDCouplingFieldDouble::convertToIntField;
%newobject MEDCoupling::MEDCouplingFieldDouble::convertToFloatField;
%newobject MEDCoupling::MEDCouplingFieldDouble::doublyContractedProduct;
%newobject MEDCoupling::MEDCouplingFieldDouble::determinant;
%newobject MEDCoupling::MEDCouplingFieldDouble::eigenValues;
%newobject MEDCoupling::MEDCouplingFieldDouble::eigenVectors;
%newobject MEDCoupling::MEDCouplingFieldDouble::inverse;
%newobject MEDCoupling::MEDCouplingFieldDouble::trace;
%newobject MEDCoupling::MEDCouplingFieldDouble::deviator;
%newobject MEDCoupling::MEDCouplingFieldDouble::magnitude;
%newobject MEDCoupling::MEDCouplingFieldDouble::maxPerTuple;
%newobject MEDCoupling::MEDCouplingFieldDouble::keepSelectedComponents;
%newobject MEDCoupling::MEDCouplingFieldDouble::extractSlice3D;
%newobject MEDCoupling::MEDCouplingFieldDouble::DotFields;
%newobject MEDCoupling::MEDCouplingFieldDouble::dot;
%newobject MEDCoupling::MEDCouplingFieldDouble::CrossProductFields;
%newobject MEDCoupling::MEDCouplingFieldDouble::crossProduct;
%newobject MEDCoupling::MEDCouplingFieldDouble::MaxFields;
%newobject MEDCoupling::MEDCouplingFieldDouble::max;
%newobject MEDCoupling::MEDCouplingFieldDouble::MinFields;
%newobject MEDCoupling::MEDCouplingFieldDouble::AddFields;
%newobject MEDCoupling::MEDCouplingFieldDouble::SubstractFields;
%newobject MEDCoupling::MEDCouplingFieldDouble::MultiplyFields;
%newobject MEDCoupling::MEDCouplingFieldDouble::DivideFields;
%newobject MEDCoupling::MEDCouplingFieldDouble::min;
%newobject MEDCoupling::MEDCouplingFieldDouble::negate;
%newobject MEDCoupling::MEDCouplingFieldDouble::findIdsInRange;
%newobject MEDCoupling::MEDCouplingFieldDouble::buildSubPart;
%newobject MEDCoupling::MEDCouplingFieldDouble::buildSubPartRange;
%newobject MEDCoupling::MEDCouplingFieldDouble::voronoize;
%newobject MEDCoupling::MEDCouplingFieldDouble::convertQuadraticCellsToLinear;
%newobject MEDCoupling::MEDCouplingFieldDouble::__getitem__;
%newobject MEDCoupling::MEDCouplingFieldDouble::__neg__;
%newobject MEDCoupling::MEDCouplingFieldDouble::__add__;
%newobject MEDCoupling::MEDCouplingFieldDouble::__sub__;
%newobject MEDCoupling::MEDCouplingFieldDouble::__mul__;
%newobject MEDCoupling::MEDCouplingFieldDouble::__div__;
%newobject MEDCoupling::MEDCouplingFieldDouble::__pow__;
%newobject MEDCoupling::MEDCouplingFieldDouble::__radd__;
%newobject MEDCoupling::MEDCouplingFieldDouble::__rsub__;
%newobject MEDCoupling::MEDCouplingFieldDouble::__rmul__;
%newobject MEDCoupling::MEDCouplingFieldDouble::__rdiv__;
%newobject MEDCoupling::MEDCouplingFieldDouble::clone;
%newobject MEDCoupling::MEDCouplingFieldDouble::cloneWithMesh;
%newobject MEDCoupling::MEDCouplingFieldDouble::deepCopy;
%newobject MEDCoupling::MEDCouplingFieldDouble::buildNewTimeReprFromThis;
%newobject MEDCoupling::MEDCouplingFieldDouble::nodeToCellDiscretization;
%newobject MEDCoupling::MEDCouplingFieldDouble::cellToNodeDiscretization;
%newobject MEDCoupling::MEDCouplingFieldDouble::getValueOnMulti;
%newobject MEDCoupling::MEDCouplingFieldDouble::computeVectorFieldCyl;
%newobject MEDCoupling::MEDCouplingFieldInt32::New;
%newobject MEDCoupling::MEDCouplingFieldInt32::convertToDblField;
%newobject MEDCoupling::MEDCouplingFieldInt32::getArray;
%newobject MEDCoupling::MEDCouplingFieldInt32::deepCopy;
%newobject MEDCoupling::MEDCouplingFieldInt32::clone;
%newobject MEDCoupling::MEDCouplingFieldInt32::cloneWithMesh;
%newobject MEDCoupling::MEDCouplingFieldInt32::buildSubPart;
%newobject MEDCoupling::MEDCouplingFieldInt32::buildSubPartRange;
%newobject MEDCoupling::MEDCouplingFieldInt32::__getitem__;
%newobject MEDCoupling::MEDCouplingFieldInt64::New;
%newobject MEDCoupling::MEDCouplingFieldInt64::convertToDblField;
%newobject MEDCoupling::MEDCouplingFieldInt64::getArray;
%newobject MEDCoupling::MEDCouplingFieldInt64::deepCopy;
%newobject MEDCoupling::MEDCouplingFieldInt64::clone;
%newobject MEDCoupling::MEDCouplingFieldInt64::cloneWithMesh;
%newobject MEDCoupling::MEDCouplingFieldInt64::buildSubPart;
%newobject MEDCoupling::MEDCouplingFieldInt64::buildSubPartRange;
%newobject MEDCoupling::MEDCouplingFieldInt64::__getitem__;
%newobject MEDCoupling::MEDCouplingFieldFloat::New;
%newobject MEDCoupling::MEDCouplingFieldFloat::convertToDblField;
%newobject MEDCoupling::MEDCouplingFieldFloat::getArray;
%newobject MEDCoupling::MEDCouplingFieldFloat::deepCopy;
%newobject MEDCoupling::MEDCouplingFieldFloat::clone;
%newobject MEDCoupling::MEDCouplingFieldFloat::cloneWithMesh;
%newobject MEDCoupling::MEDCouplingFieldFloat::buildSubPart;
%newobject MEDCoupling::MEDCouplingFieldFloat::buildSubPartRange;
%newobject MEDCoupling::MEDCouplingFieldFloat::__getitem__;
%newobject MEDCoupling::MEDCouplingFieldTemplate::New;
%newobject MEDCoupling::MEDCouplingMesh::deepCopy;
%newobject MEDCoupling::MEDCouplingMesh::clone;
%newobject MEDCoupling::MEDCouplingMesh::checkDeepEquivalOnSameNodesWith;
%newobject MEDCoupling::MEDCouplingMesh::checkTypeConsistencyAndContig;
%newobject MEDCoupling::MEDCouplingMesh::computeNbOfNodesPerCell;
%newobject MEDCoupling::MEDCouplingMesh::computeNbOfFacesPerCell;
%newobject MEDCoupling::MEDCouplingMesh::computeEffectiveNbOfNodesPerCell;
%newobject MEDCoupling::MEDCouplingMesh::buildPartRange;
%newobject MEDCoupling::MEDCouplingMesh::giveCellsWithType;
%newobject MEDCoupling::MEDCouplingMesh::getCoordinatesAndOwner;
%newobject MEDCoupling::MEDCouplingMesh::computeMeshCenterOfMass;
%newobject MEDCoupling::MEDCouplingMesh::computeCellCenterOfMass;
%newobject MEDCoupling::MEDCouplingMesh::computeIsoBarycenterOfNodesPerCell;
%newobject MEDCoupling::MEDCouplingMesh::buildOrthogonalField;
%newobject MEDCoupling::MEDCouplingMesh::getCellIdsFullyIncludedInNodeIds;
%newobject MEDCoupling::MEDCouplingMesh::mergeMyselfWith;
%newobject MEDCoupling::MEDCouplingMesh::fillFromAnalytic;
%newobject MEDCoupling::MEDCouplingMesh::fillFromAnalyticCompo;
%newobject MEDCoupling::MEDCouplingMesh::fillFromAnalyticNamedCompo;
%newobject MEDCoupling::MEDCouplingMesh::getMeasureField;
%newobject MEDCoupling::MEDCouplingMesh::simplexize;
%newobject MEDCoupling::MEDCouplingMesh::buildUnstructured;
%newobject MEDCoupling::MEDCouplingMesh::MergeMeshes;
%newobject MEDCoupling::MEDCouplingMesh::getDirectAccessOfCoordsArrIfInStructure;
%newobject MEDCoupling::MEDCouplingPointSet::zipCoordsTraducer;
%newobject MEDCoupling::MEDCouplingPointSet::getCellsInBoundingBox;
%newobject MEDCoupling::MEDCouplingPointSet::findBoundaryNodes;
%newobject MEDCoupling::MEDCouplingPointSet::buildBoundaryMesh;
%newobject MEDCoupling::MEDCouplingPointSet::MergeNodesArray;
%newobject MEDCoupling::MEDCouplingPointSet::buildPartOfMySelfSlice;
%newobject MEDCoupling::MEDCouplingPointSet::BuildInstanceFromMeshType;
%newobject MEDCoupling::MEDCouplingPointSet::zipConnectivityTraducer;
%newobject MEDCoupling::MEDCouplingPointSet::mergeMyselfWithOnSameCoords;
%newobject MEDCoupling::MEDCouplingPointSet::fillCellIdsToKeepFromNodeIds;
%newobject MEDCoupling::MEDCouplingPointSet::getCellIdsLyingOnNodes;
%newobject MEDCoupling::MEDCouplingPointSet::deepCopyConnectivityOnly;
%newobject MEDCoupling::MEDCouplingPointSet::getBoundingBoxForBBTree;
%newobject MEDCoupling::MEDCouplingPointSet::computeFetchedNodeIds;
%newobject MEDCoupling::MEDCouplingPointSet::ComputeNbOfInteractionsWithSrcCells;
%newobject MEDCoupling::MEDCouplingPointSet::computeDiameterField;
%newobject MEDCoupling::MEDCouplingPointSet::__getitem__;
%newobject MEDCoupling::MEDCouplingUMesh::New;
%newobject MEDCoupling::MEDCouplingUMesh::getNodalConnectivity;
%newobject MEDCoupling::MEDCouplingUMesh::getNodalConnectivityIndex;
%newobject MEDCoupling::MEDCouplingUMesh::__iter__;
%newobject MEDCoupling::MEDCouplingUMesh::cellsByType;
%newobject MEDCoupling::MEDCouplingUMesh::buildDescendingConnectivity;
%newobject MEDCoupling::MEDCouplingUMesh::buildDescendingConnectivity2;
%newobject MEDCoupling::MEDCouplingUMesh::explode3DMeshTo1D;
%newobject MEDCoupling::MEDCouplingUMesh::explodeMeshIntoMicroEdges;
%newobject MEDCoupling::MEDCouplingUMesh::buildExtrudedMesh;
%newobject MEDCoupling::MEDCouplingUMesh::buildSpreadZonesWithPoly;
%newobject MEDCoupling::MEDCouplingUMesh::MergeUMeshes;
%newobject MEDCoupling::MEDCouplingUMesh::MergeUMeshesOnSameCoords;
%newobject MEDCoupling::MEDCouplingUMesh::ComputeSpreadZoneGradually;
%newobject MEDCoupling::MEDCouplingUMesh::ComputeSpreadZoneGraduallyFromSeed;
%newobject MEDCoupling::MEDCouplingUMesh::buildNewNumberingFromCommNodesFrmt;
%newobject MEDCoupling::MEDCouplingUMesh::conformize2D;
%newobject MEDCoupling::MEDCouplingUMesh::conformize3D;
%newobject MEDCoupling::MEDCouplingUMesh::colinearize2D;
%newobject MEDCoupling::MEDCouplingUMesh::colinearizeKeepingConform2D;
%newobject MEDCoupling::MEDCouplingUMesh::rearrange2ConsecutiveCellTypes;
%newobject MEDCoupling::MEDCouplingUMesh::sortCellsInMEDFileFrmt;
%newobject MEDCoupling::MEDCouplingUMesh::getRenumArrForMEDFileFrmt;
%newobject MEDCoupling::MEDCouplingUMesh::convertCellArrayPerGeoType;
%newobject MEDCoupling::MEDCouplingUMesh::getRenumArrForConsecutiveCellTypesSpec;
%newobject MEDCoupling::MEDCouplingUMesh::findNodesToDuplicate;
%newobject MEDCoupling::MEDCouplingUMesh::buildDirectionVectorField;
%newobject MEDCoupling::MEDCouplingUMesh::convertLinearCellsToQuadratic;
%newobject MEDCoupling::MEDCouplingUMesh::convertToQuadraticBasedOnSeg3;
%newobject MEDCoupling::MEDCouplingUMesh::extrudeConnectivity;
%newobject MEDCoupling::MEDCouplingUMesh::getEdgeRatioField;
%newobject MEDCoupling::MEDCouplingUMesh::getAspectRatioField;
%newobject MEDCoupling::MEDCouplingUMesh::getWarpField;
%newobject MEDCoupling::MEDCouplingUMesh::getSkewField;
%newobject MEDCoupling::MEDCouplingUMesh::getPartBarycenterAndOwner;
%newobject MEDCoupling::MEDCouplingUMesh::computePlaneEquationOf3DFaces;
%newobject MEDCoupling::MEDCouplingUMesh::getPartMeasureField;
%newobject MEDCoupling::MEDCouplingUMesh::buildPartOrthogonalField;
%newobject MEDCoupling::MEDCouplingUMesh::keepCellIdsByType;
%newobject MEDCoupling::MEDCouplingUMesh::Build0DMeshFromCoords;
%newobject MEDCoupling::MEDCouplingUMesh::Build1DMeshFromCoords;
%newobject MEDCoupling::MEDCouplingUMesh::findAndCorrectBadOriented3DExtrudedCells;
%newobject MEDCoupling::MEDCouplingUMesh::findAndCorrectBadOriented3DCells;
%newobject MEDCoupling::MEDCouplingUMesh::convertIntoSingleGeoTypeMesh;
%newobject MEDCoupling::MEDCouplingUMesh::convertNodalConnectivityToStaticGeoTypeMesh;
%newobject MEDCoupling::MEDCouplingUMesh::findCellIdsOnBoundary;
%newobject MEDCoupling::MEDCouplingUMesh::computeSkin;
%newobject MEDCoupling::MEDCouplingUMesh::buildSetInstanceFromThis;
%newobject MEDCoupling::MEDCouplingUMesh::getCellIdsCrossingPlane;
%newobject MEDCoupling::MEDCouplingUMesh::convexEnvelop2D;
%newobject MEDCoupling::MEDCouplingUMesh::ComputeRangesFromTypeDistribution;
%newobject MEDCoupling::MEDCouplingUMesh::buildUnionOf2DMesh;
%newobject MEDCoupling::MEDCouplingUMesh::buildUnionOf3DMesh;
%newobject MEDCoupling::MEDCouplingUMesh::generateGraph;
%newobject MEDCoupling::MEDCouplingUMesh::orderConsecutiveCells1D;
%newobject MEDCoupling::MEDCouplingUMesh::clipSingle3DCellByPlane;
%newobject MEDCoupling::MEDCouplingUMesh::getBoundingBoxForBBTreeFast;
%newobject MEDCoupling::MEDCouplingUMesh::getBoundingBoxForBBTree2DQuadratic;
%newobject MEDCoupling::MEDCouplingUMesh::getBoundingBoxForBBTree1DQuadratic;
%newobject MEDCoupling::MEDCouplingUMesh::convertDegeneratedCellsAndRemoveFlatOnes;
%newobject MEDCoupling::MEDCouplingUMeshCellByTypeEntry::__iter__;
%newobject MEDCoupling::MEDCouplingUMeshCellEntry::__iter__;
%newobject MEDCoupling::MEDCoupling1GTUMesh::New;
%newobject MEDCoupling::MEDCoupling1GTUMesh::getNodalConnectivity;
%newobject MEDCoupling::MEDCoupling1GTUMesh::AggregateOnSameCoordsToUMesh;
%newobject MEDCoupling::MEDCoupling1SGTUMesh::New;
%newobject MEDCoupling::MEDCoupling1SGTUMesh::buildSetInstanceFromThis;
%newobject MEDCoupling::MEDCoupling1SGTUMesh::computeDualMesh;
%newobject MEDCoupling::MEDCoupling1SGTUMesh::explodeEachHexa8To6Quad4;
%newobject MEDCoupling::MEDCoupling1SGTUMesh::computeTriangleHeight;
%newobject MEDCoupling::MEDCoupling1SGTUMesh::sortHexa8EachOther;
%newobject MEDCoupling::MEDCoupling1SGTUMesh::Merge1SGTUMeshes;
%newobject MEDCoupling::MEDCoupling1SGTUMesh::Merge1SGTUMeshesOnSameCoords;
%newobject MEDCoupling::MEDCoupling1DGTUMesh::New;
%newobject MEDCoupling::MEDCoupling1DGTUMesh::getNodalConnectivityIndex;
%newobject MEDCoupling::MEDCoupling1DGTUMesh::buildSetInstanceFromThis;
%newobject MEDCoupling::MEDCoupling1DGTUMesh::Merge1DGTUMeshes;
%newobject MEDCoupling::MEDCoupling1DGTUMesh::Merge1DGTUMeshesOnSameCoords;
%newobject MEDCoupling::MEDCouplingMappedExtrudedMesh::New;
%newobject MEDCoupling::MEDCouplingMappedExtrudedMesh::build3DUnstructuredMesh;
%newobject MEDCoupling::MEDCouplingStructuredMesh::buildStructuredSubPart;
%newobject MEDCoupling::MEDCouplingStructuredMesh::build1SGTUnstructured;
%newobject MEDCoupling::MEDCouplingStructuredMesh::build1SGTSubLevelMesh;
%newobject MEDCoupling::MEDCouplingStructuredMesh::BuildExplicitIdsFrom;
%newobject MEDCoupling::MEDCouplingStructuredMesh::ExtractFieldOfDoubleFrom;
%newobject MEDCoupling::MEDCouplingStructuredMesh::Build1GTNodalConnectivity;
%newobject MEDCoupling::MEDCouplingStructuredMesh::Build1GTNodalConnectivityOfSubLevelMesh;
%newobject MEDCoupling::MEDCouplingStructuredMesh::ComputeCornersGhost;
%newobject MEDCoupling::MEDCouplingCMesh::New;
%newobject MEDCoupling::MEDCouplingCMesh::getCoordsAt;
%newobject MEDCoupling::MEDCouplingCMesh::buildCurveLinear;
%newobject MEDCoupling::MEDCouplingIMesh::New;
%newobject MEDCoupling::MEDCouplingIMesh::asSingleCell;
%newobject MEDCoupling::MEDCouplingIMesh::buildWithGhost;
%newobject MEDCoupling::MEDCouplingIMesh::convertToCartesian;
%newobject MEDCoupling::MEDCouplingCurveLinearMesh::New;
%newobject MEDCoupling::MEDCouplingCurveLinearMesh::getCoords;
%newobject MEDCoupling::MEDCouplingMultiFields::New;
%newobject MEDCoupling::MEDCouplingMultiFields::deepCopy;
%newobject MEDCoupling::MEDCouplingFieldOverTime::New;
%newobject MEDCoupling::MEDCouplingCartesianAMRPatchGen::getMesh;
%newobject MEDCoupling::MEDCouplingCartesianAMRPatchGen::__getitem__;
%newobject MEDCoupling::MEDCouplingCartesianAMRMeshGen::deepCopy;
%newobject MEDCoupling::MEDCouplingCartesianAMRMeshGen::buildUnstructured;
%newobject MEDCoupling::MEDCouplingCartesianAMRMeshGen::extractGhostFrom;
%newobject MEDCoupling::MEDCouplingCartesianAMRMeshGen::buildMeshFromPatchEnvelop;
%newobject MEDCoupling::MEDCouplingCartesianAMRMeshGen::buildMeshOfDirectChildrenOnly;
%newobject MEDCoupling::MEDCouplingCartesianAMRMeshGen::getImageMesh;
%newobject MEDCoupling::MEDCouplingCartesianAMRMeshGen::getGodFather;
%newobject MEDCoupling::MEDCouplingCartesianAMRMeshGen::getFather;
%newobject MEDCoupling::MEDCouplingCartesianAMRMeshGen::getPatch;
%newobject MEDCoupling::MEDCouplingCartesianAMRMeshGen::createCellFieldOnPatch;
%newobject MEDCoupling::MEDCouplingCartesianAMRMeshGen::findPatchesInTheNeighborhoodOf;
%newobject MEDCoupling::MEDCouplingCartesianAMRMeshGen::getPatchAtPosition;
%newobject MEDCoupling::MEDCouplingCartesianAMRMeshGen::getMeshAtPosition;
%newobject MEDCoupling::MEDCouplingCartesianAMRMeshGen::__getitem__;
%newobject MEDCoupling::MEDCouplingCartesianAMRMesh::New;
%newobject MEDCoupling::MEDCouplingDataForGodFather::getMyGodFather;
%newobject MEDCoupling::MEDCouplingAMRAttribute::New;
%newobject MEDCoupling::MEDCouplingAMRAttribute::deepCopy;
%newobject MEDCoupling::MEDCouplingAMRAttribute::deepCpyWithoutGodFather;
%newobject MEDCoupling::MEDCouplingAMRAttribute::getFieldOn;
%newobject MEDCoupling::MEDCouplingAMRAttribute::projectTo;
%newobject MEDCoupling::MEDCouplingAMRAttribute::buildCellFieldOnRecurseWithoutOverlapWithoutGhost;
%newobject MEDCoupling::MEDCouplingAMRAttribute::buildCellFieldOnWithGhost;
%newobject MEDCoupling::MEDCouplingAMRAttribute::buildCellFieldOnWithoutGhost;
%newobject MEDCoupling::DenseMatrix::New;
%newobject MEDCoupling::DenseMatrix::deepCopy;
%newobject MEDCoupling::DenseMatrix::shallowCpy;
%newobject MEDCoupling::DenseMatrix::getData;
%newobject MEDCoupling::DenseMatrix::matVecMult;
%newobject MEDCoupling::DenseMatrix::MatVecMult;
%newobject MEDCoupling::DenseMatrix::__add__;
%newobject MEDCoupling::DenseMatrix::__sub__;
%newobject MEDCoupling::DenseMatrix::__mul__;
%newobject MEDCoupling::MEDCouplingGaussLocalization::localizePtsInRefCooForEachCell;
%newobject MEDCoupling::MEDCouplingGaussLocalization::buildRefCell;
%newobject MEDCoupling::MEDCouplingGaussLocalization::getShapeFunctionValues;
%newobject MEDCoupling::MEDCouplingGaussLocalization::getDerivativeOfShapeFunctionValues;
%newobject MEDCoupling::MEDCouplingGaussLocalization::GetDefaultReferenceCoordinatesOf;
%newobject MEDCoupling::MEDCouplingSkyLineArray::BuildFromPolyhedronConn;
%newobject MEDCoupling::MEDCouplingSkyLineArray::getSuperIndexArray;
%newobject MEDCoupling::MEDCouplingSkyLineArray::getIndexArray;
%newobject MEDCoupling::MEDCouplingSkyLineArray::getValuesArray;
%newobject MEDCoupling::MEDCouplingSkyLineArray::groupPacks;
%newobject MEDCoupling::MEDCouplingSkyLineArray::uniqueNotSortedByPack;
%newobject MEDCoupling::MEDCouplingSkyLineArray::AggregatePacks;
%newobject MEDCoupling::MEDCouplingSkyLineArray::deepCopy;

%feature("unref") MEDCouplingPointSet "$this->decrRef();"
%feature("unref") MEDCouplingMesh "$this->decrRef();"
%feature("unref") MEDCouplingUMesh "$this->decrRef();"
%feature("unref") MEDCoupling1GTUMesh "$this->decrRef();"
%feature("unref") MEDCoupling1SGTUMesh "$this->decrRef();"
%feature("unref") MEDCoupling1DGTUMesh "$this->decrRef();"
%feature("unref") MEDCouplingMappedExtrudedMesh "$this->decrRef();"
%feature("unref") MEDCouplingCMesh "$this->decrRef();"
%feature("unref") MEDCouplingIMesh "$this->decrRef();"
%feature("unref") MEDCouplingCurveLinearMesh "$this->decrRef();"
%feature("unref") MEDCouplingField "$this->decrRef();"
%feature("unref") MEDCouplingFieldDiscretizationP0 "$this->decrRef();"
%feature("unref") MEDCouplingFieldDiscretizationP1 "$this->decrRef();"
%feature("unref") MEDCouplingFieldDiscretizationGauss "$this->decrRef();"
%feature("unref") MEDCouplingFieldDiscretizationGaussNE "$this->decrRef();"
%feature("unref") MEDCouplingFieldDiscretizationKriging "$this->decrRef();"
%feature("unref") MEDCouplingFieldDiscretizationOnNodesFE "$this->decrRef();"
%feature("unref") MEDCouplingFieldDouble "$this->decrRef();"
%feature("unref") MEDCouplingFieldFloat "$this->decrRef();"
%feature("unref") MEDCouplingFieldInt32 "$this->decrRef();"
%feature("unref") MEDCouplingFieldInt64 "$this->decrRef();"
%feature("unref") MEDCouplingMultiFields "$this->decrRef();"
%feature("unref") MEDCouplingFieldTemplate "$this->decrRef();"
%feature("unref") MEDCouplingMultiFields "$this->decrRef();"
%feature("unref") MEDCouplingCartesianAMRMeshGen "$this->decrRef();"
%feature("unref") MEDCouplingCartesianAMRMesh "$this->decrRef();"
%feature("unref") MEDCouplingCartesianAMRMeshSub "$this->decrRef();"
%feature("unref") MEDCouplingCartesianAMRPatchGen "$this->decrRef();"
%feature("unref") MEDCouplingCartesianAMRPatchGF "$this->decrRef();"
%feature("unref") MEDCouplingCartesianAMRPatch "$this->decrRef();"
%feature("unref") MEDCouplingDataForGodFather "$this->decrRef();"
%feature("unref") MEDCouplingAMRAttribute "$this->decrRef();"
%feature("unref") DenseMatrix "$this->decrRef();"
%feature("unref") MEDCouplingSkyLineArray "$this->decrRef();"

%rename(assign) *::operator=;
%ignore MEDCoupling::MEDCouplingGaussLocalization::pushTinySerializationIntInfo;
%ignore MEDCoupling::MEDCouplingGaussLocalization::pushTinySerializationDblInfo;
%ignore MEDCoupling::MEDCouplingGaussLocalization::fillWithValues;
%ignore MEDCoupling::MEDCouplingGaussLocalization::buildNewInstanceFromTinyInfo;

%nodefaultctor;

// ABN: Instruct SWIG that INTERP_KERNEL::Exception is an exception class and that it should inherit Exception
// on the Python side. Must be put BEFORE the %rename clause:
%exceptionclass INTERP_KERNEL::Exception;
%rename (InterpKernelException) INTERP_KERNEL::Exception;

%include "MEDCouplingRefCountObject.i"
%include "MEDCouplingMemArray.i"

%{
  void initializeMe()
  {// AGY : here initialization of C++ traits in MEDCouplingDataArrayTypemaps.i for code factorization. Awful, I know, but no other solutions.
    SWIGTITraits<double>::TI=SWIGTYPE_p_MEDCoupling__DataArrayDouble;
    SWIGTITraits<float>::TI=SWIGTYPE_p_MEDCoupling__DataArrayFloat;
    SWIGTITraits<Int32>::TI=SWIGTYPE_p_MEDCoupling__DataArrayInt32;
    SWIGTITraits<Int64>::TI=SWIGTYPE_p_MEDCoupling__DataArrayInt64;
    SWIGTITraits<double>::TI_TUPLE=SWIGTYPE_p_MEDCoupling__DataArrayDoubleTuple;
    SWIGTITraits<float>::TI_TUPLE=SWIGTYPE_p_MEDCoupling__DataArrayFloatTuple;
    SWIGTITraits<Int32>::TI_TUPLE=SWIGTYPE_p_MEDCoupling__DataArrayInt32Tuple;
    SWIGTITraits<Int64>::TI_TUPLE=SWIGTYPE_p_MEDCoupling__DataArrayInt64Tuple;
  }
%}

%inline
{
  PyObject *med2vtk_cell_types()
  {
    Py_ssize_t sz(sizeof(MEDCOUPLING2VTKTYPETRADUCER)/sizeof(decltype(MEDCOUPLING2VTKTYPETRADUCER[0])));
    PyObject *ret(PyList_New(sz));
    for(Py_ssize_t i=0;i<sz;i++)
      {
        mcIdType elt = MEDCOUPLING2VTKTYPETRADUCER[i]!=MEDCOUPLING2VTKTYPETRADUCER_NONE ? MEDCOUPLING2VTKTYPETRADUCER[i] : -1;
        PyList_SetItem(ret,i,PyInt_FromLong(elt));
      }
    return ret;
  }

  PyObject *vtk2med_cell_types()
  {
    Py_ssize_t sz(sizeof(MEDCOUPLING2VTKTYPETRADUCER)/sizeof(decltype(MEDCOUPLING2VTKTYPETRADUCER[0])));
    auto maxElt(*std::max_element(MEDCOUPLING2VTKTYPETRADUCER,MEDCOUPLING2VTKTYPETRADUCER+sz,[](unsigned char a, unsigned char b) { if(b==MEDCOUPLING2VTKTYPETRADUCER_NONE) return false; else return a<b; } ));
    auto szOut(maxElt+1);
    std::vector< mcIdType > retCpp(szOut,-1);
    mcIdType id(0);
    for(const unsigned char *it=MEDCOUPLING2VTKTYPETRADUCER;it!=MEDCOUPLING2VTKTYPETRADUCER+sz;it++,id++)
      {
        if(*it!=MEDCOUPLING2VTKTYPETRADUCER_NONE)
          retCpp[*it]=id;
      }
    //
    PyObject *ret(PyList_New(szOut));
    id = 0;
    for(auto it=retCpp.begin();it!=retCpp.end();it++,id++)
      PyList_SetItem(ret,id,PyInt_FromLong(*it));
    return ret;
  }

  PyObject *AllGeometricTypes()
  {
    Py_ssize_t sz(MEDCouplingUMesh::N_MEDMEM_ORDER);
    PyObject *ret(PyList_New(sz));
    for(Py_ssize_t i=0;i<sz;i++)
      PyList_SetItem(ret,i,PyInt_FromLong(MEDCouplingUMesh::MEDMEM_ORDER[i]));
    return ret;
  }
}

namespace INTERP_KERNEL
{
  /*!
   * \class BoxSplittingOptions
   * Class defining the options for box splitting used for AMR algorithm like creation of patches following a criterion.
   */
  class BoxSplittingOptions
  {
  public:
    BoxSplittingOptions();
    void init();
    double getEfficiencyGoal() const;
    void setEfficiencyGoal(double efficiency);
    double getEfficiencyThreshold() const;
    void setEfficiencyThreshold(double efficiencyThreshold);
    int getMinimumPatchLength() const;
    void setMinimumPatchLength(int minPatchLength);
    int getMaximumPatchLength() const;
    void setMaximumPatchLength(int maxPatchLength);
    int getMaximumNbOfCellsInPatch() const;
    void setMaximumNbOfCellsInPatch(int maxNbCellsInPatch);
    void copyOptions(const BoxSplittingOptions & other);
    std::string printOptions() const;
    %extend
    {
      std::string __str__() const
      {
        return self->printOptions();
      }
    }
  };
}

namespace MEDCoupling
{
  typedef enum
    {
      ON_CELLS = 0,
      ON_NODES = 1,
      ON_GAUSS_PT = 2,
      ON_GAUSS_NE = 3,
      ON_NODES_KR = 4,
      ON_NODES_FE = 5
    } TypeOfField;

  typedef enum
    {
      NO_TIME = 4,
      ONE_TIME = 5,
      LINEAR_TIME = 6,
      CONST_ON_TIME_INTERVAL = 7
    } TypeOfTimeDiscretization;

  typedef enum
    {
      UNSTRUCTURED = 5,
      CARTESIAN = 7,
      EXTRUDED = 8,
      CURVE_LINEAR = 9,
      SINGLE_STATIC_GEO_TYPE_UNSTRUCTURED = 10,
      SINGLE_DYNAMIC_GEO_TYPE_UNSTRUCTURED = 11,
      IMAGE_GRID = 12
    } MEDCouplingMeshType;

  class DataArrayInt32;
  class DataArrayInt64;
  class DataArrayDouble;
  class MEDCouplingUMesh;
  class MEDCouplingCMesh;
  class MEDCouplingFieldDouble;

  %extend RefCountObject
  {
    std::string getHiddenCppPointer() const
    {
      std::ostringstream oss; oss << "C++ Pointer address is : " << self;
      return oss.str();
    }

    // Hack to allow retrieving of underlying C++ pointer whatever the situation
    // This allows for example to mix different types of Python binding (SWIG and PyBind for example)
    long long getHiddenCppPointerAsLongLong() const
    {
      return (long long) self;
    }

  }

  %extend MEDCouplingGaussLocalization
  {
    std::string __str__() const
    {
      return self->getStringRepr();
    }

    std::string __repr__() const
    {
      std::ostringstream oss; oss << "MEDCouplingGaussLocalization C++ instance at " << self << "." << std::endl;
      oss << self->getStringRepr();
      return oss.str();
    }
  }

  //== MEDCouplingMesh

  class MEDCouplingMesh : public RefCountObject, public TimeLabel
  {
  public:
    void setName(const std::string& name);
    std::string getName() const;
    void setDescription(const std::string& descr);
    std::string getDescription() const;
    void setTime(double val, int iteration, int order);
    void setTimeUnit(const std::string& unit);
    std::string getTimeUnit() const;
    virtual MEDCouplingMeshType getType() const;
    bool isStructured() const;
    virtual MEDCouplingMesh *deepCopy() const;
    virtual MEDCouplingMesh *clone(bool recDeepCpy) const;
    virtual bool isEqual(const MEDCouplingMesh *other, double prec) const;
    virtual bool isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const;
    virtual void checkFastEquivalWith(const MEDCouplingMesh *other, double prec) const;
    virtual void copyTinyStringsFrom(const MEDCouplingMesh *other);
    virtual void copyTinyInfoFrom(const MEDCouplingMesh *other);
    virtual void checkConsistencyLight() const;
    virtual void checkConsistency(double eps=1e-12) const;
    virtual int getNumberOfCells() const;
    virtual int getNumberOfNodes() const;
    virtual int getSpaceDimension() const;
    virtual int getMeshDimension() const;
    virtual DataArrayDouble *getCoordinatesAndOwner() const;
    virtual DataArrayDouble *computeCellCenterOfMass() const;
    virtual DataArrayDouble *computeIsoBarycenterOfNodesPerCell() const;
    virtual DataArrayIdType *giveCellsWithType(INTERP_KERNEL::NormalizedCellType type) const;
    virtual DataArrayIdType *computeNbOfNodesPerCell() const;
    virtual DataArrayIdType *computeNbOfFacesPerCell() const;
    virtual DataArrayIdType *computeEffectiveNbOfNodesPerCell() const;
    virtual MEDCouplingMesh *buildPartRange(int beginCellIds, int endCellIds, int stepCellIds) const;
    virtual int getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType type) const;
    virtual INTERP_KERNEL::NormalizedCellType getTypeOfCell(int cellId) const;
    virtual std::string simpleRepr() const;
    virtual std::string advancedRepr() const;
    std::string writeVTK(const std::string& fileName, bool isBinary=true) const;
    virtual std::string getVTKFileExtension() const;
    std::string getVTKFileNameOf(const std::string& fileName) const;
    // tools
    virtual MEDCouplingFieldDouble *getMeasureField(bool isAbs) const;
    virtual MEDCouplingFieldDouble *getMeasureFieldOnNode(bool isAbs) const;
    virtual MEDCouplingFieldDouble *fillFromAnalytic(TypeOfField t, int nbOfComp, const std::string& func) const;
    virtual MEDCouplingFieldDouble *fillFromAnalyticCompo(TypeOfField t, int nbOfComp, const std::string& func) const;
    virtual MEDCouplingFieldDouble *fillFromAnalyticNamedCompo(TypeOfField t, int nbOfComp, const std::vector<std::string>& varsOrder, const std::string& func) const;
    virtual MEDCouplingFieldDouble *buildOrthogonalField() const;
    virtual MEDCouplingUMesh *buildUnstructured() const;
    virtual MEDCouplingMesh *mergeMyselfWith(const MEDCouplingMesh *other) const;
    virtual bool areCompatibleForMerge(const MEDCouplingMesh *other) const;
    virtual DataArrayIdType *simplexize(int policy);
    virtual void unserialization(const std::vector<double>& tinyInfoD, const std::vector<mcIdType>& tinyInfo, const DataArrayIdType *a1, DataArrayDouble *a2, const std::vector<std::string>& littleStrings);
    static MEDCouplingMesh *MergeMeshes(const MEDCouplingMesh *mesh1, const MEDCouplingMesh *mesh2);
    static bool IsStaticGeometricType(INTERP_KERNEL::NormalizedCellType type);
    static bool IsLinearGeometricType(INTERP_KERNEL::NormalizedCellType type);
    static INTERP_KERNEL::NormalizedCellType GetCorrespondingPolyType(INTERP_KERNEL::NormalizedCellType type);
    static int GetNumberOfNodesOfGeometricType(INTERP_KERNEL::NormalizedCellType type);
    static int GetDimensionOfGeometricType(INTERP_KERNEL::NormalizedCellType type);
    static const char *GetReprOfGeometricType(INTERP_KERNEL::NormalizedCellType type);
    %extend
       {
         std::string __str__() const
         {
           return self->simpleRepr();
         }

          DataArrayDouble *computeMeshCenterOfMass() const
          {
            MCAuto<DataArrayDouble> ret(self->computeMeshCenterOfMass());
            return ret.retn();
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

         DataArrayDouble *getDirectAccessOfCoordsArrIfInStructure() const
         {
           const DataArrayDouble *ret(self->getDirectAccessOfCoordsArrIfInStructure());
           DataArrayDouble *ret2(const_cast<DataArrayDouble *>(ret));
           if(ret2)
             ret2->incrRef();
           return ret2;
         }

         mcIdType getCellContainingPoint(PyObject *p, double eps) const
         {
           double val;
           DataArrayDouble *a;
           DataArrayDoubleTuple *aa;
           std::vector<double> bb;
           mcIdType sw;
           int spaceDim=self->getSpaceDimension();
           const char msg[]="Python wrap of MEDCouplingMesh::getCellContainingPoint : ";
           const double *pos=convertObjToPossibleCpp5_Safe(p,sw,val,a,aa,bb,msg,1,spaceDim,true);
           return self->getCellContainingPoint(pos,eps);
         }

         PyObject *getCellsContainingPoints(PyObject *p, int nbOfPoints, double eps) const
         {
           double val;
           DataArrayDouble *a;
           DataArrayDoubleTuple *aa;
           std::vector<double> bb;
           mcIdType sw;
           int spaceDim=self->getSpaceDimension();
           const char msg[]="Python wrap of MEDCouplingMesh::getCellsContainingPoint : ";
           const double *pos=convertObjToPossibleCpp5_Safe(p,sw,val,a,aa,bb,msg,nbOfPoints,spaceDim,true);
           MCAuto<DataArrayIdType> elts,eltsIndex;
           self->getCellsContainingPoints(pos,nbOfPoints,eps,elts,eltsIndex);
           PyObject *ret=PyTuple_New(2);
           PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(elts.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
           PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(eltsIndex.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
           return ret;
         }

         PyObject *getCellsContainingPointsLinearPartOnlyOnNonDynType(PyObject *p, int nbOfPoints, double eps) const
         {
           double val;
           DataArrayDouble *a;
           DataArrayDoubleTuple *aa;
           std::vector<double> bb;
           mcIdType sw;
           int spaceDim=self->getSpaceDimension();
           const char msg[]="Python wrap of MEDCouplingMesh::getCellsContainingPointsLinearPartOnlyOnNonDynType : ";
           const double *pos=convertObjToPossibleCpp5_Safe(p,sw,val,a,aa,bb,msg,nbOfPoints,spaceDim,true);
           MCAuto<DataArrayIdType> elts,eltsIndex;
           self->getCellsContainingPointsLinearPartOnlyOnNonDynType(pos,nbOfPoints,eps,elts,eltsIndex);
           PyObject *ret=PyTuple_New(2);
           PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(elts.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
           PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(eltsIndex.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
           return ret;
         }

         PyObject *getCellsContainingPoints(PyObject *p, double eps) const
         {
           auto getCellsContainingPointsFunc=[self](const double *a, int b,double c, MCAuto<DataArrayIdType>& d, MCAuto<DataArrayIdType>& e) { self->getCellsContainingPoints(a,b,c,d,e); };
           return Mesh_getCellsContainingPointsLike(p,eps,self,getCellsContainingPointsFunc);
         }

         PyObject *getCellsContainingPointsLinearPartOnlyOnNonDynType(PyObject *p, double eps) const
         {
           auto getCellsContainingPointsFunc=[self](const double *a, int b,double c, MCAuto<DataArrayIdType>& d, MCAuto<DataArrayIdType>& e) { self->getCellsContainingPointsLinearPartOnlyOnNonDynType(a,b,c,d,e); };
           return Mesh_getCellsContainingPointsLike(p,eps,self,getCellsContainingPointsFunc);
         }

         PyObject *getCellsContainingPoint(PyObject *p, double eps) const
         {
           double val;
           DataArrayDouble *a;
           DataArrayDoubleTuple *aa;
           std::vector<double> bb;
           mcIdType sw;
           int spaceDim=self->getSpaceDimension();
           const char msg[]="Python wrap of MEDCouplingUMesh::getCellsContainingPoint : ";
           const double *pos=convertObjToPossibleCpp5_Safe(p,sw,val,a,aa,bb,msg,1,spaceDim,true);
           std::vector<mcIdType> elts;
           self->getCellsContainingPoint(pos,eps,elts);
           DataArrayIdType *ret=DataArrayIdType::New();
           ret->alloc((int)elts.size(),1);
           std::copy(elts.begin(),elts.end(),ret->getPointer());
           return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 );
         }

         virtual PyObject *getReverseNodalConnectivity() const
         {
           MCAuto<DataArrayIdType> d0=DataArrayIdType::New();
           MCAuto<DataArrayIdType> d1=DataArrayIdType::New();
           self->getReverseNodalConnectivity(d0,d1);
           PyObject *ret=PyTuple_New(2);
           PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(d0.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
           PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(d1.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
           return ret;
         }

         void renumberCells(PyObject *li, bool check=true)
         {
           mcIdType sw,sz(-1);
           mcIdType v0; std::vector<mcIdType> v1;
           const mcIdType *ids(convertIntStarLikePyObjToCppIntStar(li,sw,sz,v0,v1));
           self->renumberCells(ids,check);
         }

         PyObject *checkGeoEquivalWith(const MEDCouplingMesh *other, int levOfCheck, double prec) const
         {
           DataArrayIdType *cellCor, *nodeCor;
           self->checkGeoEquivalWith(other,levOfCheck,prec,cellCor,nodeCor);
           PyObject *res = PyList_New(2);
           PyList_SetItem(res,0,SWIG_NewPointerObj(SWIG_as_voidptr(cellCor),SWIGTITraits<mcIdType>::TI, cellCor?SWIG_POINTER_OWN | 0:0 ));
           PyList_SetItem(res,1,SWIG_NewPointerObj(SWIG_as_voidptr(nodeCor),SWIGTITraits<mcIdType>::TI, nodeCor?SWIG_POINTER_OWN | 0:0 ));
           return res;
         }

         PyObject *checkDeepEquivalWith(const MEDCouplingMesh *other, int cellCompPol, double prec) const
         {
           DataArrayIdType *cellCor=0,*nodeCor=0;
           self->checkDeepEquivalWith(other,cellCompPol,prec,cellCor,nodeCor);
           PyObject *res = PyList_New(2);
           PyList_SetItem(res,0,SWIG_NewPointerObj(SWIG_as_voidptr(cellCor),SWIGTITraits<mcIdType>::TI, cellCor?SWIG_POINTER_OWN | 0:0 ));
           PyList_SetItem(res,1,SWIG_NewPointerObj(SWIG_as_voidptr(nodeCor),SWIGTITraits<mcIdType>::TI, nodeCor?SWIG_POINTER_OWN | 0:0 ));
           return res;
         }

         DataArrayIdType *checkDeepEquivalOnSameNodesWith(const MEDCouplingMesh *other, int cellCompPol, double prec) const
         {
           DataArrayIdType *cellCor=0;
           self->checkDeepEquivalOnSameNodesWith(other,cellCompPol,prec,cellCor);
           return cellCor;
         }

         DataArrayIdType *getCellIdsFullyIncludedInNodeIds(PyObject *li) const
         {
           void *da=0;
           int res1=SWIG_ConvertPtr(li,&da,SWIGTITraits<mcIdType>::TI, 0 |  0 );
           if (!SWIG_IsOK(res1))
             {
               mcIdType size;
               INTERP_KERNEL::AutoPtr<mcIdType> tmp=convertPyToNewIntArr2(li,&size);
               return self->getCellIdsFullyIncludedInNodeIds(tmp,((const mcIdType *)tmp)+size);
             }
           else
             {
               DataArrayIdType *da2=reinterpret_cast< DataArrayIdType * >(da);
               if(!da2)
                 throw INTERP_KERNEL::Exception("Not null DataArrayIdType instance expected !");
               da2->checkAllocated();
               return self->getCellIdsFullyIncludedInNodeIds(da2->getConstPointer(),da2->getConstPointer()+da2->getNbOfElems());
             }
         }
         PyObject *getNodeIdsOfCell(int cellId) const
         {
           std::vector<mcIdType> conn;
           self->getNodeIdsOfCell(cellId,conn);
           return convertIntArrToPyList2(conn);
         }

         PyObject *getCoordinatesOfNode(mcIdType nodeId) const
         {
           std::vector<double> coo;
           self->getCoordinatesOfNode(nodeId,coo);
           return convertDblArrToPyList2(coo);
         }

         void scale(PyObject *point, double factor)
         {
           double val;
           DataArrayDouble *a;
           DataArrayDoubleTuple *aa;
           std::vector<double> bb;
           mcIdType sw;
           int spaceDim=self->getSpaceDimension();
           const char msg[]="Python wrap of MEDCouplingPointSet::scale : ";
           const double *pointPtr=convertObjToPossibleCpp5_Safe(point,sw,val,a,aa,bb,msg,1,spaceDim,true);
           self->scale(pointPtr,factor);
         }

         PyObject *getBoundingBox() const
         {
           int spaceDim=self->getSpaceDimension();
           INTERP_KERNEL::AutoPtr<double> tmp=new double[2*spaceDim];
           self->getBoundingBox(tmp);
           PyObject *ret=convertDblArrToPyListOfTuple<double>(tmp,2,spaceDim);
           return ret;
         }

         PyObject *isEqualIfNotWhy(const MEDCouplingMesh *other, double prec) const
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

         PyObject *buildPart(PyObject *li) const
         {
           mcIdType szArr,sw,iTypppArr;
           std::vector<mcIdType> stdvecTyyppArr;
           const mcIdType *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
           MEDCouplingMesh *ret=self->buildPart(tmp,tmp+szArr);
           if(sw==3)//DataArrayIdType
             {
               void *argp; SWIG_ConvertPtr(li,&argp,SWIGTITraits<mcIdType>::TI,0|0);
               DataArrayIdType *argpt=reinterpret_cast< MEDCoupling::DataArrayIdType * >(argp);
               std::string name=argpt->getName();
               if(!name.empty())
                 ret->setName(name.c_str());
             }
           return convertMesh(ret, SWIG_POINTER_OWN | 0 );
         }

         PyObject *buildPartAndReduceNodes(PyObject *li) const
         {
           mcIdType szArr,sw,iTypppArr;
           std::vector<mcIdType> stdvecTyyppArr;
           DataArrayIdType *arr=0;
           const mcIdType *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
           MEDCouplingMesh *ret=self->buildPartAndReduceNodes(tmp,tmp+szArr,arr);
           if(sw==3)//DataArrayIdType
             {
               void *argp; SWIG_ConvertPtr(li,&argp,SWIGTITraits<mcIdType>::TI,0|0);
               DataArrayIdType *argpt=reinterpret_cast< MEDCoupling::DataArrayIdType * >(argp);
               std::string name=argpt->getName();
               if(!name.empty())
                 ret->setName(name.c_str());
             }
           //
           PyObject *res = PyList_New(2);
           PyObject *obj0=convertMesh(ret, SWIG_POINTER_OWN | 0 );
           PyObject *obj1=SWIG_NewPointerObj(SWIG_as_voidptr(arr),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 );
           PyList_SetItem(res,0,obj0);
           PyList_SetItem(res,1,obj1);
           return res;
         }

         PyObject *buildPartRangeAndReduceNodes(mcIdType beginCellIds, mcIdType endCellIds, mcIdType stepCellIds) const
         {
           mcIdType a,b,c;
           DataArrayIdType *arr=0;
           MEDCouplingMesh *ret=self->buildPartRangeAndReduceNodes(beginCellIds,endCellIds,stepCellIds,a,b,c,arr);
           PyObject *res = PyTuple_New(2);
           PyObject *obj0=convertMesh(ret, SWIG_POINTER_OWN | 0 );
           PyObject *obj1=0;
           if(arr)
             obj1=SWIG_NewPointerObj(SWIG_as_voidptr(arr),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 );
           else
             obj1=PySlice_New(PyInt_FromLong(a),PyInt_FromLong(b),PyInt_FromLong(b));
           PyTuple_SetItem(res,0,obj0);
           PyTuple_SetItem(res,1,obj1);
           return res;
         }

        PyObject *getDistributionOfTypes() const
        {
          std::vector<mcIdType> vals=self->getDistributionOfTypes();
          if(vals.size()%3!=0)
            throw INTERP_KERNEL::Exception("Internal Error detected in wrap python ! code returned by MEDCouplingMesh::getDistributionOfTypes is not so that %3==0 !");
          PyObject *ret=PyList_New((mcIdType)vals.size()/3);
          for(std::size_t j=0;j<vals.size()/3;j++)
             {
               PyObject *ret1=PyList_New(3);
               PyList_SetItem(ret1,0,PyInt_FromLong(vals[3*j]));
               PyList_SetItem(ret1,1,PyInt_FromLong(vals[3*j+1]));
               PyList_SetItem(ret1,2,PyInt_FromLong(vals[3*j+2]));
               PyList_SetItem(ret,j,ret1);
             }
          return ret;
        }

        DataArrayIdType *checkTypeConsistencyAndContig(PyObject *li, PyObject *li2) const
        {
          std::vector<mcIdType> code;
          std::vector<const DataArrayIdType *> idsPerType;
          convertFromPyObjVectorOfObj<const MEDCoupling::DataArrayIdType *>(li2,SWIGTITraits<mcIdType>::TI,"DataArrayIdType",idsPerType);
          convertPyToNewIntArr4(li,1,3,code);
          return self->checkTypeConsistencyAndContig(code,idsPerType);
        }

        PyObject *splitProfilePerType(const DataArrayIdType *profile, bool smartPflKiller=true) const
        {
          std::vector<mcIdType> code;
          std::vector<DataArrayIdType *> idsInPflPerType;
          std::vector<DataArrayIdType *> idsPerType;
          self->splitProfilePerType(profile,code,idsInPflPerType,idsPerType,smartPflKiller);
          PyObject *ret=PyTuple_New(3);
          //
          if(code.size()%3!=0)
            throw INTERP_KERNEL::Exception("Internal Error detected in wrap python ! code returned by MEDCouplingMesh::splitProfilePerType is not so that %3==0 !");
          PyObject *ret0=PyList_New((mcIdType)code.size()/3);
          for(std::size_t j=0;j<code.size()/3;j++)
             {
               PyObject *ret00=PyList_New(3);
               PyList_SetItem(ret00,0,PyInt_FromLong(code[3*j]));
               PyList_SetItem(ret00,1,PyInt_FromLong(code[3*j+1]));
               PyList_SetItem(ret00,2,PyInt_FromLong(code[3*j+2]));
               PyList_SetItem(ret0,j,ret00);
             }
          PyTuple_SetItem(ret,0,ret0);
          //
          PyObject *ret1=PyList_New(idsInPflPerType.size());
          for(std::size_t j=0;j<idsInPflPerType.size();j++)
            PyList_SetItem(ret1,j,SWIG_NewPointerObj(SWIG_as_voidptr(idsInPflPerType[j]),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
          PyTuple_SetItem(ret,1,ret1);
          std::size_t n=idsPerType.size();
          PyObject *ret2=PyList_New(n);
          for(std::size_t i=0;i<n;i++)
            PyList_SetItem(ret2,i,SWIG_NewPointerObj(SWIG_as_voidptr(idsPerType[i]),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
          PyTuple_SetItem(ret,2,ret2);
          return ret;
        }

        void translate(PyObject *vector)
        {
          double val;
          DataArrayDouble *a;
          DataArrayDoubleTuple *aa;
          std::vector<double> bb;
          mcIdType sw;
          int spaceDim=self->getSpaceDimension();
          const char msg[]="Python wrap of MEDCouplingPointSet::translate : ";
          const double *vectorPtr=convertObjToPossibleCpp5_Safe(vector,sw,val,a,aa,bb,msg,1,spaceDim,true);
          self->translate(vectorPtr);
        }

         void rotate(PyObject *center, double alpha)
         {
           const char msg[]="Python wrap of MEDCouplingPointSet::rotate : ";
           double val;
           DataArrayDouble *a;
           DataArrayDoubleTuple *aa;
           std::vector<double> bb;
           mcIdType sw;
           int spaceDim=self->getSpaceDimension();
           const double *centerPtr=convertObjToPossibleCpp5_Safe(center,sw,val,a,aa,bb,msg,1,spaceDim,true);
           self->rotate(centerPtr,0,alpha);
         }

         void rotate(PyObject *center, PyObject *vector, double alpha)
         {
           const char msg[]="Python wrap of MEDCouplingPointSet::rotate : ";
           double val,val2;
           DataArrayDouble *a,*a2;
           DataArrayDoubleTuple *aa,*aa2;
           std::vector<double> bb,bb2;
           mcIdType sw;
           int spaceDim=self->getSpaceDimension();
           const double *centerPtr=convertObjToPossibleCpp5_Safe(center,sw,val,a,aa,bb,msg,1,spaceDim,true);
           const double *vectorPtr=convertObjToPossibleCpp5_Safe(vector,sw,val2,a2,aa2,bb2,msg,1,spaceDim,false);//vectorPtr can be null in case of space dim 2
           self->rotate(centerPtr,vectorPtr,alpha);
         }

         PyObject *getAllGeoTypes() const
         {
           std::set<INTERP_KERNEL::NormalizedCellType> result=self->getAllGeoTypes();
           std::set<INTERP_KERNEL::NormalizedCellType>::const_iterator iL=result.begin();
           PyObject *res=PyList_New(result.size());
           for(int i=0;iL!=result.end(); i++, iL++)
             PyList_SetItem(res,i,PyInt_FromLong(*iL));
           return res;
         }

         virtual PyObject *getTinySerializationInformation() const
         {
           std::vector<double> a0;
           std::vector<mcIdType> a1;
           std::vector<std::string> a2;
           self->getTinySerializationInformation(a0,a1,a2);
           PyObject *ret(PyTuple_New(3));
           PyTuple_SetItem(ret,0,convertDblArrToPyList2(a0));
           PyTuple_SetItem(ret,1,convertIntArrToPyList2(a1));
           std::size_t sz(a2.size());
           PyObject *ret2(PyList_New(sz));
           {
             for(std::size_t i=0;i<sz;i++)
               PyList_SetItem(ret2,i,PyString_FromString(a2[i].c_str()));
           }
           PyTuple_SetItem(ret,2,ret2);
           return ret;
         }

         virtual PyObject *serialize() const
         {
           DataArrayIdType *a0Tmp(0);
           DataArrayDouble *a1Tmp(0);
           self->serialize(a0Tmp,a1Tmp);
           PyObject *ret(PyTuple_New(2));
           PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(a0Tmp),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
           PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(a1Tmp),SWIGTYPE_p_MEDCoupling__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
           return ret;
         }

         void resizeForUnserialization(const std::vector<mcIdType>& tinyInfo, DataArrayIdType *a1, DataArrayDouble *a2) const
         {
           std::vector<std::string> littleStrings;
           self->resizeForUnserialization(tinyInfo,a1,a2,littleStrings);
         }

         PyObject *__getstate__() const
         {
           PyObject *ret0(MEDCoupling_MEDCouplingMesh_getTinySerializationInformation(self));
           PyObject *ret1(MEDCoupling_MEDCouplingMesh_serialize(self));
           PyObject *ret(PyTuple_New(2));
           PyTuple_SetItem(ret,0,ret0);
           PyTuple_SetItem(ret,1,ret1);
           return ret;
         }

         void __setstate__(PyObject *inp)
         {
           static const char MSG[]="MEDCouplingMesh.__setstate__ : expected input is a tuple of size 2 !";
           if(!PyTuple_Check(inp))
             throw INTERP_KERNEL::Exception(MSG);
           std::size_t sz(PyTuple_Size(inp));
           if(sz!=2)
             throw INTERP_KERNEL::Exception(MSG);
           PyObject *elt0(PyTuple_GetItem(inp,0));
           PyObject *elt1(PyTuple_GetItem(inp,1));
           std::vector<double> a0;
           std::vector<mcIdType> a1;
           std::vector<std::string> a2;
           DataArrayIdType *b0(0);
           DataArrayDouble *b1(0);
           {
             if(!PyTuple_Check(elt0) && PyTuple_Size(elt0)!=3)
               throw INTERP_KERNEL::Exception(MSG);
             PyObject *a0py(PyTuple_GetItem(elt0,0)),*a1py(PyTuple_GetItem(elt0,1)),*a2py(PyTuple_GetItem(elt0,2));
             mcIdType tmp(-1);
             fillArrayWithPyListDbl3(a0py,tmp,a0);
             convertPyToNewIntArr3(a1py,a1);
             fillStringVector(a2py,a2);
           }
           {
             if(!PyTuple_Check(elt1) && PyTuple_Size(elt1)!=2)
               throw INTERP_KERNEL::Exception(MSG);
             PyObject *b0py(PyTuple_GetItem(elt1,0)),*b1py(PyTuple_GetItem(elt1,1));
             void *argp(0);
             int status(SWIG_ConvertPtr(b0py,&argp,SWIGTITraits<mcIdType>::TI,0|0));
             if(!SWIG_IsOK(status))
               throw INTERP_KERNEL::Exception(MSG);
             b0=reinterpret_cast<DataArrayIdType *>(argp);
             status=SWIG_ConvertPtr(b1py,&argp,SWIGTYPE_p_MEDCoupling__DataArrayDouble,0|0);
             if(!SWIG_IsOK(status))
               throw INTERP_KERNEL::Exception(MSG);
             b1=reinterpret_cast<DataArrayDouble *>(argp);
           }
           // useless here to call resizeForUnserialization because arrays are well resized.
           self->unserialization(a0,a1,b0,b1,a2);
         }

         static MEDCouplingMesh *MergeMeshes(PyObject *li)
         {
            std::vector<const MEDCoupling::MEDCouplingMesh *> tmp;
            convertFromPyObjVectorOfObj<const MEDCoupling::MEDCouplingMesh *>(li,SWIGTYPE_p_MEDCoupling__MEDCouplingMesh,"MEDCouplingMesh",tmp);
            return MEDCouplingMesh::MergeMeshes(tmp);
         }
       }
  };
}

//== MEDCouplingMesh End

%include "NormalizedGeometricTypes"
%include "MEDCouplingNatureOfFieldEnum"
//
namespace MEDCoupling
{
  class MEDCouplingNatureOfField
  {
  public:
    static const char *GetRepr(NatureOfField nat);
    static std::string GetReprNoThrow(NatureOfField nat);
    static std::string GetAllPossibilitiesStr();
  };
}

// the MEDCouplingTimeDiscretization classes are not swigged : in case the file can help
// include "MEDCouplingTimeDiscretization.i"

namespace MEDCoupling
{
  class MEDCouplingGaussLocalization
  {
  public:
    MEDCouplingGaussLocalization(INTERP_KERNEL::NormalizedCellType type, const std::vector<double>& refCoo,
                                 const std::vector<double>& gsCoo, const std::vector<double>& w);
    MEDCouplingGaussLocalization(INTERP_KERNEL::NormalizedCellType typ);
    INTERP_KERNEL::NormalizedCellType getType() const;
    void setType(INTERP_KERNEL::NormalizedCellType typ);
    int getNumberOfGaussPt() const;
    int getDimension() const;
    int getNumberOfPtsInRefCell() const;
    std::string getStringRepr() const;
    void checkConsistencyLight() const;
    bool isEqual(const MEDCouplingGaussLocalization& other, double eps) const;
    //
    const std::vector<double>& getRefCoords() const;
    double getRefCoord(int ptIdInCell, int comp) const;
    const std::vector<double>& getGaussCoords() const;
    double getGaussCoord(int gaussPtIdInCell, int comp) const;
    const std::vector<double>& getWeights() const;
    double getWeight(int gaussPtIdInCell) const;
    void setRefCoord(int ptIdInCell, int comp, double newVal);
    void setGaussCoord(int gaussPtIdInCell, int comp, double newVal);
    void setWeight(int gaussPtIdInCell, double newVal);
    void setRefCoords(const std::vector<double>& refCoo);
    void setGaussCoords(const std::vector<double>& gsCoo);
    void setWeights(const std::vector<double>& w);
    //
    static bool AreAlmostEqual(const std::vector<double>& v1, const std::vector<double>& v2, double eps);
    //
    %extend
    {
      DataArrayDouble *localizePtsInRefCooForEachCell(const DataArrayDouble *ptsInRefCoo, const MEDCouplingUMesh *mesh) const
      {
        MCAuto<DataArrayDouble> ret(self->localizePtsInRefCooForEachCell(ptsInRefCoo,mesh));
        return ret.retn();
      }

      MEDCouplingUMesh *buildRefCell() const
      {
        MCAuto<MEDCouplingUMesh> ret(self->buildRefCell());
        return ret.retn();
      }

      DataArrayDouble *getShapeFunctionValues() const
      {
        MCAuto<DataArrayDouble> ret(self->getShapeFunctionValues());
        return ret.retn();
      }

      DataArrayDouble *getDerivativeOfShapeFunctionValues() const
      {
        MCAuto<DataArrayDouble> ret(self->getDerivativeOfShapeFunctionValues());
        return ret.retn();
      }

      static DataArrayDouble *GetDefaultReferenceCoordinatesOf(INTERP_KERNEL::NormalizedCellType type)
      {
        MCAuto<DataArrayDouble> ret(MEDCouplingGaussLocalization::GetDefaultReferenceCoordinatesOf(type));
        return ret.retn();
      }
    }
  };

  class MEDCouplingSkyLineArray
  {
  public:
    static MEDCouplingSkyLineArray *BuildFromPolyhedronConn( const DataArrayIdType* c, const DataArrayIdType* cI );

    void set( DataArrayIdType* index, DataArrayIdType* value );
    void set3( DataArrayIdType* superIndex, DataArrayIdType* index, DataArrayIdType* value );

    int getSuperNumberOf() const;
    int getNumberOf() const;
    int getLength() const;

    void deletePack(const int i, const int j);

    void deleteSimplePack(const int i);
    void deleteSimplePacks(const DataArrayIdType* idx);

    MEDCouplingSkyLineArray *groupPacks(const DataArrayIdType *indexedPacks) const;
    MEDCouplingSkyLineArray *uniqueNotSortedByPack() const;

    MEDCouplingSkyLineArray *deepCopy() const;

    %extend
    {
      MEDCouplingSkyLineArray()
      {
        return MEDCouplingSkyLineArray::New();
      }

      MEDCouplingSkyLineArray( const std::vector<mcIdType>& index, const std::vector<mcIdType>& value)
      {
        return MEDCouplingSkyLineArray::New(index, value);
      }

      MEDCouplingSkyLineArray( DataArrayIdType* index, DataArrayIdType* value )
      {
        return MEDCouplingSkyLineArray::New(index, value);
      }

      MEDCouplingSkyLineArray( const MEDCouplingSkyLineArray & other )
      {
        return MEDCouplingSkyLineArray::New(other);
      }

      std::string __str__() const
      {
        return self->simpleRepr();
      }

      DataArrayIdType *getSuperIndexArray() const
      {
        DataArrayIdType *ret(self->getSuperIndexArray());
        if(ret)
          ret->incrRef();
        return ret;
      }

      DataArrayIdType *getIndexArray() const
      {
        DataArrayIdType *ret(self->getIndexArray());
        if(ret)
          ret->incrRef();
        return ret;
      }

      DataArrayIdType *getValuesArray() const
      {
        DataArrayIdType *ret(self->getValuesArray());
        if(ret)
          ret->incrRef();
        return ret;
      }

      PyObject *getSimplePackSafe(mcIdType absolutePackId) const
      {
        std::vector<mcIdType> ret;
        self->getSimplePackSafe(absolutePackId,ret);
        return convertIntArrToPyList2(ret);
      }

      PyObject *findPackIds(PyObject *superPackIndices, PyObject *pack) const
      {
          std::vector<mcIdType> vpack, vspIdx, out;

          convertPyToNewIntArr3(superPackIndices,vspIdx);
          convertPyToNewIntArr3(pack,vpack);

          self->findPackIds(vspIdx, vpack.data(), vpack.data()+vpack.size(), out);
          return convertIntArrToPyList2(out);
      }

      void pushBackPack(const mcIdType i, PyObject *pack)
        {
          std::vector<mcIdType> vpack;
          convertPyToNewIntArr3(pack,vpack);
          self->pushBackPack(i,vpack.data(), vpack.data()+vpack.size());
        }

      void replaceSimplePack(const mcIdType idx, PyObject *pack)
        {
          std::vector<mcIdType> vpack;
          convertPyToNewIntArr3(pack,vpack);
          self->replaceSimplePack(idx, vpack.data(), vpack.data()+vpack.size());
        }

      void replaceSimplePacks(const DataArrayIdType* idx, PyObject *listePacks)
        {
          std::vector<const DataArrayIdType*> packs;
          convertFromPyObjVectorOfObj<const MEDCoupling::DataArrayIdType*>(listePacks,SWIGTITraits<mcIdType>::TI,"DataArrayIdType",packs);
          self->replaceSimplePacks(idx, packs);
        }

      static MEDCouplingSkyLineArray *AggregatePacks(PyObject *sks)
      {
        std::vector<const MEDCouplingSkyLineArray *> sksCpp;
        convertFromPyObjVectorOfObj<const MEDCoupling::MEDCouplingSkyLineArray*>(sks,SWIGTYPE_p_MEDCoupling__MEDCouplingSkyLineArray,"MEDCouplingSkyLineArray",sksCpp);
        return MEDCoupling::MEDCouplingSkyLineArray::AggregatePacks(sksCpp);
      }

      void replacePack(const mcIdType superIdx, const mcIdType idx, PyObject *pack)
        {
          std::vector<mcIdType> vpack;
          convertPyToNewIntArr3(pack,vpack);
          self->replacePack(superIdx, idx, vpack.data(), vpack.data()+vpack.size());
        }

      PyObject *convertToPolyhedronConn() const
         {
           MCAuto<DataArrayIdType> d0=DataArrayIdType::New();
           MCAuto<DataArrayIdType> d1=DataArrayIdType::New();
           self->convertToPolyhedronConn(d0,d1);
           PyObject *ret=PyTuple_New(2);
           PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(d0.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
           PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(d1.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
           return ret;
         }

      PyObject *thresholdPerPack(mcIdType threshold) const
      {
        MCAuto<MEDCouplingSkyLineArray> left, right;
        self->thresholdPerPack(threshold,left,right);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(left.retn()),SWIGTYPE_p_MEDCoupling__MEDCouplingSkyLineArray, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(right.retn()),SWIGTYPE_p_MEDCoupling__MEDCouplingSkyLineArray, SWIG_POINTER_OWN | 0 ));
        return ret;
      }
    }
  };
}

%include "MEDCouplingFieldDiscretization.i"

//== MEDCouplingPointSet

namespace MEDCoupling
{
  class MEDCouplingPointSet : public MEDCoupling::MEDCouplingMesh
    {
    public:
      void setCoords(const DataArrayDouble *coords);
      DataArrayDouble *getCoordinatesAndOwner() const;
      bool areCoordsEqual(const MEDCouplingPointSet& other, double prec) const;
      void zipCoords();
      double getCaracteristicDimension() const;
      void recenterForMaxPrecision(double eps);
      void changeSpaceDimension(int newSpaceDim, double dftVal=0.);
      void tryToShareSameCoords(const MEDCouplingPointSet& other, double epsilon);
      virtual void shallowCopyConnectivityFrom(const MEDCouplingPointSet *other);
      virtual MEDCouplingPointSet *buildPartOfMySelfSlice(int start, int end, int step) const;
      virtual void tryToShareSameCoordsPermute(const MEDCouplingPointSet& other, double epsilon);
      static DataArrayDouble *MergeNodesArray(const MEDCouplingPointSet *m1, const MEDCouplingPointSet *m2);
      static MEDCouplingPointSet *BuildInstanceFromMeshType(MEDCouplingMeshType type);
      static DataArrayIdType *ComputeNbOfInteractionsWithSrcCells(const MEDCouplingPointSet *srcMesh, const MEDCouplingPointSet *trgMesh, double eps);
      virtual DataArrayIdType *computeFetchedNodeIds() const;
      virtual int getNumberOfNodesInCell(int cellId) const;
      virtual MEDCouplingPointSet *buildBoundaryMesh(bool keepCoords) const;
      virtual DataArrayIdType *getCellsInBoundingBox(const INTERP_KERNEL::DirectedBoundingBox& bbox, double eps);
      virtual DataArrayIdType *zipCoordsTraducer();
      virtual DataArrayIdType *findBoundaryNodes() const;
      virtual DataArrayIdType *zipConnectivityTraducer(int compType, int startCellId=0);
      virtual MEDCouplingPointSet *mergeMyselfWithOnSameCoords(const MEDCouplingPointSet *other) const;
      virtual void checkFullyDefined() const;
      virtual bool isEmptyMesh(const std::vector<mcIdType>& tinyInfo) const;
      virtual MEDCouplingPointSet *deepCopyConnectivityOnly() const;
      virtual DataArrayDouble *getBoundingBoxForBBTree(double arcDetEps=1e-12) const;
      virtual void renumberNodesWithOffsetInConn(int offset);
      virtual bool areAllNodesFetched() const;
      virtual MEDCouplingFieldDouble *computeDiameterField() const;
      virtual void invertOrientationOfAllCells();
      %extend
         {
           std::string __str__() const
           {
             return self->simpleRepr();
           }

           PyObject *buildNewNumberingFromCommonNodesFormat(const DataArrayIdType *comm, const DataArrayIdType *commIndex) const
           {
             mcIdType newNbOfNodes;
             DataArrayIdType *ret0=self->buildNewNumberingFromCommonNodesFormat(comm,commIndex,newNbOfNodes);
             PyObject *res = PyList_New(2);
             PyList_SetItem(res,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
             PyList_SetItem(res,1,PyInt_FromLong(newNbOfNodes));
             return res;
           }

           PyObject *findCommonNodes(double prec, mcIdType limitTupleId=-1) const
           {
             DataArrayIdType *comm, *commIndex;
             self->findCommonNodes(prec,limitTupleId,comm,commIndex);
             PyObject *res = PyList_New(2);
             PyList_SetItem(res,0,SWIG_NewPointerObj(SWIG_as_voidptr(comm),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
             PyList_SetItem(res,1,SWIG_NewPointerObj(SWIG_as_voidptr(commIndex),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
             return res;
           }

           PyObject *getCoords()
           {
             DataArrayDouble *ret1=self->getCoords();
             if (ret1)
                ret1->incrRef();
             return SWIG_NewPointerObj((void*)ret1,SWIGTYPE_p_MEDCoupling__DataArrayDouble,SWIG_POINTER_OWN | 0);
           }

           PyObject *buildPartOfMySelf(PyObject *li, bool keepCoords=true) const
           {
             mcIdType szArr,sw,iTypppArr;
             std::vector<mcIdType> stdvecTyyppArr;
             const mcIdType *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
             MEDCouplingPointSet *ret=self->buildPartOfMySelf(tmp,tmp+szArr,keepCoords);
             if(sw==3)//DataArrayIdType
               {
                 void *argp; SWIG_ConvertPtr(li,&argp,SWIGTITraits<mcIdType>::TI,0|0);
                 DataArrayIdType *argpt=reinterpret_cast< MEDCoupling::DataArrayIdType * >(argp);
                 std::string name=argpt->getName();
                 if(!name.empty())
                   ret->setName(name.c_str());
               }
             return convertMesh(ret, SWIG_POINTER_OWN | 0 );
           }

           PyObject *buildPartOfMySelfNode(PyObject *li, bool fullyIn) const
           {
             mcIdType szArr,sw,iTypppArr;
             std::vector<mcIdType> stdvecTyyppArr;
             const mcIdType *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
             MEDCouplingPointSet *ret=self->buildPartOfMySelfNode(tmp,tmp+szArr,fullyIn);
             if(sw==3)//DataArrayIdType
               {
                 void *argp; SWIG_ConvertPtr(li,&argp,SWIGTITraits<mcIdType>::TI,0|0);
                 DataArrayIdType *argpt=reinterpret_cast< MEDCoupling::DataArrayIdType * >(argp);
                 std::string name=argpt->getName();
                 if(!name.empty())
                   ret->setName(name.c_str());
               }
             return convertMesh(ret, SWIG_POINTER_OWN | 0 );
           }

           virtual PyObject *buildPartOfMySelfKeepCoords(PyObject *li) const
           {
             mcIdType szArr,sw,iTypppArr;
             std::vector<mcIdType> stdvecTyyppArr;
             const mcIdType *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
             MEDCouplingPointSet *ret=self->buildPartOfMySelfKeepCoords(tmp,tmp+szArr);
             if(sw==3)//DataArrayIdType
               {
                 void *argp; SWIG_ConvertPtr(li,&argp,SWIGTITraits<mcIdType>::TI,0|0);
                 DataArrayIdType *argpt=reinterpret_cast< MEDCoupling::DataArrayIdType * >(argp);
                 std::string name=argpt->getName();
                 if(!name.empty())
                   ret->setName(name.c_str());
               }
             return convertMesh(ret, SWIG_POINTER_OWN | 0 );
           }

           virtual PyObject *buildPartOfMySelfKeepCoordsSlice(mcIdType start, mcIdType end, mcIdType step) const
           {
             MEDCouplingPointSet *ret=self->buildPartOfMySelfKeepCoordsSlice(start,end,step);
             return convertMesh(ret, SWIG_POINTER_OWN | 0 );
           }

           PyObject *buildFacePartOfMySelfNode(PyObject *li, bool fullyIn) const
           {
             mcIdType szArr,sw,iTypppArr;
             std::vector<mcIdType> stdvecTyyppArr;
             const mcIdType *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
             MEDCouplingPointSet *ret=self->buildFacePartOfMySelfNode(tmp,tmp+szArr,fullyIn);
             if(sw==3)//DataArrayIdType
               {
                 void *argp; SWIG_ConvertPtr(li,&argp,SWIGTITraits<mcIdType>::TI,0|0);
                 DataArrayIdType *argpt=reinterpret_cast< MEDCoupling::DataArrayIdType * >(argp);
                 std::string name=argpt->getName();
                 if(!name.empty())
                   ret->setName(name.c_str());
               }
             return convertMesh(ret, SWIG_POINTER_OWN | 0 );
           }

           void renumberNodes(PyObject *li, mcIdType newNbOfNodes)
           {
             mcIdType szArr,sw,iTypppArr;
             std::vector<mcIdType> stdvecTyyppArr;
             const mcIdType *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
             self->renumberNodes(tmp,newNbOfNodes);
           }

           void renumberNodesCenter(PyObject *li, mcIdType newNbOfNodes)
           {
             mcIdType szArr,sw,iTypppArr;
             std::vector<mcIdType> stdvecTyyppArr;
             const mcIdType *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
             self->renumberNodesCenter(tmp,newNbOfNodes);
           }

           PyObject *findNodesOnLine(PyObject *pt, PyObject *vec, double eps) const
             {
               int spaceDim=self->getSpaceDimension();
               double val,val2;
               DataArrayDouble *a,*a2;
               DataArrayDoubleTuple *aa,*aa2;
               std::vector<double> bb,bb2;
               mcIdType sw;
               const char msg[]="Python wrap of MEDCouplingPointSet::findNodesOnLine : 1st parameter for point.";
               const char msg2[]="Python wrap of MEDCouplingPointSet::findNodesOnLine : 2nd parameter for vector.";
               const double *p=convertObjToPossibleCpp5_Safe(pt,sw,val,a,aa,bb,msg,1,spaceDim,true);
               const double *v=convertObjToPossibleCpp5_Safe(vec,sw,val2,a2,aa2,bb2,msg2,1,spaceDim,true);
               std::vector<mcIdType> nodes;
               self->findNodesOnLine(p,v,eps,nodes);
               DataArrayIdType *ret=DataArrayIdType::New();
               ret->alloc(nodes.size(),1);
               std::copy(nodes.begin(),nodes.end(),ret->getPointer());
               return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 );
             }
           PyObject *findNodesOnPlane(PyObject *pt, PyObject *vec, double eps) const
             {
               int spaceDim=self->getSpaceDimension();
               double val,val2;
               DataArrayDouble *a,*a2;
               DataArrayDoubleTuple *aa,*aa2;
               std::vector<double> bb,bb2;
               mcIdType sw;
               const char msg[]="Python wrap of MEDCouplingPointSet::findNodesOnPlane : 1st parameter for point.";
               const char msg2[]="Python wrap of MEDCouplingPointSet::findNodesOnPlane : 2nd parameter for vector.";
               const double *p=convertObjToPossibleCpp5_Safe(pt,sw,val,a,aa,bb,msg,1,spaceDim,true);
               const double *v=convertObjToPossibleCpp5_Safe(vec,sw,val2,a2,aa2,bb2,msg2,1,spaceDim,true);
               std::vector<mcIdType> nodes;
               self->findNodesOnPlane(p,v,eps,nodes);
               DataArrayIdType *ret=DataArrayIdType::New();
               ret->alloc(nodes.size(),1);
               std::copy(nodes.begin(),nodes.end(),ret->getPointer());
               return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 );
             }

           PyObject *getNodeIdsNearPoint(PyObject *pt, double eps) const
           {
             double val;
             DataArrayDouble *a;
             DataArrayDoubleTuple *aa;
             std::vector<double> bb;
             mcIdType sw;
             int spaceDim=self->getSpaceDimension();
             const char msg[]="Python wrap of MEDCouplingPointSet::getNodeIdsNearPoint : ";
             const double *pos=convertObjToPossibleCpp5_Safe(pt,sw,val,a,aa,bb,msg,1,spaceDim,true);
             DataArrayIdType *ret=self->getNodeIdsNearPoint(pos,eps);
             return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 );
           }

           PyObject *getNodeIdsNearPoints(PyObject *pt, mcIdType nbOfPoints, double eps) const
           {
             DataArrayIdType *c=0,*cI=0;
             //
             double val;
             DataArrayDouble *a;
             DataArrayDoubleTuple *aa;
             std::vector<double> bb;
             mcIdType sw;
             int spaceDim=self->getSpaceDimension();
             const char msg[]="Python wrap of MEDCouplingPointSet::getNodeIdsNearPoints : ";
             const double *pos=convertObjToPossibleCpp5_Safe(pt,sw,val,a,aa,bb,msg,nbOfPoints,spaceDim,true);
             self->getNodeIdsNearPoints(pos,nbOfPoints,eps,c,cI);
             PyObject *ret=PyTuple_New(2);
             PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(c),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
             PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(cI),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
             return ret;
           }

           PyObject *getNodeIdsNearPoints(PyObject *pt, double eps) const
           {
             DataArrayIdType *c=0,*cI=0;
             int spaceDim=self->getSpaceDimension();
             double val;
             DataArrayDouble *a;
             DataArrayDoubleTuple *aa;
             std::vector<double> bb;
             mcIdType sw;
             mcIdType nbOfTuples=-1;
             const double *ptPtr=convertObjToPossibleCpp5_Safe2(pt,sw,val,a,aa,bb,"Python wrap of MEDCouplingUMesh::getNodeIdsNearPoints",spaceDim,true,nbOfTuples);
             self->getNodeIdsNearPoints(ptPtr,nbOfTuples,eps,c,cI);
             //
             PyObject *ret=PyTuple_New(2);
             PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(c),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
             PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(cI),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
             return ret;
           }

           PyObject *getCellsInBoundingBox(PyObject *bbox, double eps) const
           {
             double val;
             DataArrayDouble *a;
             DataArrayDoubleTuple *aa;
             std::vector<double> bb;
             mcIdType sw;
             int spaceDim=self->getSpaceDimension();
             const char msg[]="Python wrap of MEDCouplingPointSet::getCellsInBoundingBox : ";
             const double *tmp=convertObjToPossibleCpp5_Safe(bbox,sw,val,a,aa,bb,msg,spaceDim,2,true);
             //
             DataArrayIdType *elems=self->getCellsInBoundingBox(tmp,eps);
             return SWIG_NewPointerObj(SWIG_as_voidptr(elems),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 );
           }

           void duplicateNodesInCoords(PyObject *li)
           {
             mcIdType sw;
             mcIdType singleVal;
             std::vector<mcIdType> multiVal;
             std::pair<mcIdType, std::pair<mcIdType,mcIdType> > slic;
             MEDCoupling::DataArrayIdType *daIntTyypp=0;
             convertIntStarOrSliceLikePyObjToCpp(li,self->getNumberOfNodes(),sw,singleVal,multiVal,slic,daIntTyypp);
             switch(sw)
               {
               case 1:
                 return self->duplicateNodesInCoords(&singleVal,&singleVal+1);
               case 2:
                 return self->duplicateNodesInCoords(&multiVal[0],&multiVal[0]+multiVal.size());
               case 4:
                 return self->duplicateNodesInCoords(daIntTyypp->begin(),daIntTyypp->end());
               default:
                 throw INTERP_KERNEL::Exception("MEDCouplingPointSet::duplicateNodesInCoords : unrecognized type entered, expected list of int, tuple of int or DataArrayIdType !");
               }
           }

           virtual PyObject *findCommonCells(int compType, mcIdType startCellId=0) const
           {
             DataArrayIdType *v0(nullptr),*v1(nullptr);
             self->findCommonCells(compType,startCellId,v0,v1);
             PyObject *res = PyList_New(2);
             PyList_SetItem(res,0,SWIG_NewPointerObj(SWIG_as_voidptr(v0),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
             PyList_SetItem(res,1,SWIG_NewPointerObj(SWIG_as_voidptr(v1),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
             return res;
           }


           virtual void renumberNodesInConn(PyObject *li)
           {
             void *da(nullptr);
             {
               int res1(SWIG_ConvertPtr(li,&da,SWIGTYPE_p_MEDCoupling__MapII, 0 |  0 ));
               if(SWIG_IsOK(res1))
                 {
                   MapII *da2(reinterpret_cast<MapII *>(da));
                   self->renumberNodesInConn(da2->data());
                   return ;
                 }
             }
             int res1(SWIG_ConvertPtr(li,&da,SWIGTITraits<mcIdType>::TI, 0 | 0 ));
             if (!SWIG_IsOK(res1))
               {
                 mcIdType size;
                 INTERP_KERNEL::AutoPtr<mcIdType> tmp=convertPyToNewIntArr2(li,&size);
                 self->renumberNodesInConn(tmp);
               }
             else
               {
                 DataArrayIdType *da2(reinterpret_cast< DataArrayIdType * >(da));
                 if(!da2)
                   throw INTERP_KERNEL::Exception("Not null DataArrayIdType instance expected !");
                 da2->checkAllocated();
                 self->renumberNodesInConn(da2->getConstPointer());
               }
           }

           virtual PyObject *getNodeIdsInUse() const
           {
             mcIdType ret1=-1;
             DataArrayIdType *ret0=self->getNodeIdsInUse(ret1);
             PyObject *ret=PyTuple_New(2);
             PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
             PyTuple_SetItem(ret,1,PyInt_FromLong(ret1));
             return ret;
           }

           virtual DataArrayIdType *fillCellIdsToKeepFromNodeIds(PyObject *li, bool fullyIn) const
           {
             DataArrayIdType *ret(nullptr);
             //
             mcIdType szArr,sw,iTypppArr;
             std::vector<mcIdType> stdvecTyyppArr;
             const mcIdType *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
             self->fillCellIdsToKeepFromNodeIds(tmp,tmp+szArr,fullyIn,ret);
             return ret;
           }

           virtual PyObject *mergeNodes(double precision)
           {
             bool ret1;
             mcIdType ret2;
             DataArrayIdType *ret0=self->mergeNodes(precision,ret1,ret2);
             PyObject *res = PyList_New(3);
             PyList_SetItem(res,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
             PyList_SetItem(res,1,SWIG_From_bool(ret1));
             PyList_SetItem(res,2,PyInt_FromLong(ret2));
             return res;
           }

           virtual PyObject *mergeNodesCenter(double precision)
           {
             bool ret1;
             mcIdType ret2;
             DataArrayIdType *ret0=self->mergeNodesCenter(precision,ret1,ret2);
             PyObject *res = PyList_New(3);
             PyList_SetItem(res,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
             PyList_SetItem(res,1,SWIG_From_bool(ret1));
             PyList_SetItem(res,2,PyInt_FromLong(ret2));
             return res;
           }

           DataArrayIdType *getCellIdsLyingOnNodes(PyObject *li, bool fullyIn) const
           {
             void *da=0;
             int res1=SWIG_ConvertPtr(li,&da,SWIGTITraits<mcIdType>::TI, 0 |  0 );
             if (!SWIG_IsOK(res1))
               {
                 mcIdType size;
                 INTERP_KERNEL::AutoPtr<mcIdType> tmp=convertPyToNewIntArr2(li,&size);
                 return self->getCellIdsLyingOnNodes(tmp,((const mcIdType *)tmp)+size,fullyIn);
               }
             else
               {
                 DataArrayIdType *da2=reinterpret_cast< DataArrayIdType * >(da);
                 if(!da2)
                   throw INTERP_KERNEL::Exception("Not null DataArrayIdType instance expected !");
                 da2->checkAllocated();
                 return self->getCellIdsLyingOnNodes(da2->getConstPointer(),da2->getConstPointer()+da2->getNbOfElems(),fullyIn);
               }
           }

           MEDCouplingPointSet *__getitem__(PyObject *listOrDataArrI)
           {
             mcIdType sw;
             mcIdType singleVal;
             std::vector<mcIdType> multiVal;
             std::pair<mcIdType, std::pair<mcIdType,mcIdType> > slic;
             MEDCoupling::DataArrayIdType *daIntTyypp=0;
             mcIdType nbc=self->getNumberOfCells();
             convertIntStarOrSliceLikePyObjToCpp(listOrDataArrI,nbc,sw,singleVal,multiVal,slic,daIntTyypp);
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
                           mcIdType tmp=nbc+singleVal;
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
                   return static_cast<MEDCouplingPointSet *>(self->buildPartOfMySelf(&multiVal[0],&multiVal[0]+multiVal.size(),true));
                 }
               case 3:
                 {
                   return self->buildPartOfMySelfSlice(slic.first,slic.second.first,slic.second.second,true);
                 }
               case 4:
                 {
                   if(!daIntTyypp)
                     throw INTERP_KERNEL::Exception("MEDCouplingUMesh::__getitem__ : null instance has been given in input !");
                   daIntTyypp->checkAllocated();
                   return self->buildPartOfMySelf(daIntTyypp->begin(),daIntTyypp->end(),true);
                 }
               default:
                 throw INTERP_KERNEL::Exception("MEDCouplingUMesh::__getitem__ : unrecognized type in input ! Possibilities are : int, list or tuple of int DataArrayIdType instance !");
               }
           }

           static void Rotate2DAlg(PyObject *center, double angle, mcIdType nbNodes, PyObject *coords)
           {
             mcIdType sz;
             INTERP_KERNEL::AutoCPtr<double> c=convertPyToNewDblArr2(center,&sz);
             INTERP_KERNEL::AutoCPtr<double> coo=convertPyToNewDblArr2(coords,&sz);
             MEDCoupling::DataArrayDouble::Rotate2DAlg(c,angle,nbNodes,coo,coo);
             for(mcIdType i=0;i<sz;i++)
               PyList_SetItem(coords,i,PyFloat_FromDouble(coo[i]));
           }

           static void Rotate2DAlg(PyObject *center, double angle, PyObject *coords)
           {
             mcIdType sz;
             INTERP_KERNEL::AutoCPtr<double> c=convertPyToNewDblArr2(center,&sz);
             mcIdType sw,nbNodes=0;
             double val0;  MEDCoupling::DataArrayDouble *val1=0; MEDCoupling::DataArrayDoubleTuple *val2=0;
             std::vector<double> val3;
             const double *coo=convertObjToPossibleCpp5_Safe2(coords,sw,val0,val1,val2,val3,
                                                            "Rotate2DAlg",2,true,nbNodes);
             if(sw!=2 && sw!=3)
               throw INTERP_KERNEL::Exception("Invalid call to MEDCouplingPointSet::Rotate2DAlg : try another overload method !");
             MEDCoupling::DataArrayDouble::Rotate2DAlg(c,angle,nbNodes,coo,const_cast<double *>(coo));
           }

           static void Rotate3DAlg(PyObject *center, PyObject *vect, double angle, mcIdType nbNodes, PyObject *coords)
           {
             mcIdType sz,sz2;
             INTERP_KERNEL::AutoCPtr<double> c=convertPyToNewDblArr2(center,&sz);
             INTERP_KERNEL::AutoCPtr<double> coo=convertPyToNewDblArr2(coords,&sz);
             INTERP_KERNEL::AutoCPtr<double> v=convertPyToNewDblArr2(vect,&sz2);
             MEDCoupling::DataArrayDouble::Rotate3DAlg(c,v,angle,nbNodes,coo,coo);
             for(mcIdType i=0;i<sz;i++)
               PyList_SetItem(coords,i,PyFloat_FromDouble(coo[i]));
           }

           static void Rotate3DAlg(PyObject *center, PyObject *vect, double angle, PyObject *coords)
           {
             mcIdType sz,sz2;
             INTERP_KERNEL::AutoCPtr<double> c=convertPyToNewDblArr2(center,&sz);
             mcIdType sw,nbNodes=0;
             double val0;  MEDCoupling::DataArrayDouble *val1=0; MEDCoupling::DataArrayDoubleTuple *val2=0;
             std::vector<double> val3;
             const double *coo=convertObjToPossibleCpp5_Safe2(coords,sw,val0,val1,val2,val3,
                                                            "Rotate3DAlg",3,true,nbNodes);
             if(sw!=2 && sw!=3)
               throw INTERP_KERNEL::Exception("Invalid call to MEDCouplingPointSet::Rotate3DAlg : try another overload method !");
             INTERP_KERNEL::AutoCPtr<double> v=convertPyToNewDblArr2(vect,&sz2);
             MEDCoupling::DataArrayDouble::Rotate3DAlg(c,v,angle,nbNodes,coo,const_cast<double *>(coo));
           }
         }
    };

  //== MEDCouplingPointSet End

  class MEDCouplingUMeshCell
  {
  public:
    INTERP_KERNEL::NormalizedCellType getType() const;
    %extend
      {
        std::string __str__() const
        {
          return self->repr();
        }

        PyObject *getAllConn() const
        {
          mcIdType ret2;
          const mcIdType *r=self->getAllConn(ret2);
          PyObject *ret=PyTuple_New(ret2);
          for(mcIdType i=0;i<ret2;i++)
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
            return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_MEDCoupling__MEDCouplingUMeshCell,0|0);
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
            return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_MEDCoupling__MEDCouplingUMeshCellEntry,SWIG_POINTER_OWN | 0);
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

  //== MEDCouplingUMesh

  class MEDCouplingUMesh : public MEDCoupling::MEDCouplingPointSet
  {
  public:
    static MEDCouplingUMesh *New();
    static MEDCouplingUMesh *New(const char *meshName, int meshDim);
    void checkConsistencyLight() const;
    void checkGeomConsistency(double eps=1e-12) const;
    void setMeshDimension(int meshDim);
    void allocateCells(int nbOfCells=0);
    void finishInsertingCells();
    MEDCouplingUMeshCellByTypeEntry *cellsByType();
    void setConnectivity(DataArrayIdType *conn, DataArrayIdType *connIndex, bool isComputingTypes=true);
    INTERP_KERNEL::NormalizedCellType getTypeOfCell(int cellId) const;
    void setPartOfMySelfSlice(int start, int end, int step, const MEDCouplingUMesh& otherOnSameCoordsThanThis);
    int getNodalConnectivityArrayLen() const;
    void computeTypes();
    std::string reprConnectivityOfThis() const;
    MEDCouplingUMesh *buildSetInstanceFromThis(int spaceDim) const;
    //tools
    DataArrayIdType *conformize2D(double eps);
    DataArrayIdType *conformize3D(double eps);
    DataArrayIdType *colinearize2D(double eps);
    DataArrayIdType *colinearizeKeepingConform2D(double eps);
    void shiftNodeNumbersInConn(int delta);
    std::vector<bool> getQuadraticStatus() const;
    DataArrayIdType *findCellIdsOnBoundary() const;
    MEDCouplingUMesh *computeSkin() const;
    bool checkConsecutiveCellTypes() const;
    bool checkConsecutiveCellTypesForMEDFileFrmt() const;
    DataArrayIdType *rearrange2ConsecutiveCellTypes();
    DataArrayIdType *sortCellsInMEDFileFrmt();
    DataArrayIdType *getRenumArrForMEDFileFrmt() const;
    DataArrayIdType *convertCellArrayPerGeoType(const DataArrayIdType *da) const;
    MEDCouplingUMesh *buildDescendingConnectivity(DataArrayIdType *desc, DataArrayIdType *descIndx, DataArrayIdType *revDesc, DataArrayIdType *revDescIndx) const;
    MEDCouplingUMesh *buildDescendingConnectivity2(DataArrayIdType *desc, DataArrayIdType *descIndx, DataArrayIdType *revDesc, DataArrayIdType *revDescIndx) const;
    MEDCouplingUMesh *explode3DMeshTo1D(DataArrayIdType *desc, DataArrayIdType *descIndx, DataArrayIdType *revDesc, DataArrayIdType *revDescIndx) const;
    MEDCouplingUMesh *explodeMeshIntoMicroEdges(DataArrayIdType *desc, DataArrayIdType *descIndx, DataArrayIdType *revDesc, DataArrayIdType *revDescIndx) const;
    void orientCorrectlyPolyhedrons();
    bool isPresenceOfQuadratic() const;
    bool isFullyQuadratic() const;
    MEDCouplingFieldDouble *buildDirectionVectorField() const;
    bool isContiguous1D() const;
    void tessellate2D(double eps);
    void convertQuadraticCellsToLinear();
    DataArrayIdType *convertLinearCellsToQuadratic(int conversionType=0);
    void convertDegeneratedCells();
    DataArrayIdType *convertDegeneratedCellsAndRemoveFlatOnes();
    bool removeDegenerated1DCells();
    bool areOnlySimplexCells() const;
    MEDCouplingFieldDouble *getEdgeRatioField() const;
    MEDCouplingFieldDouble *getAspectRatioField() const;
    MEDCouplingFieldDouble *getWarpField() const;
    MEDCouplingFieldDouble *getSkewField() const;
    DataArrayDouble *computePlaneEquationOf3DFaces() const;
    DataArrayIdType *convexEnvelop2D();
    std::string cppRepr() const;
    DataArrayIdType *findAndCorrectBadOriented3DExtrudedCells();
    DataArrayIdType *findAndCorrectBadOriented3DCells();
    MEDCoupling::MEDCoupling1GTUMesh *convertIntoSingleGeoTypeMesh() const;
    MEDCouplingSkyLineArray *generateGraph() const;
    DataArrayIdType *convertNodalConnectivityToStaticGeoTypeMesh() const;
    DataArrayIdType *buildUnionOf2DMesh() const;
    DataArrayIdType *buildUnionOf3DMesh() const;
    DataArrayIdType *orderConsecutiveCells1D() const;
    DataArrayDouble *getBoundingBoxForBBTreeFast() const;
    DataArrayDouble *getBoundingBoxForBBTree2DQuadratic(double arcDetEps=1e-12) const;
    DataArrayDouble *getBoundingBoxForBBTree1DQuadratic(double arcDetEps=1e-12) const;
    void changeOrientationOfCells();
    void orientCorrectly2DCells(const MEDCouplingUMesh *refFaces);
    void orientCorrectly3DCells();
    DataArrayDouble *computeCellCenterOfMassWithPrecision(double eps);
    int split2DCells(const DataArrayIdType *desc, const DataArrayIdType *descI, const DataArrayIdType *subNodesInSeg, const DataArrayIdType *subNodesInSegI, const DataArrayIdType *midOpt=0, const DataArrayIdType *midOptI=0);
    static MEDCouplingUMesh *Build0DMeshFromCoords(DataArrayDouble *da);
    static MEDCouplingUMesh *MergeUMeshes(const MEDCouplingUMesh *mesh1, const MEDCouplingUMesh *mesh2);
    static MEDCouplingUMesh *MergeUMeshesOnSameCoords(const MEDCouplingUMesh *mesh1, const MEDCouplingUMesh *mesh2);
    static DataArrayIdType *ComputeSpreadZoneGradually(const DataArrayIdType *arrIn, const DataArrayIdType *arrIndxIn);
    static DataArrayIdType *ComputeRangesFromTypeDistribution(const std::vector<mcIdType>& code);
    %extend {
      MEDCouplingUMesh()
      {
        return MEDCouplingUMesh::New();
      }

      MEDCouplingUMesh(const char *meshName, int meshDim)
      {
        return MEDCouplingUMesh::New(meshName,meshDim);
      }

      std::string __str__() const
      {
        return self->simpleRepr();
      }

      std::string __repr__() const
      {
        std::ostringstream oss;
        self->reprQuickOverview(oss);
        return oss.str();
      }

      MEDCouplingUMeshCellIterator *__iter__()
      {
        return self->cellIterator();
      }

      static MEDCouplingUMesh *Build1DMeshFromCoords(DataArrayDouble *da)
      {
        MCAuto<MEDCouplingUMesh> ret(MEDCouplingUMesh::Build1DMeshFromCoords(da));
        return ret.retn();
      }

      MEDCouplingUMesh *convertToQuadraticBasedOnSeg3(MEDCoupling1SGTUMesh *seg3) const
      {
        MCAuto<MEDCouplingUMesh> ret( self->convertToQuadraticBasedOnSeg3(seg3) );
        return ret.retn();
      }

      MEDCouplingUMesh *extrudeConnectivity(mcIdType nbOfCellsToExtrude) const
      {
        MCAuto<MEDCouplingUMesh> ret( self->extrudeConnectivity(nbOfCellsToExtrude) );
        return ret.retn();
      }

      PyObject *getAllGeoTypesSorted() const
      {
        std::vector<INTERP_KERNEL::NormalizedCellType> result=self->getAllGeoTypesSorted();
        std::vector<INTERP_KERNEL::NormalizedCellType>::const_iterator iL=result.begin();
        PyObject *res=PyList_New(result.size());
        for(int i=0;iL!=result.end(); i++, iL++)
          PyList_SetItem(res,i,PyInt_FromLong(*iL));
        return res;
      }

      void setPartOfMySelf(PyObject *li, const MEDCouplingUMesh& otherOnSameCoordsThanThis)
      {
        mcIdType sw;
        mcIdType singleVal;
        std::vector<mcIdType> multiVal;
        std::pair<mcIdType, std::pair<mcIdType,mcIdType> > slic;
        MEDCoupling::DataArrayIdType *daIntTyypp=0;
        mcIdType nbc=self->getNumberOfCells();
        convertIntStarOrSliceLikePyObjToCpp(li,nbc,sw,singleVal,multiVal,slic,daIntTyypp);
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
                      mcIdType tmp=nbc+singleVal;
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
            throw INTERP_KERNEL::Exception("MEDCouplingUMesh::setPartOfMySelf : unrecognized type in input ! Possibilities are : int, list or tuple of int DataArrayIdType instance !");
          }
      }

      void __setitem__(PyObject *li, const MEDCouplingUMesh& otherOnSameCoordsThanThis)
      {
        mcIdType sw;
        mcIdType singleVal;
        std::vector<mcIdType> multiVal;
        std::pair<mcIdType, std::pair<mcIdType,mcIdType> > slic;
        MEDCoupling::DataArrayIdType *daIntTyypp=0;
        mcIdType nbc=self->getNumberOfCells();
        convertIntStarOrSliceLikePyObjToCpp(li,nbc,sw,singleVal,multiVal,slic,daIntTyypp);
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
                      mcIdType tmp=nbc+singleVal;
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
              self->setPartOfMySelfSlice(slic.first,slic.second.first,slic.second.second,otherOnSameCoordsThanThis);
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
            throw INTERP_KERNEL::Exception("MEDCouplingUMesh::__setitem__ : unrecognized type in input ! Possibilities are : int, list or tuple of int, slice, DataArrayIdType instance !");
          }
      }

      void insertNextCell(INTERP_KERNEL::NormalizedCellType type, mcIdType size, PyObject *li)
      {
        mcIdType szArr,sw,iTypppArr;
        std::vector<mcIdType> stdvecTyyppArr;
        const mcIdType *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
        if(size>szArr)
          {
            std::ostringstream oss; oss << "Wrap of MEDCouplingUMesh::insertNextCell : request of connectivity with length " << size << " whereas the length of input is " << szArr << " !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
        self->insertNextCell(type,size,tmp);
      }

      void insertNextCell(INTERP_KERNEL::NormalizedCellType type, PyObject *li)
      {
        mcIdType szArr,sw,iTypppArr;
        std::vector<mcIdType> stdvecTyyppArr;
        const mcIdType *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
        self->insertNextCell(type,szArr,tmp);
      }

      DataArrayIdType *getNodalConnectivity()
      {
        DataArrayIdType *ret=self->getNodalConnectivity();
        if(ret)
          ret->incrRef();
        return ret;
      }
      DataArrayIdType *getNodalConnectivityIndex()
      {
        DataArrayIdType *ret=self->getNodalConnectivityIndex();
        if(ret)
          ret->incrRef();
        return ret;
      }

      static PyObject *ComputeSpreadZoneGraduallyFromSeed(PyObject *seed, const DataArrayIdType *arrIn, const DataArrayIdType *arrIndxIn, mcIdType nbOfDepthPeeling=-1)
      {
        mcIdType szArr,sw,iTypppArr;
        std::vector<mcIdType> stdvecTyyppArr;
        const mcIdType *seedPtr=convertIntStarLikePyObjToCppIntStar(seed,sw,szArr,iTypppArr,stdvecTyyppArr);
        mcIdType nbOfDepthPeelingPerformed=0;
        DataArrayIdType *ret0=MEDCouplingUMesh::ComputeSpreadZoneGraduallyFromSeed(seedPtr,seedPtr+szArr,arrIn,arrIndxIn,nbOfDepthPeeling,nbOfDepthPeelingPerformed);
        PyObject *res=PyTuple_New(2);
        PyTuple_SetItem(res,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(res,1,PyInt_FromLong(nbOfDepthPeelingPerformed));
        return res;
      }

      static PyObject *FindCommonCellsAlg(int compType, mcIdType startCellId, const DataArrayIdType *nodal, const DataArrayIdType *nodalI, const DataArrayIdType *revNodal, const DataArrayIdType *revNodalI)
      {
        DataArrayIdType *v0=0,*v1=0;
        MEDCouplingUMesh::FindCommonCellsAlg(compType,startCellId,nodal,nodalI,revNodal,revNodalI,v0,v1);
        PyObject *res = PyList_New(2);
        PyList_SetItem(res,0,SWIG_NewPointerObj(SWIG_as_voidptr(v0),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyList_SetItem(res,1,SWIG_NewPointerObj(SWIG_as_voidptr(v1),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return res;
      }

      PyObject *distanceToPoint(PyObject *point) const
      {
        double val;
        DataArrayDouble *a;
        DataArrayDoubleTuple *aa;
        std::vector<double> bb;
        mcIdType sw;
        int nbOfCompo=self->getSpaceDimension();
        const double *pt=convertObjToPossibleCpp5_Safe(point,sw,val,a,aa,bb,"Python wrap of MEDCouplingUMesh::distanceToPoint",1,nbOfCompo,true);
        //
        mcIdType cellId=-1;
        double ret0=self->distanceToPoint(pt,pt+nbOfCompo,cellId);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,PyFloat_FromDouble(ret0));
        PyTuple_SetItem(ret,1,PyInt_FromLong(cellId));
        return ret;
      }

      PyObject *distanceToPoints(const DataArrayDouble *pts) const
      {
        DataArrayIdType *ret1=0;
        DataArrayDouble *ret0=self->distanceToPoints(pts,ret1);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_MEDCoupling__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *tetrahedrize(int policy)
      {
        mcIdType ret2(-1);
        DataArrayIdType *ret1(0);
        MEDCoupling1SGTUMesh *ret0(self->tetrahedrize(policy,ret1,ret2));
        PyObject *ret=PyTuple_New(3);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_MEDCoupling__MEDCoupling1SGTUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,PyInt_FromLong(ret2));
        return ret;
      }

      PyObject *checkButterflyCells(double eps=1e-12)
      {
        std::vector<mcIdType> cells;
        self->checkButterflyCells(cells,eps);
        DataArrayIdType *ret=DataArrayIdType::New();
        ret->alloc(cells.size(),1);
        std::copy(cells.begin(),cells.end(),ret->getPointer());
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 );
      }

      PyObject *splitByType() const
      {
        std::vector<MEDCouplingUMesh *> ms=self->splitByType();
        std::size_t sz=ms.size();
        PyObject *ret = PyList_New(sz);
        for(std::size_t i=0;i<sz;i++)
          PyList_SetItem(ret,i,SWIG_NewPointerObj(SWIG_as_voidptr(ms[i]),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *partitionBySpreadZone() const
      {
        std::vector<DataArrayIdType *> retCpp=self->partitionBySpreadZone();
        std::size_t sz=retCpp.size();
        PyObject *ret=PyList_New(sz);
        for(std::size_t i=0;i<sz;i++)
          PyList_SetItem(ret,i,SWIG_NewPointerObj(SWIG_as_voidptr(retCpp[i]),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      static PyObject *PartitionBySpreadZone(const DataArrayIdType *arrIn, const DataArrayIdType *arrIndxIn)
      {
        std::vector<DataArrayIdType *> retCpp(MEDCouplingUMesh::PartitionBySpreadZone(arrIn,arrIndxIn));
        std::size_t sz=retCpp.size();
        PyObject *ret=PyList_New(sz);
        for(std::size_t i=0;i<sz;i++)
          PyList_SetItem(ret,i,SWIG_NewPointerObj(SWIG_as_voidptr(retCpp[i]),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *keepSpecifiedCells(INTERP_KERNEL::NormalizedCellType type, PyObject *ids) const
      {
        mcIdType size;
        INTERP_KERNEL::AutoPtr<mcIdType> tmp=convertPyToNewIntArr2(ids,&size);
        MEDCouplingUMesh *ret=self->keepSpecifiedCells(type,tmp,tmp+size);
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 );
      }

      bool checkConsecutiveCellTypesAndOrder(PyObject *li) const
      {
        mcIdType sz;
        INTERP_KERNEL::AutoPtr<INTERP_KERNEL::NormalizedCellType> order=convertPyToNewIntArr2<INTERP_KERNEL::NormalizedCellType>(li,&sz);
        bool ret=self->checkConsecutiveCellTypesAndOrder(order,order+sz);
        return ret;
      }

      DataArrayIdType *getRenumArrForConsecutiveCellTypesSpec(PyObject *li) const
      {
        mcIdType sz;
        INTERP_KERNEL::AutoPtr<INTERP_KERNEL::NormalizedCellType> order=convertPyToNewIntArr2<INTERP_KERNEL::NormalizedCellType>(li,&sz);
        DataArrayIdType *ret=self->getRenumArrForConsecutiveCellTypesSpec(order,(INTERP_KERNEL::NormalizedCellType *)order+sz);
        return ret;
      }

      DataArrayIdType *findNodesToDuplicate(const MEDCouplingUMesh& otherDimM1OnSameCoords) const
      {
        DataArrayIdType *ret=self->findNodesToDuplicate(otherDimM1OnSameCoords);
        return ret;
      }

      PyObject *findCellsToRenumber(const MEDCouplingUMesh& otherDimM1OnSameCoords, const DataArrayIdType *dupNodes) const
      {
        DataArrayIdType *tmp0=0,*tmp1=0;
        self->findCellsToRenumber(otherDimM1OnSameCoords,dupNodes->begin(), dupNodes->end(), tmp0,tmp1);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(tmp0),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(tmp1),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *findCellIdsLyingOn(const MEDCouplingUMesh& otherDimM1OnSameCoords) const
      {
        DataArrayIdType *tmp0=0,*tmp1=0;
        self->findCellIdsLyingOn(otherDimM1OnSameCoords,tmp0,tmp1);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(tmp0),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(tmp1),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      void duplicateNodes(PyObject *li)
      {
        mcIdType sw;
        mcIdType singleVal;
        std::vector<mcIdType> multiVal;
        std::pair<mcIdType, std::pair<mcIdType,mcIdType> > slic;
        MEDCoupling::DataArrayIdType *daIntTyypp=0;
        convertIntStarOrSliceLikePyObjToCpp(li,self->getNumberOfNodes(),sw,singleVal,multiVal,slic,daIntTyypp);
        switch(sw)
          {
          case 1:
            return self->duplicateNodes(&singleVal,&singleVal+1);
          case 2:
            return self->duplicateNodes(&multiVal[0],&multiVal[0]+multiVal.size());
          case 4:
            return self->duplicateNodes(daIntTyypp->begin(),daIntTyypp->end());
          default:
            throw INTERP_KERNEL::Exception("MEDCouplingUMesh::duplicateNodes : unrecognized type entered, expected list of int, tuple of int or DataArrayIdType !");
          }
      }

      void duplicateNodesInConn(PyObject *li, mcIdType offset)
      {
        mcIdType sw;
        mcIdType singleVal;
        std::vector<mcIdType> multiVal;
        std::pair<mcIdType, std::pair<mcIdType,mcIdType> > slic;
        MEDCoupling::DataArrayIdType *daIntTyypp=0;
        convertIntStarOrSliceLikePyObjToCpp(li,self->getNumberOfNodes(),sw,singleVal,multiVal,slic,daIntTyypp);
        switch(sw)
          {
          case 1:
            return self->duplicateNodesInConn(&singleVal,&singleVal+1,offset);
          case 2:
            return self->duplicateNodesInConn(&multiVal[0],&multiVal[0]+multiVal.size(),offset);
          case 4:
            return self->duplicateNodesInConn(daIntTyypp->begin(),daIntTyypp->end(),offset);
          default:
            throw INTERP_KERNEL::Exception("MEDCouplingUMesh::duplicateNodesInConn : unrecognized type entered, expected list of int, tuple of int or DataArrayIdType !");
          }
      }

      void attractSeg3MidPtsAroundNodes(double ratio, PyObject *nodeIds)
      {
        mcIdType szArr,sw,iTypppArr;
        std::vector<mcIdType> stdvecTyyppArr;
        const mcIdType *nodeIdsPtr(convertIntStarLikePyObjToCppIntStar(nodeIds,sw,szArr,iTypppArr,stdvecTyyppArr));
        self->attractSeg3MidPtsAroundNodes(ratio,nodeIdsPtr,nodeIdsPtr+szArr);
      }

      PyObject *getLevArrPerCellTypes(PyObject *li) const
      {
        mcIdType sz;
        INTERP_KERNEL::AutoPtr<INTERP_KERNEL::NormalizedCellType> order=convertPyToNewIntArr2<INTERP_KERNEL::NormalizedCellType>(li,&sz);
        DataArrayIdType *tmp0,*tmp1=0;
        tmp0=self->getLevArrPerCellTypes(order,(INTERP_KERNEL::NormalizedCellType *)order+sz,tmp1);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(tmp0),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(tmp1),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *convertNodalConnectivityToDynamicGeoTypeMesh() const
      {
        DataArrayIdType *ret0=0,*ret1=0;
        self->convertNodalConnectivityToDynamicGeoTypeMesh(ret0,ret1);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      static PyObject *AggregateSortedByTypeMeshesOnSameCoords(PyObject *ms)
      {
        std::vector<const MEDCoupling::MEDCouplingUMesh *> meshes;
        convertFromPyObjVectorOfObj<const MEDCoupling::MEDCouplingUMesh *>(ms,SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh,"MEDCouplingUMesh",meshes);
        DataArrayIdType *ret1=0,*ret2=0;
        MEDCouplingUMesh *ret0=MEDCouplingUMesh::AggregateSortedByTypeMeshesOnSameCoords(meshes,ret1,ret2);
        PyObject *ret=PyTuple_New(3);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(ret2),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      static PyObject *MergeUMeshesOnSameCoords(PyObject *ms)
      {
        std::vector<const MEDCoupling::MEDCouplingUMesh *> meshes;
        convertFromPyObjVectorOfObj<const MEDCoupling::MEDCouplingUMesh *>(ms,SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh,"MEDCouplingUMesh",meshes);
        MEDCouplingUMesh *ret=MEDCouplingUMesh::MergeUMeshesOnSameCoords(meshes);
        return convertMesh(ret, SWIG_POINTER_OWN | 0 );
      }

      static PyObject *FuseUMeshesOnSameCoords(PyObject *ms, int compType)
      {
        std::size_t sz;
        std::vector<const MEDCouplingUMesh *> meshes;
        convertFromPyObjVectorOfObj<const MEDCoupling::MEDCouplingUMesh *>(ms,SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh,"MEDCouplingUMesh",meshes);
        std::vector<DataArrayIdType *> corr;
        MEDCouplingUMesh *um=MEDCouplingUMesh::FuseUMeshesOnSameCoords(meshes,compType,corr);
        sz=corr.size();
        PyObject *ret1=PyList_New(sz);
        for(std::size_t i=0;i<sz;i++)
          PyList_SetItem(ret1,i,SWIG_NewPointerObj(SWIG_as_voidptr(corr[i]),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyObject *ret=PyList_New(2);
        PyList_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(um),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyList_SetItem(ret,1,ret1);
        return ret;
      }

      static void PutUMeshesOnSameAggregatedCoords(PyObject *ms)
      {
        std::vector<MEDCouplingUMesh *> meshes;
        convertFromPyObjVectorOfObj<MEDCoupling::MEDCouplingUMesh *>(ms,SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh,"MEDCouplingUMesh",meshes);
        MEDCouplingUMesh::PutUMeshesOnSameAggregatedCoords(meshes);
      }

      static void MergeNodesOnUMeshesSharingSameCoords(PyObject *ms, double eps)
      {
        std::vector<MEDCouplingUMesh *> meshes;
        convertFromPyObjVectorOfObj<MEDCoupling::MEDCouplingUMesh *>(ms,SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh,"MEDCouplingUMesh",meshes);
        MEDCouplingUMesh::MergeNodesOnUMeshesSharingSameCoords(meshes,eps);
      }

      PyObject *are2DCellsNotCorrectlyOriented(PyObject *vec, bool polyOnly) const
      {
        double val;
        DataArrayDouble *a;
        DataArrayDoubleTuple *aa;
        std::vector<double> bb;
        mcIdType sw;
        int spaceDim=self->getSpaceDimension();
        const char msg[]="Python wrap of MEDCouplingUMesh::are2DCellsNotCorrectlyOriented : ";
        const double *v=convertObjToPossibleCpp5_Safe(vec,sw,val,a,aa,bb,msg,1,spaceDim,true);
        //
        std::vector<mcIdType> cells;
        self->are2DCellsNotCorrectlyOriented(v,polyOnly,cells);
        DataArrayIdType *ret=DataArrayIdType::New();
        ret->alloc(cells.size(),1);
        std::copy(cells.begin(),cells.end(),ret->getPointer());
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 );
      }

      void orientCorrectly2DCells(PyObject *vec, bool polyOnly)
      {
        double val;
        DataArrayDouble *a;
        DataArrayDoubleTuple *aa;
        std::vector<double> bb;
        mcIdType sw;
        int spaceDim=self->getSpaceDimension();
        const char msg[]="Python wrap of MEDCouplingUMesh::orientCorrectly2DCells : ";
        const double *v=convertObjToPossibleCpp5_Safe(vec,sw,val,a,aa,bb,msg,1,spaceDim,true);
        self->orientCorrectly2DCells(v,polyOnly);
      }

      PyObject *arePolyhedronsNotCorrectlyOriented() const
      {
        std::vector<mcIdType> cells;
        self->arePolyhedronsNotCorrectlyOriented(cells);
        DataArrayIdType *ret=DataArrayIdType::New();
        ret->alloc(cells.size(),1);
        std::copy(cells.begin(),cells.end(),ret->getPointer());
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 );
      }

      PyObject *getFastAveragePlaneOfThis() const
      {
        double vec[3];
        double pos[3];
        self->getFastAveragePlaneOfThis(vec,pos);
        double vals[6];
        std::copy(vec,vec+3,vals);
        std::copy(pos,pos+3,vals+3);
        return convertDblArrToPyListOfTuple<double>(vals,3,2);
      }

      static MEDCouplingUMesh *MergeUMeshes(PyObject *li)
      {
        std::vector<const MEDCoupling::MEDCouplingUMesh *> tmp;
        convertFromPyObjVectorOfObj<const MEDCoupling::MEDCouplingUMesh *>(li,SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh,"MEDCouplingUMesh",tmp);
        return MEDCouplingUMesh::MergeUMeshes(tmp);
      }

      PyObject *areCellsIncludedIn(const MEDCouplingUMesh *other, int compType) const
      {
        DataArrayIdType *ret1;
        bool ret0=self->areCellsIncludedIn(other,compType,ret1);
        PyObject *ret=PyTuple_New(2);
        PyObject *ret0Py=ret0?Py_True:Py_False;
        Py_XINCREF(ret0Py);
        PyTuple_SetItem(ret,0,ret0Py);
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *areCellsIncludedInPolicy7(const MEDCouplingUMesh *other) const
      {
        DataArrayIdType *ret1;
        bool ret0=self->areCellsIncludedInPolicy7(other,ret1);
        PyObject *ret=PyTuple_New(2);
        PyObject *ret0Py=ret0?Py_True:Py_False;
        Py_XINCREF(ret0Py);
        PyTuple_SetItem(ret,0,ret0Py);
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *explode3DMeshTo1D() const
      {
        MCAuto<DataArrayIdType> d0=DataArrayIdType::New();
        MCAuto<DataArrayIdType> d1=DataArrayIdType::New();
        MCAuto<DataArrayIdType> d2=DataArrayIdType::New();
        MCAuto<DataArrayIdType> d3=DataArrayIdType::New();
        MEDCouplingUMesh *m=self->explode3DMeshTo1D(d0,d1,d2,d3);
        PyObject *ret=PyTuple_New(5);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(m),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(d0.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(d1.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,3,SWIG_NewPointerObj(SWIG_as_voidptr(d2.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,4,SWIG_NewPointerObj(SWIG_as_voidptr(d3.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *explodeMeshTo(int targetDeltaLevel) const
      {
        MCAuto<DataArrayIdType> desc,descIndx,revDesc,revDescIndx;
        MCAuto<MEDCouplingUMesh> m=self->explodeMeshTo(targetDeltaLevel,desc,descIndx,revDesc,revDescIndx);
        PyObject *ret=PyTuple_New(5);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(m.retn()),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(desc.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(descIndx.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,3,SWIG_NewPointerObj(SWIG_as_voidptr(revDesc.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,4,SWIG_NewPointerObj(SWIG_as_voidptr(revDescIndx.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *explodeIntoEdges() const
      {
        MCAuto<DataArrayIdType> desc,descIndex,revDesc,revDescIndx;
        MCAuto<MEDCouplingUMesh> m(self->explodeIntoEdges(desc,descIndex,revDesc,revDescIndx));
        PyObject *ret=PyTuple_New(5);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(m.retn()),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(desc.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(descIndex.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,3,SWIG_NewPointerObj(SWIG_as_voidptr(revDesc.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,4,SWIG_NewPointerObj(SWIG_as_voidptr(revDescIndx.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *explodeMeshIntoMicroEdges() const
      {
        MCAuto<DataArrayIdType> d0=DataArrayIdType::New();
        MCAuto<DataArrayIdType> d1=DataArrayIdType::New();
        MCAuto<DataArrayIdType> d2=DataArrayIdType::New();
        MCAuto<DataArrayIdType> d3=DataArrayIdType::New();
        MEDCouplingUMesh *m=self->explodeMeshIntoMicroEdges(d0,d1,d2,d3);
        PyObject *ret=PyTuple_New(5);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(m),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(d0.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(d1.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,3,SWIG_NewPointerObj(SWIG_as_voidptr(d2.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,4,SWIG_NewPointerObj(SWIG_as_voidptr(d3.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *buildDescendingConnectivity() const
      {
        MCAuto<DataArrayIdType> d0=DataArrayIdType::New();
        MCAuto<DataArrayIdType> d1=DataArrayIdType::New();
        MCAuto<DataArrayIdType> d2=DataArrayIdType::New();
        MCAuto<DataArrayIdType> d3=DataArrayIdType::New();
        MEDCouplingUMesh *m=self->buildDescendingConnectivity(d0,d1,d2,d3);
        PyObject *ret=PyTuple_New(5);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(m),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(d0.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(d1.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,3,SWIG_NewPointerObj(SWIG_as_voidptr(d2.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,4,SWIG_NewPointerObj(SWIG_as_voidptr(d3.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *buildDescendingConnectivity2() const
      {
        MCAuto<DataArrayIdType> d0=DataArrayIdType::New();
        MCAuto<DataArrayIdType> d1=DataArrayIdType::New();
        MCAuto<DataArrayIdType> d2=DataArrayIdType::New();
        MCAuto<DataArrayIdType> d3=DataArrayIdType::New();
        MEDCouplingUMesh *m=self->buildDescendingConnectivity2(d0,d1,d2,d3);
        PyObject *ret=PyTuple_New(5);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(m),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(d0.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(d1.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,3,SWIG_NewPointerObj(SWIG_as_voidptr(d2.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,4,SWIG_NewPointerObj(SWIG_as_voidptr(d3.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *computeNeighborsOfCells() const
      {
        DataArrayIdType *neighbors=0,*neighborsIdx=0;
        self->computeNeighborsOfCells(neighbors,neighborsIdx);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(neighbors),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(neighborsIdx),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *computeNeighborsOfNodes() const
      {
        DataArrayIdType *neighbors=0,*neighborsIdx=0;
        self->computeNeighborsOfNodes(neighbors,neighborsIdx);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(neighbors),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(neighborsIdx),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *computeEnlargedNeighborsOfNodes() const
      {
        MCAuto<DataArrayIdType> neighbors,neighborsIdx;
        self->computeEnlargedNeighborsOfNodes(neighbors,neighborsIdx);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(neighbors.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(neighborsIdx.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *computeCellNeighborhoodFromNodesOne(const DataArrayIdType *nodeNeigh, const DataArrayIdType *nodeNeighI) const
      {
        MCAuto<DataArrayIdType> cellNeigh,cellNeighIndex;
        self->computeCellNeighborhoodFromNodesOne(nodeNeigh,nodeNeighI,cellNeigh,cellNeighIndex);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(cellNeigh.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(cellNeighIndex.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      static PyObject *ComputeNeighborsOfCellsAdv(const DataArrayIdType *desc, const DataArrayIdType *descI, const DataArrayIdType *revDesc, const DataArrayIdType *revDescI)
      {
        DataArrayIdType *neighbors=0,*neighborsIdx=0;
        MEDCouplingUMesh::ComputeNeighborsOfCellsAdv(desc,descI,revDesc,revDescI,neighbors,neighborsIdx);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(neighbors),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(neighborsIdx),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *emulateMEDMEMBDC(const MEDCouplingUMesh *nM1LevMesh)
      {
        MCAuto<DataArrayIdType> d0=DataArrayIdType::New();
        MCAuto<DataArrayIdType> d1=DataArrayIdType::New();
        DataArrayIdType *d2,*d3,*d4,*dd5;
        MEDCouplingUMesh *mOut=self->emulateMEDMEMBDC(nM1LevMesh,d0,d1,d2,d3,d4,dd5);
        PyObject *ret=PyTuple_New(7);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(mOut),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(d0.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(d1.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,3,SWIG_NewPointerObj(SWIG_as_voidptr(d2),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,4,SWIG_NewPointerObj(SWIG_as_voidptr(d3),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,5,SWIG_NewPointerObj(SWIG_as_voidptr(d4),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,6,SWIG_NewPointerObj(SWIG_as_voidptr(dd5),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      DataArrayDouble *getPartBarycenterAndOwner(DataArrayIdType *da) const
      {
        if(!da)
          throw INTERP_KERNEL::Exception("Not null DataArrayIdType instance expected !");
        da->checkAllocated();
        return self->getPartBarycenterAndOwner(da->getConstPointer(),da->getConstPointer()+da->getNbOfElems());
      }

      DataArrayDouble *getPartMeasureField(bool isAbs, DataArrayIdType *da) const
      {
        if(!da)
          throw INTERP_KERNEL::Exception("Not null DataArrayIdType instance expected !");
        da->checkAllocated();
        return self->getPartMeasureField(isAbs,da->getConstPointer(),da->getConstPointer()+da->getNbOfElems());
      }

      MEDCouplingFieldDouble *buildPartOrthogonalField(DataArrayIdType *da) const
      {
        if(!da)
          throw INTERP_KERNEL::Exception("Not null DataArrayIdType instance expected !");
        da->checkAllocated();
        return self->buildPartOrthogonalField(da->getConstPointer(),da->getConstPointer()+da->getNbOfElems());
      }

      PyObject *getTypesOfPart(DataArrayIdType *da) const
      {
        if(!da)
          throw INTERP_KERNEL::Exception("Not null DataArrayIdType instance expected !");
        da->checkAllocated();
        std::set<INTERP_KERNEL::NormalizedCellType> result=self->getTypesOfPart(da->getConstPointer(),da->getConstPointer()+da->getNbOfElems());
        std::set<INTERP_KERNEL::NormalizedCellType>::const_iterator iL=result.begin();
        PyObject *res = PyList_New(result.size());
        for (int i=0;iL!=result.end(); i++, iL++)
          PyList_SetItem(res,i,PyInt_FromLong(*iL));
        return res;
      }

      DataArrayIdType *keepCellIdsByType(INTERP_KERNEL::NormalizedCellType type, DataArrayIdType *da) const
      {
        if(!da)
          throw INTERP_KERNEL::Exception("Not null DataArrayIdType instance expected !");
        da->checkAllocated();
        DataArrayIdType *ret=self->keepCellIdsByType(type,da->getConstPointer(),da->getConstPointer()+da->getNbOfElems());
        ret->setName(da->getName().c_str());
        return ret;
      }

      static PyObject *Intersect2DMeshes(const MEDCouplingUMesh *m1, const MEDCouplingUMesh *m2, double eps)
      {
        DataArrayIdType *cellNb1=0,*cellNb2=0;
        MEDCouplingUMesh *mret=MEDCouplingUMesh::Intersect2DMeshes(m1,m2,eps,cellNb1,cellNb2);
        PyObject *ret=PyTuple_New(3);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(mret),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(cellNb1),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(cellNb2),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      static PyObject *Intersect2DMeshWith1DLine(const MEDCouplingUMesh *mesh2D, const MEDCouplingUMesh *mesh1D, double eps)
      {
        MEDCouplingUMesh *splitMesh2D(0),*splitMesh1D(0);
        DataArrayIdType *cellIdInMesh2D(0),*cellIdInMesh1D(0);
        MEDCouplingUMesh::Intersect2DMeshWith1DLine(mesh2D,mesh1D,eps,splitMesh2D,splitMesh1D,cellIdInMesh2D,cellIdInMesh1D);
        PyObject *ret(PyTuple_New(4));
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(splitMesh2D),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(splitMesh1D),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(cellIdInMesh2D),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,3,SWIG_NewPointerObj(SWIG_as_voidptr(cellIdInMesh1D),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *buildSlice3D(PyObject *origin, PyObject *vec, double eps) const
      {
        int spaceDim=self->getSpaceDimension();
        if(spaceDim!=3)
          throw INTERP_KERNEL::Exception("Python wrap of MEDCouplingUMesh::buildSlice3D : works only for spaceDim 3 !");
        double val,val2;
        DataArrayDouble *a,*a2;
        DataArrayDoubleTuple *aa,*aa2;
        std::vector<double> bb,bb2;
        mcIdType sw;
        const char msg[]="Python wrap of MEDCouplingUMesh::buildSlice3D : 1st parameter for origin.";
        const char msg2[]="Python wrap of MEDCouplingUMesh::buildSlice3D : 2nd parameter for vector.";
        const double *orig=convertObjToPossibleCpp5_Safe(origin,sw,val,a,aa,bb,msg,1,spaceDim,true);
        const double *vect=convertObjToPossibleCpp5_Safe(vec,sw,val2,a2,aa2,bb2,msg2,1,spaceDim,true);
        //
        DataArrayIdType *cellIds=0;
        MEDCouplingUMesh *ret0=self->buildSlice3D(orig,vect,eps,cellIds);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(cellIds),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *buildSlice3DSurf(PyObject *origin, PyObject *vec, double eps) const
      {
        int spaceDim=self->getSpaceDimension();
        if(spaceDim!=3)
          throw INTERP_KERNEL::Exception("Python wrap of MEDCouplingUMesh::buildSlice3DSurf : works only for spaceDim 3 !");
        double val,val2;
        DataArrayDouble *a,*a2;
        DataArrayDoubleTuple *aa,*aa2;
        std::vector<double> bb,bb2;
        mcIdType sw;
        const char msg[]="Python wrap of MEDCouplingUMesh::buildSlice3DSurf : 1st parameter for origin.";
        const char msg2[]="Python wrap of MEDCouplingUMesh::buildSlice3DSurf : 2nd parameter for vector.";
        const double *orig=convertObjToPossibleCpp5_Safe(origin,sw,val,a,aa,bb,msg,1,spaceDim,true);
        const double *vect=convertObjToPossibleCpp5_Safe(vec,sw,val2,a2,aa2,bb2,msg2,1,spaceDim,true);
        //
        DataArrayIdType *cellIds=0;
        MEDCouplingUMesh *ret0=self->buildSlice3DSurf(orig,vect,eps,cellIds);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(cellIds),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      MEDCouplingUMesh *clipSingle3DCellByPlane(PyObject *origin, PyObject *vec, double eps) const
      {
        double val,val2;
        DataArrayDouble *a,*a2;
        DataArrayDoubleTuple *aa,*aa2;
        std::vector<double> bb,bb2;
        mcIdType sw;
        const char msg[]="Python wrap of MEDCouplingUMesh::clipSingle3DCellByPlane : 1st parameter for origin.";
        const char msg2[]="Python wrap of MEDCouplingUMesh::clipSingle3DCellByPlane : 2nd parameter for vector.";
        const double *orig=convertObjToPossibleCpp5_Safe(origin,sw,val,a,aa,bb,msg,1,3,true);
        const double *vect=convertObjToPossibleCpp5_Safe(vec,sw,val2,a2,aa2,bb2,msg2,1,3,true);
        MCAuto<MEDCouplingUMesh> ret(self->clipSingle3DCellByPlane(orig,vect,eps));
        return ret.retn();
      }

      DataArrayIdType *getCellIdsCrossingPlane(PyObject *origin, PyObject *vec, double eps) const
      {
        int spaceDim=self->getSpaceDimension();
        if(spaceDim!=3)
          throw INTERP_KERNEL::Exception("Python wrap of MEDCouplingUMesh::getCellIdsCrossingPlane : works only for spaceDim 3 !");
        double val,val2;
        DataArrayDouble *a,*a2;
        DataArrayDoubleTuple *aa,*aa2;
        std::vector<double> bb,bb2;
        mcIdType sw;
        const char msg[]="Python wrap of MEDCouplingUMesh::getCellIdsCrossingPlane : 1st parameter for origin.";
        const char msg2[]="Python wrap of MEDCouplingUMesh::getCellIdsCrossingPlane : 2nd parameter for vector.";
        const double *orig=convertObjToPossibleCpp5_Safe(origin,sw,val,a,aa,bb,msg,1,spaceDim,true);
        const double *vect=convertObjToPossibleCpp5_Safe(vec,sw,val2,a2,aa2,bb2,msg2,1,spaceDim,true);
        return self->getCellIdsCrossingPlane(orig,vect,eps);
      }

      void convertToPolyTypes(PyObject *li)
      {
        mcIdType sw;
        mcIdType pos1;
        std::vector<mcIdType> pos2;
        DataArrayIdType *pos3=0;
        DataArrayIdTypeTuple *pos4=0;
        convertIntStarLikePyObjToCpp(li,sw,pos1,pos2,pos3,pos4);
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
    void convertExtrudedPolyhedra();
    bool unPolyze();
    void simplifyPolyhedra(double eps);
    void colinearizeEdges(double eps);
    MEDCouplingUMesh *buildSpreadZonesWithPoly() const;
    MEDCouplingUMesh *buildExtrudedMesh(const MEDCouplingUMesh *mesh1D, int policy);
  };

  //== MEDCouplingUMesh End

  //== MEDCouplingMappedExtrudedMesh

  class MEDCouplingMappedExtrudedMesh : public MEDCoupling::MEDCouplingMesh
  {
  public:
    static MEDCouplingMappedExtrudedMesh *New(const MEDCouplingUMesh *mesh3D, const MEDCouplingUMesh *mesh2D, int cell2DId);
    static MEDCouplingMappedExtrudedMesh *New(const MEDCouplingCMesh *mesh3D);
    MEDCouplingUMesh *build3DUnstructuredMesh() const;
    int get2DCellIdForExtrusion() const;
    %extend {
      MEDCouplingMappedExtrudedMesh(const MEDCouplingUMesh *mesh3D, const MEDCouplingUMesh *mesh2D, mcIdType cell2DId)
      {
        return MEDCouplingMappedExtrudedMesh::New(mesh3D,mesh2D,cell2DId);
      }

      MEDCouplingMappedExtrudedMesh(const MEDCouplingCMesh *mesh3D)
      {
        return MEDCouplingMappedExtrudedMesh::New(mesh3D);
      }

      MEDCouplingMappedExtrudedMesh()
      {
        return MEDCouplingMappedExtrudedMesh::New();
      }

      std::string __str__() const
      {
        return self->simpleRepr();
      }

      std::string __repr__() const
      {
        std::ostringstream oss;
        self->reprQuickOverview(oss);
        return oss.str();
      }

      PyObject *getMesh2D() const
      {
        MEDCouplingUMesh *ret=self->getMesh2D();
        if(ret)
          ret->incrRef();
        return convertMesh(ret, SWIG_POINTER_OWN | 0 );
      }
      PyObject *getMesh1D() const
      {
        MEDCouplingUMesh *ret=self->getMesh1D();
        if(ret)
          ret->incrRef();
        return convertMesh(ret, SWIG_POINTER_OWN | 0 );
      }
      PyObject *getMesh3DIds() const
      {
        DataArrayIdType *ret=self->getMesh3DIds();
        if(ret)
          ret->incrRef();
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 );
      }
    }
  };

  //== MEDCouplingMappedExtrudedMesh End

  class MEDCoupling1GTUMesh : public MEDCoupling::MEDCouplingPointSet
  {
  public:
    static MEDCoupling1GTUMesh *New(const std::string& name, INTERP_KERNEL::NormalizedCellType type);
    static MEDCoupling1GTUMesh *New(const MEDCouplingUMesh *m);
    INTERP_KERNEL::NormalizedCellType getCellModelEnum() const;
    int getNodalConnectivityLength() const;
    virtual void allocateCells(int nbOfCells=0);
    virtual void checkConsistencyOfConnectivity() const;
    %extend
    {
      virtual void insertNextCell(PyObject *li)
      {
        mcIdType szArr,sw,iTypppArr;
        std::vector<mcIdType> stdvecTyyppArr;
        const mcIdType *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
        self->insertNextCell(tmp,tmp+szArr);
      }

      virtual DataArrayIdType *getNodalConnectivity() const
      {
        DataArrayIdType *ret=self->getNodalConnectivity();
        if(ret) ret->incrRef();
        return ret;
      }

      static MEDCouplingUMesh *AggregateOnSameCoordsToUMesh(PyObject *li)
      {
        std::vector< const MEDCoupling1GTUMesh *> parts;
        convertFromPyObjVectorOfObj<const MEDCoupling::MEDCoupling1GTUMesh *>(li,SWIGTYPE_p_MEDCoupling__MEDCoupling1GTUMesh,"MEDCoupling1GTUMesh",parts);
        return MEDCoupling1GTUMesh::AggregateOnSameCoordsToUMesh(parts);
      }
    }
  };

  //== MEDCoupling1SGTUMesh

  class MEDCoupling1SGTUMesh : public MEDCoupling::MEDCoupling1GTUMesh
  {
  public:
    static MEDCoupling1SGTUMesh *New(const std::string& name, INTERP_KERNEL::NormalizedCellType type);
    static MEDCoupling1SGTUMesh *New(const MEDCouplingUMesh *m);
    void setNodalConnectivity(DataArrayIdType *nodalConn);
    int getNumberOfNodesPerCell() const;
    static MEDCoupling1SGTUMesh *Merge1SGTUMeshes(const MEDCoupling1SGTUMesh *mesh1, const MEDCoupling1SGTUMesh *mesh2);
    MEDCoupling1SGTUMesh *buildSetInstanceFromThis(int spaceDim) const;
    MEDCoupling1GTUMesh *computeDualMesh() const;
    MEDCoupling1SGTUMesh *explodeEachHexa8To6Quad4() const;
    DataArrayIdType *sortHexa8EachOther();
    %extend
    {
      MEDCoupling1SGTUMesh()
      {
        return MEDCoupling1SGTUMesh::New();
      }

      MEDCoupling1SGTUMesh(const std::string& name, INTERP_KERNEL::NormalizedCellType type)
      {
        return MEDCoupling1SGTUMesh::New(name,type);
      }

      MEDCoupling1SGTUMesh(const MEDCouplingUMesh *m)
      {
        return MEDCoupling1SGTUMesh::New(m);
      }

      std::string __str__() const
      {
        return self->simpleRepr();
      }

      std::string __repr__() const
      {
        std::ostringstream oss;
        self->reprQuickOverview(oss);
        return oss.str();
      }

      PyObject *structurizeMe(double eps=1e-12) const
      {
        DataArrayIdType *cellPerm(0),*nodePerm(0);
        MEDCouplingCMesh *retCpp(self->structurizeMe(cellPerm,nodePerm,eps));
        PyObject *ret(PyTuple_New(3));
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(retCpp),SWIGTYPE_p_MEDCoupling__MEDCouplingCMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(cellPerm),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(nodePerm),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      DataArrayDouble *MEDCoupling1SGTUMesh::computeTriangleHeight() const
      {
        MCAuto<DataArrayDouble> ret = self->computeTriangleHeight();
        return ret.retn();
      }

      static MEDCoupling1SGTUMesh *Merge1SGTUMeshes(PyObject *li)
      {
        std::vector<const MEDCoupling::MEDCoupling1SGTUMesh *> tmp;
        convertFromPyObjVectorOfObj<const MEDCoupling::MEDCoupling1SGTUMesh *>(li,SWIGTYPE_p_MEDCoupling__MEDCoupling1SGTUMesh,"MEDCoupling1SGTUMesh",tmp);
        return MEDCoupling1SGTUMesh::Merge1SGTUMeshes(tmp);
      }

      static MEDCoupling1SGTUMesh *Merge1SGTUMeshesOnSameCoords(PyObject *li)
      {
        std::vector<const MEDCoupling::MEDCoupling1SGTUMesh *> tmp;
        convertFromPyObjVectorOfObj<const MEDCoupling::MEDCoupling1SGTUMesh *>(li,SWIGTYPE_p_MEDCoupling__MEDCoupling1SGTUMesh,"MEDCoupling1SGTUMesh",tmp);
        return MEDCoupling1SGTUMesh::Merge1SGTUMeshesOnSameCoords(tmp);
      }
    }
  };

  //== MEDCoupling1SGTUMesh End

  //== MEDCoupling1DGTUMesh

  class MEDCoupling1DGTUMesh : public MEDCoupling::MEDCoupling1GTUMesh
  {
  public:
    static MEDCoupling1DGTUMesh *New(const std::string& name, INTERP_KERNEL::NormalizedCellType type);
    static MEDCoupling1DGTUMesh *New(const MEDCouplingUMesh *m);
    void setNodalConnectivity(DataArrayIdType *nodalConn, DataArrayIdType *nodalConnIndex);
    MEDCoupling1DGTUMesh *buildSetInstanceFromThis(int spaceDim) const;
    bool isPacked() const;
    %extend
    {
      MEDCoupling1DGTUMesh()
      {
        return MEDCoupling1DGTUMesh::New();
      }
      MEDCoupling1DGTUMesh(const std::string& name, INTERP_KERNEL::NormalizedCellType type)
      {
        return MEDCoupling1DGTUMesh::New(name,type);
      }

      MEDCoupling1DGTUMesh(const MEDCouplingUMesh *m)
      {
        return MEDCoupling1DGTUMesh::New(m);
      }

      std::string __str__() const
      {
        return self->simpleRepr();
      }

      std::string __repr__() const
      {
        std::ostringstream oss;
        self->reprQuickOverview(oss);
        return oss.str();
      }

      DataArrayIdType *getNodalConnectivityIndex() const
      {
        DataArrayIdType *ret=self->getNodalConnectivityIndex();
        if(ret) ret->incrRef();
        return ret;
      }

      PyObject *retrievePackedNodalConnectivity() const
      {
        DataArrayIdType *ret1=0,*ret2=0;
        bool ret0=self->retrievePackedNodalConnectivity(ret1,ret2);
        PyObject *ret0Py=ret0?Py_True:Py_False;
        Py_XINCREF(ret0Py);
        PyObject *ret=PyTuple_New(3);
        PyTuple_SetItem(ret,0,ret0Py);
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(ret2),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *copyWithNodalConnectivityPacked() const
      {
        bool ret1;
        MEDCoupling1DGTUMesh *ret0=self->copyWithNodalConnectivityPacked(ret1);
        PyObject *ret=PyTuple_New(2);
        PyObject *ret1Py=ret1?Py_True:Py_False; Py_XINCREF(ret1Py);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_MEDCoupling__MEDCoupling1DGTUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,ret1Py);
        return ret;
      }

      static MEDCoupling1DGTUMesh *Merge1DGTUMeshes(PyObject *li)
      {
        std::vector<const MEDCoupling::MEDCoupling1DGTUMesh *> tmp;
        convertFromPyObjVectorOfObj<const MEDCoupling::MEDCoupling1DGTUMesh *>(li,SWIGTYPE_p_MEDCoupling__MEDCoupling1DGTUMesh,"MEDCoupling1DGTUMesh",tmp);
        return MEDCoupling1DGTUMesh::Merge1DGTUMeshes(tmp);
      }

      static MEDCoupling1DGTUMesh *Merge1DGTUMeshesOnSameCoords(PyObject *li)
      {
        std::vector<const MEDCoupling::MEDCoupling1DGTUMesh *> tmp;
        convertFromPyObjVectorOfObj<const MEDCoupling::MEDCoupling1DGTUMesh *>(li,SWIGTYPE_p_MEDCoupling__MEDCoupling1DGTUMesh,"MEDCoupling1DGTUMesh",tmp);
        return MEDCoupling1DGTUMesh::Merge1DGTUMeshesOnSameCoords(tmp);
      }

      static DataArrayIdType *AggregateNodalConnAndShiftNodeIds(PyObject *li, const std::vector<mcIdType>& offsetInNodeIdsPerElt)
      {
        std::vector<const MEDCoupling::DataArrayIdType *> tmp;
        convertFromPyObjVectorOfObj<const MEDCoupling::DataArrayIdType *>(li,SWIGTITraits<mcIdType>::TI,"DataArrayIdType",tmp);
        return MEDCoupling1DGTUMesh::AggregateNodalConnAndShiftNodeIds(tmp,offsetInNodeIdsPerElt);
      }
    }
  };

  //== MEDCoupling1DGTUMeshEnd

  class MEDCouplingStructuredMesh : public MEDCoupling::MEDCouplingMesh
  {
  public:
    mcIdType getCellIdFromPos(mcIdType i, mcIdType j, mcIdType k) const;
    mcIdType getNodeIdFromPos(mcIdType i, mcIdType j, mcIdType k) const;
    mcIdType getNumberOfCellsOfSubLevelMesh() const;
    int getSpaceDimensionOnNodeStruct() const;
    double computeSquareness() const;
    virtual std::vector<mcIdType> getNodeGridStructure() const;
    std::vector<mcIdType> getCellGridStructure() const;
    MEDCoupling1SGTUMesh *build1SGTUnstructured() const;
    std::vector<mcIdType> getLocationFromCellId(mcIdType cellId) const;
    std::vector<mcIdType> getLocationFromNodeId(mcIdType cellId) const;
    static INTERP_KERNEL::NormalizedCellType GetGeoTypeGivenMeshDimension(int meshDim);
    MEDCoupling1SGTUMesh *build1SGTSubLevelMesh() const;
    static mcIdType DeduceNumberOfGivenStructure(const std::vector<mcIdType>& st);
    static DataArrayIdType *ComputeCornersGhost(const std::vector<mcIdType>& st, mcIdType ghostLev);
    static std::vector<mcIdType> GetSplitVectFromStruct(const std::vector<mcIdType>& strct);
    %extend
    {
      virtual MEDCouplingStructuredMesh *buildStructuredSubPart(PyObject *cellPart) const
      {
        mcIdType tmpp1=-1,tmpp2=-1;
        std::vector<mcIdType> tmp=fillArrayWithPyListInt2(cellPart,tmpp1,tmpp2);
        std::vector< std::pair<mcIdType,mcIdType> > inp;
        if(tmpp2==2)
          {
            inp.resize(tmpp1);
            for(mcIdType i=0;i<tmpp1;i++)
              { inp[i].first=tmp[2*i]; inp[i].second=tmp[2*i+1]; }
          }
        else if(tmpp2==1)
          {
            if(tmpp1%2!=0)
              throw INTERP_KERNEL::Exception("Wrap of MEDCouplingStructuredMesh.buildStructuredSubPart : invalid input size ! Must be even size !");
            inp.resize(tmpp1/2);
            for(mcIdType i=0;i<tmpp1/2;i++)
              { inp[i].first=tmp[2*i]; inp[i].second=tmp[2*i+1]; }
          }
        else
          throw INTERP_KERNEL::Exception("Wrap of MEDCouplingStructuredMesh.buildStructuredSubPart : invalid input size !");
        return self->buildStructuredSubPart(inp);
      }

      static DataArrayIdType *BuildExplicitIdsFrom(PyObject *st, PyObject *part)
      {
        std::vector< std::pair<mcIdType,mcIdType> > inp;
        convertPyToVectorPairInt(part,inp);
        //
        mcIdType szArr,sw,iTypppArr;
        std::vector<mcIdType> stdvecTyyppArr;
        const mcIdType *tmp4=convertIntStarLikePyObjToCppIntStar(st,sw,szArr,iTypppArr,stdvecTyyppArr);
        std::vector<mcIdType> tmp5(tmp4,tmp4+szArr);
        //
        return MEDCouplingStructuredMesh::BuildExplicitIdsFrom(tmp5,inp);
      }

      static void MultiplyPartOf(const std::vector<mcIdType>& st, PyObject *part, double factor, DataArrayDouble *da)
      {
        std::vector< std::pair<mcIdType,mcIdType> > inp;
        convertPyToVectorPairInt(part,inp);
        MEDCouplingStructuredMesh::MultiplyPartOf(st,inp,factor,da);
      }

      static void MultiplyPartOfByGhost(const std::vector<mcIdType>& st, PyObject *part, mcIdType ghostSize, double factor, DataArrayDouble *da)
      {
        std::vector< std::pair<mcIdType,mcIdType> > inp;
        convertPyToVectorPairInt(part,inp);
        MEDCouplingStructuredMesh::MultiplyPartOfByGhost(st,inp,ghostSize,factor,da);
      }

      static PyObject *PutInGhostFormat(mcIdType ghostSize, const std::vector<mcIdType>& st, PyObject *part)
      {
        std::vector< std::pair<mcIdType,mcIdType> > inp;
        convertPyToVectorPairInt(part,inp);
        std::vector<mcIdType> stWithGhost;
        std::vector< std::pair<mcIdType,mcIdType> > partWithGhost;
        MEDCouplingStructuredMesh::PutInGhostFormat(ghostSize,st,inp,stWithGhost,partWithGhost);
        PyObject *ret(PyTuple_New(2));
        PyTuple_SetItem(ret,0,convertIntArrToPyList2(stWithGhost));
        PyTuple_SetItem(ret,1,convertFromVectorPairInt(partWithGhost));
        return ret;
      }

      static DataArrayDouble *ExtractFieldOfDoubleFrom(const std::vector<mcIdType>& st, const DataArrayDouble *fieldOfDbl, PyObject *partCompactFormat)
      {
        std::vector< std::pair<mcIdType,mcIdType> > inp;
        convertPyToVectorPairInt(partCompactFormat,inp);
        return MEDCouplingStructuredMesh::ExtractFieldOfDoubleFrom(st,fieldOfDbl,inp);
      }

      static void AssignPartOfFieldOfDoubleUsing(const std::vector<mcIdType>& st, DataArrayDouble *fieldOfDbl, PyObject *partCompactFormat, const DataArrayDouble *other)
      {
        std::vector< std::pair<mcIdType,mcIdType> > inp;
        convertPyToVectorPairInt(partCompactFormat,inp);
        MEDCouplingStructuredMesh::AssignPartOfFieldOfDoubleUsing(st,fieldOfDbl,inp,other);
      }

      static mcIdType DeduceNumberOfGivenRangeInCompactFrmt(PyObject *part)
      {
        std::vector< std::pair<mcIdType,mcIdType> > inp;
        convertPyToVectorPairInt(part,inp);
        return MEDCouplingStructuredMesh::DeduceNumberOfGivenRangeInCompactFrmt(inp);
      }

      static DataArrayIdType *Build1GTNodalConnectivity(PyObject *li)
      {
        mcIdType szArr,sw,iTypppArr;
        std::vector<mcIdType> stdvecTyyppArr;
        const mcIdType *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
        return MEDCouplingStructuredMesh::Build1GTNodalConnectivity(tmp,tmp+szArr);
      }

      static DataArrayIdType *Build1GTNodalConnectivityOfSubLevelMesh(PyObject *li)
      {
        mcIdType szArr,sw,iTypppArr;
        std::vector<mcIdType> stdvecTyyppArr;
        const mcIdType *tmp(convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr));
        return MEDCouplingStructuredMesh::Build1GTNodalConnectivityOfSubLevelMesh(tmp,tmp+szArr);
      }

      static std::vector<mcIdType> GetDimensionsFromCompactFrmt(PyObject *partCompactFormat)
      {
        std::vector< std::pair<mcIdType,mcIdType> > inp;
        convertPyToVectorPairInt(partCompactFormat,inp);
        return MEDCouplingStructuredMesh::GetDimensionsFromCompactFrmt(inp);
      }

      static PyObject *GetCompactFrmtFromDimensions(const std::vector<mcIdType>& dims)
      {
        std::vector< std::pair<mcIdType,mcIdType> > ret(MEDCouplingStructuredMesh::GetCompactFrmtFromDimensions(dims));
        PyObject *retPy=PyList_New(ret.size());
        for(std::size_t i=0;i<ret.size();i++)
          {
            PyObject *tmp=PyTuple_New(2);
            PyTuple_SetItem(tmp,0,PyInt_FromLong(ret[i].first));
            PyTuple_SetItem(tmp,1,PyInt_FromLong(ret[i].second));
            PyList_SetItem(retPy,i,tmp);
          }
        return retPy;
      }

      static PyObject *IntersectRanges(PyObject *r1, PyObject *r2)
      {
        std::vector< std::pair<mcIdType,mcIdType> > r1Cpp,r2Cpp;
        convertPyToVectorPairInt(r1,r1Cpp);
        convertPyToVectorPairInt(r2,r2Cpp);
        std::vector< std::pair<mcIdType,mcIdType> > ret(MEDCouplingStructuredMesh::IntersectRanges(r1Cpp,r2Cpp));
        PyObject *retPy=PyList_New(ret.size());
        for(std::size_t i=0;i<ret.size();i++)
          {
            PyObject *tmp=PyTuple_New(2);
            PyTuple_SetItem(tmp,0,PyInt_FromLong(ret[i].first));
            PyTuple_SetItem(tmp,1,PyInt_FromLong(ret[i].second));
            PyList_SetItem(retPy,i,tmp);
          }
        return retPy;
      }

      static bool AreRangesIntersect(PyObject *r1, PyObject *r2)
      {
        std::vector< std::pair<mcIdType,mcIdType> > r1Cpp,r2Cpp;
        convertPyToVectorPairInt(r1,r1Cpp);
        convertPyToVectorPairInt(r2,r2Cpp);
        return MEDCouplingStructuredMesh::AreRangesIntersect(r1Cpp,r2Cpp);
      }

      static PyObject *IsPartStructured(PyObject *li, PyObject *st)
      {
        mcIdType szArr,sw,iTypppArr;
        std::vector<mcIdType> stdvecTyyppArr;
        const mcIdType *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
        mcIdType szArr2,sw2,iTypppArr2;
        std::vector<mcIdType> stdvecTyyppArr2;
        const mcIdType *tmp2=convertIntStarLikePyObjToCppIntStar(st,sw2,szArr2,iTypppArr2,stdvecTyyppArr2);
        std::vector<mcIdType> tmp3(tmp2,tmp2+szArr2);
        std::vector< std::pair<mcIdType,mcIdType> > partCompactFormat;
        bool ret0=MEDCouplingStructuredMesh::IsPartStructured(tmp,tmp+szArr,tmp3,partCompactFormat);
        PyObject *ret=PyTuple_New(2);
        PyObject *ret0Py=ret0?Py_True:Py_False; Py_XINCREF(ret0Py);
        PyTuple_SetItem(ret,0,ret0Py);
        PyObject *ret1Py=PyList_New(partCompactFormat.size());
        for(std::size_t i=0;i<partCompactFormat.size();i++)
          {
            PyObject *tmp4=PyTuple_New(2);
            PyTuple_SetItem(tmp4,0,PyInt_FromLong(partCompactFormat[i].first));
            PyTuple_SetItem(tmp4,1,PyInt_FromLong(partCompactFormat[i].second));
            PyList_SetItem(ret1Py,i,tmp4);
          }
        PyTuple_SetItem(ret,1,ret1Py);
        return ret;
      }

      static PyObject *ChangeReferenceFromGlobalOfCompactFrmt(PyObject *bigInAbs, PyObject *partOfBigInAbs, bool check=true)
      {
        std::vector< std::pair<mcIdType,mcIdType> > param0,param1,ret;
        convertPyToVectorPairInt(bigInAbs,param0);
        convertPyToVectorPairInt(partOfBigInAbs,param1);
        MEDCouplingStructuredMesh::ChangeReferenceFromGlobalOfCompactFrmt(param0,param1,ret,check);
        PyObject *retPy(PyList_New(ret.size()));
        for(std::size_t i=0;i<ret.size();i++)
          {
            PyObject *tmp(PyTuple_New(2));
            PyTuple_SetItem(tmp,0,PyInt_FromLong(ret[i].first));
            PyTuple_SetItem(tmp,1,PyInt_FromLong(ret[i].second));
            PyList_SetItem(retPy,i,tmp);
          }
        return retPy;
      }

      static PyObject *TranslateCompactFrmt(PyObject *part, const std::vector<mcIdType>& translation)
      {
        std::vector< std::pair<mcIdType,mcIdType> > param0;
        convertPyToVectorPairInt(part,param0);
        std::vector< std::pair<mcIdType,mcIdType> > ret(MEDCouplingStructuredMesh::TranslateCompactFrmt(param0,translation));
        PyObject *retPy(PyList_New(ret.size()));
        for(std::size_t i=0;i<ret.size();i++)
          {
            PyObject *tmp(PyTuple_New(2));
            PyTuple_SetItem(tmp,0,PyInt_FromLong(ret[i].first));
            PyTuple_SetItem(tmp,1,PyInt_FromLong(ret[i].second));
            PyList_SetItem(retPy,i,tmp);
          }
        return retPy;
      }

      static std::vector<mcIdType> FindTranslationFrom(PyObject *startingFrom, PyObject *goingTo)
      {
        std::vector< std::pair<mcIdType,mcIdType> > param0,param1;
        convertPyToVectorPairInt(startingFrom,param0);
        convertPyToVectorPairInt(goingTo,param1);
        return  MEDCouplingStructuredMesh::FindTranslationFrom(param0,param1);
      }

      static PyObject *ChangeReferenceToGlobalOfCompactFrmt(PyObject *bigInAbs, PyObject *partOfBigRelativeToBig, bool check=true)
      {
        std::vector< std::pair<mcIdType,mcIdType> > param0,param1,ret;
        convertPyToVectorPairInt(bigInAbs,param0);
        convertPyToVectorPairInt(partOfBigRelativeToBig,param1);
        MEDCouplingStructuredMesh::ChangeReferenceToGlobalOfCompactFrmt(param0,param1,ret,check);
        PyObject *retPy(PyList_New(ret.size()));
        for(std::size_t i=0;i<ret.size();i++)
          {
            PyObject *tmp(PyTuple_New(2));
            PyTuple_SetItem(tmp,0,PyInt_FromLong(ret[i].first));
            PyTuple_SetItem(tmp,1,PyInt_FromLong(ret[i].second));
            PyList_SetItem(retPy,i,tmp);
          }
        return retPy;
      }
    }
  };

  class MEDCouplingCurveLinearMesh;

  //== MEDCouplingCMesh

  class MEDCouplingCMesh : public MEDCoupling::MEDCouplingStructuredMesh
  {
  public:
    static MEDCouplingCMesh *New();
    static MEDCouplingCMesh *New(const std::string& meshName);
    void setCoords(const DataArrayDouble *coordsX,
                   const DataArrayDouble *coordsY=0,
                   const DataArrayDouble *coordsZ=0);
    void setCoordsAt(int i, const DataArrayDouble *arr);
    MEDCouplingCurveLinearMesh *buildCurveLinear() const;
    %extend {
      MEDCouplingCMesh()
      {
        return MEDCouplingCMesh::New();
      }
      MEDCouplingCMesh(const std::string& meshName)
      {
        return MEDCouplingCMesh::New(meshName);
      }
      std::string __str__() const
      {
        return self->simpleRepr();
      }
      std::string __repr__() const
      {
        std::ostringstream oss;
        self->reprQuickOverview(oss);
        return oss.str();
      }
      DataArrayDouble *getCoordsAt(int i)
      {
        DataArrayDouble *ret=self->getCoordsAt(i);
        if(ret)
          ret->incrRef();
        return ret;
      }
    }
  };

  //== MEDCouplingCMesh End

  //== MEDCouplingCurveLinearMesh

  class MEDCouplingCurveLinearMesh : public MEDCoupling::MEDCouplingStructuredMesh
  {
  public:
    static MEDCouplingCurveLinearMesh *New();
    static MEDCouplingCurveLinearMesh *New(const std::string& meshName);
    void setCoords(const DataArrayDouble *coords);
    %extend {
      MEDCouplingCurveLinearMesh()
      {
        return MEDCouplingCurveLinearMesh::New();
      }
      MEDCouplingCurveLinearMesh(const std::string& meshName)
      {
        return MEDCouplingCurveLinearMesh::New(meshName);
      }
      std::string __str__() const
      {
        return self->simpleRepr();
      }
      std::string __repr__() const
      {
        std::ostringstream oss;
        self->reprQuickOverview(oss);
        return oss.str();
      }
      DataArrayDouble *getCoords()
      {
        DataArrayDouble *ret=self->getCoords();
        if(ret)
          ret->incrRef();
        return ret;
      }
      void setNodeGridStructure(PyObject *gridStruct)
      {
        mcIdType szArr,sw,iTypppArr;
        std::vector<mcIdType> stdvecTyyppArr;
        const mcIdType *tmp=convertIntStarLikePyObjToCppIntStar(gridStruct,sw,szArr,iTypppArr,stdvecTyyppArr);
        self->setNodeGridStructure(tmp,tmp+szArr);
      }
    }
  };

  //== MEDCouplingCurveLinearMesh End

  //== MEDCouplingIMesh

  class MEDCouplingIMesh : public MEDCoupling::MEDCouplingStructuredMesh
  {
  public:
    static MEDCouplingIMesh *New();
    //
    void setSpaceDimension(int spaceDim);
    std::vector<mcIdType> getNodeStruct() const;
    std::vector<double> getOrigin() const;
    std::vector<double> getDXYZ() const;
    void setAxisUnit(const std::string& unitName);
    std::string getAxisUnit() const;
    double getMeasureOfAnyCell() const;
    MEDCouplingCMesh *convertToCartesian() const;
    void refineWithFactor(const std::vector<mcIdType>& factors);
    MEDCouplingIMesh *asSingleCell() const;
    MEDCouplingIMesh *buildWithGhost(mcIdType ghostLev) const;
    %extend
    {
      MEDCouplingIMesh()
      {
        return MEDCouplingIMesh::New();
      }
      static MEDCouplingIMesh *New(const std::string& meshName, int spaceDim, PyObject *nodeStrct, PyObject *origin, PyObject *dxyz)
      {
        static const char msg0[]="MEDCouplingIMesh::New : error on 'origin' parameter !";
        static const char msg1[]="MEDCouplingIMesh::New : error on 'dxyz' parameter !";
        const mcIdType *nodeStrctPtr(0);
        const double *originPtr(0),*dxyzPtr(0);
        mcIdType sw,sz,val0;
        std::vector<mcIdType> bb0;
        nodeStrctPtr=convertIntStarLikePyObjToCppIntStar(nodeStrct,sw,sz,val0,bb0);
        //
        double val,val2;
        std::vector<double> bb,bb2;
        mcIdType sz1,sz2;
        originPtr=convertObjToPossibleCpp5_SingleCompo(origin,sw,val,bb,msg0,false,sz1);
        dxyzPtr=convertObjToPossibleCpp5_SingleCompo(dxyz,sw,val2,bb2,msg1,false,sz2);
        //
        return MEDCouplingIMesh::New(meshName,spaceDim,nodeStrctPtr,nodeStrctPtr+sz,originPtr,originPtr+sz1,dxyzPtr,dxyzPtr+sz2);
      }

      MEDCouplingIMesh(const std::string& meshName, int spaceDim, PyObject *nodeStrct, PyObject *origin, PyObject *dxyz)
      {
        return MEDCoupling_MEDCouplingIMesh_New__SWIG_1(meshName,spaceDim,nodeStrct,origin,dxyz);
      }

      void setNodeStruct(PyObject *nodeStrct)
      {
        mcIdType sw,sz,val0;
        std::vector<mcIdType> bb0;
        const mcIdType *nodeStrctPtr(convertIntStarLikePyObjToCppIntStar(nodeStrct,sw,sz,val0,bb0));
        self->setNodeStruct(nodeStrctPtr,nodeStrctPtr+sz);
      }

      void setOrigin(PyObject *origin)
      {
        static const char msg[]="MEDCouplingIMesh::setOrigin : invalid input 'origin' parameter ! integer, float, list/tuple of float, DataArrayDouble or DataArrayDoubleTuple supported !";
        double val;
        std::vector<double> bb;
        mcIdType sw,nbTuples;
        const double *originPtr(convertObjToPossibleCpp5_SingleCompo(origin,sw,val,bb,msg,false,nbTuples));
        self->setOrigin(originPtr,originPtr+nbTuples);
      }

      void setDXYZ(PyObject *dxyz)
      {
        static const char msg[]="MEDCouplingIMesh::setDXYZ : invalid input 'dxyz' parameter ! integer, float, list/tuple of float, DataArrayDouble or DataArrayDoubleTuple supported !";
        double val;
        std::vector<double> bb;
        mcIdType sw,nbTuples;
        const double *originPtr(convertObjToPossibleCpp5_SingleCompo(dxyz,sw,val,bb,msg,false,nbTuples));
        self->setDXYZ(originPtr,originPtr+nbTuples);
      }

      static void CondenseFineToCoarse(const std::vector<mcIdType>& coarseSt, const DataArrayDouble *fineDA, PyObject *fineLocInCoarse, const std::vector<mcIdType>& facts, DataArrayDouble *coarseDA)
      {
        std::vector< std::pair<mcIdType,mcIdType> > inp;
        convertPyToVectorPairInt(fineLocInCoarse,inp);
        MEDCouplingIMesh::CondenseFineToCoarse(coarseSt,fineDA,inp,facts,coarseDA);
      }

      static void CondenseFineToCoarseGhost(const std::vector<mcIdType>& coarseSt, const DataArrayDouble *fineDA, PyObject *fineLocInCoarse, const std::vector<mcIdType>& facts, DataArrayDouble *coarseDA, mcIdType ghostSize)
      {
        std::vector< std::pair<mcIdType,mcIdType> > inp;
        convertPyToVectorPairInt(fineLocInCoarse,inp);
        MEDCouplingIMesh::CondenseFineToCoarseGhost(coarseSt,fineDA,inp,facts,coarseDA,ghostSize);
      }

      static void SpreadCoarseToFine(const DataArrayDouble *coarseDA, const std::vector<mcIdType>& coarseSt, DataArrayDouble *fineDA, PyObject *fineLocInCoarse, const std::vector<mcIdType>& facts)
      {
        std::vector< std::pair<mcIdType,mcIdType> > inp;
        convertPyToVectorPairInt(fineLocInCoarse,inp);
        MEDCouplingIMesh::SpreadCoarseToFine(coarseDA,coarseSt,fineDA,inp,facts);
      }

      static void SpreadCoarseToFineGhost(const DataArrayDouble *coarseDA, const std::vector<mcIdType>& coarseSt, DataArrayDouble *fineDA, PyObject *fineLocInCoarse, const std::vector<mcIdType>& facts, mcIdType ghostSize)
      {
        std::vector< std::pair<mcIdType,mcIdType> > inp;
        convertPyToVectorPairInt(fineLocInCoarse,inp);
        MEDCouplingIMesh::SpreadCoarseToFineGhost(coarseDA,coarseSt,fineDA,inp,facts,ghostSize);
      }

      static void SpreadCoarseToFineGhostZone(const DataArrayDouble *coarseDA, const std::vector<mcIdType>& coarseSt, DataArrayDouble *fineDA, PyObject *fineLocInCoarse, const std::vector<mcIdType>& facts, mcIdType ghostSize)
      {
        std::vector< std::pair<mcIdType,mcIdType> > inp;
        convertPyToVectorPairInt(fineLocInCoarse,inp);
        MEDCouplingIMesh::SpreadCoarseToFineGhostZone(coarseDA,coarseSt,fineDA,inp,facts,ghostSize);
      }

      std::string __str__() const
      {
        return self->simpleRepr();
      }
      std::string __repr__() const
      {
        std::ostringstream oss;
        self->reprQuickOverview(oss);
        return oss.str();
      }
    }
  };

  //== MEDCouplingIMesh End

}

namespace MEDCoupling
{
  class MEDCouplingField : public MEDCoupling::RefCountObject, public MEDCoupling::TimeLabel
  {
  public:
    virtual void checkConsistencyLight() const;
    virtual bool areCompatibleForMerge(const MEDCouplingField *other) const;
    bool areStrictlyCompatible(const MEDCouplingField *other) const;
    bool areStrictlyCompatibleForMulDiv(const MEDCouplingField *other) const;
    virtual void copyTinyStringsFrom(const MEDCouplingField *other);
    void setMesh(const MEDCoupling::MEDCouplingMesh *mesh);
    void setName(const char *name);
    std::string getDescription() const;
    void setDescription(const char *desc);
    std::string getName() const;
    TypeOfField getTypeOfField() const;
    NatureOfField getNature() const;
    virtual void setNature(NatureOfField nat);
    DataArrayDouble *getLocalizationOfDiscr() const;
    MEDCouplingFieldDouble *buildMeasureField(bool isAbs) const;
    mcIdType getNumberOfTuplesExpected() const;
    mcIdType getNumberOfMeshPlacesExpected() const;
    void setGaussLocalizationOnType(INTERP_KERNEL::NormalizedCellType type, const std::vector<double>& refCoo,
                                    const std::vector<double>& gsCoo, const std::vector<double>& wg);
    void clearGaussLocalizations();
    MEDCouplingGaussLocalization& getGaussLocalization(int locId);
    mcIdType getNbOfGaussLocalization() const;
    mcIdType getGaussLocalizationIdOfOneCell(mcIdType cellId) const;
    const MEDCouplingGaussLocalization& getGaussLocalization(int locId) const;
    mcIdType getGaussLocalizationIdOfOneType(INTERP_KERNEL::NormalizedCellType type) const;
    void setDiscretization(MEDCouplingFieldDiscretization *newDisc);
    %extend {
      PyObject *getMesh() const
      {
        MEDCouplingMesh *ret1=const_cast<MEDCouplingMesh *>(self->getMesh());
        if(ret1)
          ret1->incrRef();
        return convertMesh(ret1,SWIG_POINTER_OWN | 0 );
      }

      PyObject *getDiscretization()
      {
        MEDCouplingFieldDiscretization *ret=self->getDiscretization();
        if(ret)
          ret->incrRef();
        return convertFieldDiscretization(ret,SWIG_POINTER_OWN | 0 );
      }

      PyObject *getGaussLocalizationIdsOfOneType(INTERP_KERNEL::NormalizedCellType type) const
      {
        std::set<mcIdType> ret=self->getGaussLocalizationIdsOfOneType(type);
        return convertIntArrToPyList3(ret);
      }

      PyObject *buildSubMeshData(PyObject *li) const
      {
        DataArrayIdType *ret1=0;
        MEDCouplingMesh *ret0=0;
        void *da=0;
        int res1=SWIG_ConvertPtr(li,&da,SWIGTITraits<mcIdType>::TI, 0 |  0 );
        if (!SWIG_IsOK(res1))
          {
            mcIdType size;
            INTERP_KERNEL::AutoPtr<mcIdType> tmp=convertPyToNewIntArr2(li,&size);
            ret0=self->buildSubMeshData(tmp,tmp+size,ret1);
          }
        else
          {
            DataArrayIdType *da2=reinterpret_cast< DataArrayIdType * >(da);
            if(!da2)
              throw INTERP_KERNEL::Exception("Not null DataArrayIdType instance expected !");
            da2->checkAllocated();
            ret0=self->buildSubMeshData(da2->getConstPointer(),da2->getConstPointer()+da2->getNbOfElems(),ret1);
          }
        PyObject *res = PyList_New(2);
        PyList_SetItem(res,0,convertMesh(ret0, SWIG_POINTER_OWN | 0 ));
        PyList_SetItem(res,1,SWIG_NewPointerObj((void*)ret1,SWIGTITraits<mcIdType>::TI,SWIG_POINTER_OWN | 0));
        return res;
      }

      PyObject *buildSubMeshDataRange(mcIdType begin, mcIdType end, mcIdType step) const
      {
        DataArrayIdType *ret1=0;
        mcIdType bb,ee,ss;
        MEDCouplingMesh *ret0=self->buildSubMeshDataRange(begin,end,step,bb,ee,ss,ret1);
        PyObject *res=PyTuple_New(2);
        PyTuple_SetItem(res,0,convertMesh(ret0, SWIG_POINTER_OWN | 0 ));
        if(ret1)
          PyTuple_SetItem(res,1,SWIG_NewPointerObj((void*)ret1,SWIGTITraits<mcIdType>::TI,SWIG_POINTER_OWN | 0));
        else
          {
            PyObject *res1=PySlice_New(PyInt_FromLong(bb),PyInt_FromLong(ee),PyInt_FromLong(ss));
            PyTuple_SetItem(res,1,res1);
          }
        return res;
      }

      DataArrayIdType *computeTupleIdsToSelectFromCellIds(PyObject *cellIds) const
      {
        mcIdType sw,sz(-1);
        mcIdType v0; std::vector<mcIdType> v1;
        const mcIdType *cellIdsBg(convertIntStarLikePyObjToCppIntStar(cellIds,sw,sz,v0,v1));
        return self->computeTupleIdsToSelectFromCellIds(cellIdsBg,cellIdsBg+sz);
      }

      void setGaussLocalizationOnCells(PyObject *li, const std::vector<double>& refCoo,
                                       const std::vector<double>& gsCoo, const std::vector<double>& wg)
      {
        void *da=0;
        int res1=SWIG_ConvertPtr(li,&da,SWIGTITraits<mcIdType>::TI, 0 |  0 );
        if (!SWIG_IsOK(res1))
          {
            mcIdType size;
            INTERP_KERNEL::AutoPtr<mcIdType> tmp=convertPyToNewIntArr2(li,&size);
            self->setGaussLocalizationOnCells(tmp,((mcIdType *)tmp)+size,refCoo,gsCoo,wg);
          }
        else
          {
            DataArrayIdType *da2=reinterpret_cast< DataArrayIdType * >(da);
            if(!da2)
              throw INTERP_KERNEL::Exception("Not null DataArrayIdType instance expected !");
            da2->checkAllocated();
            self->setGaussLocalizationOnCells(da2->getConstPointer(),da2->getConstPointer()+da2->getNbOfElems(),refCoo,gsCoo,wg);
          }
      }

      PyObject *getCellIdsHavingGaussLocalization(int locId) const
      {
        std::vector<mcIdType> tmp;
        self->getCellIdsHavingGaussLocalization(locId,tmp);
        DataArrayIdType *ret=DataArrayIdType::New();
        ret->alloc((mcIdType)tmp.size(),1);
        std::copy(tmp.begin(),tmp.end(),ret->getPointer());
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 );
      }

      mcIdType getNumberOfTuplesExpectedRegardingCode(PyObject *code, PyObject *idsPerType) const
      {
        std::vector<mcIdType> inp0;
        convertPyToNewIntArr4(code,1,3,inp0);
        std::vector<const DataArrayIdType *> inp1;
        convertFromPyObjVectorOfObj<const MEDCoupling::DataArrayIdType *>(idsPerType,SWIGTITraits<mcIdType>::TI,"DataArrayIdType",inp1);
        return self->getNumberOfTuplesExpectedRegardingCode(inp0,inp1);
      }
    }
  };

  class MEDCouplingFieldTemplate : public MEDCoupling::MEDCouplingField
  {
  public:
    static MEDCouplingFieldTemplate *New(const MEDCouplingFieldDouble& f);
    static MEDCouplingFieldTemplate *New(const MEDCouplingFieldFloat& f);
    static MEDCouplingFieldTemplate *New(const MEDCouplingFieldInt32& f);
    static MEDCouplingFieldTemplate *New(const MEDCouplingFieldInt64& f);
    static MEDCouplingFieldTemplate *New(TypeOfField type);
    std::string simpleRepr() const;
    std::string advancedRepr() const;
    bool isEqual(const MEDCouplingFieldTemplate *other, double meshPrec) const;
    bool isEqualWithoutConsideringStr(const MEDCouplingFieldTemplate *other, double meshPrec) const;
    %extend
       {
         MEDCouplingFieldTemplate(const MEDCouplingFieldDouble& f)
         {
           return MEDCouplingFieldTemplate::New(f);
         }

         MEDCouplingFieldTemplate(const MEDCouplingFieldFloat& f)
         {
           return MEDCouplingFieldTemplate::New(f);
         }

         MEDCouplingFieldTemplate(const MEDCouplingFieldInt32& f)
         {
           return MEDCouplingFieldTemplate::New(f);
         }

         MEDCouplingFieldTemplate(const MEDCouplingFieldInt64& f)
         {
           return MEDCouplingFieldTemplate::New(f);
         }

         MEDCouplingFieldTemplate(TypeOfField type)
         {
           return MEDCouplingFieldTemplate::New(type);
         }

         std::string __str__() const
         {
           return self->simpleRepr();
         }

         std::string __repr__() const
         {
           std::ostringstream oss;
           self->reprQuickOverview(oss);
           return oss.str();
         }

         PyObject *isEqualIfNotWhy(const MEDCouplingFieldTemplate *other, double meshPrec) const
         {
           std::string ret1;
           bool ret0=self->isEqualIfNotWhy(other,meshPrec,ret1);
           PyObject *ret=PyTuple_New(2);
           PyObject *ret0Py=ret0?Py_True:Py_False;
           Py_XINCREF(ret0Py);
           PyTuple_SetItem(ret,0,ret0Py);
           PyTuple_SetItem(ret,1,PyString_FromString(ret1.c_str()));
           return ret;
         }
       }
  };

  template<class T>
 class MEDCouplingFieldT : public MEDCoupling::MEDCouplingField
  {
  public:
    TypeOfTimeDiscretization getTimeDiscretization() const;
  protected:
    MEDCouplingFieldT();
    ~MEDCouplingFieldT();
  };

  %template(MEDCouplingFieldTdouble) MEDCoupling::MEDCouplingFieldT<double>;
  %template(MEDCouplingFieldTfloat) MEDCoupling::MEDCouplingFieldT<float>;
  %template(MEDCouplingFieldTint) MEDCoupling::MEDCouplingFieldT<int>;

  class MEDCouplingFieldInt32;
  class MEDCouplingFieldInt64;
  class MEDCouplingFieldFloat;

  class MEDCouplingFieldDouble : public MEDCouplingFieldT<double>
  {
  public:
    static MEDCouplingFieldDouble *New(TypeOfField type, TypeOfTimeDiscretization td=ONE_TIME);
    static MEDCouplingFieldDouble *New(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td=ONE_TIME);
    bool isEqual(const MEDCouplingFieldDouble *other, double meshPrec, double valsPrec) const;
    bool isEqualWithoutConsideringStr(const MEDCouplingFieldDouble *other, double meshPrec, double valsPrec) const;
    void setTimeUnit(const std::string& unit);
    std::string getTimeUnit() const;
    void synchronizeTimeWithSupport();
    void copyTinyAttrFrom(const MEDCouplingFieldDouble *other);
    void copyAllTinyAttrFrom(const MEDCouplingFieldDouble *other);
    std::string simpleRepr() const;
    std::string advancedRepr() const;
    std::string  writeVTK(const std::string& fileName, bool isBinary=true) const;
    MEDCouplingFieldInt32 *convertToIntField() const;
    MEDCouplingFieldFloat *convertToFloatField() const;
    MEDCouplingFieldDouble *clone(bool recDeepCpy) const;
    MEDCouplingFieldDouble *cloneWithMesh(bool recDeepCpy) const;
    MEDCouplingFieldDouble *deepCopy() const;
    MEDCouplingFieldDouble *buildNewTimeReprFromThis(TypeOfTimeDiscretization td, bool deepCopy) const;
    MEDCouplingFieldDouble *nodeToCellDiscretization() const;
    MEDCouplingFieldDouble *cellToNodeDiscretization() const;
    double getIJ(int tupleId, int compoId) const;
    double getIJK(int cellId, int nodeIdInCell, int compoId) const;
    void synchronizeTimeWithMesh();
    void setArray(DataArrayDouble *array);
    void setEndArray(DataArrayDouble *array);
    void setTime(double val, int iteration, int order);
    void setStartTime(double val, int iteration, int order);
    void setEndTime(double val, int iteration, int order);
    void applyLin(double a, double b, int compoId);
    void applyLin(double a, double b);
    int getNumberOfComponents() const;
    int getNumberOfTuples() const;
    int getNumberOfValues() const;
    void setTimeTolerance(double val);
    double getTimeTolerance() const;
    void setIteration(int it);
    void setEndIteration(int it);
    void setOrder(int order);
    void setEndOrder(int order);
    void setTimeValue(double val);
    void setEndTimeValue(double val);
    void changeUnderlyingMesh(const MEDCouplingMesh *other, int levOfCheck, double precOnMesh, double eps=1e-15);
    void substractInPlaceDM(const MEDCouplingFieldDouble *f, int levOfCheck, double precOnMesh, double eps=1e-15);
    bool mergeNodes(double eps, double epsOnVals=1e-15);
    bool mergeNodesCenter(double eps, double epsOnVals=1e-15);
    bool zipCoords(double epsOnVals=1e-15);
    bool zipConnectivity(int compType,double epsOnVals=1e-15);
    bool simplexize(int policy);
    MEDCouplingFieldDouble *doublyContractedProduct() const;
    MEDCouplingFieldDouble *determinant() const;
    MEDCouplingFieldDouble *eigenValues() const;
    MEDCouplingFieldDouble *eigenVectors() const;
    MEDCouplingFieldDouble *inverse() const;
    MEDCouplingFieldDouble *trace() const;
    MEDCouplingFieldDouble *deviator() const;
    MEDCouplingFieldDouble *magnitude() const;
    MEDCouplingFieldDouble *maxPerTuple() const;
    void changeNbOfComponents(std::size_t newNbOfComp, double dftValue=0.);
    void sortPerTuple(bool asc);
    MEDCouplingFieldDouble &operator=(double value);
    void fillFromAnalytic(int nbOfComp, const std::string& func);
    void fillFromAnalyticCompo(int nbOfComp, const std::string& func);
    void fillFromAnalyticNamedCompo(int nbOfComp, const std::vector<std::string>& varsOrder, const std::string& func);
    void applyFunc(int nbOfComp, const std::string& func);
    void applyFuncCompo(int nbOfComp, const std::string& func);
    void applyFuncNamedCompo(int nbOfComp, const std::vector<std::string>& varsOrder, const std::string& func);
    void applyFunc(int nbOfComp, double val);
    void applyFunc(const std::string& func);
    void applyFuncFast32(const std::string& func);
    void applyFuncFast64(const std::string& func);
    double accumulate(int compId) const;
    double getMaxValue() const;
    double getMinValue() const;
    double getAverageValue() const;
    double norm2() const;
    //do not put a default value to isWAbs because confusion in python with overloaded getWeightedAverageValue method
    double getWeightedAverageValue(int compId, bool isWAbs) const;
    double integral(int compId, bool isWAbs) const;
    double normL1(int compId) const;
    double normL2(int compId) const;
    double normMax(int compId) const;
    DataArrayIdType *findIdsInRange(double vmin, double vmax) const;
    MEDCouplingFieldDouble *buildSubPartRange(int begin, int end, int step) const;
    static MEDCouplingFieldDouble *MergeFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    static MEDCouplingFieldDouble *MeldFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    static MEDCouplingFieldDouble *DotFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    MEDCouplingFieldDouble *dot(const MEDCouplingFieldDouble& other) const;
    static MEDCouplingFieldDouble *CrossProductFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    MEDCouplingFieldDouble *crossProduct(const MEDCouplingFieldDouble& other) const;
    static MEDCouplingFieldDouble *MaxFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    MEDCouplingFieldDouble *max(const MEDCouplingFieldDouble& other) const;
    static MEDCouplingFieldDouble *MinFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    static MEDCouplingFieldDouble *AddFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    static MEDCouplingFieldDouble *SubstractFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    static MEDCouplingFieldDouble *MultiplyFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    static MEDCouplingFieldDouble *DivideFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
    MEDCouplingFieldDouble *min(const MEDCouplingFieldDouble& other) const;
    MEDCouplingFieldDouble *negate() const;
    %extend {
      MEDCouplingFieldDouble(TypeOfField type, TypeOfTimeDiscretization td=ONE_TIME)
      {
        return MEDCouplingFieldDouble::New(type,td);
      }

      MEDCouplingFieldDouble(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td=ONE_TIME)
      {
        return MEDCouplingFieldDouble::New(ft,td);
      }

      std::string __str__() const
      {
        return self->simpleRepr();
      }

      std::string __repr__() const
      {
        std::ostringstream oss;
        self->reprQuickOverview(oss);
        return oss.str();
      }

      PyObject *isEqualIfNotWhy(const MEDCouplingFieldDouble *other, double meshPrec, double valsPrec) const
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

      MEDCouplingFieldDouble *voronoize(double eps) const
      {
        MCAuto<MEDCouplingFieldDouble> ret(self->voronoize(eps));
        return ret.retn();
      }

      MEDCouplingFieldDouble *convertQuadraticCellsToLinear() const
      {
        MCAuto<MEDCouplingFieldDouble> ret(self->convertQuadraticCellsToLinear());
        return ret.retn();
      }

      MEDCouplingFieldDouble *computeVectorFieldCyl(PyObject *center, PyObject *vector) const
      {
        const char msg[]="Python wrap of MEDCouplingFieldDouble::computeVectorFieldCyl : ";
        double val,val2;
        DataArrayDouble *a,*a2;
        DataArrayDoubleTuple *aa,*aa2;
        std::vector<double> bb,bb2;
        mcIdType sw;
        const double *centerPtr=convertObjToPossibleCpp5_Safe(center,sw,val,a,aa,bb,msg,1,3,true);
        const double *vectorPtr=convertObjToPossibleCpp5_Safe(vector,sw,val2,a2,aa2,bb2,msg,1,3,true);
        return self->computeVectorFieldCyl(centerPtr,vectorPtr);
      }

      DataArrayDouble *getArray()
      {
        DataArrayDouble *ret=self->getArray();
        if(ret)
          ret->incrRef();
        return ret;
      }

      PyObject *getArrays() const
      {
        std::vector<DataArrayDouble *> arrs=self->getArrays();
        for(std::vector<DataArrayDouble *>::iterator it=arrs.begin();it!=arrs.end();it++)
          if(*it)
            (*it)->incrRef();
        std::size_t sz=arrs.size();
        PyObject *ret=PyTuple_New(sz);
        for(std::size_t i=0;i<sz;i++)
          {
            if(arrs[i])
              PyTuple_SetItem(ret,i,SWIG_NewPointerObj(SWIG_as_voidptr(arrs[i]),SWIGTYPE_p_MEDCoupling__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
            else
              PyTuple_SetItem(ret,i,SWIG_NewPointerObj(SWIG_as_voidptr(0),SWIGTYPE_p_MEDCoupling__DataArrayDouble, 0 | 0 ));
          }
        return ret;
      }

      void setArrays(PyObject *ls)
      {
        std::vector<const DataArrayDouble *> tmp;
        convertFromPyObjVectorOfObj<const DataArrayDouble *>(ls,SWIGTYPE_p_MEDCoupling__DataArrayDouble,"DataArrayDouble",tmp);
        std::size_t sz=tmp.size();
        std::vector<DataArrayDouble *> arrs(sz);
        for(std::size_t i=0;i<sz;i++)
          arrs[i]=const_cast<DataArrayDouble *>(tmp[i]);
        self->setArrays(arrs);
      }

      DataArrayDouble *getEndArray()
      {
        DataArrayDouble *ret=self->getEndArray();
        if(ret)
          ret->incrRef();
        return ret;
      }

      PyObject *getValueOn(PyObject *sl) const
      {
        double val;
        DataArrayDouble *a;
        DataArrayDoubleTuple *aa;
        std::vector<double> bb;
        mcIdType sw;
        const MEDCouplingMesh *mesh=self->getMesh();
        if(!mesh)
          throw INTERP_KERNEL::Exception("Python wrap of MEDCouplingFieldDouble::getValueOn : no underlying mesh !");
        int spaceDim=mesh->getSpaceDimension();
        const char msg[]="Python wrap of MEDCouplingFieldDouble::getValueOn : ";
        const double *spaceLoc=convertObjToPossibleCpp5_Safe(sl,sw,val,a,aa,bb,msg,1,spaceDim,true);
        //
        mcIdType sz=ToIdType(self->getNumberOfComponents());
        INTERP_KERNEL::AutoPtr<double> res=new double[sz];
        self->getValueOn(spaceLoc,res);
        return convertDblArrToPyList<double>(res,sz);
      }

       PyObject *getValueOnPos(mcIdType i, mcIdType j, mcIdType k) const
       {
         mcIdType sz=ToIdType(self->getNumberOfComponents());
         INTERP_KERNEL::AutoPtr<double> res=new double[sz];
         self->getValueOnPos(i,j,k,res);
         return convertDblArrToPyList<double>(res,sz);
       }

      DataArrayDouble *getValueOnMulti(PyObject *locs) const
      {
        const MEDCouplingMesh *mesh(self->getMesh());
        if(!mesh)
          throw INTERP_KERNEL::Exception("Python wrap MEDCouplingFieldDouble::getValueOnMulti : lying on a null mesh !");
        //
        mcIdType sw,nbPts;
        double v0; MEDCoupling::DataArrayDouble *v1(0); MEDCoupling::DataArrayDoubleTuple *v2(0); std::vector<double> v3;
        const double *inp=convertObjToPossibleCpp5_Safe2(locs,sw,v0,v1,v2,v3,"wrap of MEDCouplingFieldDouble::getValueOnMulti",
                                                         mesh->getSpaceDimension(),true,nbPts);
        return self->getValueOnMulti(inp,(int)nbPts);
      }

      PyObject *getValueOn(PyObject *sl, double time) const
      {
        double val;
        DataArrayDouble *a;
        DataArrayDoubleTuple *aa;
        std::vector<double> bb;
        mcIdType sw;
        const MEDCouplingMesh *mesh=self->getMesh();
        if(!mesh)
          throw INTERP_KERNEL::Exception("Python wrap of MEDCouplingFieldDouble::getValueOn : no underlying mesh !");
        int spaceDim=mesh->getSpaceDimension();
        const char msg[]="Python wrap of MEDCouplingFieldDouble::getValueOn : ";
        const double *spaceLoc=convertObjToPossibleCpp5_Safe(sl,sw,val,a,aa,bb,msg,1,spaceDim,true);
        //
        //
        mcIdType sz=ToIdType(self->getNumberOfComponents());
        INTERP_KERNEL::AutoPtr<double> res=new double[sz];
        self->getValueOn(spaceLoc,time,res);
        return convertDblArrToPyList<double>(res,sz);
      }

      void setValues(PyObject *li, PyObject *nbOfTuples=0, PyObject *nbOfComp=0)
      {
        if(self->getArray()!=0)
          MEDCoupling_DataArrayDouble_setValues__SWIG_0(self->getArray(),li,nbOfTuples,nbOfComp);
        else
          {
            MCAuto<DataArrayDouble> arr=DataArrayDouble::New();
            MEDCoupling_DataArrayDouble_setValues__SWIG_0(arr,li,nbOfTuples,nbOfComp);
            self->setArray(arr);
          }
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
        mcIdType sz=ToIdType(self->getNumberOfComponents());
        INTERP_KERNEL::AutoPtr<double> tmp=new double[sz];
        self->accumulate(tmp);
        return convertDblArrToPyList<double>(tmp,sz);
      }
      PyObject *integral(bool isWAbs) const
      {
        mcIdType sz=ToIdType(self->getNumberOfComponents());
        INTERP_KERNEL::AutoPtr<double> tmp=new double[sz];
        self->integral(isWAbs,tmp);
        return convertDblArrToPyList<double>(tmp,sz);
      }
      PyObject *getWeightedAverageValue(bool isWAbs=true) const
      {
        mcIdType sz=ToIdType(self->getNumberOfComponents());
        INTERP_KERNEL::AutoPtr<double> tmp=new double[sz];
        self->getWeightedAverageValue(tmp,isWAbs);
        return convertDblArrToPyList<double>(tmp,sz);
      }
      PyObject *normL1() const
      {
        mcIdType sz=ToIdType(self->getNumberOfComponents());
        INTERP_KERNEL::AutoPtr<double> tmp=new double[sz];
        self->normL1(tmp);
        return convertDblArrToPyList<double>(tmp,sz);
      }
      PyObject *normL2() const
      {
        mcIdType sz=ToIdType(self->getNumberOfComponents());
        INTERP_KERNEL::AutoPtr<double> tmp=new double[sz];
        self->normL2(tmp);
        return convertDblArrToPyList<double>(tmp,sz);
      }
      PyObject *normMax() const
      {
        mcIdType sz=ToIdType(self->getNumberOfComponents());
        INTERP_KERNEL::AutoPtr<double> tmp=new double[sz];
        self->normMax(tmp);
        return convertDblArrToPyList<double>(tmp,sz);
      }
      void renumberCells(PyObject *li, bool check=true)
      {
        mcIdType szArr,sw,iTypppArr;
        std::vector<mcIdType> stdvecTyyppArr;
        const mcIdType *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
        self->renumberCells(tmp,check);
      }

      void renumberCellsWithoutMesh(PyObject *li, bool check=true)
      {
        mcIdType szArr,sw,iTypppArr;
        std::vector<mcIdType> stdvecTyyppArr;
        const mcIdType *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
        self->renumberCellsWithoutMesh(tmp,check);
      }

      void renumberNodes(PyObject *li, double eps=1e-15)
      {
        mcIdType szArr,sw,iTypppArr;
        std::vector<mcIdType> stdvecTyyppArr;
        const mcIdType *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
        self->renumberNodes(tmp,eps);
      }

      void renumberNodesWithoutMesh(PyObject *li, mcIdType newNbOfNodes, double eps=1e-15)
      {
        mcIdType szArr,sw,iTypppArr;
        std::vector<mcIdType> stdvecTyyppArr;
        const mcIdType *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
        self->renumberNodesWithoutMesh(tmp,newNbOfNodes,eps);
      }

      MEDCouplingFieldDouble *buildSubPart(PyObject *li) const
      {
        return fieldT_buildSubPart(self,li);
      }

      MEDCouplingFieldDouble *__getitem__(PyObject *li) const
      {
        return fieldT__getitem__(self,li);
      }

      PyObject *getMaxValue2() const
      {
        DataArrayIdType *tmp;
        double r1=self->getMaxValue2(tmp);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,PyFloat_FromDouble(r1));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(tmp),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *getMinValue2() const
      {
        DataArrayIdType *tmp;
        double r1=self->getMinValue2(tmp);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,PyFloat_FromDouble(r1));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(tmp),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      MEDCouplingFieldDouble *keepSelectedComponents(PyObject *li) const
      {
        std::vector<std::size_t> tmp;
        convertPyToNewIntArr3(li,tmp);
        return self->keepSelectedComponents(tmp);
      }

      void setSelectedComponents(const MEDCouplingFieldDouble *f, PyObject *li)
      {
        std::vector<std::size_t> tmp;
        convertPyToNewIntArr3(li,tmp);
        self->setSelectedComponents(f,tmp);
      }

      MEDCouplingFieldDouble *extractSlice3D(PyObject *origin, PyObject *vec, double eps) const
      {
        double val,val2;
        DataArrayDouble *a,*a2;
        DataArrayDoubleTuple *aa,*aa2;
        std::vector<double> bb,bb2;
        mcIdType sw;
        int spaceDim=3;
        const char msg[]="Python wrap of MEDCouplingFieldDouble::extractSlice3D : 1st parameter for origin.";
        const char msg2[]="Python wrap of MEDCouplingFieldDouble::extractSlice3D : 2nd parameter for vector.";
        const double *orig=convertObjToPossibleCpp5_Safe(origin,sw,val,a,aa,bb,msg,1,spaceDim,true);
        const double *vect=convertObjToPossibleCpp5_Safe(vec,sw,val2,a2,aa2,bb2,msg2,1,spaceDim,true);
        //
        return self->extractSlice3D(orig,vect,eps);
      }

      MEDCouplingFieldDouble *__add__(PyObject *obj)
      {
        return MEDCoupling_MEDCouplingFieldDouble___add__Impl(self,obj);
      }

      MEDCouplingFieldDouble *__radd__(PyObject *obj)
      {
        return MEDCoupling_MEDCouplingFieldDouble___radd__Impl(self,obj);
      }

      MEDCouplingFieldDouble *__sub__(PyObject *obj)
      {
        const char msg[]="Unexpected situation in MEDCouplingFieldDouble.__sub__ ! Expecting a not null MEDCouplingFieldDouble or DataArrayDouble or DataArrayDoubleTuple instance, or a list of double, or a double.";
        const char msg2[]="in MEDCouplingFieldDouble.__sub__ : self field has no Array of values set !";
        void *argp;
        //
        if(SWIG_IsOK(SWIG_ConvertPtr(obj,&argp,SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDouble,0|0)))
          {
            MEDCouplingFieldDouble *other=reinterpret_cast< MEDCoupling::MEDCouplingFieldDouble * >(argp);
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
        mcIdType sw;
        convertDoubleStarLikePyObjToCpp_2(obj,sw,val,a,aa,bb);
        switch(sw)
          {
          case 1:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MCAuto<DataArrayDouble> ret=self->getArray()->deepCopy();
              ret->applyLin(1.,-val);
              MCAuto<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(ret);
              return ret2.retn();
            }
          case 2:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MCAuto<DataArrayDouble> ret=DataArrayDouble::Substract(self->getArray(),a);
              MCAuto<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(ret);
              return ret2.retn();
            }
          case 3:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MCAuto<DataArrayDouble> aaa=aa->buildDADouble(1,self->getNumberOfComponents());
              MCAuto<DataArrayDouble> ret=DataArrayDouble::Substract(self->getArray(),aaa);
              MCAuto<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(ret);
              return ret2.retn();
            }
          case 4:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MCAuto<DataArrayDouble> aaa=DataArrayDouble::New(); aaa->useArray(&bb[0],false,DeallocType::CPP_DEALLOC,1,bb.size());
              MCAuto<DataArrayDouble> ret=DataArrayDouble::Substract(self->getArray(),aaa);
              MCAuto<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(ret);
              return ret2.retn();
            }
          default:
            { throw INTERP_KERNEL::Exception(msg); }
          }
      }

      MEDCouplingFieldDouble *__rsub__(PyObject *obj)
      {
        return MEDCoupling_MEDCouplingFieldDouble___rsub__Impl(self,obj);
      }

      MEDCouplingFieldDouble *__mul__(PyObject *obj)
      {
        return MEDCoupling_MEDCouplingFieldDouble___mul__Impl(self,obj);
      }

      MEDCouplingFieldDouble *__rmul__(PyObject *obj)
      {
        return MEDCoupling_MEDCouplingFieldDouble___rmul__Impl(self,obj);
      }

      MEDCouplingFieldDouble *__div__(PyObject *obj)
      {
        const char msg[]="Unexpected situation in MEDCouplingFieldDouble.__div__ ! Expecting a not null MEDCouplingFieldDouble or DataArrayDouble or DataArrayDoubleTuple instance, or a list of double, or a double.";
        const char msg2[]="in MEDCouplingFieldDouble.__div__ : self field has no Array of values set !";
        void *argp;
        //
        if(SWIG_IsOK(SWIG_ConvertPtr(obj,&argp,SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDouble,0|0)))
          {
            MEDCouplingFieldDouble *other=reinterpret_cast< MEDCoupling::MEDCouplingFieldDouble * >(argp);
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
        mcIdType sw;
        convertDoubleStarLikePyObjToCpp_2(obj,sw,val,a,aa,bb);
        switch(sw)
          {
          case 1:
            {
              if(val==0.)
                throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble.__div__ : trying to divide by zero !");
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MCAuto<DataArrayDouble> ret=self->getArray()->deepCopy();
              ret->applyLin(1./val,0);
              MCAuto<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(ret);
              return ret2.retn();
            }
          case 2:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MCAuto<DataArrayDouble> ret=DataArrayDouble::Divide(self->getArray(),a);
              MCAuto<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(ret);
              return ret2.retn();
            }
          case 3:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MCAuto<DataArrayDouble> aaa=aa->buildDADouble(1,self->getNumberOfComponents());
              MCAuto<DataArrayDouble> ret=DataArrayDouble::Divide(self->getArray(),aaa);
              MCAuto<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(ret);
              return ret2.retn();
            }
          case 4:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MCAuto<DataArrayDouble> aaa=DataArrayDouble::New(); aaa->useArray(&bb[0],false,DeallocType::CPP_DEALLOC,1,bb.size());
              MCAuto<DataArrayDouble> ret=DataArrayDouble::Divide(self->getArray(),aaa);
              MCAuto<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(ret);
              return ret2.retn();
            }
          default:
            { throw INTERP_KERNEL::Exception(msg); }
          }
      }

      MEDCouplingFieldDouble *__rdiv__(PyObject *obj)
      {
        return MEDCoupling_MEDCouplingFieldDouble___rdiv__Impl(self,obj);
      }

      MEDCouplingFieldDouble *__pow__(PyObject *obj)
      {
        const char msg[]="Unexpected situation in MEDCouplingFieldDouble.__pow__ ! Expecting a not null MEDCouplingFieldDouble or DataArrayDouble or DataArrayDoubleTuple instance, or a list of double, or a double.";
        const char msg2[]="in MEDCouplingFieldDouble.__pow__ : self field has no Array of values set !";
        void *argp;
        //
        if(SWIG_IsOK(SWIG_ConvertPtr(obj,&argp,SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDouble,0|0)))
          {
            MEDCouplingFieldDouble *other=reinterpret_cast< MEDCoupling::MEDCouplingFieldDouble * >(argp);
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
        mcIdType sw;
        convertDoubleStarLikePyObjToCpp_2(obj,sw,val,a,aa,bb);
        switch(sw)
          {
          case 1:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MCAuto<DataArrayDouble> ret=self->getArray()->deepCopy();
              ret->applyPow(val);
              MCAuto<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(ret);
              return ret2.retn();
            }
          case 2:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MCAuto<DataArrayDouble> ret=DataArrayDouble::Pow(self->getArray(),a);
              MCAuto<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(ret);
              return ret2.retn();
            }
          case 3:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MCAuto<DataArrayDouble> aaa=aa->buildDADouble(1,self->getNumberOfComponents());
              MCAuto<DataArrayDouble> ret=DataArrayDouble::Pow(self->getArray(),aaa);
              MCAuto<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(ret);
              return ret2.retn();
            }
          case 4:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MCAuto<DataArrayDouble> aaa=DataArrayDouble::New(); aaa->useArray(&bb[0],false,DeallocType::CPP_DEALLOC,1,bb.size());
              MCAuto<DataArrayDouble> ret=DataArrayDouble::Pow(self->getArray(),aaa);
              MCAuto<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(ret);
              return ret2.retn();
            }
          default:
            { throw INTERP_KERNEL::Exception(msg); }
          }
      }

      MEDCouplingFieldDouble *__neg__() const
      {
        return self->negate();
      }

      PyObject *___iadd___(PyObject *trueSelf, PyObject *obj)
      {
        const char msg[]="Unexpected situation in MEDCouplingFieldDouble.__iadd__ ! Expecting a not null MEDCouplingFieldDouble or DataArrayDouble or DataArrayDoubleTuple instance, or a list of double, or a double.";
        const char msg2[]="in MEDCouplingFieldDouble.__iadd__ : self field has no Array of values set !";
        void *argp;
        //
        if(SWIG_IsOK(SWIG_ConvertPtr(obj,&argp,SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDouble,0|0)))
          {
            MEDCouplingFieldDouble *other=reinterpret_cast< MEDCoupling::MEDCouplingFieldDouble * >(argp);
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
        mcIdType sw;
        convertDoubleStarLikePyObjToCpp_2(obj,sw,val,a,aa,bb);
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
              MCAuto<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(a);
              *self+=*ret2;
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 3:
            {
              MCAuto<DataArrayDouble> aaa=aa->buildDADouble(1,self->getNumberOfComponents());
              MCAuto<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(aaa);
              *self+=*ret2;
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 4:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MCAuto<DataArrayDouble> aaa=DataArrayDouble::New(); aaa->useArray(&bb[0],false,DeallocType::CPP_DEALLOC,1,bb.size());
              self->getArray()->addEqual(aaa);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          default:
            { throw INTERP_KERNEL::Exception(msg); }
          }
      }

      PyObject *___isub___(PyObject *trueSelf, PyObject *obj)
      {
        const char msg[]="Unexpected situation in MEDCouplingFieldDouble.__isub__ ! Expecting a not null MEDCouplingFieldDouble or DataArrayDouble or DataArrayDoubleTuple instance, or a list of double, or a double.";
        const char msg2[]="in MEDCouplingFieldDouble.__isub__ : self field has no Array of values set !";
        void *argp;
        //
        if(SWIG_IsOK(SWIG_ConvertPtr(obj,&argp,SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDouble,0|0)))
          {
            MEDCouplingFieldDouble *other=reinterpret_cast< MEDCoupling::MEDCouplingFieldDouble * >(argp);
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
        mcIdType sw;
        convertDoubleStarLikePyObjToCpp_2(obj,sw,val,a,aa,bb);
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
              MCAuto<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(a);
              *self-=*ret2;
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 3:
            {
              MCAuto<DataArrayDouble> aaa=aa->buildDADouble(1,self->getNumberOfComponents());
              MCAuto<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(aaa);
              *self-=*ret2;
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 4:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MCAuto<DataArrayDouble> aaa=DataArrayDouble::New(); aaa->useArray(&bb[0],false,DeallocType::CPP_DEALLOC,1,bb.size());
              self->getArray()->substractEqual(aaa);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          default:
            { throw INTERP_KERNEL::Exception(msg); }
          }
      }

      PyObject *___imul___(PyObject *trueSelf, PyObject *obj)
      {
        const char msg[]="Unexpected situation in MEDCouplingFieldDouble.__imul__ ! Expecting a not null MEDCouplingFieldDouble or DataArrayDouble or DataArrayDoubleTuple instance, or a list of double, or a double.";
        const char msg2[]="in MEDCouplingFieldDouble.__imul__ : self field has no Array of values set !";
        void *argp;
        //
        if(SWIG_IsOK(SWIG_ConvertPtr(obj,&argp,SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDouble,0|0)))
          {
            MEDCouplingFieldDouble *other=reinterpret_cast< MEDCoupling::MEDCouplingFieldDouble * >(argp);
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
        mcIdType sw;
        convertDoubleStarLikePyObjToCpp_2(obj,sw,val,a,aa,bb);
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
              MCAuto<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(a);
              *self*=*ret2;
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 3:
            {
              MCAuto<DataArrayDouble> aaa=aa->buildDADouble(1,self->getNumberOfComponents());
              MCAuto<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(aaa);
              *self*=*ret2;
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 4:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MCAuto<DataArrayDouble> aaa=DataArrayDouble::New(); aaa->useArray(&bb[0],false,DeallocType::CPP_DEALLOC,1,(int)bb.size());
              self->getArray()->multiplyEqual(aaa);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          default:
            { throw INTERP_KERNEL::Exception(msg); }
          }
      }

      PyObject *___idiv___(PyObject *trueSelf, PyObject *obj)
      {
        const char msg[]="Unexpected situation in MEDCouplingFieldDouble.__idiv__ ! Expecting a not null MEDCouplingFieldDouble or DataArrayDouble or DataArrayDoubleTuple instance, or a list of double, or a double.";
        const char msg2[]="in MEDCouplingFieldDouble.__idiv__ : self field has no Array of values set !";
        void *argp;
        //
        if(SWIG_IsOK(SWIG_ConvertPtr(obj,&argp,SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDouble,0|0)))
          {
            MEDCouplingFieldDouble *other=reinterpret_cast< MEDCoupling::MEDCouplingFieldDouble * >(argp);
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
        mcIdType sw;
        convertDoubleStarLikePyObjToCpp_2(obj,sw,val,a,aa,bb);
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
              MCAuto<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(a);
              *self/=*ret2;
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 3:
            {
              MCAuto<DataArrayDouble> aaa=aa->buildDADouble(1,self->getNumberOfComponents());
              MCAuto<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(aaa);
              *self/=*ret2;
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 4:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MCAuto<DataArrayDouble> aaa=DataArrayDouble::New(); aaa->useArray(&bb[0],false,DeallocType::CPP_DEALLOC,1,bb.size());
              self->getArray()->divideEqual(aaa);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          default:
            { throw INTERP_KERNEL::Exception(msg); }
          }
      }

      PyObject *___ipow___(PyObject *trueSelf, PyObject *obj)
      {
        const char msg[]="Unexpected situation in MEDCouplingFieldDouble.__ipow__ ! Expecting a not null MEDCouplingFieldDouble or DataArrayDouble or DataArrayDoubleTuple instance, or a list of double, or a double.";
        const char msg2[]="in MEDCouplingFieldDouble.__ipow__ : self field has no Array of values set !";
        void *argp;
        //
        if(SWIG_IsOK(SWIG_ConvertPtr(obj,&argp,SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDouble,0|0)))
          {
            MEDCouplingFieldDouble *other=reinterpret_cast< MEDCoupling::MEDCouplingFieldDouble * >(argp);
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
        mcIdType sw;
        convertDoubleStarLikePyObjToCpp_2(obj,sw,val,a,aa,bb);
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
              MCAuto<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(a);
              *self^=*ret2;
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 3:
            {
              MCAuto<DataArrayDouble> aaa=aa->buildDADouble(1,self->getNumberOfComponents());
              MCAuto<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(aaa);
              *self^=*ret2;
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 4:
            {
              if(!self->getArray())
                throw INTERP_KERNEL::Exception(msg2);
              MCAuto<DataArrayDouble> aaa=DataArrayDouble::New(); aaa->useArray(&bb[0],false,DeallocType::CPP_DEALLOC,1,bb.size());
              self->getArray()->powEqual(aaa);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          default:
            { throw INTERP_KERNEL::Exception(msg); }
          }
      }

      static MEDCouplingFieldDouble *MergeFields(PyObject *li)
      {
        std::vector<const MEDCouplingFieldDouble *> tmp;
        convertFromPyObjVectorOfObj<const MEDCoupling::MEDCouplingFieldDouble *>(li,SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDouble,"MEDCouplingFieldDouble",tmp);
        return MEDCouplingFieldDouble::MergeFields(tmp);
      }

      static std::string WriteVTK(const char *fileName, PyObject *li, bool isBinary=true)
      {
        std::vector<const MEDCouplingFieldDouble *> tmp;
        convertFromPyObjVectorOfObj<const MEDCoupling::MEDCouplingFieldDouble *>(li,SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDouble,"MEDCouplingFieldDouble",tmp);
        return MEDCouplingFieldDouble::WriteVTK(fileName,tmp,isBinary);
      }

      PyObject *getTinySerializationInformation() const
      {
        return field_getTinySerializationInformation<MEDCouplingFieldDouble>(self);
      }

      PyObject *serialize() const
      {
        return field_serialize<double>(self);
      }

      PyObject *__getstate__() const
      {
        return field__getstate__<MEDCouplingFieldDouble>(self,MEDCoupling_MEDCouplingFieldDouble_getTinySerializationInformation,MEDCoupling_MEDCouplingFieldDouble_serialize);
      }

      void __setstate__(PyObject *inp)
      {
        field__setstate__<double>(self,inp);
      }
    }
  };

  class MEDCouplingMultiFields : public RefCountObject, public TimeLabel
  {
  public:
    int getNumberOfFields() const;
    MEDCouplingMultiFields *deepCopy() const;
    virtual std::string simpleRepr() const;
    virtual std::string advancedRepr() const;
    virtual bool isEqual(const MEDCouplingMultiFields *other, double meshPrec, double valsPrec) const;
    virtual bool isEqualWithoutConsideringStr(const MEDCouplingMultiFields *other, double meshPrec, double valsPrec) const;
    virtual void checkConsistencyLight() const;
    %extend
       {
         std::string __str__() const
         {
           return self->simpleRepr();
         }
         static MEDCouplingMultiFields *New(PyObject *li)
         {
           std::vector<const MEDCoupling::MEDCouplingFieldDouble *> tmp;
           convertFromPyObjVectorOfObj<const MEDCoupling::MEDCouplingFieldDouble *>(li,SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDouble,"MEDCouplingFieldDouble",tmp);
           std::size_t sz=tmp.size();
           std::vector<MEDCouplingFieldDouble *> fs(sz);
           for(std::size_t i=0;i<sz;i++)
             fs[i]=const_cast<MEDCouplingFieldDouble *>(tmp[i]);
           return MEDCouplingMultiFields::New(fs);
         }
         MEDCouplingMultiFields(PyObject *li)
         {
           std::vector<const MEDCoupling::MEDCouplingFieldDouble *> tmp;
           convertFromPyObjVectorOfObj<const MEDCoupling::MEDCouplingFieldDouble *>(li,SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDouble,"MEDCouplingFieldDouble",tmp);
           std::size_t sz=tmp.size();
           std::vector<MEDCouplingFieldDouble *> fs(sz);
           for(std::size_t i=0;i<sz;i++)
             fs[i]=const_cast<MEDCouplingFieldDouble *>(tmp[i]);
           return MEDCouplingMultiFields::New(fs);
         }
         PyObject *getFields() const
         {
           std::vector<const MEDCouplingFieldDouble *> fields=self->getFields();
           std::size_t sz=fields.size();
           PyObject *res = PyList_New(sz);
           for(std::size_t i=0;i<sz;i++)
             {
               if(fields[i])
                 {
                   fields[i]->incrRef();
                   PyList_SetItem(res,i,SWIG_NewPointerObj(SWIG_as_voidptr(fields[i]),SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDouble, SWIG_POINTER_OWN | 0 ));
                 }
               else
                 {
                   PyList_SetItem(res,i,SWIG_NewPointerObj(SWIG_as_voidptr(0),SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDouble, 0 ));
                 }
             }
           return res;
         }
         PyObject *getFieldAtPos(int id) const
         {
           const MEDCouplingFieldDouble *ret=self->getFieldAtPos(id);
           if(ret)
             {
               ret->incrRef();
               return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDouble, SWIG_POINTER_OWN | 0 );
             }
           else
             return SWIG_NewPointerObj(SWIG_as_voidptr(0),SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDouble, 0 );
         }
         PyObject *getMeshes() const
         {
           std::vector<MEDCouplingMesh *> ms=self->getMeshes();
           std::size_t sz=ms.size();
           PyObject *res = PyList_New(sz);
           for(std::size_t i=0;i<sz;i++)
             {
               if(ms[i])
                 {
                   ms[i]->incrRef();
                   PyList_SetItem(res,i,convertMesh(ms[i], SWIG_POINTER_OWN | 0 ));
                 }
               else
                 {
                   PyList_SetItem(res,i,SWIG_NewPointerObj(SWIG_as_voidptr(0),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, 0 ));
                 }
             }
           return res;
         }
         PyObject *getDifferentMeshes() const
         {
           std::vector<int> refs;
           std::vector<MEDCouplingMesh *> ms=self->getDifferentMeshes(refs);
           std::size_t sz=ms.size();
           PyObject *res = PyList_New(sz);
           for(std::size_t i=0;i<sz;i++)
             {
               if(ms[i])
                 {
                   ms[i]->incrRef();
                   PyList_SetItem(res,i,convertMesh(ms[i], SWIG_POINTER_OWN | 0 ));
                 }
               else
                 {
                   PyList_SetItem(res,i,SWIG_NewPointerObj(SWIG_as_voidptr(0),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, 0 ));
                 }
             }
           //
           PyObject *ret=PyTuple_New(2);
           PyTuple_SetItem(ret,0,res);
           PyTuple_SetItem(ret,1,convertIntArrToPyList2(refs));
           return ret;
         }
         PyObject *getArrays() const
         {
           std::vector<DataArrayDouble *> ms=self->getArrays();
           std::size_t sz=ms.size();
           PyObject *res = PyList_New(sz);
           for(std::size_t i=0;i<sz;i++)
             {
               if(ms[i])
                 {
                   ms[i]->incrRef();
                   PyList_SetItem(res,i,SWIG_NewPointerObj(SWIG_as_voidptr(ms[i]),SWIGTYPE_p_MEDCoupling__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
                 }
               else
                 {
                   PyList_SetItem(res,i,SWIG_NewPointerObj(SWIG_as_voidptr(0),SWIGTYPE_p_MEDCoupling__DataArrayDouble, 0 ));
                 }
             }
           return res;
         }
         PyObject *getDifferentArrays() const
         {
           std::vector< std::vector<int> > refs;
           std::vector<DataArrayDouble *> ms=self->getDifferentArrays(refs);
           std::size_t sz=ms.size();
           PyObject *res = PyList_New(sz);
           PyObject *res2 = PyList_New(sz);
           for(std::size_t i=0;i<sz;i++)
             {
               if(ms[i])
                 {
                   ms[i]->incrRef();
                   PyList_SetItem(res,i,SWIG_NewPointerObj(SWIG_as_voidptr(ms[i]),SWIGTYPE_p_MEDCoupling__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
                 }
               else
                 {
                   PyList_SetItem(res,i,SWIG_NewPointerObj(SWIG_as_voidptr(0),SWIGTYPE_p_MEDCoupling__DataArrayDouble, 0 ));
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

  class MEDCouplingFieldInt32 : public MEDCouplingFieldT<int>
  {
  public:
    static MEDCouplingFieldInt32 *New(TypeOfField type, TypeOfTimeDiscretization td=ONE_TIME);
    static MEDCouplingFieldInt32 *New(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td=ONE_TIME);
    bool isEqual(const MEDCouplingFieldInt32 *other, double meshPrec, int valsPrec) const;
    bool isEqualWithoutConsideringStr(const MEDCouplingFieldInt32 *other, double meshPrec, int valsPrec) const;
    void setTimeUnit(const std::string& unit);
    std::string getTimeUnit() const;
    void setTime(double val, int iteration, int order);
    void setArray(DataArrayInt32 *array);
    MEDCouplingFieldInt32 *deepCopy() const;
    MEDCouplingFieldInt32 *clone(bool recDeepCpy) const;
    MEDCouplingFieldInt32 *cloneWithMesh(bool recDeepCpy) const;
    MEDCouplingFieldDouble *convertToDblField() const;
    MEDCouplingFieldInt32 *buildSubPartRange(int begin, int end, int step) const;
    %extend {
      MEDCouplingFieldInt32(TypeOfField type, TypeOfTimeDiscretization td=ONE_TIME)
      {
        return MEDCouplingFieldInt32::New(type,td);
      }

      MEDCouplingFieldInt32(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td=ONE_TIME)
      {
        return MEDCouplingFieldInt32::New(ft,td);
      }

      PyObject *isEqualIfNotWhy(const MEDCouplingFieldInt32 *other, double meshPrec, int valsPrec) const
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

      std::string __str__() const
      {
        return self->simpleRepr();
      }

      std::string __repr__() const
      {
        std::ostringstream oss;
        self->reprQuickOverview(oss);
        return oss.str();
      }

      MEDCouplingFieldInt32 *buildSubPart(PyObject *li) const
      {
        return fieldT_buildSubPart(self,li);
      }

      MEDCouplingFieldInt32 *__getitem__(PyObject *li) const
      {
        return fieldT__getitem__(self,li);
      }

      DataArrayInt32 *getArray()
      {
        DataArrayInt32 *ret=self->getArray();
        if(ret)
          ret->incrRef();
        return ret;
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

      PyObject *getTinySerializationInformation() const
      {
        return field_getTinySerializationInformation<MEDCouplingFieldInt32>(self);
      }

      PyObject *serialize() const
      {
        return field_serialize<int>(self);
      }

      PyObject *__getstate__() const
      {
        return field__getstate__<MEDCouplingFieldInt32>(self,MEDCoupling_MEDCouplingFieldInt32_getTinySerializationInformation,MEDCoupling_MEDCouplingFieldInt32_serialize);
      }

      void __setstate__(PyObject *inp)
      {
        field__setstate__<int>(self,inp);
      }
    }
  };

  class MEDCouplingFieldInt64 : public MEDCouplingFieldT<int>
  {
  public:
    static MEDCouplingFieldInt64 *New(TypeOfField type, TypeOfTimeDiscretization td=ONE_TIME);
    static MEDCouplingFieldInt64 *New(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td=ONE_TIME);
    bool isEqual(const MEDCouplingFieldInt64 *other, double meshPrec, int valsPrec) const;
    bool isEqualWithoutConsideringStr(const MEDCouplingFieldInt64 *other, double meshPrec, int valsPrec) const;
    void setTimeUnit(const std::string& unit);
    std::string getTimeUnit() const;
    void setTime(double val, int iteration, int order);
    void setArray(DataArrayInt64 *array);
    MEDCouplingFieldInt64 *deepCopy() const;
    MEDCouplingFieldInt64 *clone(bool recDeepCpy) const;
    MEDCouplingFieldInt64 *cloneWithMesh(bool recDeepCpy) const;
    MEDCouplingFieldDouble *convertToDblField() const;
    MEDCouplingFieldInt64 *buildSubPartRange(int begin, int end, int step) const;
    %extend {
      MEDCouplingFieldInt64(TypeOfField type, TypeOfTimeDiscretization td=ONE_TIME)
      {
        return MEDCouplingFieldInt64::New(type,td);
      }

      MEDCouplingFieldInt64(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td=ONE_TIME)
      {
        return MEDCouplingFieldInt64::New(ft,td);
      }

      PyObject *isEqualIfNotWhy(const MEDCouplingFieldInt64 *other, double meshPrec, int valsPrec) const
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

      std::string __str__() const
      {
        return self->simpleRepr();
      }

      std::string __repr__() const
      {
        std::ostringstream oss;
        self->reprQuickOverview(oss);
        return oss.str();
      }

      MEDCouplingFieldInt64 *buildSubPart(PyObject *li) const
      {
        return fieldT_buildSubPart(self,li);
      }

      MEDCouplingFieldInt64 *__getitem__(PyObject *li) const
      {
        return fieldT__getitem__(self,li);
      }

      DataArrayInt64 *getArray()
      {
        DataArrayInt64 *ret=self->getArray();
        if(ret)
          ret->incrRef();
        return ret;
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

      PyObject *getTinySerializationInformation() const
      {
        return field_getTinySerializationInformation<MEDCouplingFieldInt64>(self);
      }

      PyObject *serialize() const
      {
        return field_serialize<Int64>(self);
      }

      PyObject *__getstate__() const
      {
        return field__getstate__<MEDCouplingFieldInt64>(self,MEDCoupling_MEDCouplingFieldInt64_getTinySerializationInformation,MEDCoupling_MEDCouplingFieldInt64_serialize);
      }

      void __setstate__(PyObject *inp)
      {
        field__setstate__<Int64>(self,inp);
      }
    }
  };

  class MEDCouplingFieldFloat : public MEDCouplingFieldT<float>
  {
  public:
    static MEDCouplingFieldFloat *New(TypeOfField type, TypeOfTimeDiscretization td=ONE_TIME);
    static MEDCouplingFieldFloat *New(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td=ONE_TIME);
    bool isEqual(const MEDCouplingFieldFloat *other, double meshPrec, float valsPrec) const;
    bool isEqualWithoutConsideringStr(const MEDCouplingFieldFloat *other, double meshPrec, float valsPrec) const;
    void setTimeUnit(const std::string& unit);
    std::string getTimeUnit() const;
    void setTime(double val, int iteration, int order);
    void setArray(DataArrayFloat *array);
    MEDCouplingFieldFloat *deepCopy() const;
    MEDCouplingFieldFloat *clone(bool recDeepCpy) const;
    MEDCouplingFieldFloat *cloneWithMesh(bool recDeepCpy) const;
    MEDCouplingFieldDouble *convertToDblField() const;
    MEDCouplingFieldFloat *buildSubPartRange(int begin, int end, int step) const;
    %extend {
      MEDCouplingFieldFloat(TypeOfField type, TypeOfTimeDiscretization td=ONE_TIME)
      {
        return MEDCouplingFieldFloat::New(type,td);
      }

      MEDCouplingFieldFloat(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td=ONE_TIME)
      {
        return MEDCouplingFieldFloat::New(ft,td);
      }

      PyObject *isEqualIfNotWhy(const MEDCouplingFieldFloat *other, double meshPrec, float valsPrec) const
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

      std::string __str__() const
      {
        return self->simpleRepr();
      }

      std::string __repr__() const
      {
        std::ostringstream oss;
        self->reprQuickOverview(oss);
        return oss.str();
      }

      MEDCouplingFieldFloat *buildSubPart(PyObject *li) const
      {
        return fieldT_buildSubPart(self,li);
      }

      MEDCouplingFieldFloat *__getitem__(PyObject *li) const
      {
        return fieldT__getitem__(self,li);
      }

      DataArrayFloat *getArray()
      {
        DataArrayFloat *ret=self->getArray();
        if(ret)
          ret->incrRef();
        return ret;
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

      PyObject *getTinySerializationInformation() const
      {
        return field_getTinySerializationInformation<MEDCouplingFieldFloat>(self);
      }

      PyObject *serialize() const
      {
        return field_serialize<float>(self);
      }

      PyObject *__getstate__() const
      {
        return field__getstate__<MEDCouplingFieldFloat>(self,MEDCoupling_MEDCouplingFieldFloat_getTinySerializationInformation,MEDCoupling_MEDCouplingFieldFloat_serialize);
      }

      void __setstate__(PyObject *inp)
      {
        field__setstate__<float>(self,inp);
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

        PyObject *getIdsOnTimeRight(double tm) const
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

        PyObject *getIdsOnTimeLeft(double tm) const
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
    double getTimeTolerance() const;
    MEDCouplingDefinitionTime getDefinitionTimeZone() const;

    %extend
      {
        MEDCouplingFieldOverTime(PyObject *li)
          {
            std::vector<const MEDCoupling::MEDCouplingFieldDouble *> tmp;
            convertFromPyObjVectorOfObj<const MEDCoupling::MEDCouplingFieldDouble *>(li,SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDouble,"MEDCouplingFieldDouble",tmp);
            std::size_t sz=tmp.size();
            std::vector<MEDCouplingFieldDouble *> fs(sz);
            for(std::size_t i=0;i<sz;i++)
              fs[i]=const_cast<MEDCouplingFieldDouble *>(tmp[i]);
            return MEDCouplingFieldOverTime::New(fs);
          }
        std::string __str__() const
          {
            return self->simpleRepr();
          }
        static MEDCouplingFieldOverTime *New(PyObject *li)
          {
            std::vector<const MEDCoupling::MEDCouplingFieldDouble *> tmp;
            convertFromPyObjVectorOfObj<const MEDCoupling::MEDCouplingFieldDouble *>(li,SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDouble,"MEDCouplingFieldDouble",tmp);
            std::size_t sz=tmp.size();
            std::vector<MEDCouplingFieldDouble *> fs(sz);
            for(std::size_t i=0;i<sz;i++)
              fs[i]=const_cast<MEDCouplingFieldDouble *>(tmp[i]);
            return MEDCouplingFieldOverTime::New(fs);
          }
      }
  };

  class MEDCouplingCartesianAMRMesh;

  class MEDCouplingCartesianAMRPatchGen : public RefCountObject
  {
  public:
    int getNumberOfCellsRecursiveWithOverlap() const;
    int getNumberOfCellsRecursiveWithoutOverlap() const;
    int getMaxNumberOfLevelsRelativeToThis() const;
    %extend
    {
      MEDCouplingCartesianAMRMeshGen *getMesh() const
      {
        MEDCouplingCartesianAMRMeshGen *ret(const_cast<MEDCouplingCartesianAMRMeshGen *>(self->getMesh()));
        if(ret)
          ret->incrRef();
        return ret;
      }
    }
  };

  class MEDCouplingCartesianAMRPatch : public MEDCouplingCartesianAMRPatchGen
  {
  public:
    int getNumberOfOverlapedCellsForFather() const;
    bool isInMyNeighborhood(const MEDCouplingCartesianAMRPatch *other, int ghostLev) const;
    std::vector<mcIdType> computeCellGridSt() const;
    %extend
    {
      PyObject *getBLTRRange() const
      {
        const std::vector< std::pair<mcIdType,mcIdType> >& ret(self->getBLTRRange());
        return convertFromVectorPairInt(ret);
      }

      PyObject *getBLTRRangeRelativeToGF() const
      {
        std::vector< std::pair<mcIdType,mcIdType> > ret(self->getBLTRRangeRelativeToGF());
        return convertFromVectorPairInt(ret);
      }

      void addPatch(PyObject *bottomLeftTopRight, const std::vector<mcIdType>& factors)
      {
        std::vector< std::pair<mcIdType,mcIdType> > inp;
        convertPyToVectorPairInt(bottomLeftTopRight,inp);
        self->addPatch(inp,factors);
      }

      MEDCouplingCartesianAMRPatch *__getitem__(mcIdType patchId) const
      {
        const MEDCouplingCartesianAMRMeshGen *mesh(self->getMesh());
        if(!mesh)
          throw INTERP_KERNEL::Exception("wrap MEDCouplingCartesianAMRPatchGen.__getitem__ : no underlying mesh !");
        if(patchId==mesh->getNumberOfPatches())
          {
            std::ostringstream oss;
            oss << "Requesting for patchId " << patchId << " having only " << mesh->getNumberOfPatches() << " patches !";
            PyErr_SetString(PyExc_StopIteration,oss.str().c_str());
            return 0;
          }
        MEDCouplingCartesianAMRPatch *ret(const_cast<MEDCouplingCartesianAMRPatch *>(mesh->getPatch(patchId)));
        if(ret)
          ret->incrRef();
        return ret;
      }

      void __delitem__(mcIdType patchId)
      {
        MEDCouplingCartesianAMRMeshGen *mesh(const_cast<MEDCouplingCartesianAMRMeshGen *>(self->getMesh()));
        if(!mesh)
          throw INTERP_KERNEL::Exception("wrap MEDCouplingCartesianAMRPatch.__delitem__ : no underlying mesh !");
        mesh->removePatch(patchId);
      }

      mcIdType __len__() const
      {
        const MEDCouplingCartesianAMRMeshGen *mesh(self->getMesh());
        if(!mesh)
          throw INTERP_KERNEL::Exception("wrap MEDCouplingCartesianAMRPatch.__len__ : no underlying mesh !");
        return mesh->getNumberOfPatches();
      }
    }
  };

  class MEDCouplingCartesianAMRPatchGF : public MEDCouplingCartesianAMRPatchGen
  {
  };

  class MEDCouplingCartesianAMRMeshGen : public RefCountObject, public TimeLabel
  {
  public:
    mcIdType getAbsoluteLevel() const;
    mcIdType getAbsoluteLevelRelativeTo(const MEDCouplingCartesianAMRMeshGen *ref) const;
    std::vector<mcIdType> getPositionRelativeTo(const MEDCouplingCartesianAMRMeshGen *ref) const;
    int getSpaceDimension() const;
    const std::vector<mcIdType>& getFactors() const;
    void setFactors(const std::vector<mcIdType>& newFactors);
    mcIdType getMaxNumberOfLevelsRelativeToThis() const;
    mcIdType getNumberOfCellsAtCurrentLevel() const;
    mcIdType getNumberOfCellsAtCurrentLevelGhost(mcIdType ghostLev) const;
    mcIdType getNumberOfCellsRecursiveWithOverlap() const;
    mcIdType getNumberOfCellsRecursiveWithoutOverlap() const;
    bool isPatchInNeighborhoodOf(mcIdType patchId1, mcIdType patchId2, mcIdType ghostLev) const;
   virtual void detachFromFather();
    //
    mcIdType getNumberOfPatches() const;
    mcIdType getPatchIdFromChildMesh(const MEDCouplingCartesianAMRMeshGen *mesh) const;
    MEDCouplingUMesh *buildUnstructured() const;
    DataArrayDouble *extractGhostFrom(mcIdType ghostSz, const DataArrayDouble *arr) const;
    std::vector<mcIdType> getPatchIdsInTheNeighborhoodOf(mcIdType patchId, mcIdType ghostLev) const;
    MEDCoupling1SGTUMesh *buildMeshFromPatchEnvelop() const;
    MEDCoupling1SGTUMesh *buildMeshOfDirectChildrenOnly() const;
    void removeAllPatches();
    void removePatch(mcIdType patchId);
    void createPatchesFromCriterion(const INTERP_KERNEL::BoxSplittingOptions& bso, const DataArrayByte *criterion, const std::vector<mcIdType>& factors);
    void createPatchesFromCriterion(const INTERP_KERNEL::BoxSplittingOptions& bso, const DataArrayDouble *criterion, const std::vector<mcIdType>& factors, double eps);
    DataArrayDouble *createCellFieldOnPatch(mcIdType patchId, const DataArrayDouble *cellFieldOnThis) const;
    void fillCellFieldOnPatch(mcIdType patchId, const DataArrayDouble *cellFieldOnThis, DataArrayDouble *cellFieldOnPatch, bool isConservative=true) const;
    void fillCellFieldOnPatchGhost(mcIdType patchId, const DataArrayDouble *cellFieldOnThis, DataArrayDouble *cellFieldOnPatch, mcIdType ghostLev, bool isConservative=true) const;
    void fillCellFieldOnPatchOnlyOnGhostZone(mcIdType patchId, const DataArrayDouble *cellFieldOnThis, DataArrayDouble *cellFieldOnPatch, mcIdType ghostLev) const;
    void fillCellFieldOnPatchOnlyOnGhostZoneWith(mcIdType ghostLev, const MEDCouplingCartesianAMRPatch *patchToBeModified, const MEDCouplingCartesianAMRPatch *neighborPatch, DataArrayDouble *cellFieldOnPatch, const DataArrayDouble *cellFieldNeighbor) const;
    void fillCellFieldComingFromPatch(mcIdType patchId, const DataArrayDouble *cellFieldOnPatch, DataArrayDouble *cellFieldOnThis, bool isConservative=true) const;
    void fillCellFieldComingFromPatchGhost(mcIdType patchId, const DataArrayDouble *cellFieldOnPatch, DataArrayDouble *cellFieldOnThis, mcIdType ghostLev, bool isConservative=true) const;
    DataArrayIdType *findPatchesInTheNeighborhoodOf(mcIdType patchId, mcIdType ghostLev) const;
    std::string buildPythonDumpOfThis() const;
    %extend
    {
      void addPatch(PyObject *bottomLeftTopRight, const std::vector<mcIdType>& factors)
      {
        std::vector< std::pair<mcIdType,mcIdType> > inp;
        convertPyToVectorPairInt(bottomLeftTopRight,inp);
        self->addPatch(inp,factors);
      }

      PyObject *getPatches() const
      {
        std::vector< const MEDCouplingCartesianAMRPatch *> ps(self->getPatches());
        std::size_t sz(ps.size());
        PyObject *ret = PyList_New(sz);
        for(std::size_t i=0;i<sz;i++)
          {
            MEDCouplingCartesianAMRPatch *elt(const_cast<MEDCouplingCartesianAMRPatch *>(ps[i]));
            if(elt)
              elt->incrRef();
            PyList_SetItem(ret,i,convertCartesianAMRPatch(elt, SWIG_POINTER_OWN | 0 ));
          }
        return ret;
      }

      // agy : don't know why typemap fails here ??? let it in the extend section
      PyObject *deepCopy(MEDCouplingCartesianAMRMeshGen *father) const
      {
        return convertCartesianAMRMesh(self->deepCopy(father), SWIG_POINTER_OWN | 0 );
      }

      MEDCouplingCartesianAMRPatch *getPatchAtPosition(const std::vector<mcIdType>& pos) const
      {
        const MEDCouplingCartesianAMRPatch *ret(self->getPatchAtPosition(pos));
        MEDCouplingCartesianAMRPatch *ret2(const_cast<MEDCouplingCartesianAMRPatch *>(ret));
        if(ret2)
          ret2->incrRef();
        return ret2;
      }

      MEDCouplingCartesianAMRMeshGen *getMeshAtPosition(const std::vector<mcIdType>& pos) const
      {
        const MEDCouplingCartesianAMRMeshGen *ret(self->getMeshAtPosition(pos));
        MEDCouplingCartesianAMRMeshGen *ret2(const_cast<MEDCouplingCartesianAMRMeshGen *>(ret));
        if(ret2)
          ret2->incrRef();
        return ret2;
      }

      virtual PyObject *positionRelativeToGodFather() const
      {
        std::vector<mcIdType> out1;
        std::vector< std::pair<mcIdType,mcIdType> > out0(self->positionRelativeToGodFather(out1));
        PyObject *ret(PyTuple_New(2));
        PyTuple_SetItem(ret,0,convertFromVectorPairInt(out0));
        PyTuple_SetItem(ret,1,convertIntArrToPyList2(out1));
        return ret;
      }

      virtual PyObject *retrieveGridsAt(mcIdType absoluteLev) const
      {
        std::vector<MEDCouplingCartesianAMRPatchGen *> ps(self->retrieveGridsAt(absoluteLev));
        std::size_t sz(ps.size());
        PyObject *ret = PyList_New(sz);
        for(std::size_t i=0;i<sz;i++)
          PyList_SetItem(ret,i,convertCartesianAMRPatch(ps[i], SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      MEDCouplingFieldDouble *buildCellFieldOnRecurseWithoutOverlapWithoutGhost(mcIdType ghostSz, PyObject *recurseArrs) const
      {
        std::vector<const DataArrayDouble *> inp;
        convertFromPyObjVectorOfObj<const MEDCoupling::DataArrayDouble *>(recurseArrs,SWIGTYPE_p_MEDCoupling__DataArrayDouble,"DataArrayDouble",inp);
        return self->buildCellFieldOnRecurseWithoutOverlapWithoutGhost(ghostSz,inp);
      }

      virtual MEDCouplingCartesianAMRMeshGen *getFather() const
      {
        MEDCouplingCartesianAMRMeshGen *ret(const_cast<MEDCouplingCartesianAMRMeshGen *>(self->getFather()));
        if(ret)
          ret->incrRef();
        return ret;
      }

      virtual MEDCouplingCartesianAMRMeshGen *getGodFather() const
      {
        MEDCouplingCartesianAMRMeshGen *ret(const_cast<MEDCouplingCartesianAMRMeshGen *>(self->getGodFather()));
        if(ret)
          ret->incrRef();
        return ret;
      }

      MEDCouplingCartesianAMRPatch *getPatch(mcIdType patchId) const
      {
        MEDCouplingCartesianAMRPatch *ret(const_cast<MEDCouplingCartesianAMRPatch *>(self->getPatch(patchId)));
        if(ret)
          ret->incrRef();
        return ret;
      }

      MEDCouplingIMesh *getImageMesh() const
      {
        const MEDCouplingIMesh *ret(self->getImageMesh());
        if(ret)
          ret->incrRef();
        return const_cast<MEDCouplingIMesh *>(ret);
      }

      MEDCouplingCartesianAMRPatch *__getitem__(mcIdType patchId) const
      {
        if(patchId==self->getNumberOfPatches())
          {
            std::ostringstream oss;
            oss << "Requesting for patchId " << patchId << " having only " << self->getNumberOfPatches() << " patches !";
            PyErr_SetString(PyExc_StopIteration,oss.str().c_str());
            return 0;
          }
        MEDCouplingCartesianAMRPatch *ret(const_cast<MEDCouplingCartesianAMRPatch *>(self->getPatch(patchId)));
        if(ret)
          ret->incrRef();
        return ret;
      }

      void fillCellFieldOnPatchGhostAdv(mcIdType patchId, const DataArrayDouble *cellFieldOnThis, mcIdType ghostLev, PyObject *arrsOnPatches, bool isConservative=true) const
      {
        std::vector<const MEDCoupling::DataArrayDouble *> arrsOnPatches2;
        convertFromPyObjVectorOfObj<const MEDCoupling::DataArrayDouble *>(arrsOnPatches,SWIGTYPE_p_MEDCoupling__DataArrayDouble,"DataArrayDouble",arrsOnPatches2);
        self->fillCellFieldOnPatchGhostAdv(patchId,cellFieldOnThis,ghostLev,arrsOnPatches2,isConservative);
      }

      void fillCellFieldOnPatchOnlyGhostAdv(mcIdType patchId, mcIdType ghostLev, PyObject *arrsOnPatches) const
      {
        std::vector<const MEDCoupling::DataArrayDouble *> arrsOnPatches2;
        convertFromPyObjVectorOfObj<const MEDCoupling::DataArrayDouble *>(arrsOnPatches,SWIGTYPE_p_MEDCoupling__DataArrayDouble,"DataArrayDouble",arrsOnPatches2);
        self->fillCellFieldOnPatchOnlyGhostAdv(patchId,ghostLev,arrsOnPatches2);
      }

      void __delitem__(mcIdType patchId)
      {
        self->removePatch(patchId);
      }

      mcIdType __len__() const
      {
        return self->getNumberOfPatches();
      }
    }
  };

  class MEDCouplingCartesianAMRMeshSub : public MEDCouplingCartesianAMRMeshGen
  {
  };

  class MEDCouplingCartesianAMRMesh : public MEDCouplingCartesianAMRMeshGen
  {
  public:
    static MEDCouplingCartesianAMRMesh *New(MEDCouplingIMesh *mesh);
    %extend
    {
      static MEDCouplingCartesianAMRMesh *New(const std::string& meshName, int spaceDim, PyObject *nodeStrct, PyObject *origin, PyObject *dxyz)
      {
        static const char msg0[]="MEDCouplingCartesianAMRMesh::New : error on 'origin' parameter !";
        static const char msg1[]="MEDCouplingCartesianAMRMesh::New : error on 'dxyz' parameter !";
        const mcIdType *nodeStrctPtr(0);
        const double *originPtr(0),*dxyzPtr(0);
        mcIdType sw,sz,val0;
        std::vector<mcIdType> bb0;
        nodeStrctPtr=convertIntStarLikePyObjToCppIntStar(nodeStrct,sw,sz,val0,bb0);
        //
        double val,val2;
        std::vector<double> bb,bb2;
        mcIdType sz1,sz2;
        originPtr=convertObjToPossibleCpp5_SingleCompo(origin,sw,val,bb,msg0,false,sz1);
        dxyzPtr=convertObjToPossibleCpp5_SingleCompo(dxyz,sw,val2,bb2,msg1,false,sz2);
        //
        return MEDCouplingCartesianAMRMesh::New(meshName,spaceDim,nodeStrctPtr,nodeStrctPtr+sz,originPtr,originPtr+sz1,dxyzPtr,dxyzPtr+sz2);
      }

      void createPatchesFromCriterionML(PyObject *bso, const DataArrayDouble *criterion, PyObject *factors, double eps)
      {
        std::vector<const INTERP_KERNEL::BoxSplittingOptions *> inp0;
        convertFromPyObjVectorOfObj<const INTERP_KERNEL::BoxSplittingOptions *>(bso,SWIGTYPE_p_INTERP_KERNEL__BoxSplittingOptions,"BoxSplittingOptions",inp0);
        std::vector< std::vector<mcIdType> > inp2;
        convertPyToVectorOfVectorOfInt(factors,inp2);
        self->createPatchesFromCriterionML(inp0,criterion,inp2,eps);
      }

      MEDCouplingCartesianAMRMesh(const std::string& meshName, int spaceDim, PyObject *nodeStrct, PyObject *origin, PyObject *dxyz)
      {
        return MEDCoupling_MEDCouplingCartesianAMRMesh_New__SWIG_1(meshName,spaceDim,nodeStrct,origin,dxyz);
      }

      MEDCouplingCartesianAMRMesh(MEDCouplingIMesh *mesh)
      {
        return MEDCouplingCartesianAMRMesh::New(mesh);
      }
    }
  };

  class MEDCouplingDataForGodFather : public RefCountObject
  {
  public:
    virtual void synchronizeFineToCoarse();
    virtual void synchronizeFineToCoarseBetween(mcIdType fromLev, mcIdType toLev);
    virtual void synchronizeCoarseToFine();
    virtual void synchronizeCoarseToFineBetween(mcIdType fromLev, mcIdType toLev);
    virtual void synchronizeAllGhostZones();
    virtual void synchronizeAllGhostZonesOfDirectChidrenOf(const MEDCouplingCartesianAMRMeshGen *mesh);
    virtual void synchronizeAllGhostZonesAtASpecifiedLevel(mcIdType level);
    virtual void synchronizeAllGhostZonesAtASpecifiedLevelUsingOnlyFather(mcIdType level);
    virtual void alloc();
    virtual void dealloc();
    %extend
    {
      MEDCouplingCartesianAMRMesh *getMyGodFather()
      {
        MEDCouplingCartesianAMRMesh *ret(self->getMyGodFather());
        if(ret)
          ret->incrRef();
        return ret;
      }
    }
  };

  class MEDCouplingAMRAttribute : public MEDCouplingDataForGodFather, public TimeLabel
  {
  public:
    mcIdType getNumberOfLevels() const;
    MEDCouplingAMRAttribute *deepCopy() const;
    MEDCouplingAMRAttribute *deepCpyWithoutGodFather() const;
    MEDCouplingFieldDouble *buildCellFieldOnRecurseWithoutOverlapWithoutGhost(MEDCouplingCartesianAMRMeshGen *mesh, const std::string& fieldName) const;
    MEDCouplingFieldDouble *buildCellFieldOnWithGhost(MEDCouplingCartesianAMRMeshGen *mesh, const std::string& fieldName) const;
    MEDCouplingFieldDouble *buildCellFieldOnWithoutGhost(MEDCouplingCartesianAMRMeshGen *mesh, const std::string& fieldName) const;
    bool changeGodFather(MEDCouplingCartesianAMRMesh *gf);
    MEDCouplingAMRAttribute *projectTo(MEDCouplingCartesianAMRMesh *targetGF) const;
    std::string writeVTHB(const std::string& fileName) const;
    %extend
    {
      static MEDCouplingAMRAttribute *New(MEDCouplingCartesianAMRMesh *gf, PyObject *fieldNames, mcIdType ghostLev)
      {
        std::vector< std::pair<std::string,int> > fieldNamesCpp0;
        std::vector< std::pair<std::string, std::vector<std::string> > > fieldNamesCpp1;
        MEDCouplingAMRAttribute *ret(0);
        try
          {
            convertPyToVectorPairStringInt(fieldNames,fieldNamesCpp0);
            ret=MEDCouplingAMRAttribute::New(gf,fieldNamesCpp0,ghostLev);
          }
        catch(INTERP_KERNEL::Exception&)
          {
            convertPyToVectorPairStringVecString(fieldNames,fieldNamesCpp1);
            ret=MEDCouplingAMRAttribute::New(gf,fieldNamesCpp1,ghostLev);
          }
        return ret;
      }

      MEDCouplingAMRAttribute(MEDCouplingCartesianAMRMesh *gf, PyObject *fieldNames, mcIdType ghostLev)
      {
        return MEDCoupling_MEDCouplingAMRAttribute_New(gf,fieldNames,ghostLev);
      }

      DataArrayDouble *getFieldOn(MEDCouplingCartesianAMRMeshGen *mesh, const std::string& fieldName) const
      {
        const DataArrayDouble *ret(self->getFieldOn(mesh,fieldName));
        DataArrayDouble *ret2(const_cast<DataArrayDouble *>(ret));
        if(ret2)
          ret2->incrRef();
        return ret2;
      }

      void spillInfoOnComponents(PyObject *compNames)
      {
        std::vector< std::vector<std::string> > compNamesCpp;
        convertPyToVectorOfVectorOfString(compNames,compNamesCpp);
        self->spillInfoOnComponents(compNamesCpp);
      }

      void spillNatures(PyObject *nfs)
      {
        std::vector<mcIdType> inp0;
        if(!fillIntVector(nfs,inp0))
          throw INTERP_KERNEL::Exception("wrap of MEDCouplingAMRAttribute::spillNatures : vector of NatureOfField enum expected !");
        std::size_t sz(inp0.size());
        std::vector<NatureOfField> inp00(sz);
        for(std::size_t i=0;i<sz;i++)
          inp00[i]=(NatureOfField)inp0[i];
        self->spillNatures(inp00);
      }

      PyObject *retrieveFieldsOn(MEDCouplingCartesianAMRMeshGen *mesh) const
      {
        std::vector<DataArrayDouble *> ret(self->retrieveFieldsOn(mesh));
        std::size_t sz(ret.size());
        PyObject *retPy(PyList_New(sz));
        for(std::size_t i=0;i<sz;i++)
          PyList_SetItem(retPy,i,SWIG_NewPointerObj(SWIG_as_voidptr(ret[i]),SWIGTYPE_p_MEDCoupling__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
        return retPy;
      }
    }
  };

  class DenseMatrix : public RefCountObject, public TimeLabel
  {
  public:
    static DenseMatrix *New(mcIdType nbRows, mcIdType nbCols);
    static DenseMatrix *New(DataArrayDouble *array, mcIdType nbRows, mcIdType nbCols);
    DenseMatrix *deepCopy() const;
    DenseMatrix *shallowCpy() const;
    //
    mcIdType getNumberOfRows() const;
    mcIdType getNumberOfCols() const;
    mcIdType getNbOfElems() const;
    void reBuild(DataArrayDouble *array, mcIdType nbRows=-1, mcIdType nbCols=-1);
    void reShape(mcIdType nbRows, mcIdType nbCols);
    void transpose();
    //
    bool isEqual(const DenseMatrix& other, double eps) const;
    DataArrayDouble *matVecMult(const DataArrayDouble *vec) const;
    static DataArrayDouble *MatVecMult(const DenseMatrix *mat, const DataArrayDouble *vec);
    %extend
    {
      DenseMatrix(mcIdType nbRows, mcIdType nbCols)
      {
        return DenseMatrix::New(nbRows,nbCols);
      }

      DenseMatrix(DataArrayDouble *array, mcIdType nbRows, mcIdType nbCols)
      {
        return DenseMatrix::New(array,nbRows,nbCols);
      }

      PyObject *isEqualIfNotWhy(const DenseMatrix& other, double eps) const
      {
        std::string ret1;
        bool ret0=self->isEqualIfNotWhy(other,eps,ret1);
        PyObject *ret=PyTuple_New(2);
        PyObject *ret0Py=ret0?Py_True:Py_False;
        Py_XINCREF(ret0Py);
        PyTuple_SetItem(ret,0,ret0Py);
        PyTuple_SetItem(ret,1,PyString_FromString(ret1.c_str()));
        return ret;
      }

      DataArrayDouble *getData()
      {
        DataArrayDouble *ret(self->getData());
        if(ret)
          ret->incrRef();
        return ret;
      }

      DenseMatrix *__add__(const DenseMatrix *other)
      {
        return MEDCoupling::DenseMatrix::Add(self,other);
      }

      DenseMatrix *__sub__(const DenseMatrix *other)
      {
        return MEDCoupling::DenseMatrix::Substract(self,other);
      }

      DenseMatrix *__mul__(const DenseMatrix *other)
      {
        return MEDCoupling::DenseMatrix::Multiply(self,other);
      }

      DenseMatrix *__mul__(const DataArrayDouble *other)
      {
        return MEDCoupling::DenseMatrix::Multiply(self,other);
      }

      PyObject *___iadd___(PyObject *trueSelf, const DenseMatrix *other)
      {
        self->addEqual(other);
        Py_XINCREF(trueSelf);
        return trueSelf;
      }

      PyObject *___isub___(PyObject *trueSelf, const DenseMatrix *other)
      {
        self->substractEqual(other);
        Py_XINCREF(trueSelf);
        return trueSelf;
      }
#ifdef WITH_NUMPY
      PyObject *toNumPyMatrix() // not const. It is not a bug !
      {
        PyObject *obj(ToNumPyArrayUnderground<DataArrayDouble,double>(self->getData(),NPY_DOUBLE,"DataArrayDouble",self->getNumberOfRows(),self->getNumberOfCols()));
        return obj;
      }
#endif
    }
  };
}

%pythoncode %{
def MEDCouplingUMeshReduce(self):
    return MEDCouplingStdReduceFunct,(MEDCouplingUMesh,((),(self.__getstate__()),))
def MEDCouplingCMeshReduce(self):
    return MEDCouplingStdReduceFunct,(MEDCouplingCMesh,((),(self.__getstate__()),))
def MEDCouplingIMeshReduce(self):
    return MEDCouplingStdReduceFunct,(MEDCouplingIMesh,((),(self.__getstate__()),))
def MEDCouplingMappedExtrudedMeshReduce(self):
    return MEDCouplingStdReduceFunct,(MEDCouplingMappedExtrudedMesh,((),(self.__getstate__()),))
def MEDCouplingCurveLinearMeshReduce(self):
    return MEDCouplingStdReduceFunct,(MEDCouplingCurveLinearMesh,((),(self.__getstate__()),))
def MEDCoupling1SGTUMeshReduce(self):
    return MEDCouplingStdReduceFunct,(MEDCoupling1SGTUMesh,((),(self.__getstate__()),))
def MEDCoupling1DGTUMeshReduce(self):
    return MEDCouplingStdReduceFunct,(MEDCoupling1DGTUMesh,((),(self.__getstate__()),))
def MEDCouplingFieldDoubleReduce(self):
    self.checkConsistencyLight()
    d=(self.getTypeOfField(),self.getTimeDiscretization())
    return MEDCouplingStdReduceFunct,(MEDCouplingFieldDouble,(d,(self.__getstate__()),))
def MEDCouplingFieldInt32Reduce(self):
    self.checkConsistencyLight()
    d=(self.getTypeOfField(),self.getTimeDiscretization())
    return MEDCouplingStdReduceFunct,(MEDCouplingFieldInt32,(d,(self.__getstate__()),))
def MEDCouplingFieldInt64Reduce(self):
    self.checkConsistencyLight()
    d=(self.getTypeOfField(),self.getTimeDiscretization())
    return MEDCouplingStdReduceFunct,(MEDCouplingFieldInt64,(d,(self.__getstate__()),))
def MEDCouplingFieldFloatReduce(self):
    self.checkConsistencyLight()
    d=(self.getTypeOfField(),self.getTimeDiscretization())
    return MEDCouplingStdReduceFunct,(MEDCouplingFieldFloat,(d,(self.__getstate__()),))
def MEDCouplingFTReduceFunct(cls,params):
    a,b=params
    ret=object.__new__(cls)
    ret.__init__(*a)
    return ret

def MEDCouplingFieldTemplateReduce(self):
    ret = MEDCouplingFieldDouble(self)
    nbTuples = self.getNumberOfTuplesExpected()
    arr = DataArrayDouble(nbTuples) ; arr[:] = 0.
    ret.setArray(arr)
    return MEDCouplingFTReduceFunct,(MEDCouplingFieldTemplate,((ret,),()))
#
# Forwarding DataArrayInt functions to MEDCouplingUMesh:
#
MEDCouplingUMesh.ExtractFromIndexedArrays           = DataArrayInt.ExtractFromIndexedArrays
MEDCouplingUMesh.ExtractFromIndexedArraysSlice      = DataArrayInt.ExtractFromIndexedArraysSlice
MEDCouplingUMesh.SetPartOfIndexedArrays             = DataArrayInt.SetPartOfIndexedArrays
MEDCouplingUMesh.SetPartOfIndexedArraysSameIdx      = DataArrayInt.SetPartOfIndexedArraysSameIdx
MEDCouplingUMesh.RemoveIdsFromIndexedArrays         = DataArrayInt.RemoveIdsFromIndexedArrays
MEDCouplingFieldInt = MEDCouplingFieldInt32

if MEDCouplingUse64BitIDs():
  MEDCouplingFieldID = MEDCouplingFieldInt64
else:
  MEDCouplingFieldID = MEDCouplingFieldInt32

%}

%pythoncode %{
import os
__filename=os.environ.get('PYTHONSTARTUP')
if __filename and os.path.isfile(__filename):
    with open(__filename) as __fp:
        exec(__fp.read())
%}
