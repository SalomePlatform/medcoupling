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
// Author : Anthony Geay (CEA/DEN)

%module MEDCoupling

#ifdef WITH_DOCSTRINGS
%include MEDCoupling_doc.i
#endif

%include std_vector.i
%include std_string.i

%{
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingMappedExtrudedMesh.hxx"
#include "MEDCouplingCMesh.hxx"
#include "MEDCouplingIMesh.hxx"
#include "MEDCouplingCurveLinearMesh.hxx"
#include "MEDCoupling1GTUMesh.hxx"
#include "MEDCouplingField.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingFieldInt.hxx"
#include "MEDCouplingFieldFloat.hxx"
#include "MEDCouplingFieldTemplate.hxx"
#include "MEDCouplingGaussLocalization.hxx"
#include "MCAuto.hxx"
#include "MEDCouplingMultiFields.hxx"
#include "MEDCouplingFieldOverTime.hxx"
#include "MEDCouplingDefinitionTime.hxx"
#include "MEDCouplingFieldDiscretization.hxx"
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

%template(ivec) std::vector<int>;
%template(dvec) std::vector<double>;
%template(svec) std::vector<std::string>;

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
%newobject MEDCoupling::MEDCouplingFieldInt::New;
%newobject MEDCoupling::MEDCouplingFieldInt::convertToDblField;
%newobject MEDCoupling::MEDCouplingFieldInt::getArray;
%newobject MEDCoupling::MEDCouplingFieldInt::deepCopy;
%newobject MEDCoupling::MEDCouplingFieldInt::clone;
%newobject MEDCoupling::MEDCouplingFieldInt::cloneWithMesh;
%newobject MEDCoupling::MEDCouplingFieldInt::buildSubPart;
%newobject MEDCoupling::MEDCouplingFieldInt::buildSubPartRange;
%newobject MEDCoupling::MEDCouplingFieldInt::__getitem__;
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
%newobject MEDCoupling::MEDCouplingUMesh::rearrange2ConsecutiveCellTypes;
%newobject MEDCoupling::MEDCouplingUMesh::sortCellsInMEDFileFrmt;
%newobject MEDCoupling::MEDCouplingUMesh::getRenumArrForMEDFileFrmt;
%newobject MEDCoupling::MEDCouplingUMesh::convertCellArrayPerGeoType;
%newobject MEDCoupling::MEDCouplingUMesh::getRenumArrForConsecutiveCellTypesSpec;
%newobject MEDCoupling::MEDCouplingUMesh::buildDirectionVectorField;
%newobject MEDCoupling::MEDCouplingUMesh::convertLinearCellsToQuadratic;
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
%newobject MEDCoupling::MEDCouplingUMeshCellByTypeEntry::__iter__;
%newobject MEDCoupling::MEDCouplingUMeshCellEntry::__iter__;
%newobject MEDCoupling::MEDCoupling1GTUMesh::New;
%newobject MEDCoupling::MEDCoupling1GTUMesh::getNodalConnectivity;
%newobject MEDCoupling::MEDCoupling1GTUMesh::AggregateOnSameCoordsToUMesh;
%newobject MEDCoupling::MEDCoupling1SGTUMesh::New;
%newobject MEDCoupling::MEDCoupling1SGTUMesh::buildSetInstanceFromThis;
%newobject MEDCoupling::MEDCoupling1SGTUMesh::computeDualMesh;
%newobject MEDCoupling::MEDCoupling1SGTUMesh::explodeEachHexa8To6Quad4;
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
%newobject MEDCoupling::MEDCouplingSkyLineArray::BuildFromPolyhedronConn;
%newobject MEDCoupling::MEDCouplingSkyLineArray::getSuperIndexArray;
%newobject MEDCoupling::MEDCouplingSkyLineArray::getIndexArray;
%newobject MEDCoupling::MEDCouplingSkyLineArray::getValuesArray;

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
%feature("unref") MEDCouplingFieldDouble "$this->decrRef();"
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

%rename (InterpKernelException) INTERP_KERNEL::Exception;

%include "MEDCouplingRefCountObject.i"
%include "MEDCouplingMemArray.i"

%{
  void initializeMe()
  {// AGY : here initialization of C++ traits in MEDCouplingDataArrayTypemaps.i for code factorization. Awful, I know, but no other solutions.
    SWIGTITraits<double>::TI=SWIGTYPE_p_MEDCoupling__DataArrayDouble;
    SWIGTITraits<float>::TI=SWIGTYPE_p_MEDCoupling__DataArrayFloat;
    SWIGTITraits<int>::TI=SWIGTYPE_p_MEDCoupling__DataArrayInt;
    SWIGTITraits<double>::TI_TUPLE=SWIGTYPE_p_MEDCoupling__DataArrayDoubleTuple;
    SWIGTITraits<float>::TI_TUPLE=SWIGTYPE_p_MEDCoupling__DataArrayFloatTuple;
    SWIGTITraits<int>::TI_TUPLE=SWIGTYPE_p_MEDCoupling__DataArrayIntTuple;
  }
%}

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
    void init() throw(INTERP_KERNEL::Exception);
    double getEfficiencyGoal() const throw(INTERP_KERNEL::Exception);
    void setEfficiencyGoal(double efficiency) throw(INTERP_KERNEL::Exception);
    double getEfficiencyThreshold() const throw(INTERP_KERNEL::Exception);
    void setEfficiencyThreshold(double efficiencyThreshold) throw(INTERP_KERNEL::Exception);
    int getMinimumPatchLength() const throw(INTERP_KERNEL::Exception);
    void setMinimumPatchLength(int minPatchLength) throw(INTERP_KERNEL::Exception);
    int getMaximumPatchLength() const throw(INTERP_KERNEL::Exception);
    void setMaximumPatchLength(int maxPatchLength) throw(INTERP_KERNEL::Exception);
    int getMaximumNbOfCellsInPatch() const throw(INTERP_KERNEL::Exception);
    void setMaximumNbOfCellsInPatch(int maxNbCellsInPatch) throw(INTERP_KERNEL::Exception);
    void copyOptions(const BoxSplittingOptions & other) throw(INTERP_KERNEL::Exception);
    std::string printOptions() const throw(INTERP_KERNEL::Exception);
    %extend
    {
      std::string __str__() const throw(INTERP_KERNEL::Exception)
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
      ON_NODES_KR = 4
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

  class DataArrayInt;
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
    virtual MEDCouplingMeshType getType() const throw(INTERP_KERNEL::Exception);
    bool isStructured() const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingMesh *deepCopy() const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingMesh *clone(bool recDeepCpy) const throw(INTERP_KERNEL::Exception);
    virtual bool isEqual(const MEDCouplingMesh *other, double prec) const throw(INTERP_KERNEL::Exception);
    virtual bool isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const throw(INTERP_KERNEL::Exception);
    virtual void checkFastEquivalWith(const MEDCouplingMesh *other, double prec) const throw(INTERP_KERNEL::Exception);
    virtual void copyTinyStringsFrom(const MEDCouplingMesh *other) throw(INTERP_KERNEL::Exception);
    virtual void copyTinyInfoFrom(const MEDCouplingMesh *other) throw(INTERP_KERNEL::Exception);
    virtual void checkConsistencyLight() const throw(INTERP_KERNEL::Exception);
    virtual void checkConsistency(double eps=1e-12) const throw(INTERP_KERNEL::Exception);
    virtual int getNumberOfCells() const throw(INTERP_KERNEL::Exception);
    virtual int getNumberOfNodes() const throw(INTERP_KERNEL::Exception);
    virtual int getSpaceDimension() const throw(INTERP_KERNEL::Exception);
    virtual int getMeshDimension() const throw(INTERP_KERNEL::Exception);
    virtual DataArrayDouble *getCoordinatesAndOwner() const throw(INTERP_KERNEL::Exception);
    virtual DataArrayDouble *computeCellCenterOfMass() const throw(INTERP_KERNEL::Exception);
    virtual DataArrayDouble *computeIsoBarycenterOfNodesPerCell() const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *giveCellsWithType(INTERP_KERNEL::NormalizedCellType type) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *computeNbOfNodesPerCell() const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *computeNbOfFacesPerCell() const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *computeEffectiveNbOfNodesPerCell() const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingMesh *buildPartRange(int beginCellIds, int endCellIds, int stepCellIds) const throw(INTERP_KERNEL::Exception);
    virtual int getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType type) const throw(INTERP_KERNEL::Exception);
    virtual INTERP_KERNEL::NormalizedCellType getTypeOfCell(int cellId) const throw(INTERP_KERNEL::Exception);
    virtual std::string simpleRepr() const throw(INTERP_KERNEL::Exception);
    virtual std::string advancedRepr() const throw(INTERP_KERNEL::Exception);
    std::string writeVTK(const std::string& fileName, bool isBinary=true) const throw(INTERP_KERNEL::Exception);
    virtual std::string getVTKFileExtension() const;
    std::string getVTKFileNameOf(const std::string& fileName) const;
    // tools
    virtual MEDCouplingFieldDouble *getMeasureField(bool isAbs) const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingFieldDouble *getMeasureFieldOnNode(bool isAbs) const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingFieldDouble *fillFromAnalytic(TypeOfField t, int nbOfComp, const std::string& func) const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingFieldDouble *fillFromAnalyticCompo(TypeOfField t, int nbOfComp, const std::string& func) const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingFieldDouble *fillFromAnalyticNamedCompo(TypeOfField t, int nbOfComp, const std::vector<std::string>& varsOrder, const std::string& func) const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingFieldDouble *buildOrthogonalField() const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingUMesh *buildUnstructured() const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingMesh *mergeMyselfWith(const MEDCouplingMesh *other) const throw(INTERP_KERNEL::Exception);
    virtual bool areCompatibleForMerge(const MEDCouplingMesh *other) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *simplexize(int policy) throw(INTERP_KERNEL::Exception);
    virtual void unserialization(const std::vector<double>& tinyInfoD, const std::vector<int>& tinyInfo, const DataArrayInt *a1, DataArrayDouble *a2, const std::vector<std::string>& littleStrings) throw(INTERP_KERNEL::Exception);
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

         DataArrayDouble *getDirectAccessOfCoordsArrIfInStructure() const throw(INTERP_KERNEL::Exception)
         {
           const DataArrayDouble *ret(self->getDirectAccessOfCoordsArrIfInStructure());
           DataArrayDouble *ret2(const_cast<DataArrayDouble *>(ret));
           if(ret2)
             ret2->incrRef();
           return ret2;
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
           MCAuto<DataArrayInt> elts,eltsIndex;
           self->getCellsContainingPoints(pos,nbOfPoints,eps,elts,eltsIndex);
           PyObject *ret=PyTuple_New(2);
           PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(elts.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
           PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(eltsIndex.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
           return ret;
         }

         PyObject *getCellsContainingPoints(PyObject *p, double eps) const throw(INTERP_KERNEL::Exception)
         {
           MCAuto<DataArrayInt> elts,eltsIndex;
           int spaceDim=self->getSpaceDimension();
           void *da=0;
           int res1=SWIG_ConvertPtr(p,&da,SWIGTYPE_p_MEDCoupling__DataArrayDouble, 0 |  0 );
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
           PyObject *ret=PyTuple_New(2);
           PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(elts.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
           PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(eltsIndex.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
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
           return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 );
         }
         
         virtual PyObject *getReverseNodalConnectivity() const throw(INTERP_KERNEL::Exception)
         {
           MCAuto<DataArrayInt> d0=DataArrayInt::New();
           MCAuto<DataArrayInt> d1=DataArrayInt::New();
           self->getReverseNodalConnectivity(d0,d1);
           PyObject *ret=PyTuple_New(2);
           PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(d0.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
           PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(d1.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
           return ret;
         }
         
         void renumberCells(PyObject *li, bool check=true) throw(INTERP_KERNEL::Exception)
         {
           int sw,sz(-1);
           int v0; std::vector<int> v1;
           const int *ids(convertIntStarLikePyObjToCppIntStar(li,sw,sz,v0,v1));
           self->renumberCells(ids,check);
         }

         PyObject *checkGeoEquivalWith(const MEDCouplingMesh *other, int levOfCheck, double prec) const throw(INTERP_KERNEL::Exception)
         {
           DataArrayInt *cellCor, *nodeCor;
           self->checkGeoEquivalWith(other,levOfCheck,prec,cellCor,nodeCor);
           PyObject *res = PyList_New(2);
           PyList_SetItem(res,0,SWIG_NewPointerObj(SWIG_as_voidptr(cellCor),SWIGTYPE_p_MEDCoupling__DataArrayInt, cellCor?SWIG_POINTER_OWN | 0:0 ));
           PyList_SetItem(res,1,SWIG_NewPointerObj(SWIG_as_voidptr(nodeCor),SWIGTYPE_p_MEDCoupling__DataArrayInt, nodeCor?SWIG_POINTER_OWN | 0:0 ));
           return res;
         }

         PyObject *checkDeepEquivalWith(const MEDCouplingMesh *other, int cellCompPol, double prec) const throw(INTERP_KERNEL::Exception)
         {
           DataArrayInt *cellCor=0,*nodeCor=0;
           self->checkDeepEquivalWith(other,cellCompPol,prec,cellCor,nodeCor);
           PyObject *res = PyList_New(2);
           PyList_SetItem(res,0,SWIG_NewPointerObj(SWIG_as_voidptr(cellCor),SWIGTYPE_p_MEDCoupling__DataArrayInt, cellCor?SWIG_POINTER_OWN | 0:0 ));
           PyList_SetItem(res,1,SWIG_NewPointerObj(SWIG_as_voidptr(nodeCor),SWIGTYPE_p_MEDCoupling__DataArrayInt, nodeCor?SWIG_POINTER_OWN | 0:0 ));
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
           int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_MEDCoupling__DataArrayInt, 0 |  0 );
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
           PyObject *ret=convertDblArrToPyListOfTuple<double>(tmp,2,spaceDim);
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
           const int *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
           MEDCouplingMesh *ret=self->buildPart(tmp,tmp+szArr);
           if(sw==3)//DataArrayInt
             { 
               void *argp; SWIG_ConvertPtr(li,&argp,SWIGTYPE_p_MEDCoupling__DataArrayInt,0|0);
               DataArrayInt *argpt=reinterpret_cast< MEDCoupling::DataArrayInt * >(argp);
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
           const int *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
           MEDCouplingMesh *ret=self->buildPartAndReduceNodes(tmp,tmp+szArr,arr);
           if(sw==3)//DataArrayInt
             { 
               void *argp; SWIG_ConvertPtr(li,&argp,SWIGTYPE_p_MEDCoupling__DataArrayInt,0|0);
               DataArrayInt *argpt=reinterpret_cast< MEDCoupling::DataArrayInt * >(argp);
               std::string name=argpt->getName();
               if(!name.empty())
                 ret->setName(name.c_str());
             }
           //
           PyObject *res = PyList_New(2);
           PyObject *obj0=convertMesh(ret, SWIG_POINTER_OWN | 0 );
           PyObject *obj1=SWIG_NewPointerObj(SWIG_as_voidptr(arr),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 );
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
             obj1=SWIG_NewPointerObj(SWIG_as_voidptr(arr),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 );
           else
             obj1=PySlice_New(PyInt_FromLong(a),PyInt_FromLong(b),PyInt_FromLong(b));
           PyTuple_SetItem(res,0,obj0);
           PyTuple_SetItem(res,1,obj1);
           return res;
         }

        PyObject *getDistributionOfTypes() const throw(INTERP_KERNEL::Exception)
        {
          std::vector<int> vals=self->getDistributionOfTypes();
          if(vals.size()%3!=0)
            throw INTERP_KERNEL::Exception("Internal Error detected in wrap python ! code returned by MEDCouplingMesh::getDistributionOfTypes is not so that %3==0 !");
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
          convertFromPyObjVectorOfObj<const MEDCoupling::DataArrayInt *>(li2,SWIGTYPE_p_MEDCoupling__DataArrayInt,"DataArrayInt",idsPerType);
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
          //
          if(code.size()%3!=0)
            throw INTERP_KERNEL::Exception("Internal Error detected in wrap python ! code returned by MEDCouplingMesh::splitProfilePerType is not so that %3==0 !");
          PyObject *ret0=PyList_New((int)code.size()/3);
          for(int j=0;j<(int)code.size()/3;j++)
             {
               PyObject *ret00=PyList_New(3);
               PyList_SetItem(ret00,0,SWIG_From_int(code[3*j]));
               PyList_SetItem(ret00,1,SWIG_From_int(code[3*j+1]));
               PyList_SetItem(ret00,2,SWIG_From_int(code[3*j+2]));
               PyList_SetItem(ret0,j,ret00);
             }
          PyTuple_SetItem(ret,0,ret0);
          //
          PyObject *ret1=PyList_New(idsInPflPerType.size());
          for(std::size_t j=0;j<idsInPflPerType.size();j++)
            PyList_SetItem(ret1,j,SWIG_NewPointerObj(SWIG_as_voidptr(idsInPflPerType[j]),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
          PyTuple_SetItem(ret,1,ret1);
          int n=idsPerType.size();
          PyObject *ret2=PyList_New(n);
          for(int i=0;i<n;i++)
            PyList_SetItem(ret2,i,SWIG_NewPointerObj(SWIG_as_voidptr(idsPerType[i]),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
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

         virtual PyObject *getTinySerializationInformation() const throw(INTERP_KERNEL::Exception)
         {
           std::vector<double> a0;
           std::vector<int> a1;
           std::vector<std::string> a2;
           self->getTinySerializationInformation(a0,a1,a2);
           PyObject *ret(PyTuple_New(3));
           PyTuple_SetItem(ret,0,convertDblArrToPyList2(a0));
           PyTuple_SetItem(ret,1,convertIntArrToPyList2(a1));
           int sz(a2.size());
           PyObject *ret2(PyList_New(sz));
           {
             for(int i=0;i<sz;i++)
               PyList_SetItem(ret2,i,PyString_FromString(a2[i].c_str()));
           }
           PyTuple_SetItem(ret,2,ret2);
           return ret;
         }

         virtual PyObject *serialize() const throw(INTERP_KERNEL::Exception)
         {
           DataArrayInt *a0Tmp(0);
           DataArrayDouble *a1Tmp(0);
           self->serialize(a0Tmp,a1Tmp);
           PyObject *ret(PyTuple_New(2));
           PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(a0Tmp),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
           PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(a1Tmp),SWIGTYPE_p_MEDCoupling__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
           return ret;
         }

         void resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2) const throw(INTERP_KERNEL::Exception)
         {
           std::vector<std::string> littleStrings;
           self->resizeForUnserialization(tinyInfo,a1,a2,littleStrings);
         }
         
         PyObject *__getnewargs__() throw(INTERP_KERNEL::Exception)
         {// put an empty dict in input to say to __new__ to call __init__...
           PyObject *ret(PyTuple_New(1));
           PyObject *ret0(PyDict_New());
           PyTuple_SetItem(ret,0,ret0);
           return ret;
         }
         
         PyObject *__getstate__() const throw(INTERP_KERNEL::Exception)
         {
           PyObject *ret0(MEDCoupling_MEDCouplingMesh_getTinySerializationInformation(self));
           PyObject *ret1(MEDCoupling_MEDCouplingMesh_serialize(self));
           PyObject *ret(PyTuple_New(2));
           PyTuple_SetItem(ret,0,ret0);
           PyTuple_SetItem(ret,1,ret1);
           return ret;
         }

         void __setstate__(PyObject *inp) throw(INTERP_KERNEL::Exception)
         {
           static const char MSG[]="MEDCouplingMesh.__setstate__ : expected input is a tuple of size 2 !";
           if(!PyTuple_Check(inp))
             throw INTERP_KERNEL::Exception(MSG);
           int sz(PyTuple_Size(inp));
           if(sz!=2)
             throw INTERP_KERNEL::Exception(MSG);
           PyObject *elt0(PyTuple_GetItem(inp,0));
           PyObject *elt1(PyTuple_GetItem(inp,1));
           std::vector<double> a0;
           std::vector<int> a1;
           std::vector<std::string> a2;
           DataArrayInt *b0(0);
           DataArrayDouble *b1(0);
           {
             if(!PyTuple_Check(elt0) && PyTuple_Size(elt0)!=3)
               throw INTERP_KERNEL::Exception(MSG);
             PyObject *a0py(PyTuple_GetItem(elt0,0)),*a1py(PyTuple_GetItem(elt0,1)),*a2py(PyTuple_GetItem(elt0,2));
             int tmp(-1);
             fillArrayWithPyListDbl3(a0py,tmp,a0);
             convertPyToNewIntArr3(a1py,a1);
             fillStringVector(a2py,a2);
           }
           {
             if(!PyTuple_Check(elt1) && PyTuple_Size(elt1)!=2)
               throw INTERP_KERNEL::Exception(MSG);
             PyObject *b0py(PyTuple_GetItem(elt1,0)),*b1py(PyTuple_GetItem(elt1,1));
             void *argp(0);
             int status(SWIG_ConvertPtr(b0py,&argp,SWIGTYPE_p_MEDCoupling__DataArrayInt,0|0));
             if(!SWIG_IsOK(status))
               throw INTERP_KERNEL::Exception(MSG);
             b0=reinterpret_cast<DataArrayInt *>(argp);
             status=SWIG_ConvertPtr(b1py,&argp,SWIGTYPE_p_MEDCoupling__DataArrayDouble,0|0);
             if(!SWIG_IsOK(status))
               throw INTERP_KERNEL::Exception(MSG);
             b1=reinterpret_cast<DataArrayDouble *>(argp);
           }
           // useless here to call resizeForUnserialization because arrays are well resized.
           self->unserialization(a0,a1,b0,b1,a2);
         }
         
         static MEDCouplingMesh *MergeMeshes(PyObject *li) throw(INTERP_KERNEL::Exception)
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
    static const char *GetRepr(NatureOfField nat) throw(INTERP_KERNEL::Exception);
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
                                 const std::vector<double>& gsCoo, const std::vector<double>& w) throw(INTERP_KERNEL::Exception);
    MEDCouplingGaussLocalization(INTERP_KERNEL::NormalizedCellType typ) throw(INTERP_KERNEL::Exception);
    INTERP_KERNEL::NormalizedCellType getType() const throw(INTERP_KERNEL::Exception);
    void setType(INTERP_KERNEL::NormalizedCellType typ) throw(INTERP_KERNEL::Exception);
    int getNumberOfGaussPt() const throw(INTERP_KERNEL::Exception);
    int getDimension() const throw(INTERP_KERNEL::Exception);
    int getNumberOfPtsInRefCell() const throw(INTERP_KERNEL::Exception);
    std::string getStringRepr() const throw(INTERP_KERNEL::Exception);
    void checkConsistencyLight() const throw(INTERP_KERNEL::Exception);
    bool isEqual(const MEDCouplingGaussLocalization& other, double eps) const throw(INTERP_KERNEL::Exception);
    //
    const std::vector<double>& getRefCoords() const throw(INTERP_KERNEL::Exception);
    double getRefCoord(int ptIdInCell, int comp) const throw(INTERP_KERNEL::Exception);
    const std::vector<double>& getGaussCoords() const throw(INTERP_KERNEL::Exception);
    double getGaussCoord(int gaussPtIdInCell, int comp) const throw(INTERP_KERNEL::Exception);
    const std::vector<double>& getWeights() const throw(INTERP_KERNEL::Exception);
    double getWeight(int gaussPtIdInCell, double newVal) const throw(INTERP_KERNEL::Exception);
    void setRefCoord(int ptIdInCell, int comp, double newVal) throw(INTERP_KERNEL::Exception);
    void setGaussCoord(int gaussPtIdInCell, int comp, double newVal) throw(INTERP_KERNEL::Exception);
    void setWeight(int gaussPtIdInCell, double newVal) throw(INTERP_KERNEL::Exception);
    void setRefCoords(const std::vector<double>& refCoo) throw(INTERP_KERNEL::Exception);
    void setGaussCoords(const std::vector<double>& gsCoo) throw(INTERP_KERNEL::Exception);
    void setWeights(const std::vector<double>& w) throw(INTERP_KERNEL::Exception);
    //
    static bool AreAlmostEqual(const std::vector<double>& v1, const std::vector<double>& v2, double eps);
    //
    %extend 
    {
      DataArrayDouble *localizePtsInRefCooForEachCell(const DataArrayDouble *ptsInRefCoo, const MEDCouplingUMesh *mesh) const throw(INTERP_KERNEL::Exception)
      {
        MCAuto<DataArrayDouble> ret(self->localizePtsInRefCooForEachCell(ptsInRefCoo,mesh));
        return ret.retn();
      }

      MEDCouplingUMesh *buildRefCell() const throw(INTERP_KERNEL::Exception)
      {
        MCAuto<MEDCouplingUMesh> ret(self->buildRefCell());
        return ret.retn();
      }
    }
  };

  class MEDCouplingSkyLineArray
  {
  public:  
    static MEDCouplingSkyLineArray *BuildFromPolyhedronConn( const DataArrayInt* c, const DataArrayInt* cI ) throw(INTERP_KERNEL::Exception);
  
    void set( DataArrayInt* index, DataArrayInt* value );
    void set3( DataArrayInt* superIndex, DataArrayInt* index, DataArrayInt* value );
    
    int getSuperNumberOf() const;
    int getNumberOf() const;
    int getLength() const;
    
    void deletePack(const int i, const int j) throw(INTERP_KERNEL::Exception);
    
    %extend 
    {
      MEDCouplingSkyLineArray() throw(INTERP_KERNEL::Exception)
      {
        return MEDCouplingSkyLineArray::New();
      }

      MEDCouplingSkyLineArray( const std::vector<int>& index, const std::vector<int>& value) throw(INTERP_KERNEL::Exception)
      {
        return MEDCouplingSkyLineArray::New(index, value);
      }

      MEDCouplingSkyLineArray( DataArrayInt* index, DataArrayInt* value ) throw(INTERP_KERNEL::Exception)
      {
        return MEDCouplingSkyLineArray::New(index, value);
      }

      MEDCouplingSkyLineArray( const MEDCouplingSkyLineArray & other ) throw(INTERP_KERNEL::Exception)
      {
        return MEDCouplingSkyLineArray::New(other);
      }

      std::string __str__() const throw(INTERP_KERNEL::Exception)
      {
        return self->simpleRepr();
      }
      
      DataArrayInt *getSuperIndexArray() const
      {
        DataArrayInt *ret(self->getSuperIndexArray());
        if(ret)
          ret->incrRef();
        return ret;
      }
      
      DataArrayInt *getIndexArray() const
      {
        DataArrayInt *ret(self->getIndexArray());
        if(ret)
          ret->incrRef();
        return ret;
      }
      
      DataArrayInt *getValuesArray() const
      {
        DataArrayInt *ret(self->getValuesArray());
        if(ret)
          ret->incrRef();
        return ret;
      }
     
      PyObject *getSimplePackSafe(int absolutePackId) const throw(INTERP_KERNEL::Exception)
      {
        std::vector<int> ret;
        self->getSimplePackSafe(absolutePackId,ret);
        return convertIntArrToPyList2(ret);
      }

      PyObject *findPackIds(PyObject *superPackIndices, PyObject *pack) const throw(INTERP_KERNEL::Exception)
      {
          std::vector<int> vpack, vspIdx, out;
          
          convertPyToNewIntArr3(superPackIndices,vspIdx);
          convertPyToNewIntArr3(pack,vpack);
          
          self->findPackIds(vspIdx, vpack.data(), vpack.data()+vpack.size(), out);
          return convertIntArrToPyList2(out);
      }
      
      void pushBackPack(const int i, PyObject *pack) throw(INTERP_KERNEL::Exception)
        {
          std::vector<int> vpack;
          convertPyToNewIntArr3(pack,vpack);
          self->pushBackPack(i,vpack.data(), vpack.data()+vpack.size());
        }
        
      void replaceSimplePack(const int idx, PyObject *pack) throw(INTERP_KERNEL::Exception)
        {
          std::vector<int> vpack;
          convertPyToNewIntArr3(pack,vpack);
          self->replaceSimplePack(idx, vpack.data(), vpack.data()+vpack.size());
        }
        
      void replacePack(const int superIdx, const int idx, PyObject *pack) throw(INTERP_KERNEL::Exception)
        {
          std::vector<int> vpack;
          convertPyToNewIntArr3(pack,vpack);
          self->replacePack(superIdx, idx, vpack.data(), vpack.data()+vpack.size());
        }

      PyObject *convertToPolyhedronConn() const throw(INTERP_KERNEL::Exception)
         {
           MCAuto<DataArrayInt> d0=DataArrayInt::New();
           MCAuto<DataArrayInt> d1=DataArrayInt::New();
           self->convertToPolyhedronConn(d0,d1);
           PyObject *ret=PyTuple_New(2);
           PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(d0.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
           PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(d1.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
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
      void setCoords(const DataArrayDouble *coords) throw(INTERP_KERNEL::Exception);
      DataArrayDouble *getCoordinatesAndOwner() const throw(INTERP_KERNEL::Exception);
      bool areCoordsEqual(const MEDCouplingPointSet& other, double prec) const throw(INTERP_KERNEL::Exception);
      void zipCoords() throw(INTERP_KERNEL::Exception);
      double getCaracteristicDimension() const throw(INTERP_KERNEL::Exception);
      void recenterForMaxPrecision(double eps) throw(INTERP_KERNEL::Exception);
      void changeSpaceDimension(int newSpaceDim, double dftVal=0.) throw(INTERP_KERNEL::Exception);
      void tryToShareSameCoords(const MEDCouplingPointSet& other, double epsilon) throw(INTERP_KERNEL::Exception);
      virtual void shallowCopyConnectivityFrom(const MEDCouplingPointSet *other) throw(INTERP_KERNEL::Exception);
      virtual MEDCouplingPointSet *buildPartOfMySelfSlice(int start, int end, int step) const throw(INTERP_KERNEL::Exception);
      virtual void tryToShareSameCoordsPermute(const MEDCouplingPointSet& other, double epsilon) throw(INTERP_KERNEL::Exception);
      static DataArrayDouble *MergeNodesArray(const MEDCouplingPointSet *m1, const MEDCouplingPointSet *m2) throw(INTERP_KERNEL::Exception);
      static MEDCouplingPointSet *BuildInstanceFromMeshType(MEDCouplingMeshType type) throw(INTERP_KERNEL::Exception);
      static DataArrayInt *ComputeNbOfInteractionsWithSrcCells(const MEDCouplingPointSet *srcMesh, const MEDCouplingPointSet *trgMesh, double eps) throw(INTERP_KERNEL::Exception);
      virtual DataArrayInt *computeFetchedNodeIds() const throw(INTERP_KERNEL::Exception);
      virtual int getNumberOfNodesInCell(int cellId) const throw(INTERP_KERNEL::Exception);
      virtual MEDCouplingPointSet *buildBoundaryMesh(bool keepCoords) const throw(INTERP_KERNEL::Exception);
      virtual DataArrayInt *getCellsInBoundingBox(const INTERP_KERNEL::DirectedBoundingBox& bbox, double eps) throw(INTERP_KERNEL::Exception);
      virtual DataArrayInt *zipCoordsTraducer() throw(INTERP_KERNEL::Exception);
      virtual DataArrayInt *findBoundaryNodes() const;
      virtual DataArrayInt *zipConnectivityTraducer(int compType, int startCellId=0) throw(INTERP_KERNEL::Exception);
      virtual MEDCouplingPointSet *mergeMyselfWithOnSameCoords(const MEDCouplingPointSet *other) const throw(INTERP_KERNEL::Exception);
      virtual void checkFullyDefined() const throw(INTERP_KERNEL::Exception);
      virtual bool isEmptyMesh(const std::vector<int>& tinyInfo) const throw(INTERP_KERNEL::Exception);
      virtual MEDCouplingPointSet *deepCopyConnectivityOnly() const throw(INTERP_KERNEL::Exception);
      virtual DataArrayDouble *getBoundingBoxForBBTree(double arcDetEps=1e-12) const throw(INTERP_KERNEL::Exception);
      virtual void renumberNodesWithOffsetInConn(int offset) throw(INTERP_KERNEL::Exception);
      virtual bool areAllNodesFetched() const throw(INTERP_KERNEL::Exception);
      virtual MEDCouplingFieldDouble *computeDiameterField() const throw(INTERP_KERNEL::Exception);
      virtual void invertOrientationOfAllCells() throw(INTERP_KERNEL::Exception);
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
             PyList_SetItem(res,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
             PyList_SetItem(res,1,SWIG_From_int(newNbOfNodes));
             return res;
           }
           
           PyObject *findCommonNodes(double prec, int limitTupleId=-1) const throw(INTERP_KERNEL::Exception)
           {
             DataArrayInt *comm, *commIndex;
             self->findCommonNodes(prec,limitTupleId,comm,commIndex);
             PyObject *res = PyList_New(2);
             PyList_SetItem(res,0,SWIG_NewPointerObj(SWIG_as_voidptr(comm),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
             PyList_SetItem(res,1,SWIG_NewPointerObj(SWIG_as_voidptr(commIndex),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
             return res;
           }
           
           PyObject *getCoords() throw(INTERP_KERNEL::Exception)
           {
             DataArrayDouble *ret1=self->getCoords();
             if (ret1)
                ret1->incrRef();
             return SWIG_NewPointerObj((void*)ret1,SWIGTYPE_p_MEDCoupling__DataArrayDouble,SWIG_POINTER_OWN | 0);
           }
           
           PyObject *buildPartOfMySelf(PyObject *li, bool keepCoords=true) const throw(INTERP_KERNEL::Exception)
           {
             int szArr,sw,iTypppArr;
             std::vector<int> stdvecTyyppArr;
             const int *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
             MEDCouplingPointSet *ret=self->buildPartOfMySelf(tmp,tmp+szArr,keepCoords);
             if(sw==3)//DataArrayInt
               { 
                 void *argp; SWIG_ConvertPtr(li,&argp,SWIGTYPE_p_MEDCoupling__DataArrayInt,0|0);
                 DataArrayInt *argpt=reinterpret_cast< MEDCoupling::DataArrayInt * >(argp);
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
             const int *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
             MEDCouplingPointSet *ret=self->buildPartOfMySelfNode(tmp,tmp+szArr,fullyIn);
             if(sw==3)//DataArrayInt
               { 
                 void *argp; SWIG_ConvertPtr(li,&argp,SWIGTYPE_p_MEDCoupling__DataArrayInt,0|0);
                 DataArrayInt *argpt=reinterpret_cast< MEDCoupling::DataArrayInt * >(argp);
                 std::string name=argpt->getName();
                 if(!name.empty())
                   ret->setName(name.c_str());
               }
             return convertMesh(ret, SWIG_POINTER_OWN | 0 );
           }

           virtual PyObject *buildPartOfMySelfKeepCoords(PyObject *li) const throw(INTERP_KERNEL::Exception)
           {
             int szArr,sw,iTypppArr;
             std::vector<int> stdvecTyyppArr;
             const int *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
             MEDCouplingPointSet *ret=self->buildPartOfMySelfKeepCoords(tmp,tmp+szArr);
             if(sw==3)//DataArrayInt
               { 
                 void *argp; SWIG_ConvertPtr(li,&argp,SWIGTYPE_p_MEDCoupling__DataArrayInt,0|0);
                 DataArrayInt *argpt=reinterpret_cast< MEDCoupling::DataArrayInt * >(argp);
                 std::string name=argpt->getName();
                 if(!name.empty())
                   ret->setName(name.c_str());
               }
             return convertMesh(ret, SWIG_POINTER_OWN | 0 );
           }

           virtual PyObject *buildPartOfMySelfKeepCoordsSlice(int start, int end, int step) const throw(INTERP_KERNEL::Exception)
           {
             MEDCouplingPointSet *ret=self->buildPartOfMySelfKeepCoordsSlice(start,end,step);
             return convertMesh(ret, SWIG_POINTER_OWN | 0 );
           }

           PyObject *buildFacePartOfMySelfNode(PyObject *li, bool fullyIn) const throw(INTERP_KERNEL::Exception)
           {
             int szArr,sw,iTypppArr;
             std::vector<int> stdvecTyyppArr;
             const int *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
             MEDCouplingPointSet *ret=self->buildFacePartOfMySelfNode(tmp,tmp+szArr,fullyIn);
             if(sw==3)//DataArrayInt
               { 
                 void *argp; SWIG_ConvertPtr(li,&argp,SWIGTYPE_p_MEDCoupling__DataArrayInt,0|0);
                 DataArrayInt *argpt=reinterpret_cast< MEDCoupling::DataArrayInt * >(argp);
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
             const int *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
             self->renumberNodes(tmp,newNbOfNodes);
           }

           void renumberNodesCenter(PyObject *li, int newNbOfNodes) throw(INTERP_KERNEL::Exception)
           {
             int szArr,sw,iTypppArr;
             std::vector<int> stdvecTyyppArr;
             const int *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
             self->renumberNodesCenter(tmp,newNbOfNodes);
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
               return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 );
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
               return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 );
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
             return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 );
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
             PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(c),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
             PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(cI),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
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
             PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(c),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
             PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(cI),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
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
             return SWIG_NewPointerObj(SWIG_as_voidptr(elems),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 );
           }

           void duplicateNodesInCoords(PyObject *li) throw(INTERP_KERNEL::Exception)
           {
             int sw;
             int singleVal;
             std::vector<int> multiVal;
             std::pair<int, std::pair<int,int> > slic;
             MEDCoupling::DataArrayInt *daIntTyypp=0;
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
                 throw INTERP_KERNEL::Exception("MEDCouplingPointSet::duplicateNodesInCoords : unrecognized type entered, expected list of int, tuple of int or DataArrayInt !");
               }
           }

           virtual PyObject *findCommonCells(int compType, int startCellId=0) const throw(INTERP_KERNEL::Exception)
           {
             DataArrayInt *v0=0,*v1=0;
             self->findCommonCells(compType,startCellId,v0,v1);
             PyObject *res = PyList_New(2);
             PyList_SetItem(res,0,SWIG_NewPointerObj(SWIG_as_voidptr(v0),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
             PyList_SetItem(res,1,SWIG_NewPointerObj(SWIG_as_voidptr(v1),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
             return res;
           }

      
           virtual void renumberNodesInConn(PyObject *li) throw(INTERP_KERNEL::Exception)
           {
             void *da=0;
             int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_MEDCoupling__DataArrayInt, 0 | 0 );
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

           virtual PyObject *getNodeIdsInUse() const throw(INTERP_KERNEL::Exception)
           {
             int ret1=-1;
             DataArrayInt *ret0=self->getNodeIdsInUse(ret1);
             PyObject *ret=PyTuple_New(2);
             PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
             PyTuple_SetItem(ret,1,PyInt_FromLong(ret1));
             return ret;
           }

           virtual DataArrayInt *fillCellIdsToKeepFromNodeIds(PyObject *li, bool fullyIn) const
           {
             DataArrayInt *ret=0;
             //
             int szArr,sw,iTypppArr;
             std::vector<int> stdvecTyyppArr;
             const int *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
             self->fillCellIdsToKeepFromNodeIds(tmp,tmp+szArr,fullyIn,ret);
             return ret;
           }

           virtual PyObject *mergeNodes(double precision) throw(INTERP_KERNEL::Exception)
           {
             bool ret1;
             int ret2;
             DataArrayInt *ret0=self->mergeNodes(precision,ret1,ret2);
             PyObject *res = PyList_New(3);
             PyList_SetItem(res,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
             PyList_SetItem(res,1,SWIG_From_bool(ret1));
             PyList_SetItem(res,2,SWIG_From_int(ret2));
             return res;
           }
           
           virtual PyObject *mergeNodesCenter(double precision) throw(INTERP_KERNEL::Exception)
           {
             bool ret1;
             int ret2;
             DataArrayInt *ret0=self->mergeNodesCenter(precision,ret1,ret2);
             PyObject *res = PyList_New(3);
             PyList_SetItem(res,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
             PyList_SetItem(res,1,SWIG_From_bool(ret1));
             PyList_SetItem(res,2,SWIG_From_int(ret2));
             return res;
           }
           
           DataArrayInt *getCellIdsLyingOnNodes(PyObject *li, bool fullyIn) const throw(INTERP_KERNEL::Exception)
           {
             void *da=0;
             int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_MEDCoupling__DataArrayInt, 0 |  0 );
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

           MEDCouplingPointSet *__getitem__(PyObject *listOrDataArrI) throw(INTERP_KERNEL::Exception)
           {
             int sw;
             int singleVal;
             std::vector<int> multiVal;
             std::pair<int, std::pair<int,int> > slic;
             MEDCoupling::DataArrayInt *daIntTyypp=0;
             int nbc=self->getNumberOfCells();
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
                 throw INTERP_KERNEL::Exception("MEDCouplingUMesh::__getitem__ : unrecognized type in input ! Possibilities are : int, list or tuple of int DataArrayInt instance !");
               }
           }
           
           static void Rotate2DAlg(PyObject *center, double angle, int nbNodes, PyObject *coords) throw(INTERP_KERNEL::Exception)
           {
             int sz;
             INTERP_KERNEL::AutoCPtr<double> c=convertPyToNewDblArr2(center,&sz);
             INTERP_KERNEL::AutoCPtr<double> coo=convertPyToNewDblArr2(coords,&sz);
             MEDCoupling::DataArrayDouble::Rotate2DAlg(c,angle,nbNodes,coo,coo);
             for(int i=0;i<sz;i++)
               PyList_SetItem(coords,i,PyFloat_FromDouble(coo[i]));
           }
           
           static void Rotate2DAlg(PyObject *center, double angle, PyObject *coords) throw(INTERP_KERNEL::Exception)
           {
             int sz;
             INTERP_KERNEL::AutoCPtr<double> c=convertPyToNewDblArr2(center,&sz);
             int sw,nbNodes=0;
             double val0;  MEDCoupling::DataArrayDouble *val1=0; MEDCoupling::DataArrayDoubleTuple *val2=0;
             std::vector<double> val3;
             const double *coo=convertObjToPossibleCpp5_Safe2(coords,sw,val0,val1,val2,val3,
                                                            "Rotate2DAlg",2,true,nbNodes);
             if(sw!=2 && sw!=3)
               throw INTERP_KERNEL::Exception("Invalid call to MEDCouplingPointSet::Rotate2DAlg : try another overload method !");
             MEDCoupling::DataArrayDouble::Rotate2DAlg(c,angle,nbNodes,coo,const_cast<double *>(coo));
           }
           
           static void Rotate3DAlg(PyObject *center, PyObject *vect, double angle, int nbNodes, PyObject *coords) throw(INTERP_KERNEL::Exception)
           {
             int sz,sz2;
             INTERP_KERNEL::AutoCPtr<double> c=convertPyToNewDblArr2(center,&sz);
             INTERP_KERNEL::AutoCPtr<double> coo=convertPyToNewDblArr2(coords,&sz);
             INTERP_KERNEL::AutoCPtr<double> v=convertPyToNewDblArr2(vect,&sz2);
             MEDCoupling::DataArrayDouble::Rotate3DAlg(c,v,angle,nbNodes,coo,coo);
             for(int i=0;i<sz;i++)
               PyList_SetItem(coords,i,PyFloat_FromDouble(coo[i]));
           }
           
           static void Rotate3DAlg(PyObject *center, PyObject *vect, double angle, PyObject *coords) throw(INTERP_KERNEL::Exception)
           {
             int sz,sz2;
             INTERP_KERNEL::AutoCPtr<double> c=convertPyToNewDblArr2(center,&sz);
             int sw,nbNodes=0;
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
    static MEDCouplingUMesh *New() throw(INTERP_KERNEL::Exception);
    static MEDCouplingUMesh *New(const char *meshName, int meshDim) throw(INTERP_KERNEL::Exception);
    void checkConsistencyLight() const throw(INTERP_KERNEL::Exception);
    void setMeshDimension(int meshDim) throw(INTERP_KERNEL::Exception);
    void allocateCells(int nbOfCells=0) throw(INTERP_KERNEL::Exception);
    void finishInsertingCells() throw(INTERP_KERNEL::Exception);
    MEDCouplingUMeshCellByTypeEntry *cellsByType() throw(INTERP_KERNEL::Exception);
    void setConnectivity(DataArrayInt *conn, DataArrayInt *connIndex, bool isComputingTypes=true) throw(INTERP_KERNEL::Exception);
    INTERP_KERNEL::NormalizedCellType getTypeOfCell(int cellId) const throw(INTERP_KERNEL::Exception);
    void setPartOfMySelfSlice(int start, int end, int step, const MEDCouplingUMesh& otherOnSameCoordsThanThis) throw(INTERP_KERNEL::Exception);
    int getNodalConnectivityArrayLen() const throw(INTERP_KERNEL::Exception);
    void computeTypes() throw(INTERP_KERNEL::Exception);
    std::string reprConnectivityOfThis() const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *buildSetInstanceFromThis(int spaceDim) const throw(INTERP_KERNEL::Exception);
    //tools
    DataArrayInt *conformize2D(double eps) throw(INTERP_KERNEL::Exception);
    DataArrayInt *conformize3D(double eps) throw(INTERP_KERNEL::Exception);
    DataArrayInt *colinearize2D(double eps) throw(INTERP_KERNEL::Exception);
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
    MEDCouplingUMesh *buildDescendingConnectivity(DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *revDesc, DataArrayInt *revDescIndx) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *buildDescendingConnectivity2(DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *revDesc, DataArrayInt *revDescIndx) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *explode3DMeshTo1D(DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *revDesc, DataArrayInt *revDescIndx) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *explodeMeshIntoMicroEdges(DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *revDesc, DataArrayInt *revDescIndx) const throw(INTERP_KERNEL::Exception);
    void orientCorrectlyPolyhedrons() throw(INTERP_KERNEL::Exception);
    bool isPresenceOfQuadratic() const throw(INTERP_KERNEL::Exception);
    bool isFullyQuadratic() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *buildDirectionVectorField() const throw(INTERP_KERNEL::Exception);
    bool isContiguous1D() const throw(INTERP_KERNEL::Exception);
    void tessellate2D(double eps) throw(INTERP_KERNEL::Exception);
    void convertQuadraticCellsToLinear() throw(INTERP_KERNEL::Exception);
    DataArrayInt *convertLinearCellsToQuadratic(int conversionType=0) throw(INTERP_KERNEL::Exception);
    void convertDegeneratedCells() throw(INTERP_KERNEL::Exception);
    bool areOnlySimplexCells() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getEdgeRatioField() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getAspectRatioField() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getWarpField() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getSkewField() const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *computePlaneEquationOf3DFaces() const throw(INTERP_KERNEL::Exception);
    DataArrayInt *convexEnvelop2D() throw(INTERP_KERNEL::Exception);
    std::string cppRepr() const throw(INTERP_KERNEL::Exception);
    DataArrayInt *findAndCorrectBadOriented3DExtrudedCells() throw(INTERP_KERNEL::Exception);
    DataArrayInt *findAndCorrectBadOriented3DCells() throw(INTERP_KERNEL::Exception);
    MEDCoupling::MEDCoupling1GTUMesh *convertIntoSingleGeoTypeMesh() const throw(INTERP_KERNEL::Exception);
    MEDCouplingSkyLineArray *generateGraph() const throw(INTERP_KERNEL::Exception);
    DataArrayInt *convertNodalConnectivityToStaticGeoTypeMesh() const throw(INTERP_KERNEL::Exception);
    DataArrayInt *buildUnionOf2DMesh() const throw(INTERP_KERNEL::Exception);
    DataArrayInt *buildUnionOf3DMesh() const throw(INTERP_KERNEL::Exception);
    DataArrayInt *orderConsecutiveCells1D() const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getBoundingBoxForBBTreeFast() const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getBoundingBoxForBBTree2DQuadratic(double arcDetEps=1e-12) const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getBoundingBoxForBBTree1DQuadratic(double arcDetEps=1e-12) const throw(INTERP_KERNEL::Exception);
    void changeOrientationOfCells() throw(INTERP_KERNEL::Exception);
    int split2DCells(const DataArrayInt *desc, const DataArrayInt *descI, const DataArrayInt *subNodesInSeg, const DataArrayInt *subNodesInSegI, const DataArrayInt *midOpt=0, const DataArrayInt *midOptI=0) throw(INTERP_KERNEL::Exception);
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

      // serialization
      static PyObject *___new___(PyObject *cls, PyObject *args) throw(INTERP_KERNEL::Exception)
      {
        return NewMethWrapCallInitOnlyIfEmptyDictInInput(cls,args,"MEDCouplingUMesh");
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

      static MEDCouplingUMesh *Build1DMeshFromCoords(DataArrayDouble *da) throw(INTERP_KERNEL::Exception)
      {
        MCAuto<MEDCouplingUMesh> ret(MEDCouplingUMesh::Build1DMeshFromCoords(da));
        return ret.retn();
      }
      
      PyObject *getAllGeoTypesSorted() const throw(INTERP_KERNEL::Exception)
      {
        std::vector<INTERP_KERNEL::NormalizedCellType> result=self->getAllGeoTypesSorted();
        std::vector<INTERP_KERNEL::NormalizedCellType>::const_iterator iL=result.begin();
        PyObject *res=PyList_New(result.size());
        for(int i=0;iL!=result.end(); i++, iL++)
          PyList_SetItem(res,i,PyInt_FromLong(*iL));
        return res;
      }
      
      void setPartOfMySelf(PyObject *li, const MEDCouplingUMesh& otherOnSameCoordsThanThis) throw(INTERP_KERNEL::Exception)
      {
        int sw;
        int singleVal;
        std::vector<int> multiVal;
        std::pair<int, std::pair<int,int> > slic;
        MEDCoupling::DataArrayInt *daIntTyypp=0;
        int nbc=self->getNumberOfCells();
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
        MEDCoupling::DataArrayInt *daIntTyypp=0;
        int nbc=self->getNumberOfCells();
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
            throw INTERP_KERNEL::Exception("MEDCouplingUMesh::__setitem__ : unrecognized type in input ! Possibilities are : int, list or tuple of int, slice, DataArrayInt instance !");
          }
      }

      void insertNextCell(INTERP_KERNEL::NormalizedCellType type, int size, PyObject *li) throw(INTERP_KERNEL::Exception)
      {
        int szArr,sw,iTypppArr;
        std::vector<int> stdvecTyyppArr;
        const int *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
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
        const int *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
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
      
      static PyObject *ComputeSpreadZoneGraduallyFromSeed(PyObject *seed, const DataArrayInt *arrIn, const DataArrayInt *arrIndxIn, int nbOfDepthPeeling=-1) throw(INTERP_KERNEL::Exception)
      {
        int szArr,sw,iTypppArr;
        std::vector<int> stdvecTyyppArr;
        const int *seedPtr=convertIntStarLikePyObjToCppIntStar(seed,sw,szArr,iTypppArr,stdvecTyyppArr);
        int nbOfDepthPeelingPerformed=0;
        DataArrayInt *ret0=MEDCouplingUMesh::ComputeSpreadZoneGraduallyFromSeed(seedPtr,seedPtr+szArr,arrIn,arrIndxIn,nbOfDepthPeeling,nbOfDepthPeelingPerformed);
        PyObject *res=PyTuple_New(2);
        PyTuple_SetItem(res,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(res,1,PyInt_FromLong(nbOfDepthPeelingPerformed));
        return res;
      }

      static PyObject *FindCommonCellsAlg(int compType, int startCellId, const DataArrayInt *nodal, const DataArrayInt *nodalI, const DataArrayInt *revNodal, const DataArrayInt *revNodalI) throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *v0=0,*v1=0;
        MEDCouplingUMesh::FindCommonCellsAlg(compType,startCellId,nodal,nodalI,revNodal,revNodalI,v0,v1);
        PyObject *res = PyList_New(2);
        PyList_SetItem(res,0,SWIG_NewPointerObj(SWIG_as_voidptr(v0),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyList_SetItem(res,1,SWIG_NewPointerObj(SWIG_as_voidptr(v1),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
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
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_MEDCoupling__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *tetrahedrize(int policy) throw(INTERP_KERNEL::Exception)
      {
        int ret2(-1);
        DataArrayInt *ret1(0);
        MEDCoupling1SGTUMesh *ret0(self->tetrahedrize(policy,ret1,ret2));
        PyObject *ret=PyTuple_New(3);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_MEDCoupling__MEDCoupling1SGTUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,PyInt_FromLong(ret2));
        return ret;
      }
      
      PyObject *checkButterflyCells(double eps=1e-12) throw(INTERP_KERNEL::Exception)
      {
        std::vector<int> cells;
        self->checkButterflyCells(cells,eps);
        DataArrayInt *ret=DataArrayInt::New();
        ret->alloc((int)cells.size(),1);
        std::copy(cells.begin(),cells.end(),ret->getPointer());
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 );
      }

      PyObject *splitByType() const throw(INTERP_KERNEL::Exception)
      {
        std::vector<MEDCouplingUMesh *> ms=self->splitByType();
        int sz=ms.size();
        PyObject *ret = PyList_New(sz);
        for(int i=0;i<sz;i++)
          PyList_SetItem(ret,i,SWIG_NewPointerObj(SWIG_as_voidptr(ms[i]),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *partitionBySpreadZone() const throw(INTERP_KERNEL::Exception)
      {
        std::vector<DataArrayInt *> retCpp=self->partitionBySpreadZone();
        int sz=retCpp.size();
        PyObject *ret=PyList_New(sz);
        for(int i=0;i<sz;i++)
          PyList_SetItem(ret,i,SWIG_NewPointerObj(SWIG_as_voidptr(retCpp[i]),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      static PyObject *PartitionBySpreadZone(const DataArrayInt *arrIn, const DataArrayInt *arrIndxIn) throw(INTERP_KERNEL::Exception)
      {
        std::vector<DataArrayInt *> retCpp(MEDCouplingUMesh::PartitionBySpreadZone(arrIn,arrIndxIn));
        int sz=retCpp.size();
        PyObject *ret=PyList_New(sz);
        for(int i=0;i<sz;i++)
          PyList_SetItem(ret,i,SWIG_NewPointerObj(SWIG_as_voidptr(retCpp[i]),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *keepSpecifiedCells(INTERP_KERNEL::NormalizedCellType type, PyObject *ids) const throw(INTERP_KERNEL::Exception)
      {
        int size;
        INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(ids,&size);
        MEDCouplingUMesh *ret=self->keepSpecifiedCells(type,tmp,tmp+size);
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 );
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
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(tmp0),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(tmp1),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(tmp2),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *findCellIdsLyingOn(const MEDCouplingUMesh& otherDimM1OnSameCoords) const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *tmp0=0,*tmp1=0;
        self->findCellIdsLyingOn(otherDimM1OnSameCoords,tmp0,tmp1);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(tmp0),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(tmp1),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      void duplicateNodes(PyObject *li) throw(INTERP_KERNEL::Exception)
      {
        int sw;
        int singleVal;
        std::vector<int> multiVal;
        std::pair<int, std::pair<int,int> > slic;
        MEDCoupling::DataArrayInt *daIntTyypp=0;
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
            throw INTERP_KERNEL::Exception("MEDCouplingUMesh::duplicateNodes : unrecognized type entered, expected list of int, tuple of int or DataArrayInt !");
          }
      }

      void duplicateNodesInConn(PyObject *li, int offset) throw(INTERP_KERNEL::Exception)
      {
        int sw;
        int singleVal;
        std::vector<int> multiVal;
        std::pair<int, std::pair<int,int> > slic;
        MEDCoupling::DataArrayInt *daIntTyypp=0;
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
            throw INTERP_KERNEL::Exception("MEDCouplingUMesh::duplicateNodesInConn : unrecognized type entered, expected list of int, tuple of int or DataArrayInt !");
          }
      }

      PyObject *getLevArrPerCellTypes(PyObject *li) const throw(INTERP_KERNEL::Exception)
      {
        int sz;
        INTERP_KERNEL::AutoPtr<INTERP_KERNEL::NormalizedCellType> order=(INTERP_KERNEL::NormalizedCellType *)convertPyToNewIntArr2(li,&sz);
        DataArrayInt *tmp0,*tmp1=0;
        tmp0=self->getLevArrPerCellTypes(order,(INTERP_KERNEL::NormalizedCellType *)order+sz,tmp1);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(tmp0),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(tmp1),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *convertNodalConnectivityToDynamicGeoTypeMesh() const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *ret0=0,*ret1=0;
        self->convertNodalConnectivityToDynamicGeoTypeMesh(ret0,ret1);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      static PyObject *AggregateSortedByTypeMeshesOnSameCoords(PyObject *ms) throw(INTERP_KERNEL::Exception)
      {
        std::vector<const MEDCoupling::MEDCouplingUMesh *> meshes;
        convertFromPyObjVectorOfObj<const MEDCoupling::MEDCouplingUMesh *>(ms,SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh,"MEDCouplingUMesh",meshes);
        DataArrayInt *ret1=0,*ret2=0;
        MEDCouplingUMesh *ret0=MEDCouplingUMesh::AggregateSortedByTypeMeshesOnSameCoords(meshes,ret1,ret2);
        PyObject *ret=PyTuple_New(3);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(ret2),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      static PyObject *MergeUMeshesOnSameCoords(PyObject *ms) throw(INTERP_KERNEL::Exception)
      {
        std::vector<const MEDCoupling::MEDCouplingUMesh *> meshes;
        convertFromPyObjVectorOfObj<const MEDCoupling::MEDCouplingUMesh *>(ms,SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh,"MEDCouplingUMesh",meshes);
        MEDCouplingUMesh *ret=MEDCouplingUMesh::MergeUMeshesOnSameCoords(meshes);
        return convertMesh(ret, SWIG_POINTER_OWN | 0 );
      }

      static PyObject *FuseUMeshesOnSameCoords(PyObject *ms, int compType) throw(INTERP_KERNEL::Exception)
      {
        int sz;
        std::vector<const MEDCouplingUMesh *> meshes;
        convertFromPyObjVectorOfObj<const MEDCoupling::MEDCouplingUMesh *>(ms,SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh,"MEDCouplingUMesh",meshes);
        std::vector<DataArrayInt *> corr;
        MEDCouplingUMesh *um=MEDCouplingUMesh::FuseUMeshesOnSameCoords(meshes,compType,corr);
        sz=corr.size();
        PyObject *ret1=PyList_New(sz);
        for(int i=0;i<sz;i++)
          PyList_SetItem(ret1,i,SWIG_NewPointerObj(SWIG_as_voidptr(corr[i]),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyObject *ret=PyList_New(2);
        PyList_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(um),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyList_SetItem(ret,1,ret1);
        return ret;
      }

      static void PutUMeshesOnSameAggregatedCoords(PyObject *ms) throw(INTERP_KERNEL::Exception)
      {
        std::vector<MEDCouplingUMesh *> meshes;
        convertFromPyObjVectorOfObj<MEDCoupling::MEDCouplingUMesh *>(ms,SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh,"MEDCouplingUMesh",meshes);
        MEDCouplingUMesh::PutUMeshesOnSameAggregatedCoords(meshes);
      }

      static void MergeNodesOnUMeshesSharingSameCoords(PyObject *ms, double eps) throw(INTERP_KERNEL::Exception)
      {
        std::vector<MEDCouplingUMesh *> meshes;
        convertFromPyObjVectorOfObj<MEDCoupling::MEDCouplingUMesh *>(ms,SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh,"MEDCouplingUMesh",meshes);
        MEDCouplingUMesh::MergeNodesOnUMeshesSharingSameCoords(meshes,eps);
      }

      static bool RemoveIdsFromIndexedArrays(PyObject *li, DataArrayInt *arr, DataArrayInt *arrIndx, int offsetForRemoval=0) throw(INTERP_KERNEL::Exception)
      {
        int sw;
        int singleVal;
        std::vector<int> multiVal;
        std::pair<int, std::pair<int,int> > slic;
        MEDCoupling::DataArrayInt *daIntTyypp=0;
        if(!arrIndx)
          throw INTERP_KERNEL::Exception("MEDCouplingUMesh::RemoveIdsFromIndexedArrays : null pointer as arrIndex !");
        convertIntStarOrSliceLikePyObjToCpp(li,arrIndx->getNumberOfTuples()-1,sw,singleVal,multiVal,slic,daIntTyypp);
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
        MEDCoupling::DataArrayInt *daIntTyypp=0;
        if(!arrIndxIn)
          throw INTERP_KERNEL::Exception("MEDCouplingUMesh::ExtractFromIndexedArrays : null pointer as arrIndxIn !");
        convertIntStarOrSliceLikePyObjToCpp(li,arrIndxIn->getNumberOfTuples()-1,sw,singleVal,multiVal,slic,daIntTyypp);
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
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(arrOut),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(arrIndexOut),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      static PyObject *ExtractFromIndexedArraysSlice(int strt, int stp, int step, const DataArrayInt *arrIn, const DataArrayInt *arrIndxIn) throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *arrOut=0,*arrIndexOut=0;
        MEDCouplingUMesh::ExtractFromIndexedArraysSlice(strt,stp,step,arrIn,arrIndxIn,arrOut,arrIndexOut);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(arrOut),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(arrIndexOut),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      static PyObject *ExtractFromIndexedArraysSlice(PyObject *slic, const DataArrayInt *arrIn, const DataArrayInt *arrIndxIn) throw(INTERP_KERNEL::Exception)
      {
        if(!PySlice_Check(slic))
          throw INTERP_KERNEL::Exception("ExtractFromIndexedArraysSlice (wrap) : the first param is not a pyslice !");
        Py_ssize_t strt=2,stp=2,step=2;
        if(!arrIndxIn)
          throw INTERP_KERNEL::Exception("ExtractFromIndexedArraysSlice (wrap) : last array is null !");
        arrIndxIn->checkAllocated();
        if(arrIndxIn->getNumberOfComponents()!=1)
          throw INTERP_KERNEL::Exception("ExtractFromIndexedArraysSlice (wrap) : number of components of last argument must be equal to one !");
        GetIndicesOfSlice(slic,arrIndxIn->getNumberOfTuples(),&strt,&stp,&step,"ExtractFromIndexedArraysSlice (wrap) : Invalid slice regarding nb of elements !");
        DataArrayInt *arrOut=0,*arrIndexOut=0;
        MEDCouplingUMesh::ExtractFromIndexedArraysSlice(strt,stp,step,arrIn,arrIndxIn,arrOut,arrIndexOut);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(arrOut),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(arrIndexOut),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
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
        MEDCoupling::DataArrayInt *daIntTyypp=0;
        if(!arrIndxIn)
          throw INTERP_KERNEL::Exception("MEDCouplingUMesh::SetPartOfIndexedArrays : null pointer as arrIndex !");
        convertIntStarOrSliceLikePyObjToCpp(li,arrIndxIn->getNumberOfTuples()-1,sw,singleVal,multiVal,slic,daIntTyypp);
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
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(arrOut),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(arrIndexOut),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      static void SetPartOfIndexedArraysSameIdx(PyObject *li, DataArrayInt *arrIn, const DataArrayInt *arrIndxIn,
                                                const DataArrayInt *srcArr, const DataArrayInt *srcArrIndex) throw(INTERP_KERNEL::Exception)
      {
        int sw;
        int singleVal;
        std::vector<int> multiVal;
        std::pair<int, std::pair<int,int> > slic;
        MEDCoupling::DataArrayInt *daIntTyypp=0;
        if(!arrIndxIn)
          throw INTERP_KERNEL::Exception("MEDCouplingUMesh::SetPartOfIndexedArraysSameIdx : null pointer as arrIndex !");
        convertIntStarOrSliceLikePyObjToCpp(li,arrIndxIn->getNumberOfTuples()-1,sw,singleVal,multiVal,slic,daIntTyypp);
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
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 );
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
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 );
      }

      PyObject *getFastAveragePlaneOfThis() const throw(INTERP_KERNEL::Exception)
      {
        double vec[3];
        double pos[3];
        self->getFastAveragePlaneOfThis(vec,pos);
        double vals[6];
        std::copy(vec,vec+3,vals);
        std::copy(pos,pos+3,vals+3);
        return convertDblArrToPyListOfTuple<double>(vals,3,2);
      }
      
      static MEDCouplingUMesh *MergeUMeshes(PyObject *li) throw(INTERP_KERNEL::Exception)
      {
        std::vector<const MEDCoupling::MEDCouplingUMesh *> tmp;
        convertFromPyObjVectorOfObj<const MEDCoupling::MEDCouplingUMesh *>(li,SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh,"MEDCouplingUMesh",tmp);
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
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *areCellsIncludedInPolicy7(const MEDCouplingUMesh *other) const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *ret1;
        bool ret0=self->areCellsIncludedInPolicy7(other,ret1);
        PyObject *ret=PyTuple_New(2);
        PyObject *ret0Py=ret0?Py_True:Py_False;
        Py_XINCREF(ret0Py);
        PyTuple_SetItem(ret,0,ret0Py);
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *explode3DMeshTo1D() const throw(INTERP_KERNEL::Exception)
      {
        MCAuto<DataArrayInt> d0=DataArrayInt::New();
        MCAuto<DataArrayInt> d1=DataArrayInt::New();
        MCAuto<DataArrayInt> d2=DataArrayInt::New();
        MCAuto<DataArrayInt> d3=DataArrayInt::New();
        MEDCouplingUMesh *m=self->explode3DMeshTo1D(d0,d1,d2,d3);
        PyObject *ret=PyTuple_New(5);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(m),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(d0.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(d1.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,3,SWIG_NewPointerObj(SWIG_as_voidptr(d2.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,4,SWIG_NewPointerObj(SWIG_as_voidptr(d3.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *explodeIntoEdges() const throw(INTERP_KERNEL::Exception)
      {
        MCAuto<DataArrayInt> desc,descIndex,revDesc,revDescIndx;
        MCAuto<MEDCouplingUMesh> m(self->explodeIntoEdges(desc,descIndex,revDesc,revDescIndx));
        PyObject *ret=PyTuple_New(5);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(m.retn()),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(desc.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(descIndex.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,3,SWIG_NewPointerObj(SWIG_as_voidptr(revDesc.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,4,SWIG_NewPointerObj(SWIG_as_voidptr(revDescIndx.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *explodeMeshIntoMicroEdges() const throw(INTERP_KERNEL::Exception)
      {
        MCAuto<DataArrayInt> d0=DataArrayInt::New();
        MCAuto<DataArrayInt> d1=DataArrayInt::New();
        MCAuto<DataArrayInt> d2=DataArrayInt::New();
        MCAuto<DataArrayInt> d3=DataArrayInt::New();
        MEDCouplingUMesh *m=self->explodeMeshIntoMicroEdges(d0,d1,d2,d3);
        PyObject *ret=PyTuple_New(5);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(m),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(d0.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(d1.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,3,SWIG_NewPointerObj(SWIG_as_voidptr(d2.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,4,SWIG_NewPointerObj(SWIG_as_voidptr(d3.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *buildDescendingConnectivity() const throw(INTERP_KERNEL::Exception)
      {
        MCAuto<DataArrayInt> d0=DataArrayInt::New();
        MCAuto<DataArrayInt> d1=DataArrayInt::New();
        MCAuto<DataArrayInt> d2=DataArrayInt::New();
        MCAuto<DataArrayInt> d3=DataArrayInt::New();
        MEDCouplingUMesh *m=self->buildDescendingConnectivity(d0,d1,d2,d3);
        PyObject *ret=PyTuple_New(5);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(m),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(d0.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(d1.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,3,SWIG_NewPointerObj(SWIG_as_voidptr(d2.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,4,SWIG_NewPointerObj(SWIG_as_voidptr(d3.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *buildDescendingConnectivity2() const throw(INTERP_KERNEL::Exception)
      {
        MCAuto<DataArrayInt> d0=DataArrayInt::New();
        MCAuto<DataArrayInt> d1=DataArrayInt::New();
        MCAuto<DataArrayInt> d2=DataArrayInt::New();
        MCAuto<DataArrayInt> d3=DataArrayInt::New();
        MEDCouplingUMesh *m=self->buildDescendingConnectivity2(d0,d1,d2,d3);
        PyObject *ret=PyTuple_New(5);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(m),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(d0.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(d1.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,3,SWIG_NewPointerObj(SWIG_as_voidptr(d2.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,4,SWIG_NewPointerObj(SWIG_as_voidptr(d3.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }
      
      PyObject *computeNeighborsOfCells() const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *neighbors=0,*neighborsIdx=0;
        self->computeNeighborsOfCells(neighbors,neighborsIdx);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(neighbors),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(neighborsIdx),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *computeNeighborsOfNodes() const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *neighbors=0,*neighborsIdx=0;
        self->computeNeighborsOfNodes(neighbors,neighborsIdx);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(neighbors),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(neighborsIdx),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *computeEnlargedNeighborsOfNodes() const throw(INTERP_KERNEL::Exception)
      {
        MCAuto<DataArrayInt> neighbors,neighborsIdx;
        self->computeEnlargedNeighborsOfNodes(neighbors,neighborsIdx);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(neighbors.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(neighborsIdx.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }
      
      PyObject *computeCellNeighborhoodFromNodesOne(const DataArrayInt *nodeNeigh, const DataArrayInt *nodeNeighI) const throw(INTERP_KERNEL::Exception)
      {
        MCAuto<DataArrayInt> cellNeigh,cellNeighIndex;
        self->computeCellNeighborhoodFromNodesOne(nodeNeigh,nodeNeighI,cellNeigh,cellNeighIndex);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(cellNeigh.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(cellNeighIndex.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }
      
      static PyObject *ComputeNeighborsOfCellsAdv(const DataArrayInt *desc, const DataArrayInt *descI, const DataArrayInt *revDesc, const DataArrayInt *revDescI) throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *neighbors=0,*neighborsIdx=0;
        MEDCouplingUMesh::ComputeNeighborsOfCellsAdv(desc,descI,revDesc,revDescI,neighbors,neighborsIdx);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(neighbors),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(neighborsIdx),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *emulateMEDMEMBDC(const MEDCouplingUMesh *nM1LevMesh)
      {
        MCAuto<DataArrayInt> d0=DataArrayInt::New();
        MCAuto<DataArrayInt> d1=DataArrayInt::New();
        DataArrayInt *d2,*d3,*d4,*dd5;
        MEDCouplingUMesh *mOut=self->emulateMEDMEMBDC(nM1LevMesh,d0,d1,d2,d3,d4,dd5);
        PyObject *ret=PyTuple_New(7);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(mOut),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(d0.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(d1.retn()),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,3,SWIG_NewPointerObj(SWIG_as_voidptr(d2),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,4,SWIG_NewPointerObj(SWIG_as_voidptr(d3),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,5,SWIG_NewPointerObj(SWIG_as_voidptr(d4),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,6,SWIG_NewPointerObj(SWIG_as_voidptr(dd5),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
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

      static PyObject *Intersect2DMeshes(const MEDCouplingUMesh *m1, const MEDCouplingUMesh *m2, double eps) throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *cellNb1=0,*cellNb2=0;
        MEDCouplingUMesh *mret=MEDCouplingUMesh::Intersect2DMeshes(m1,m2,eps,cellNb1,cellNb2);
        PyObject *ret=PyTuple_New(3);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(mret),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(cellNb1),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(cellNb2),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      static PyObject *Intersect2DMeshWith1DLine(const MEDCouplingUMesh *mesh2D, const MEDCouplingUMesh *mesh1D, double eps) throw(INTERP_KERNEL::Exception)
      {
        MEDCouplingUMesh *splitMesh2D(0),*splitMesh1D(0);
        DataArrayInt *cellIdInMesh2D(0),*cellIdInMesh1D(0);
        MEDCouplingUMesh::Intersect2DMeshWith1DLine(mesh2D,mesh1D,eps,splitMesh2D,splitMesh1D,cellIdInMesh2D,cellIdInMesh1D);
        PyObject *ret(PyTuple_New(4));
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(splitMesh2D),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(splitMesh1D),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(cellIdInMesh2D),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,3,SWIG_NewPointerObj(SWIG_as_voidptr(cellIdInMesh1D),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
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
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(cellIds),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
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
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(cellIds),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      MEDCouplingUMesh *clipSingle3DCellByPlane(PyObject *origin, PyObject *vec, double eps) const throw(INTERP_KERNEL::Exception)
      {
        double val,val2;
        DataArrayDouble *a,*a2;
        DataArrayDoubleTuple *aa,*aa2;
        std::vector<double> bb,bb2;
        int sw;
        const char msg[]="Python wrap of MEDCouplingUMesh::clipSingle3DCellByPlane : 1st paramater for origin.";
        const char msg2[]="Python wrap of MEDCouplingUMesh::clipSingle3DCellByPlane : 2nd paramater for vector.";
        const double *orig=convertObjToPossibleCpp5_Safe(origin,sw,val,a,aa,bb,msg,1,3,true);
        const double *vect=convertObjToPossibleCpp5_Safe(vec,sw,val2,a2,aa2,bb2,msg2,1,3,true);
        MCAuto<MEDCouplingUMesh> ret(self->clipSingle3DCellByPlane(orig,vect,eps));
        return ret.retn();
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
    void convertExtrudedPolyhedra() throw(INTERP_KERNEL::Exception);
    bool unPolyze() throw(INTERP_KERNEL::Exception);
    void simplifyPolyhedra(double eps) throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *buildSpreadZonesWithPoly() const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *buildExtrudedMesh(const MEDCouplingUMesh *mesh1D, int policy) throw(INTERP_KERNEL::Exception);
  };

  //== MEDCouplingUMesh End

  //== MEDCouplingMappedExtrudedMesh

  class MEDCouplingMappedExtrudedMesh : public MEDCoupling::MEDCouplingMesh
  {
  public:
    static MEDCouplingMappedExtrudedMesh *New(const MEDCouplingUMesh *mesh3D, const MEDCouplingUMesh *mesh2D, int cell2DId) throw(INTERP_KERNEL::Exception);
    static MEDCouplingMappedExtrudedMesh *New(const MEDCouplingCMesh *mesh3D) throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *build3DUnstructuredMesh() const throw(INTERP_KERNEL::Exception);
    int get2DCellIdForExtrusion() const;
    %extend {
      MEDCouplingMappedExtrudedMesh(const MEDCouplingUMesh *mesh3D, const MEDCouplingUMesh *mesh2D, int cell2DId) throw(INTERP_KERNEL::Exception)
      {
        return MEDCouplingMappedExtrudedMesh::New(mesh3D,mesh2D,cell2DId);
      }

      MEDCouplingMappedExtrudedMesh(const MEDCouplingCMesh *mesh3D) throw(INTERP_KERNEL::Exception)
      {
        return MEDCouplingMappedExtrudedMesh::New(mesh3D);
      }

      MEDCouplingMappedExtrudedMesh()
      {
        return MEDCouplingMappedExtrudedMesh::New();
      }

      static PyObject *___new___(PyObject *cls, PyObject *args) throw(INTERP_KERNEL::Exception)
      {
        return NewMethWrapCallInitOnlyIfEmptyDictInInput(cls,args,"MEDCouplingMappedExtrudedMesh");
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
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 );
      } 
    }
  };

  //== MEDCouplingMappedExtrudedMesh End

  class MEDCoupling1GTUMesh : public MEDCoupling::MEDCouplingPointSet
  {
  public:
    static MEDCoupling1GTUMesh *New(const std::string& name, INTERP_KERNEL::NormalizedCellType type) throw(INTERP_KERNEL::Exception);
    static MEDCoupling1GTUMesh *New(const MEDCouplingUMesh *m) throw(INTERP_KERNEL::Exception);
    INTERP_KERNEL::NormalizedCellType getCellModelEnum() const throw(INTERP_KERNEL::Exception);
    int getNodalConnectivityLength() const throw(INTERP_KERNEL::Exception);
    virtual void allocateCells(int nbOfCells=0) throw(INTERP_KERNEL::Exception);
    virtual void checkConsistencyOfConnectivity() const throw(INTERP_KERNEL::Exception);
    %extend
    {
      virtual void insertNextCell(PyObject *li) throw(INTERP_KERNEL::Exception)
      {
        int szArr,sw,iTypppArr;
        std::vector<int> stdvecTyyppArr;
        const int *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
        self->insertNextCell(tmp,tmp+szArr);
      }

      virtual DataArrayInt *getNodalConnectivity() const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *ret=self->getNodalConnectivity();
        if(ret) ret->incrRef();
        return ret;
      }
      
      static MEDCouplingUMesh *AggregateOnSameCoordsToUMesh(PyObject *li) throw(INTERP_KERNEL::Exception)
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
    static MEDCoupling1SGTUMesh *New(const std::string& name, INTERP_KERNEL::NormalizedCellType type) throw(INTERP_KERNEL::Exception);
    static MEDCoupling1SGTUMesh *New(const MEDCouplingUMesh *m) throw(INTERP_KERNEL::Exception);
    void setNodalConnectivity(DataArrayInt *nodalConn) throw(INTERP_KERNEL::Exception);
    int getNumberOfNodesPerCell() const throw(INTERP_KERNEL::Exception);
    static MEDCoupling1SGTUMesh *Merge1SGTUMeshes(const MEDCoupling1SGTUMesh *mesh1, const MEDCoupling1SGTUMesh *mesh2) throw(INTERP_KERNEL::Exception);
    MEDCoupling1SGTUMesh *buildSetInstanceFromThis(int spaceDim) const throw(INTERP_KERNEL::Exception);
    MEDCoupling1GTUMesh *computeDualMesh() const throw(INTERP_KERNEL::Exception);
    MEDCoupling1SGTUMesh *explodeEachHexa8To6Quad4() const throw(INTERP_KERNEL::Exception);
    DataArrayInt *sortHexa8EachOther() throw(INTERP_KERNEL::Exception);
    %extend
    {
      MEDCoupling1SGTUMesh()
      {
        return MEDCoupling1SGTUMesh::New();
      }

      MEDCoupling1SGTUMesh(const std::string& name, INTERP_KERNEL::NormalizedCellType type) throw(INTERP_KERNEL::Exception)
      {
        return MEDCoupling1SGTUMesh::New(name,type);
      }

      MEDCoupling1SGTUMesh(const MEDCouplingUMesh *m) throw(INTERP_KERNEL::Exception)
      {
        return MEDCoupling1SGTUMesh::New(m);
      }

      static PyObject *___new___(PyObject *cls, PyObject *args) throw(INTERP_KERNEL::Exception)
      {
        return NewMethWrapCallInitOnlyIfEmptyDictInInput(cls,args,"MEDCoupling1SGTUMesh");
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

      PyObject *structurizeMe(double eps=1e-12) const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *cellPerm(0),*nodePerm(0);
        MEDCouplingCMesh *retCpp(self->structurizeMe(cellPerm,nodePerm,eps));
        PyObject *ret(PyTuple_New(3));
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(retCpp),SWIGTYPE_p_MEDCoupling__MEDCouplingCMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(cellPerm),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(nodePerm),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      static MEDCoupling1SGTUMesh *Merge1SGTUMeshes(PyObject *li) throw(INTERP_KERNEL::Exception)
      {
        std::vector<const MEDCoupling::MEDCoupling1SGTUMesh *> tmp;
        convertFromPyObjVectorOfObj<const MEDCoupling::MEDCoupling1SGTUMesh *>(li,SWIGTYPE_p_MEDCoupling__MEDCoupling1SGTUMesh,"MEDCoupling1SGTUMesh",tmp);
        return MEDCoupling1SGTUMesh::Merge1SGTUMeshes(tmp);
      }
      
      static MEDCoupling1SGTUMesh *Merge1SGTUMeshesOnSameCoords(PyObject *li) throw(INTERP_KERNEL::Exception)
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
    static MEDCoupling1DGTUMesh *New(const std::string& name, INTERP_KERNEL::NormalizedCellType type) throw(INTERP_KERNEL::Exception);
    static MEDCoupling1DGTUMesh *New(const MEDCouplingUMesh *m) throw(INTERP_KERNEL::Exception);
    void setNodalConnectivity(DataArrayInt *nodalConn, DataArrayInt *nodalConnIndex) throw(INTERP_KERNEL::Exception);
    MEDCoupling1DGTUMesh *buildSetInstanceFromThis(int spaceDim) const throw(INTERP_KERNEL::Exception);
    bool isPacked() const throw(INTERP_KERNEL::Exception);
    %extend
    {
      MEDCoupling1DGTUMesh()
      {
        return MEDCoupling1DGTUMesh::New();
      }
      MEDCoupling1DGTUMesh(const std::string& name, INTERP_KERNEL::NormalizedCellType type) throw(INTERP_KERNEL::Exception)
      {
        return MEDCoupling1DGTUMesh::New(name,type);
      }

      MEDCoupling1DGTUMesh(const MEDCouplingUMesh *m) throw(INTERP_KERNEL::Exception)
      {
        return MEDCoupling1DGTUMesh::New(m);
      }

      static PyObject *___new___(PyObject *cls, PyObject *args) throw(INTERP_KERNEL::Exception)
      {
        return NewMethWrapCallInitOnlyIfEmptyDictInInput(cls,args,"MEDCoupling1DGTUMesh");
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

      DataArrayInt *getNodalConnectivityIndex() const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *ret=self->getNodalConnectivityIndex();
        if(ret) ret->incrRef();
        return ret;
      }

      PyObject *retrievePackedNodalConnectivity() const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *ret1=0,*ret2=0;
        bool ret0=self->retrievePackedNodalConnectivity(ret1,ret2);
        PyObject *ret0Py=ret0?Py_True:Py_False;
        Py_XINCREF(ret0Py);
        PyObject *ret=PyTuple_New(3);
        PyTuple_SetItem(ret,0,ret0Py);
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(ret2),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }
      
      PyObject *copyWithNodalConnectivityPacked() const throw(INTERP_KERNEL::Exception)
      {
        bool ret1;
        MEDCoupling1DGTUMesh *ret0=self->copyWithNodalConnectivityPacked(ret1);
        PyObject *ret=PyTuple_New(2);
        PyObject *ret1Py=ret1?Py_True:Py_False; Py_XINCREF(ret1Py);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_MEDCoupling__MEDCoupling1DGTUMesh, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,ret1Py);
        return ret;
      }

      static MEDCoupling1DGTUMesh *Merge1DGTUMeshes(PyObject *li) throw(INTERP_KERNEL::Exception)
      {
        std::vector<const MEDCoupling::MEDCoupling1DGTUMesh *> tmp;
        convertFromPyObjVectorOfObj<const MEDCoupling::MEDCoupling1DGTUMesh *>(li,SWIGTYPE_p_MEDCoupling__MEDCoupling1DGTUMesh,"MEDCoupling1DGTUMesh",tmp);
        return MEDCoupling1DGTUMesh::Merge1DGTUMeshes(tmp);
      }
      
      static MEDCoupling1DGTUMesh *Merge1DGTUMeshesOnSameCoords(PyObject *li) throw(INTERP_KERNEL::Exception)
      {
        std::vector<const MEDCoupling::MEDCoupling1DGTUMesh *> tmp;
        convertFromPyObjVectorOfObj<const MEDCoupling::MEDCoupling1DGTUMesh *>(li,SWIGTYPE_p_MEDCoupling__MEDCoupling1DGTUMesh,"MEDCoupling1DGTUMesh",tmp);
        return MEDCoupling1DGTUMesh::Merge1DGTUMeshesOnSameCoords(tmp);
      }
      
      static DataArrayInt *AggregateNodalConnAndShiftNodeIds(PyObject *li, const std::vector<int>& offsetInNodeIdsPerElt) throw(INTERP_KERNEL::Exception)
      {
        std::vector<const MEDCoupling::DataArrayInt *> tmp;
        convertFromPyObjVectorOfObj<const MEDCoupling::DataArrayInt *>(li,SWIGTYPE_p_MEDCoupling__DataArrayInt,"DataArrayInt",tmp);
        return MEDCoupling1DGTUMesh::AggregateNodalConnAndShiftNodeIds(tmp,offsetInNodeIdsPerElt);
      }
    }
  };

  //== MEDCoupling1DGTUMeshEnd

  class MEDCouplingStructuredMesh : public MEDCoupling::MEDCouplingMesh
  {
  public:
    int getCellIdFromPos(int i, int j, int k) const throw(INTERP_KERNEL::Exception);
    int getNodeIdFromPos(int i, int j, int k) const throw(INTERP_KERNEL::Exception);
    int getNumberOfCellsOfSubLevelMesh() const throw(INTERP_KERNEL::Exception);
    int getSpaceDimensionOnNodeStruct() const throw(INTERP_KERNEL::Exception);
    double computeSquareness() const throw(INTERP_KERNEL::Exception);
    virtual std::vector<int> getNodeGridStructure() const throw(INTERP_KERNEL::Exception);
    std::vector<int> getCellGridStructure() const throw(INTERP_KERNEL::Exception);
    MEDCoupling1SGTUMesh *build1SGTUnstructured() const throw(INTERP_KERNEL::Exception);
    std::vector<int> getLocationFromCellId(int cellId) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getLocationFromNodeId(int cellId) const throw(INTERP_KERNEL::Exception);
    static INTERP_KERNEL::NormalizedCellType GetGeoTypeGivenMeshDimension(int meshDim) throw(INTERP_KERNEL::Exception);
    MEDCoupling1SGTUMesh *build1SGTSubLevelMesh() const throw(INTERP_KERNEL::Exception);
    static int DeduceNumberOfGivenStructure(const std::vector<int>& st) throw(INTERP_KERNEL::Exception);
    static DataArrayInt *ComputeCornersGhost(const std::vector<int>& st, int ghostLev) throw(INTERP_KERNEL::Exception);
    static std::vector<int> GetSplitVectFromStruct(const std::vector<int>& strct) throw(INTERP_KERNEL::Exception);
    %extend
    {
      virtual MEDCouplingStructuredMesh *buildStructuredSubPart(PyObject *cellPart) const throw(INTERP_KERNEL::Exception)
      {
        int tmpp1=-1,tmpp2=-1;
        std::vector<int> tmp=fillArrayWithPyListInt2(cellPart,tmpp1,tmpp2);
        std::vector< std::pair<int,int> > inp;
        if(tmpp2==2)
          {
            inp.resize(tmpp1);
            for(int i=0;i<tmpp1;i++)
              { inp[i].first=tmp[2*i]; inp[i].second=tmp[2*i+1]; }
          }
        else if(tmpp2==1)
          {
            if(tmpp1%2!=0)
              throw INTERP_KERNEL::Exception("Wrap of MEDCouplingStructuredMesh.buildStructuredSubPart : invalid input size ! Must be even size !");
            inp.resize(tmpp1/2);
            for(int i=0;i<tmpp1/2;i++)
              { inp[i].first=tmp[2*i]; inp[i].second=tmp[2*i+1]; }
          }
        else
          throw INTERP_KERNEL::Exception("Wrap of MEDCouplingStructuredMesh.buildStructuredSubPart : invalid input size !");
        return self->buildStructuredSubPart(inp);
      }

      static DataArrayInt *BuildExplicitIdsFrom(PyObject *st, PyObject *part) throw(INTERP_KERNEL::Exception)
      {
        std::vector< std::pair<int,int> > inp;
        convertPyToVectorPairInt(part,inp);
        //
        int szArr,sw,iTypppArr;
        std::vector<int> stdvecTyyppArr;
        const int *tmp4=convertIntStarLikePyObjToCppIntStar(st,sw,szArr,iTypppArr,stdvecTyyppArr);
        std::vector<int> tmp5(tmp4,tmp4+szArr);
        //
        return MEDCouplingStructuredMesh::BuildExplicitIdsFrom(tmp5,inp);
      }

      static void MultiplyPartOf(const std::vector<int>& st, PyObject *part, double factor, DataArrayDouble *da) throw(INTERP_KERNEL::Exception)
      {
        std::vector< std::pair<int,int> > inp;
        convertPyToVectorPairInt(part,inp);
        MEDCouplingStructuredMesh::MultiplyPartOf(st,inp,factor,da);
      }

      static void MultiplyPartOfByGhost(const std::vector<int>& st, PyObject *part, int ghostSize, double factor, DataArrayDouble *da) throw(INTERP_KERNEL::Exception)
      {
        std::vector< std::pair<int,int> > inp;
        convertPyToVectorPairInt(part,inp);
        MEDCouplingStructuredMesh::MultiplyPartOfByGhost(st,inp,ghostSize,factor,da);
      }

      static PyObject *PutInGhostFormat(int ghostSize, const std::vector<int>& st, PyObject *part) throw(INTERP_KERNEL::Exception)
      {
        std::vector< std::pair<int,int> > inp;
        convertPyToVectorPairInt(part,inp);
        std::vector<int> stWithGhost;
        std::vector< std::pair<int,int> > partWithGhost;
        MEDCouplingStructuredMesh::PutInGhostFormat(ghostSize,st,inp,stWithGhost,partWithGhost);
        PyObject *ret(PyTuple_New(2));
        PyTuple_SetItem(ret,0,convertIntArrToPyList2(stWithGhost));
        PyTuple_SetItem(ret,1,convertFromVectorPairInt(partWithGhost));
        return ret;
      }

      static DataArrayDouble *ExtractFieldOfDoubleFrom(const std::vector<int>& st, const DataArrayDouble *fieldOfDbl, PyObject *partCompactFormat) throw(INTERP_KERNEL::Exception)
      {
        std::vector< std::pair<int,int> > inp;
        convertPyToVectorPairInt(partCompactFormat,inp);
        return MEDCouplingStructuredMesh::ExtractFieldOfDoubleFrom(st,fieldOfDbl,inp);
      }

      static void AssignPartOfFieldOfDoubleUsing(const std::vector<int>& st, DataArrayDouble *fieldOfDbl, PyObject *partCompactFormat, const DataArrayDouble *other) throw(INTERP_KERNEL::Exception)
      {
        std::vector< std::pair<int,int> > inp;
        convertPyToVectorPairInt(partCompactFormat,inp);
        MEDCouplingStructuredMesh::AssignPartOfFieldOfDoubleUsing(st,fieldOfDbl,inp,other);
      }

      static int DeduceNumberOfGivenRangeInCompactFrmt(PyObject *part) throw(INTERP_KERNEL::Exception)
      {
        std::vector< std::pair<int,int> > inp;
        convertPyToVectorPairInt(part,inp);
        return MEDCouplingStructuredMesh::DeduceNumberOfGivenRangeInCompactFrmt(inp);
      }

      static DataArrayInt *Build1GTNodalConnectivity(PyObject *li) throw(INTERP_KERNEL::Exception)
      {
        int szArr,sw,iTypppArr;
        std::vector<int> stdvecTyyppArr;
        const int *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
        return MEDCouplingStructuredMesh::Build1GTNodalConnectivity(tmp,tmp+szArr);
      }

      static DataArrayInt *Build1GTNodalConnectivityOfSubLevelMesh(PyObject *li) throw(INTERP_KERNEL::Exception)
      {
        int szArr,sw,iTypppArr;
        std::vector<int> stdvecTyyppArr;
        const int *tmp(convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr));
        return MEDCouplingStructuredMesh::Build1GTNodalConnectivityOfSubLevelMesh(tmp,tmp+szArr);
      }

      static std::vector<int> GetDimensionsFromCompactFrmt(PyObject *partCompactFormat) throw(INTERP_KERNEL::Exception)
      {
        std::vector< std::pair<int,int> > inp;
        convertPyToVectorPairInt(partCompactFormat,inp);
        return MEDCouplingStructuredMesh::GetDimensionsFromCompactFrmt(inp);
      }

      static PyObject *GetCompactFrmtFromDimensions(const std::vector<int>& dims) throw(INTERP_KERNEL::Exception)
      {
        std::vector< std::pair<int,int> > ret(MEDCouplingStructuredMesh::GetCompactFrmtFromDimensions(dims));
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

      static PyObject *IntersectRanges(PyObject *r1, PyObject *r2) throw(INTERP_KERNEL::Exception)
      {
        std::vector< std::pair<int,int> > r1Cpp,r2Cpp;
        convertPyToVectorPairInt(r1,r1Cpp);
        convertPyToVectorPairInt(r2,r2Cpp);
        std::vector< std::pair<int,int> > ret(MEDCouplingStructuredMesh::IntersectRanges(r1Cpp,r2Cpp));
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
        std::vector< std::pair<int,int> > r1Cpp,r2Cpp;
        convertPyToVectorPairInt(r1,r1Cpp);
        convertPyToVectorPairInt(r2,r2Cpp);
        return MEDCouplingStructuredMesh::AreRangesIntersect(r1Cpp,r2Cpp);
      }

      static PyObject *IsPartStructured(PyObject *li, PyObject *st) throw(INTERP_KERNEL::Exception)
      {
        int szArr,sw,iTypppArr;
        std::vector<int> stdvecTyyppArr;
        const int *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
        int szArr2,sw2,iTypppArr2;
        std::vector<int> stdvecTyyppArr2;
        const int *tmp2=convertIntStarLikePyObjToCppIntStar(st,sw2,szArr2,iTypppArr2,stdvecTyyppArr2);
        std::vector<int> tmp3(tmp2,tmp2+szArr2);
        std::vector< std::pair<int,int> > partCompactFormat;
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

      static PyObject *ChangeReferenceFromGlobalOfCompactFrmt(PyObject *bigInAbs, PyObject *partOfBigInAbs, bool check=true) throw(INTERP_KERNEL::Exception)
      {
        std::vector< std::pair<int,int> > param0,param1,ret;
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

      static PyObject *TranslateCompactFrmt(PyObject *part, const std::vector<int>& translation) throw(INTERP_KERNEL::Exception)
      {
        std::vector< std::pair<int,int> > param0;
        convertPyToVectorPairInt(part,param0);
        std::vector< std::pair<int,int> > ret(MEDCouplingStructuredMesh::TranslateCompactFrmt(param0,translation));
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

      static std::vector<int> FindTranslationFrom(PyObject *startingFrom, PyObject *goingTo) throw(INTERP_KERNEL::Exception)
      {
        std::vector< std::pair<int,int> > param0,param1;
        convertPyToVectorPairInt(startingFrom,param0);
        convertPyToVectorPairInt(goingTo,param1);
        return  MEDCouplingStructuredMesh::FindTranslationFrom(param0,param1);
      }

      static PyObject *ChangeReferenceToGlobalOfCompactFrmt(PyObject *bigInAbs, PyObject *partOfBigRelativeToBig, bool check=true) throw(INTERP_KERNEL::Exception)
      {
        std::vector< std::pair<int,int> > param0,param1,ret;
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
    static MEDCouplingCMesh *New() throw(INTERP_KERNEL::Exception);
    static MEDCouplingCMesh *New(const std::string& meshName) throw(INTERP_KERNEL::Exception);
    void setCoords(const DataArrayDouble *coordsX,
                   const DataArrayDouble *coordsY=0,
                   const DataArrayDouble *coordsZ=0) throw(INTERP_KERNEL::Exception);
    void setCoordsAt(int i, const DataArrayDouble *arr) throw(INTERP_KERNEL::Exception);
    MEDCouplingCurveLinearMesh *buildCurveLinear() const throw(INTERP_KERNEL::Exception);
    %extend {
      MEDCouplingCMesh() throw(INTERP_KERNEL::Exception)
      {
        return MEDCouplingCMesh::New();
      }
      MEDCouplingCMesh(const std::string& meshName) throw(INTERP_KERNEL::Exception)
      {
        return MEDCouplingCMesh::New(meshName);
      }
      // serialization
      static PyObject *___new___(PyObject *cls, PyObject *args) throw(INTERP_KERNEL::Exception)
      {
        return NewMethWrapCallInitOnlyIfEmptyDictInInput(cls,args,"MEDCouplingCMesh");
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

  //== MEDCouplingCMesh End

  //== MEDCouplingCurveLinearMesh

  class MEDCouplingCurveLinearMesh : public MEDCoupling::MEDCouplingStructuredMesh
  {
  public:
    static MEDCouplingCurveLinearMesh *New() throw(INTERP_KERNEL::Exception);
    static MEDCouplingCurveLinearMesh *New(const std::string& meshName) throw(INTERP_KERNEL::Exception);
    void setCoords(const DataArrayDouble *coords) throw(INTERP_KERNEL::Exception);
    %extend {
      MEDCouplingCurveLinearMesh() throw(INTERP_KERNEL::Exception)
      {
        return MEDCouplingCurveLinearMesh::New();
      }
      MEDCouplingCurveLinearMesh(const std::string& meshName) throw(INTERP_KERNEL::Exception)
      {
        return MEDCouplingCurveLinearMesh::New(meshName);
      }
      static PyObject *___new___(PyObject *cls, PyObject *args) throw(INTERP_KERNEL::Exception)
      {
        return NewMethWrapCallInitOnlyIfEmptyDictInInput(cls,args,"MEDCouplingCurveLinearMesh");
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
        const int *tmp=convertIntStarLikePyObjToCppIntStar(gridStruct,sw,szArr,iTypppArr,stdvecTyyppArr);
        self->setNodeGridStructure(tmp,tmp+szArr);
      }
    }
  };

  //== MEDCouplingCurveLinearMesh End

  //== MEDCouplingIMesh

  class MEDCouplingIMesh : public MEDCoupling::MEDCouplingStructuredMesh
  {
  public:
    static MEDCouplingIMesh *New() throw(INTERP_KERNEL::Exception);
    //
    void setSpaceDimension(int spaceDim) throw(INTERP_KERNEL::Exception);
    std::vector<int> getNodeStruct() const throw(INTERP_KERNEL::Exception);
    std::vector<double> getOrigin() const throw(INTERP_KERNEL::Exception);
    std::vector<double> getDXYZ() const throw(INTERP_KERNEL::Exception);
    void setAxisUnit(const std::string& unitName) throw(INTERP_KERNEL::Exception);
    std::string getAxisUnit() const throw(INTERP_KERNEL::Exception);
    double getMeasureOfAnyCell() const throw(INTERP_KERNEL::Exception);
    MEDCouplingCMesh *convertToCartesian() const throw(INTERP_KERNEL::Exception);
    void refineWithFactor(const std::vector<int>& factors) throw(INTERP_KERNEL::Exception);
    MEDCouplingIMesh *asSingleCell() const throw(INTERP_KERNEL::Exception);
    MEDCouplingIMesh *buildWithGhost(int ghostLev) const throw(INTERP_KERNEL::Exception);
    %extend
    {
      MEDCouplingIMesh()
      {
        return MEDCouplingIMesh::New();
      }
      static MEDCouplingIMesh *New(const std::string& meshName, int spaceDim, PyObject *nodeStrct, PyObject *origin, PyObject *dxyz) throw(INTERP_KERNEL::Exception)
      {
        static const char msg0[]="MEDCouplingIMesh::New : error on 'origin' parameter !";
        static const char msg1[]="MEDCouplingIMesh::New : error on 'dxyz' parameter !";
        const int *nodeStrctPtr(0);
        const double *originPtr(0),*dxyzPtr(0);
        int sw,sz,val0;
        std::vector<int> bb0;
        nodeStrctPtr=convertIntStarLikePyObjToCppIntStar(nodeStrct,sw,sz,val0,bb0);
        //
        double val,val2;
        std::vector<double> bb,bb2;
        int sz1,sz2;
        originPtr=convertObjToPossibleCpp5_SingleCompo(origin,sw,val,bb,msg0,false,sz1);
        dxyzPtr=convertObjToPossibleCpp5_SingleCompo(dxyz,sw,val2,bb2,msg1,false,sz2);
        //
        return MEDCouplingIMesh::New(meshName,spaceDim,nodeStrctPtr,nodeStrctPtr+sz,originPtr,originPtr+sz1,dxyzPtr,dxyzPtr+sz2);
      }

      MEDCouplingIMesh(const std::string& meshName, int spaceDim, PyObject *nodeStrct, PyObject *origin, PyObject *dxyz) throw(INTERP_KERNEL::Exception)
      {
        return MEDCoupling_MEDCouplingIMesh_New__SWIG_1(meshName,spaceDim,nodeStrct,origin,dxyz);
      }

      static PyObject *___new___(PyObject *cls, PyObject *args) throw(INTERP_KERNEL::Exception)
      {
        return NewMethWrapCallInitOnlyIfEmptyDictInInput(cls,args,"MEDCouplingIMesh");
      }

      void setNodeStruct(PyObject *nodeStrct) throw(INTERP_KERNEL::Exception)
      {
        int sw,sz,val0;
        std::vector<int> bb0;
        const int *nodeStrctPtr(convertIntStarLikePyObjToCppIntStar(nodeStrct,sw,sz,val0,bb0));
        self->setNodeStruct(nodeStrctPtr,nodeStrctPtr+sz);
      }

      void setOrigin(PyObject *origin) throw(INTERP_KERNEL::Exception)
      {
        static const char msg[]="MEDCouplingIMesh::setOrigin : invalid input 'origin' parameter ! integer, float, list/tuple of float, DataArrayDouble or DataArrayDoubleTuple supported !";
        double val;
        DataArrayDouble *a;
        DataArrayDoubleTuple *aa;
        std::vector<double> bb;
        int sw,nbTuples;
        const double *originPtr(convertObjToPossibleCpp5_SingleCompo(origin,sw,val,bb,msg,false,nbTuples));
        self->setOrigin(originPtr,originPtr+nbTuples);
      }
      
      void setDXYZ(PyObject *dxyz) throw(INTERP_KERNEL::Exception)
      {
        static const char msg[]="MEDCouplingIMesh::setDXYZ : invalid input 'dxyz' parameter ! integer, float, list/tuple of float, DataArrayDouble or DataArrayDoubleTuple supported !";
        double val;
        DataArrayDouble *a;
        DataArrayDoubleTuple *aa;
        std::vector<double> bb;
        int sw,nbTuples;
        const double *originPtr(convertObjToPossibleCpp5_SingleCompo(dxyz,sw,val,bb,msg,false,nbTuples));
        self->setDXYZ(originPtr,originPtr+nbTuples);
      }

      static void CondenseFineToCoarse(const std::vector<int>& coarseSt, const DataArrayDouble *fineDA, PyObject *fineLocInCoarse, const std::vector<int>& facts, DataArrayDouble *coarseDA) throw(INTERP_KERNEL::Exception)
      {
        std::vector< std::pair<int,int> > inp;
        convertPyToVectorPairInt(fineLocInCoarse,inp);
        MEDCouplingIMesh::CondenseFineToCoarse(coarseSt,fineDA,inp,facts,coarseDA);
      }

      static void CondenseFineToCoarseGhost(const std::vector<int>& coarseSt, const DataArrayDouble *fineDA, PyObject *fineLocInCoarse, const std::vector<int>& facts, DataArrayDouble *coarseDA, int ghostSize) throw(INTERP_KERNEL::Exception)
      {
        std::vector< std::pair<int,int> > inp;
        convertPyToVectorPairInt(fineLocInCoarse,inp);
        MEDCouplingIMesh::CondenseFineToCoarseGhost(coarseSt,fineDA,inp,facts,coarseDA,ghostSize);
      }

      static void SpreadCoarseToFine(const DataArrayDouble *coarseDA, const std::vector<int>& coarseSt, DataArrayDouble *fineDA, PyObject *fineLocInCoarse, const std::vector<int>& facts) throw(INTERP_KERNEL::Exception)
      {
        std::vector< std::pair<int,int> > inp;
        convertPyToVectorPairInt(fineLocInCoarse,inp);
        MEDCouplingIMesh::SpreadCoarseToFine(coarseDA,coarseSt,fineDA,inp,facts);
      }

      static void SpreadCoarseToFineGhost(const DataArrayDouble *coarseDA, const std::vector<int>& coarseSt, DataArrayDouble *fineDA, PyObject *fineLocInCoarse, const std::vector<int>& facts, int ghostSize) throw(INTERP_KERNEL::Exception)
      {
        std::vector< std::pair<int,int> > inp;
        convertPyToVectorPairInt(fineLocInCoarse,inp);
        MEDCouplingIMesh::SpreadCoarseToFineGhost(coarseDA,coarseSt,fineDA,inp,facts,ghostSize);
      }

      static void SpreadCoarseToFineGhostZone(const DataArrayDouble *coarseDA, const std::vector<int>& coarseSt, DataArrayDouble *fineDA, PyObject *fineLocInCoarse, const std::vector<int>& facts, int ghostSize) throw(INTERP_KERNEL::Exception)
      {
        std::vector< std::pair<int,int> > inp;
        convertPyToVectorPairInt(fineLocInCoarse,inp);
        MEDCouplingIMesh::SpreadCoarseToFineGhostZone(coarseDA,coarseSt,fineDA,inp,facts,ghostSize);
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

  //== MEDCouplingIMesh End

}

namespace MEDCoupling
{
  class MEDCouplingField : public MEDCoupling::RefCountObject, public MEDCoupling::TimeLabel
  {
  public:
    virtual void checkConsistencyLight() const throw(INTERP_KERNEL::Exception);
    virtual bool areCompatibleForMerge(const MEDCouplingField *other) const throw(INTERP_KERNEL::Exception);
    virtual void copyTinyStringsFrom(const MEDCouplingField *other) throw(INTERP_KERNEL::Exception);
    void setMesh(const MEDCoupling::MEDCouplingMesh *mesh) throw(INTERP_KERNEL::Exception);
    void setName(const char *name) throw(INTERP_KERNEL::Exception);
    std::string getDescription() const throw(INTERP_KERNEL::Exception);
    void setDescription(const char *desc) throw(INTERP_KERNEL::Exception);
    std::string getName() const throw(INTERP_KERNEL::Exception);
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

      PyObject *buildSubMeshData(PyObject *li) const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *ret1=0;
        MEDCouplingMesh *ret0=0;
        void *da=0;
        int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_MEDCoupling__DataArrayInt, 0 |  0 );
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
        PyList_SetItem(res,1,SWIG_NewPointerObj((void*)ret1,SWIGTYPE_p_MEDCoupling__DataArrayInt,SWIG_POINTER_OWN | 0));
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
          PyTuple_SetItem(res,1,SWIG_NewPointerObj((void*)ret1,SWIGTYPE_p_MEDCoupling__DataArrayInt,SWIG_POINTER_OWN | 0));
        else
          {
            PyObject *res1=PySlice_New(PyInt_FromLong(bb),PyInt_FromLong(ee),PyInt_FromLong(ss));
            PyTuple_SetItem(res,1,res1);
          }
        return res;
      }

      DataArrayInt *computeTupleIdsToSelectFromCellIds(PyObject *cellIds) const
      {
        int sw,sz(-1);
        int v0; std::vector<int> v1;
        const int *cellIdsBg(convertIntStarLikePyObjToCppIntStar(cellIds,sw,sz,v0,v1));
        return self->computeTupleIdsToSelectFromCellIds(cellIdsBg,cellIdsBg+sz);
      }

      void setGaussLocalizationOnCells(PyObject *li, const std::vector<double>& refCoo,
                                       const std::vector<double>& gsCoo, const std::vector<double>& wg) throw(INTERP_KERNEL::Exception)
      {
        void *da=0;
        int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_MEDCoupling__DataArrayInt, 0 |  0 );
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
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 );
      }
      
      int getNumberOfTuplesExpectedRegardingCode(PyObject *code, PyObject *idsPerType) const throw(INTERP_KERNEL::Exception)
      {
        std::vector<int> inp0;
        convertPyToNewIntArr4(code,1,3,inp0);
        std::vector<const DataArrayInt *> inp1;
        convertFromPyObjVectorOfObj<const MEDCoupling::DataArrayInt *>(idsPerType,SWIGTYPE_p_MEDCoupling__DataArrayInt,"DataArrayInt",inp1);
        return self->getNumberOfTuplesExpectedRegardingCode(inp0,inp1);
      }
    }
  };
  
  class MEDCouplingFieldTemplate : public MEDCoupling::MEDCouplingField
  {
  public:
    static MEDCouplingFieldTemplate *New(const MEDCouplingFieldDouble& f) throw(INTERP_KERNEL::Exception);
    static MEDCouplingFieldTemplate *New(const MEDCouplingFieldFloat& f) throw(INTERP_KERNEL::Exception);
    static MEDCouplingFieldTemplate *New(const MEDCouplingFieldInt& f) throw(INTERP_KERNEL::Exception);
    static MEDCouplingFieldTemplate *New(TypeOfField type);
    std::string simpleRepr() const throw(INTERP_KERNEL::Exception);
    std::string advancedRepr() const throw(INTERP_KERNEL::Exception);
    bool isEqual(const MEDCouplingFieldTemplate *other, double meshPrec) const throw(INTERP_KERNEL::Exception);
    bool isEqualWithoutConsideringStr(const MEDCouplingFieldTemplate *other, double meshPrec) const throw(INTERP_KERNEL::Exception);
    %extend
       {
         MEDCouplingFieldTemplate(const MEDCouplingFieldDouble& f) throw(INTERP_KERNEL::Exception)
         {
           return MEDCouplingFieldTemplate::New(f);
         }

         MEDCouplingFieldTemplate(const MEDCouplingFieldFloat& f) throw(INTERP_KERNEL::Exception)
         {
           return MEDCouplingFieldTemplate::New(f);
         }
         
         MEDCouplingFieldTemplate(const MEDCouplingFieldInt& f) throw(INTERP_KERNEL::Exception)
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

         PyObject *isEqualIfNotWhy(const MEDCouplingFieldTemplate *other, double meshPrec) const throw(INTERP_KERNEL::Exception)
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
  
  class MEDCouplingFieldInt;
  class MEDCouplingFieldFloat;
  
  class MEDCouplingFieldDouble : public MEDCoupling::MEDCouplingField
  {
  public:
    static MEDCouplingFieldDouble *New(TypeOfField type, TypeOfTimeDiscretization td=ONE_TIME);
    static MEDCouplingFieldDouble *New(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td=ONE_TIME);
    bool isEqual(const MEDCouplingFieldDouble *other, double meshPrec, double valsPrec) const throw(INTERP_KERNEL::Exception);
    bool isEqualWithoutConsideringStr(const MEDCouplingFieldDouble *other, double meshPrec, double valsPrec) const throw(INTERP_KERNEL::Exception);
    void setTimeUnit(const std::string& unit);
    std::string getTimeUnit() const;
    void synchronizeTimeWithSupport() throw(INTERP_KERNEL::Exception);
    void copyTinyAttrFrom(const MEDCouplingFieldDouble *other) throw(INTERP_KERNEL::Exception);
    void copyAllTinyAttrFrom(const MEDCouplingFieldDouble *other) throw(INTERP_KERNEL::Exception);
    std::string simpleRepr() const throw(INTERP_KERNEL::Exception);
    std::string advancedRepr() const throw(INTERP_KERNEL::Exception);
    std::string  writeVTK(const std::string& fileName, bool isBinary=true) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldInt *convertToIntField() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldFloat *convertToFloatField() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *clone(bool recDeepCpy) const;
    MEDCouplingFieldDouble *cloneWithMesh(bool recDeepCpy) const;
    MEDCouplingFieldDouble *deepCopy() const;
    MEDCouplingFieldDouble *buildNewTimeReprFromThis(TypeOfTimeDiscretization td, bool deepCopy) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *nodeToCellDiscretization() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *cellToNodeDiscretization() const throw(INTERP_KERNEL::Exception);
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
    void applyLin(double a, double b) throw(INTERP_KERNEL::Exception);
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
    void changeUnderlyingMesh(const MEDCouplingMesh *other, int levOfCheck, double precOnMesh, double eps=1e-15) throw(INTERP_KERNEL::Exception);
    void substractInPlaceDM(const MEDCouplingFieldDouble *f, int levOfCheck, double precOnMesh, double eps=1e-15) throw(INTERP_KERNEL::Exception);
    bool mergeNodes(double eps, double epsOnVals=1e-15) throw(INTERP_KERNEL::Exception);
    bool mergeNodesCenter(double eps, double epsOnVals=1e-15) throw(INTERP_KERNEL::Exception);
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
    void fillFromAnalytic(int nbOfComp, const std::string& func) throw(INTERP_KERNEL::Exception);
    void fillFromAnalyticCompo(int nbOfComp, const std::string& func) throw(INTERP_KERNEL::Exception);
    void fillFromAnalyticNamedCompo(int nbOfComp, const std::vector<std::string>& varsOrder, const std::string& func) throw(INTERP_KERNEL::Exception);
    void applyFunc(int nbOfComp, const std::string& func) throw(INTERP_KERNEL::Exception);
    void applyFuncCompo(int nbOfComp, const std::string& func) throw(INTERP_KERNEL::Exception);
    void applyFuncNamedCompo(int nbOfComp, const std::vector<std::string>& varsOrder, const std::string& func) throw(INTERP_KERNEL::Exception);
    void applyFunc(int nbOfComp, double val) throw(INTERP_KERNEL::Exception);
    void applyFunc(const std::string& func) throw(INTERP_KERNEL::Exception);
    void applyFuncFast32(const std::string& func) throw(INTERP_KERNEL::Exception);
    void applyFuncFast64(const std::string& func) throw(INTERP_KERNEL::Exception);
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
    DataArrayInt *findIdsInRange(double vmin, double vmax) const throw(INTERP_KERNEL::Exception);
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

      PyObject *isEqualIfNotWhy(const MEDCouplingFieldDouble *other, double meshPrec, double valsPrec) const throw(INTERP_KERNEL::Exception)
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
      
      MEDCouplingFieldDouble *voronoize(double eps) const throw(INTERP_KERNEL::Exception)
      {
        MCAuto<MEDCouplingFieldDouble> ret(self->voronoize(eps));
        return ret.retn();
      }

      MEDCouplingFieldDouble *convertQuadraticCellsToLinear() const throw(INTERP_KERNEL::Exception)
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
        int sw;
        const double *centerPtr=convertObjToPossibleCpp5_Safe(center,sw,val,a,aa,bb,msg,1,3,true);
        const double *vectorPtr=convertObjToPossibleCpp5_Safe(vector,sw,val2,a2,aa2,bb2,msg,1,3,true);
        return self->computeVectorFieldCyl(centerPtr,vectorPtr);
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
              PyTuple_SetItem(ret,i,SWIG_NewPointerObj(SWIG_as_voidptr(arrs[i]),SWIGTYPE_p_MEDCoupling__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
            else
              PyTuple_SetItem(ret,i,SWIG_NewPointerObj(SWIG_as_voidptr(0),SWIGTYPE_p_MEDCoupling__DataArrayDouble, 0 | 0 ));
          }
        return ret;
      }

      void setArrays(PyObject *ls) throw(INTERP_KERNEL::Exception)
      {
        std::vector<const DataArrayDouble *> tmp;
        convertFromPyObjVectorOfObj<const DataArrayDouble *>(ls,SWIGTYPE_p_MEDCoupling__DataArrayDouble,"DataArrayDouble",tmp);
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
        return convertDblArrToPyList<double>(res,sz);
      }

       PyObject *getValueOnPos(int i, int j, int k) const throw(INTERP_KERNEL::Exception)
       {
         int sz=self->getNumberOfComponents();
         INTERP_KERNEL::AutoPtr<double> res=new double[sz];
         self->getValueOnPos(i,j,k,res);
         return convertDblArrToPyList<double>(res,sz);
       }

      DataArrayDouble *getValueOnMulti(PyObject *locs) const throw(INTERP_KERNEL::Exception)
      {
        const MEDCouplingMesh *mesh(self->getMesh());
        if(!mesh)
          throw INTERP_KERNEL::Exception("Python wrap MEDCouplingFieldDouble::getValueOnMulti : lying on a null mesh !");
        //
        int sw,nbPts;
        double v0; MEDCoupling::DataArrayDouble *v1(0); MEDCoupling::DataArrayDoubleTuple *v2(0); std::vector<double> v3;
        const double *inp=convertObjToPossibleCpp5_Safe2(locs,sw,v0,v1,v2,v3,"wrap of MEDCouplingFieldDouble::getValueOnMulti",
                                                         mesh->getSpaceDimension(),true,nbPts);
        return self->getValueOnMulti(inp,nbPts);
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
        return convertDblArrToPyList<double>(res,sz);
      }

      void setValues(PyObject *li, PyObject *nbOfTuples=0, PyObject *nbOfComp=0) throw(INTERP_KERNEL::Exception)
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
        return convertDblArrToPyList<double>(tmp,sz);
      }
      PyObject *integral(bool isWAbs) const throw(INTERP_KERNEL::Exception)
      {
        int sz=self->getNumberOfComponents();
        INTERP_KERNEL::AutoPtr<double> tmp=new double[sz];
        self->integral(isWAbs,tmp);
        return convertDblArrToPyList<double>(tmp,sz);
      }
      PyObject *getWeightedAverageValue(bool isWAbs=true) const throw(INTERP_KERNEL::Exception)
      {
        int sz=self->getNumberOfComponents();
        INTERP_KERNEL::AutoPtr<double> tmp=new double[sz];
        self->getWeightedAverageValue(tmp,isWAbs);
        return convertDblArrToPyList<double>(tmp,sz);
      }
      PyObject *normL1() const throw(INTERP_KERNEL::Exception)
      {
        int sz=self->getNumberOfComponents();
        INTERP_KERNEL::AutoPtr<double> tmp=new double[sz];
        self->normL1(tmp);
        return convertDblArrToPyList<double>(tmp,sz);
      }
      PyObject *normL2() const throw(INTERP_KERNEL::Exception)
      {
        int sz=self->getNumberOfComponents();
        INTERP_KERNEL::AutoPtr<double> tmp=new double[sz];
        self->normL2(tmp);
        return convertDblArrToPyList<double>(tmp,sz);
      }
      void renumberCells(PyObject *li, bool check=true) throw(INTERP_KERNEL::Exception)
      {
        int szArr,sw,iTypppArr;
        std::vector<int> stdvecTyyppArr;
        const int *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
        self->renumberCells(tmp,check);
      }
      
      void renumberCellsWithoutMesh(PyObject *li, bool check=true) throw(INTERP_KERNEL::Exception)
      {
        int szArr,sw,iTypppArr;
        std::vector<int> stdvecTyyppArr;
        const int *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
        self->renumberCellsWithoutMesh(tmp,check);
      }
      
      void renumberNodes(PyObject *li, double eps=1e-15) throw(INTERP_KERNEL::Exception)
      {
        int szArr,sw,iTypppArr;
        std::vector<int> stdvecTyyppArr;
        const int *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
        self->renumberNodes(tmp,eps);
      }

      void renumberNodesWithoutMesh(PyObject *li, int newNbOfNodes, double eps=1e-15) throw(INTERP_KERNEL::Exception)
      {
        int szArr,sw,iTypppArr;
        std::vector<int> stdvecTyyppArr;
        const int *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
        self->renumberNodesWithoutMesh(tmp,newNbOfNodes,eps);
      }

      MEDCouplingFieldDouble *buildSubPart(PyObject *li) const throw(INTERP_KERNEL::Exception)
      {
        return fieldT_buildSubPart(self,li);
      }

      MEDCouplingFieldDouble *__getitem__(PyObject *li) const throw(INTERP_KERNEL::Exception)
      {
        return fieldT__getitem__(self,li);
      }

      PyObject *getMaxValue2() const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *tmp;
        double r1=self->getMaxValue2(tmp);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,PyFloat_FromDouble(r1));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(tmp),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }
      
      PyObject *getMinValue2() const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *tmp;
        double r1=self->getMinValue2(tmp);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,PyFloat_FromDouble(r1));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(tmp),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
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
        return MEDCoupling_MEDCouplingFieldDouble___add__Impl(self,obj);
      }

      MEDCouplingFieldDouble *__radd__(PyObject *obj) throw(INTERP_KERNEL::Exception)
      {
        return MEDCoupling_MEDCouplingFieldDouble___radd__Impl(self,obj);
      }

      MEDCouplingFieldDouble *__sub__(PyObject *obj) throw(INTERP_KERNEL::Exception)
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
        int sw;
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
              MCAuto<DataArrayDouble> aaa=DataArrayDouble::New(); aaa->useArray(&bb[0],false,CPP_DEALLOC,1,(int)bb.size());
              MCAuto<DataArrayDouble> ret=DataArrayDouble::Substract(self->getArray(),aaa);
              MCAuto<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(ret);
              return ret2.retn();
            }
          default:
            { throw INTERP_KERNEL::Exception(msg); }
          }
      }

      MEDCouplingFieldDouble *__rsub__(PyObject *obj) throw(INTERP_KERNEL::Exception)
      {
        return MEDCoupling_MEDCouplingFieldDouble___rsub__Impl(self,obj);
      }

      MEDCouplingFieldDouble *__mul__(PyObject *obj) throw(INTERP_KERNEL::Exception)
      {
        return MEDCoupling_MEDCouplingFieldDouble___mul__Impl(self,obj);
      }

      MEDCouplingFieldDouble *__rmul__(PyObject *obj) throw(INTERP_KERNEL::Exception)
      {
        return MEDCoupling_MEDCouplingFieldDouble___rmul__Impl(self,obj);
      }

      MEDCouplingFieldDouble *__div__(PyObject *obj) throw(INTERP_KERNEL::Exception)
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
        int sw;
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
              MCAuto<DataArrayDouble> aaa=DataArrayDouble::New(); aaa->useArray(&bb[0],false,CPP_DEALLOC,1,(int)bb.size());
              MCAuto<DataArrayDouble> ret=DataArrayDouble::Divide(self->getArray(),aaa);
              MCAuto<MEDCouplingFieldDouble> ret2=self->clone(false);
              ret2->setArray(ret);
              return ret2.retn();
            }
          default:
            { throw INTERP_KERNEL::Exception(msg); }
          }
      }

      MEDCouplingFieldDouble *__rdiv__(PyObject *obj) throw(INTERP_KERNEL::Exception)
      {
        return MEDCoupling_MEDCouplingFieldDouble___rdiv__Impl(self,obj);
      }

      MEDCouplingFieldDouble *__pow__(PyObject *obj) throw(INTERP_KERNEL::Exception)
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
        int sw;
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
              MCAuto<DataArrayDouble> aaa=DataArrayDouble::New(); aaa->useArray(&bb[0],false,CPP_DEALLOC,1,(int)bb.size());
              MCAuto<DataArrayDouble> ret=DataArrayDouble::Pow(self->getArray(),aaa);
              MCAuto<MEDCouplingFieldDouble> ret2=self->clone(false);
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
        int sw;
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
              MCAuto<DataArrayDouble> aaa=DataArrayDouble::New(); aaa->useArray(&bb[0],false,CPP_DEALLOC,1,(int)bb.size());
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
        int sw;
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
              MCAuto<DataArrayDouble> aaa=DataArrayDouble::New(); aaa->useArray(&bb[0],false,CPP_DEALLOC,1,(int)bb.size());
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
        int sw;
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
              MCAuto<DataArrayDouble> aaa=DataArrayDouble::New(); aaa->useArray(&bb[0],false,CPP_DEALLOC,1,(int)bb.size());
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
        int sw;
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
              MCAuto<DataArrayDouble> aaa=DataArrayDouble::New(); aaa->useArray(&bb[0],false,CPP_DEALLOC,1,(int)bb.size());
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
        int sw;
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
              MCAuto<DataArrayDouble> aaa=DataArrayDouble::New(); aaa->useArray(&bb[0],false,CPP_DEALLOC,1,(int)bb.size());
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
        convertFromPyObjVectorOfObj<const MEDCoupling::MEDCouplingFieldDouble *>(li,SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDouble,"MEDCouplingFieldDouble",tmp);
        return MEDCouplingFieldDouble::MergeFields(tmp);
      }

      static std::string WriteVTK(const char *fileName, PyObject *li, bool isBinary=true) throw(INTERP_KERNEL::Exception)
      {
        std::vector<const MEDCouplingFieldDouble *> tmp;
        convertFromPyObjVectorOfObj<const MEDCoupling::MEDCouplingFieldDouble *>(li,SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDouble,"MEDCouplingFieldDouble",tmp);
        return MEDCouplingFieldDouble::WriteVTK(fileName,tmp,isBinary);
      }

      PyObject *getTinySerializationInformation() const throw(INTERP_KERNEL::Exception)
      {
        return field_getTinySerializationInformation<MEDCouplingFieldDouble>(self);
      }
      
      PyObject *serialize() const throw(INTERP_KERNEL::Exception)
      {
        return field_serialize<double>(self);
      }

      static PyObject *___new___(PyObject *cls, PyObject *args) throw(INTERP_KERNEL::Exception)
      {
        return NewMethWrapCallInitOnlyIfDictWithSingleEltInInputGeneral<SinglePyObjExpectToBeAListOfSz2>(cls,args,"MEDCouplingFieldDouble");
      }

      PyObject *__getnewargs__() throw(INTERP_KERNEL::Exception)
      {// put an empty dict in input to say to __new__ to call __init__...
        return field__getnewargs__<MEDCouplingFieldDouble>(self);
      }

      PyObject *__getstate__() const throw(INTERP_KERNEL::Exception)
      {
        return field__getstate__<MEDCouplingFieldDouble>(self,MEDCoupling_MEDCouplingFieldDouble_getTinySerializationInformation,MEDCoupling_MEDCouplingFieldDouble_serialize);
      }
      
      void __setstate__(PyObject *inp) throw(INTERP_KERNEL::Exception)
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
    virtual std::string simpleRepr() const throw(INTERP_KERNEL::Exception);
    virtual std::string advancedRepr() const throw(INTERP_KERNEL::Exception);
    virtual bool isEqual(const MEDCouplingMultiFields *other, double meshPrec, double valsPrec) const;
    virtual bool isEqualWithoutConsideringStr(const MEDCouplingMultiFields *other, double meshPrec, double valsPrec) const;
    virtual void checkConsistencyLight() const throw(INTERP_KERNEL::Exception);
    %extend
       {
         std::string __str__() const throw(INTERP_KERNEL::Exception)
         {
           return self->simpleRepr();
         }
         static MEDCouplingMultiFields *New(PyObject *li) throw(INTERP_KERNEL::Exception)
         {
           std::vector<const MEDCoupling::MEDCouplingFieldDouble *> tmp;
           convertFromPyObjVectorOfObj<const MEDCoupling::MEDCouplingFieldDouble *>(li,SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDouble,"MEDCouplingFieldDouble",tmp);
           int sz=tmp.size();
           std::vector<MEDCouplingFieldDouble *> fs(sz);
           for(int i=0;i<sz;i++)
             fs[i]=const_cast<MEDCouplingFieldDouble *>(tmp[i]);
           return MEDCouplingMultiFields::New(fs);
         }
         MEDCouplingMultiFields(PyObject *li) throw(INTERP_KERNEL::Exception)
         {
           std::vector<const MEDCoupling::MEDCouplingFieldDouble *> tmp;
           convertFromPyObjVectorOfObj<const MEDCoupling::MEDCouplingFieldDouble *>(li,SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDouble,"MEDCouplingFieldDouble",tmp);
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
                   PyList_SetItem(res,i,SWIG_NewPointerObj(SWIG_as_voidptr(fields[i]),SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDouble, SWIG_POINTER_OWN | 0 ));
                 }
               else
                 {
                   PyList_SetItem(res,i,SWIG_NewPointerObj(SWIG_as_voidptr(0),SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDouble, 0 ));
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
               return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDouble, SWIG_POINTER_OWN | 0 );
             }
           else
             return SWIG_NewPointerObj(SWIG_as_voidptr(0),SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDouble, 0 );
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
                   PyList_SetItem(res,i,SWIG_NewPointerObj(SWIG_as_voidptr(0),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, 0 ));
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
                   PyList_SetItem(res,i,SWIG_NewPointerObj(SWIG_as_voidptr(0),SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh, 0 ));
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
                   PyList_SetItem(res,i,SWIG_NewPointerObj(SWIG_as_voidptr(ms[i]),SWIGTYPE_p_MEDCoupling__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
                 }
               else
                 {
                   PyList_SetItem(res,i,SWIG_NewPointerObj(SWIG_as_voidptr(0),SWIGTYPE_p_MEDCoupling__DataArrayDouble, 0 ));
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

  class MEDCouplingFieldInt : public MEDCouplingField
  {
  public:
    static MEDCouplingFieldInt *New(TypeOfField type, TypeOfTimeDiscretization td=ONE_TIME);
    static MEDCouplingFieldInt *New(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td=ONE_TIME);
    bool isEqual(const MEDCouplingFieldInt *other, double meshPrec, int valsPrec) const throw(INTERP_KERNEL::Exception);
    bool isEqualWithoutConsideringStr(const MEDCouplingFieldInt *other, double meshPrec, int valsPrec) const throw(INTERP_KERNEL::Exception);
    void setTimeUnit(const std::string& unit) throw(INTERP_KERNEL::Exception);
    std::string getTimeUnit() const throw(INTERP_KERNEL::Exception);
    void setTime(double val, int iteration, int order) throw(INTERP_KERNEL::Exception);
    void setArray(DataArrayInt *array) throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldInt *deepCopy() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldInt *clone(bool recDeepCpy) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldInt *cloneWithMesh(bool recDeepCpy) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *convertToDblField() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldInt *buildSubPartRange(int begin, int end, int step) const throw(INTERP_KERNEL::Exception);
    %extend {
      MEDCouplingFieldInt(TypeOfField type, TypeOfTimeDiscretization td=ONE_TIME)
      {
        return MEDCouplingFieldInt::New(type,td);
      }

      MEDCouplingFieldInt(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td=ONE_TIME)
      {
        return MEDCouplingFieldInt::New(ft,td);
      }

      PyObject *isEqualIfNotWhy(const MEDCouplingFieldInt *other, double meshPrec, int valsPrec) const throw(INTERP_KERNEL::Exception)
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

      MEDCouplingFieldInt *buildSubPart(PyObject *li) const throw(INTERP_KERNEL::Exception)
      {
        return fieldT_buildSubPart(self,li);
      }

      MEDCouplingFieldInt *__getitem__(PyObject *li) const throw(INTERP_KERNEL::Exception)
      {
        return fieldT__getitem__(self,li);
      }

      DataArrayInt *getArray() throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *ret=self->getArray();
        if(ret)
          ret->incrRef();
        return ret;
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

      PyObject *getTinySerializationInformation() const throw(INTERP_KERNEL::Exception)
      {
        return field_getTinySerializationInformation<MEDCouplingFieldInt>(self);
      }
      
      PyObject *serialize() const throw(INTERP_KERNEL::Exception)
      {
        return field_serialize<int>(self);
      }

      static PyObject *___new___(PyObject *cls, PyObject *args) throw(INTERP_KERNEL::Exception)
      {
        return NewMethWrapCallInitOnlyIfDictWithSingleEltInInputGeneral<SinglePyObjExpectToBeAListOfSz2>(cls,args,"MEDCouplingFieldInt");
      }

      PyObject *__getnewargs__() throw(INTERP_KERNEL::Exception)
      {// put an empty dict in input to say to __new__ to call __init__...
        return field__getnewargs__<MEDCouplingFieldInt>(self);
      }

      PyObject *__getstate__() const throw(INTERP_KERNEL::Exception)
      {
        return field__getstate__<MEDCouplingFieldInt>(self,MEDCoupling_MEDCouplingFieldInt_getTinySerializationInformation,MEDCoupling_MEDCouplingFieldInt_serialize);
      }
      
      void __setstate__(PyObject *inp) throw(INTERP_KERNEL::Exception)
      {
        field__setstate__<int>(self,inp);
      }
    }
  };

  class MEDCouplingFieldFloat : public MEDCouplingField
  {
  public:
    static MEDCouplingFieldFloat *New(TypeOfField type, TypeOfTimeDiscretization td=ONE_TIME);
    static MEDCouplingFieldFloat *New(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td=ONE_TIME);
    bool isEqual(const MEDCouplingFieldFloat *other, double meshPrec, float valsPrec) const throw(INTERP_KERNEL::Exception);
    bool isEqualWithoutConsideringStr(const MEDCouplingFieldFloat *other, double meshPrec, float valsPrec) const throw(INTERP_KERNEL::Exception);
    void setTimeUnit(const std::string& unit) throw(INTERP_KERNEL::Exception);
    std::string getTimeUnit() const throw(INTERP_KERNEL::Exception);
    void setTime(double val, int iteration, int order) throw(INTERP_KERNEL::Exception);
    void setArray(DataArrayFloat *array) throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldFloat *deepCopy() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldFloat *clone(bool recDeepCpy) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldFloat *cloneWithMesh(bool recDeepCpy) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *convertToDblField() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldFloat *buildSubPartRange(int begin, int end, int step) const throw(INTERP_KERNEL::Exception);
    %extend {
      MEDCouplingFieldFloat(TypeOfField type, TypeOfTimeDiscretization td=ONE_TIME)
      {
        return MEDCouplingFieldFloat::New(type,td);
      }

      MEDCouplingFieldFloat(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td=ONE_TIME)
      {
        return MEDCouplingFieldFloat::New(ft,td);
      }

      PyObject *isEqualIfNotWhy(const MEDCouplingFieldFloat *other, double meshPrec, float valsPrec) const throw(INTERP_KERNEL::Exception)
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

      MEDCouplingFieldFloat *buildSubPart(PyObject *li) const throw(INTERP_KERNEL::Exception)
      {
        return fieldT_buildSubPart(self,li);
      }

      MEDCouplingFieldFloat *__getitem__(PyObject *li) const throw(INTERP_KERNEL::Exception)
      {
        return fieldT__getitem__(self,li);
      }

      DataArrayFloat *getArray() throw(INTERP_KERNEL::Exception)
      {
        DataArrayFloat *ret=self->getArray();
        if(ret)
          ret->incrRef();
        return ret;
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

      PyObject *getTinySerializationInformation() const throw(INTERP_KERNEL::Exception)
      {
        return field_getTinySerializationInformation<MEDCouplingFieldFloat>(self);
      }
      
      PyObject *serialize() const throw(INTERP_KERNEL::Exception)
      {
        return field_serialize<float>(self);
      }
        
      static PyObject *___new___(PyObject *cls, PyObject *args) throw(INTERP_KERNEL::Exception)
      {
        return NewMethWrapCallInitOnlyIfDictWithSingleEltInInputGeneral<SinglePyObjExpectToBeAListOfSz2>(cls,args,"MEDCouplingFieldFloat");
      }
      
      PyObject *__getnewargs__() throw(INTERP_KERNEL::Exception)
      {// put an empty dict in input to say to __new__ to call __init__...
        return field__getnewargs__<MEDCouplingFieldFloat>(self);
      }
      
      PyObject *__getstate__() const throw(INTERP_KERNEL::Exception)
      {
        return field__getstate__<MEDCouplingFieldFloat>(self,MEDCoupling_MEDCouplingFieldFloat_getTinySerializationInformation,MEDCoupling_MEDCouplingFieldFloat_serialize);
      }
        
      void __setstate__(PyObject *inp) throw(INTERP_KERNEL::Exception)
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
            std::vector<const MEDCoupling::MEDCouplingFieldDouble *> tmp;
            convertFromPyObjVectorOfObj<const MEDCoupling::MEDCouplingFieldDouble *>(li,SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDouble,"MEDCouplingFieldDouble",tmp);
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
          std::vector<const MEDCoupling::MEDCouplingFieldDouble *> tmp;
          convertFromPyObjVectorOfObj<const MEDCoupling::MEDCouplingFieldDouble *>(li,SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDouble,"MEDCouplingFieldDouble",tmp);
           int sz=tmp.size();
           std::vector<MEDCouplingFieldDouble *> fs(sz);
           for(int i=0;i<sz;i++)
             fs[i]=const_cast<MEDCouplingFieldDouble *>(tmp[i]);
           return MEDCouplingFieldOverTime::New(fs);
         }
      }
  };

  class MEDCouplingCartesianAMRMesh;
  
  class MEDCouplingCartesianAMRPatchGen : public RefCountObject
  {
  public:
    int getNumberOfCellsRecursiveWithOverlap() const throw(INTERP_KERNEL::Exception);
    int getNumberOfCellsRecursiveWithoutOverlap() const throw(INTERP_KERNEL::Exception);
    int getMaxNumberOfLevelsRelativeToThis() const throw(INTERP_KERNEL::Exception);
    %extend
    {
      MEDCouplingCartesianAMRMeshGen *getMesh() const throw(INTERP_KERNEL::Exception)
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
    int getNumberOfOverlapedCellsForFather() const throw(INTERP_KERNEL::Exception);
    bool isInMyNeighborhood(const MEDCouplingCartesianAMRPatch *other, int ghostLev) const throw(INTERP_KERNEL::Exception);
    std::vector<int> computeCellGridSt() const throw(INTERP_KERNEL::Exception);
    %extend
    {
      PyObject *getBLTRRange() const throw(INTERP_KERNEL::Exception)
      {
        const std::vector< std::pair<int,int> >& ret(self->getBLTRRange());
        return convertFromVectorPairInt(ret);
      }

      PyObject *getBLTRRangeRelativeToGF() const throw(INTERP_KERNEL::Exception)
      {
        std::vector< std::pair<int,int> > ret(self->getBLTRRangeRelativeToGF());
        return convertFromVectorPairInt(ret);
      }

      void addPatch(PyObject *bottomLeftTopRight, const std::vector<int>& factors) throw(INTERP_KERNEL::Exception)
      {
        std::vector< std::pair<int,int> > inp;
        convertPyToVectorPairInt(bottomLeftTopRight,inp);
        self->addPatch(inp,factors);
      }

      MEDCouplingCartesianAMRPatch *__getitem__(int patchId) const throw(INTERP_KERNEL::Exception)
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

      void __delitem__(int patchId) throw(INTERP_KERNEL::Exception)
      {
        MEDCouplingCartesianAMRMeshGen *mesh(const_cast<MEDCouplingCartesianAMRMeshGen *>(self->getMesh()));
        if(!mesh)
          throw INTERP_KERNEL::Exception("wrap MEDCouplingCartesianAMRPatch.__delitem__ : no underlying mesh !");
        mesh->removePatch(patchId);
      }

      int __len__() const throw(INTERP_KERNEL::Exception)
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
    int getAbsoluteLevel() const throw(INTERP_KERNEL::Exception);
    int getAbsoluteLevelRelativeTo(const MEDCouplingCartesianAMRMeshGen *ref) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getPositionRelativeTo(const MEDCouplingCartesianAMRMeshGen *ref) const throw(INTERP_KERNEL::Exception);
    int getSpaceDimension() const throw(INTERP_KERNEL::Exception);
    const std::vector<int>& getFactors() const throw(INTERP_KERNEL::Exception);
    void setFactors(const std::vector<int>& newFactors) throw(INTERP_KERNEL::Exception);
    int getMaxNumberOfLevelsRelativeToThis() const throw(INTERP_KERNEL::Exception);
    int getNumberOfCellsAtCurrentLevel() const throw(INTERP_KERNEL::Exception);
    int getNumberOfCellsAtCurrentLevelGhost(int ghostLev) const throw(INTERP_KERNEL::Exception);
    int getNumberOfCellsRecursiveWithOverlap() const throw(INTERP_KERNEL::Exception);
    int getNumberOfCellsRecursiveWithoutOverlap() const throw(INTERP_KERNEL::Exception);
    bool isPatchInNeighborhoodOf(int patchId1, int patchId2, int ghostLev) const throw(INTERP_KERNEL::Exception);
   virtual void detachFromFather() throw(INTERP_KERNEL::Exception);
    //
    int getNumberOfPatches() const throw(INTERP_KERNEL::Exception);
    int getPatchIdFromChildMesh(const MEDCouplingCartesianAMRMeshGen *mesh) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *buildUnstructured() const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *extractGhostFrom(int ghostSz, const DataArrayDouble *arr) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getPatchIdsInTheNeighborhoodOf(int patchId, int ghostLev) const throw(INTERP_KERNEL::Exception);
    MEDCoupling1SGTUMesh *buildMeshFromPatchEnvelop() const throw(INTERP_KERNEL::Exception);
    MEDCoupling1SGTUMesh *buildMeshOfDirectChildrenOnly() const throw(INTERP_KERNEL::Exception);
    void removeAllPatches() throw(INTERP_KERNEL::Exception);
    void removePatch(int patchId) throw(INTERP_KERNEL::Exception);
    void createPatchesFromCriterion(const INTERP_KERNEL::BoxSplittingOptions& bso, const DataArrayByte *criterion, const std::vector<int>& factors) throw(INTERP_KERNEL::Exception);
    void createPatchesFromCriterion(const INTERP_KERNEL::BoxSplittingOptions& bso, const DataArrayDouble *criterion, const std::vector<int>& factors, double eps) throw(INTERP_KERNEL::Exception);
    DataArrayDouble *createCellFieldOnPatch(int patchId, const DataArrayDouble *cellFieldOnThis) const throw(INTERP_KERNEL::Exception);
    void fillCellFieldOnPatch(int patchId, const DataArrayDouble *cellFieldOnThis, DataArrayDouble *cellFieldOnPatch, bool isConservative=true) const throw(INTERP_KERNEL::Exception);
    void fillCellFieldOnPatchGhost(int patchId, const DataArrayDouble *cellFieldOnThis, DataArrayDouble *cellFieldOnPatch, int ghostLev, bool isConservative=true) const throw(INTERP_KERNEL::Exception);
    void fillCellFieldOnPatchOnlyOnGhostZone(int patchId, const DataArrayDouble *cellFieldOnThis, DataArrayDouble *cellFieldOnPatch, int ghostLev) const throw(INTERP_KERNEL::Exception);
    void fillCellFieldOnPatchOnlyOnGhostZoneWith(int ghostLev, const MEDCouplingCartesianAMRPatch *patchToBeModified, const MEDCouplingCartesianAMRPatch *neighborPatch, DataArrayDouble *cellFieldOnPatch, const DataArrayDouble *cellFieldNeighbor) const;
    void fillCellFieldComingFromPatch(int patchId, const DataArrayDouble *cellFieldOnPatch, DataArrayDouble *cellFieldOnThis, bool isConservative=true) const throw(INTERP_KERNEL::Exception);
    void fillCellFieldComingFromPatchGhost(int patchId, const DataArrayDouble *cellFieldOnPatch, DataArrayDouble *cellFieldOnThis, int ghostLev, bool isConservative=true) const throw(INTERP_KERNEL::Exception);
    DataArrayInt *findPatchesInTheNeighborhoodOf(int patchId, int ghostLev) const throw(INTERP_KERNEL::Exception);
    std::string buildPythonDumpOfThis() const throw(INTERP_KERNEL::Exception);
    %extend
    {
      void addPatch(PyObject *bottomLeftTopRight, const std::vector<int>& factors) throw(INTERP_KERNEL::Exception)
      {
        std::vector< std::pair<int,int> > inp;
        convertPyToVectorPairInt(bottomLeftTopRight,inp);
        self->addPatch(inp,factors);
      }

      PyObject *getPatches() const throw(INTERP_KERNEL::Exception)
      {
        std::vector< const MEDCouplingCartesianAMRPatch *> ps(self->getPatches());
        int sz(ps.size());
        PyObject *ret = PyList_New(sz);
        for(int i=0;i<sz;i++)
          {
            MEDCouplingCartesianAMRPatch *elt(const_cast<MEDCouplingCartesianAMRPatch *>(ps[i]));
            if(elt)
              elt->incrRef();
            PyList_SetItem(ret,i,convertCartesianAMRPatch(elt, SWIG_POINTER_OWN | 0 ));
          }
        return ret;
      }

      // agy : don't know why typemap fails here ??? let it in the extend section
      PyObject *deepCopy(MEDCouplingCartesianAMRMeshGen *father) const throw(INTERP_KERNEL::Exception)
      {
        return convertCartesianAMRMesh(self->deepCopy(father), SWIG_POINTER_OWN | 0 );
      }

      MEDCouplingCartesianAMRPatch *getPatchAtPosition(const std::vector<int>& pos) const throw(INTERP_KERNEL::Exception)
      {
        const MEDCouplingCartesianAMRPatch *ret(self->getPatchAtPosition(pos));
        MEDCouplingCartesianAMRPatch *ret2(const_cast<MEDCouplingCartesianAMRPatch *>(ret));
        if(ret2)
          ret2->incrRef();
        return ret2;
      }

      MEDCouplingCartesianAMRMeshGen *getMeshAtPosition(const std::vector<int>& pos) const throw(INTERP_KERNEL::Exception)
      {
        const MEDCouplingCartesianAMRMeshGen *ret(self->getMeshAtPosition(pos));
        MEDCouplingCartesianAMRMeshGen *ret2(const_cast<MEDCouplingCartesianAMRMeshGen *>(ret));
        if(ret2)
          ret2->incrRef();
        return ret2;
      }

      virtual PyObject *positionRelativeToGodFather() const throw(INTERP_KERNEL::Exception)
      {
        std::vector<int> out1;
        std::vector< std::pair<int,int> > out0(self->positionRelativeToGodFather(out1));
        PyObject *ret(PyTuple_New(2));
        PyTuple_SetItem(ret,0,convertFromVectorPairInt(out0));
        PyTuple_SetItem(ret,1,convertIntArrToPyList2(out1));
        return ret;
      }

      virtual PyObject *retrieveGridsAt(int absoluteLev) const throw(INTERP_KERNEL::Exception)
      {
        std::vector<MEDCouplingCartesianAMRPatchGen *> ps(self->retrieveGridsAt(absoluteLev));
        int sz(ps.size());
        PyObject *ret = PyList_New(sz);
        for(int i=0;i<sz;i++)
          PyList_SetItem(ret,i,convertCartesianAMRPatch(ps[i], SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      MEDCouplingFieldDouble *buildCellFieldOnRecurseWithoutOverlapWithoutGhost(int ghostSz, PyObject *recurseArrs) const
      {
        std::vector<const DataArrayDouble *> inp;
        convertFromPyObjVectorOfObj<const MEDCoupling::DataArrayDouble *>(recurseArrs,SWIGTYPE_p_MEDCoupling__DataArrayDouble,"DataArrayDouble",inp);
        return self->buildCellFieldOnRecurseWithoutOverlapWithoutGhost(ghostSz,inp);
      }

      virtual MEDCouplingCartesianAMRMeshGen *getFather() const throw(INTERP_KERNEL::Exception)
      {
        MEDCouplingCartesianAMRMeshGen *ret(const_cast<MEDCouplingCartesianAMRMeshGen *>(self->getFather()));
        if(ret)
          ret->incrRef();
        return ret;
      }
      
      virtual MEDCouplingCartesianAMRMeshGen *getGodFather() const throw(INTERP_KERNEL::Exception)
      {
        MEDCouplingCartesianAMRMeshGen *ret(const_cast<MEDCouplingCartesianAMRMeshGen *>(self->getGodFather()));
        if(ret)
          ret->incrRef();
        return ret;
      }

      MEDCouplingCartesianAMRPatch *getPatch(int patchId) const throw(INTERP_KERNEL::Exception)
      {
        MEDCouplingCartesianAMRPatch *ret(const_cast<MEDCouplingCartesianAMRPatch *>(self->getPatch(patchId)));
        if(ret)
          ret->incrRef();
        return ret;
      }

      MEDCouplingIMesh *getImageMesh() const throw(INTERP_KERNEL::Exception)
      {
        const MEDCouplingIMesh *ret(self->getImageMesh());
        if(ret)
          ret->incrRef();
        return const_cast<MEDCouplingIMesh *>(ret);
      }

      MEDCouplingCartesianAMRPatch *__getitem__(int patchId) const throw(INTERP_KERNEL::Exception)
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

      void fillCellFieldOnPatchGhostAdv(int patchId, const DataArrayDouble *cellFieldOnThis, int ghostLev, PyObject *arrsOnPatches, bool isConservative=true) const throw(INTERP_KERNEL::Exception)
      {
        std::vector<const MEDCoupling::DataArrayDouble *> arrsOnPatches2;
        convertFromPyObjVectorOfObj<const MEDCoupling::DataArrayDouble *>(arrsOnPatches,SWIGTYPE_p_MEDCoupling__DataArrayDouble,"DataArrayDouble",arrsOnPatches2);
        self->fillCellFieldOnPatchGhostAdv(patchId,cellFieldOnThis,ghostLev,arrsOnPatches2,isConservative);
      }

      void fillCellFieldOnPatchOnlyGhostAdv(int patchId, int ghostLev, PyObject *arrsOnPatches) const
      {
        std::vector<const MEDCoupling::DataArrayDouble *> arrsOnPatches2;
        convertFromPyObjVectorOfObj<const MEDCoupling::DataArrayDouble *>(arrsOnPatches,SWIGTYPE_p_MEDCoupling__DataArrayDouble,"DataArrayDouble",arrsOnPatches2);
        self->fillCellFieldOnPatchOnlyGhostAdv(patchId,ghostLev,arrsOnPatches2);
      }

      void __delitem__(int patchId) throw(INTERP_KERNEL::Exception)
      {
        self->removePatch(patchId);
      }

      int __len__() const throw(INTERP_KERNEL::Exception)
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
    static MEDCouplingCartesianAMRMesh *New(MEDCouplingIMesh *mesh) throw(INTERP_KERNEL::Exception);
    %extend
    {
      static MEDCouplingCartesianAMRMesh *New(const std::string& meshName, int spaceDim, PyObject *nodeStrct, PyObject *origin, PyObject *dxyz) throw(INTERP_KERNEL::Exception)
      {
        static const char msg0[]="MEDCouplingCartesianAMRMesh::New : error on 'origin' parameter !";
        static const char msg1[]="MEDCouplingCartesianAMRMesh::New : error on 'dxyz' parameter !";
        const int *nodeStrctPtr(0);
        const double *originPtr(0),*dxyzPtr(0);
        int sw,sz,val0;
        std::vector<int> bb0;
        nodeStrctPtr=convertIntStarLikePyObjToCppIntStar(nodeStrct,sw,sz,val0,bb0);
        //
        double val,val2;
        std::vector<double> bb,bb2;
        int sz1,sz2;
        originPtr=convertObjToPossibleCpp5_SingleCompo(origin,sw,val,bb,msg0,false,sz1);
        dxyzPtr=convertObjToPossibleCpp5_SingleCompo(dxyz,sw,val2,bb2,msg1,false,sz2);
        //
        return MEDCouplingCartesianAMRMesh::New(meshName,spaceDim,nodeStrctPtr,nodeStrctPtr+sz,originPtr,originPtr+sz1,dxyzPtr,dxyzPtr+sz2);
      }

      void createPatchesFromCriterionML(PyObject *bso, const DataArrayDouble *criterion, PyObject *factors, double eps) throw(INTERP_KERNEL::Exception)
      {
        std::vector<const INTERP_KERNEL::BoxSplittingOptions *> inp0;
        convertFromPyObjVectorOfObj<const INTERP_KERNEL::BoxSplittingOptions *>(bso,SWIGTYPE_p_INTERP_KERNEL__BoxSplittingOptions,"BoxSplittingOptions",inp0);
        std::vector< std::vector<int> > inp2;
        convertPyToVectorOfVectorOfInt(factors,inp2);
        self->createPatchesFromCriterionML(inp0,criterion,inp2,eps);
      }

      MEDCouplingCartesianAMRMesh(const std::string& meshName, int spaceDim, PyObject *nodeStrct, PyObject *origin, PyObject *dxyz) throw(INTERP_KERNEL::Exception)
      {
        return MEDCoupling_MEDCouplingCartesianAMRMesh_New__SWIG_1(meshName,spaceDim,nodeStrct,origin,dxyz);
      }

      MEDCouplingCartesianAMRMesh(MEDCouplingIMesh *mesh) throw(INTERP_KERNEL::Exception)
      {
        return MEDCouplingCartesianAMRMesh::New(mesh);
      }
    }
  };

  class MEDCouplingDataForGodFather : public RefCountObject
  {
  public:
    virtual void synchronizeFineToCoarse() throw(INTERP_KERNEL::Exception);
    virtual void synchronizeFineToCoarseBetween(int fromLev, int toLev) throw(INTERP_KERNEL::Exception);
    virtual void synchronizeCoarseToFine() throw(INTERP_KERNEL::Exception);
    virtual void synchronizeCoarseToFineBetween(int fromLev, int toLev) throw(INTERP_KERNEL::Exception);
    virtual void synchronizeAllGhostZones() throw(INTERP_KERNEL::Exception);
    virtual void synchronizeAllGhostZonesOfDirectChidrenOf(const MEDCouplingCartesianAMRMeshGen *mesh) throw(INTERP_KERNEL::Exception);
    virtual void synchronizeAllGhostZonesAtASpecifiedLevel(int level) throw(INTERP_KERNEL::Exception);
    virtual void synchronizeAllGhostZonesAtASpecifiedLevelUsingOnlyFather(int level) throw(INTERP_KERNEL::Exception);
    virtual void alloc() throw(INTERP_KERNEL::Exception);
    virtual void dealloc() throw(INTERP_KERNEL::Exception);
    %extend
    {
      MEDCouplingCartesianAMRMesh *getMyGodFather() throw(INTERP_KERNEL::Exception)
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
    int getNumberOfLevels() const throw(INTERP_KERNEL::Exception);
    MEDCouplingAMRAttribute *deepCopy() const throw(INTERP_KERNEL::Exception);
    MEDCouplingAMRAttribute *deepCpyWithoutGodFather() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *buildCellFieldOnRecurseWithoutOverlapWithoutGhost(MEDCouplingCartesianAMRMeshGen *mesh, const std::string& fieldName) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *buildCellFieldOnWithGhost(MEDCouplingCartesianAMRMeshGen *mesh, const std::string& fieldName) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *buildCellFieldOnWithoutGhost(MEDCouplingCartesianAMRMeshGen *mesh, const std::string& fieldName) const throw(INTERP_KERNEL::Exception);
    bool changeGodFather(MEDCouplingCartesianAMRMesh *gf) throw(INTERP_KERNEL::Exception);
    MEDCouplingAMRAttribute *projectTo(MEDCouplingCartesianAMRMesh *targetGF) const throw(INTERP_KERNEL::Exception);
    std::string writeVTHB(const std::string& fileName) const throw(INTERP_KERNEL::Exception);
    %extend
    {
      static MEDCouplingAMRAttribute *New(MEDCouplingCartesianAMRMesh *gf, PyObject *fieldNames, int ghostLev) throw(INTERP_KERNEL::Exception)
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

      MEDCouplingAMRAttribute(MEDCouplingCartesianAMRMesh *gf, PyObject *fieldNames, int ghostLev) throw(INTERP_KERNEL::Exception)
      {
        return MEDCoupling_MEDCouplingAMRAttribute_New(gf,fieldNames,ghostLev);
      }

      DataArrayDouble *getFieldOn(MEDCouplingCartesianAMRMeshGen *mesh, const std::string& fieldName) const throw(INTERP_KERNEL::Exception)
      {
        const DataArrayDouble *ret(self->getFieldOn(mesh,fieldName));
        DataArrayDouble *ret2(const_cast<DataArrayDouble *>(ret));
        if(ret2)
          ret2->incrRef();
        return ret2;
      }

      void spillInfoOnComponents(PyObject *compNames) throw(INTERP_KERNEL::Exception)
      {
        std::vector< std::vector<std::string> > compNamesCpp;
        convertPyToVectorOfVectorOfString(compNames,compNamesCpp);
        self->spillInfoOnComponents(compNamesCpp);
      }

      void spillNatures(PyObject *nfs) throw(INTERP_KERNEL::Exception)
      {
        std::vector<int> inp0;
        if(!fillIntVector(nfs,inp0))
          throw INTERP_KERNEL::Exception("wrap of MEDCouplingAMRAttribute::spillNatures : vector of NatureOfField enum expected !");
        std::size_t sz(inp0.size());
        std::vector<NatureOfField> inp00(sz);
        for(std::size_t i=0;i<sz;i++)
          inp00[i]=(NatureOfField)inp0[i];
        self->spillNatures(inp00);
      }
      
      PyObject *retrieveFieldsOn(MEDCouplingCartesianAMRMeshGen *mesh) const throw(INTERP_KERNEL::Exception)
      {
        std::vector<DataArrayDouble *> ret(self->retrieveFieldsOn(mesh));
        int sz((int)ret.size());
        PyObject *retPy(PyList_New(sz));
        for(int i=0;i<sz;i++)
          PyList_SetItem(retPy,i,SWIG_NewPointerObj(SWIG_as_voidptr(ret[i]),SWIGTYPE_p_MEDCoupling__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
        return retPy;
      }
    }
  };

  class DenseMatrix : public RefCountObject, public TimeLabel
  {
  public:
    static DenseMatrix *New(int nbRows, int nbCols) throw(INTERP_KERNEL::Exception);
    static DenseMatrix *New(DataArrayDouble *array, int nbRows, int nbCols) throw(INTERP_KERNEL::Exception);
    DenseMatrix *deepCopy() const throw(INTERP_KERNEL::Exception);
    DenseMatrix *shallowCpy() const throw(INTERP_KERNEL::Exception);
    //
    int getNumberOfRows() const throw(INTERP_KERNEL::Exception);
    int getNumberOfCols() const throw(INTERP_KERNEL::Exception);
    int getNbOfElems() const throw(INTERP_KERNEL::Exception);
    void reBuild(DataArrayDouble *array, int nbRows=-1, int nbCols=-1) throw(INTERP_KERNEL::Exception);
    void reShape(int nbRows, int nbCols) throw(INTERP_KERNEL::Exception);
    void transpose() throw(INTERP_KERNEL::Exception);
    //
    bool isEqual(const DenseMatrix& other, double eps) const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *matVecMult(const DataArrayDouble *vec) const throw(INTERP_KERNEL::Exception);
    static DataArrayDouble *MatVecMult(const DenseMatrix *mat, const DataArrayDouble *vec) throw(INTERP_KERNEL::Exception);
    %extend
    {
      DenseMatrix(int nbRows, int nbCols) throw(INTERP_KERNEL::Exception)
      {
        return DenseMatrix::New(nbRows,nbCols);
      }

      DenseMatrix(DataArrayDouble *array, int nbRows, int nbCols) throw(INTERP_KERNEL::Exception)
      {
        return DenseMatrix::New(array,nbRows,nbCols);
      }

      PyObject *isEqualIfNotWhy(const DenseMatrix& other, double eps) const throw(INTERP_KERNEL::Exception)
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

      DataArrayDouble *getData() throw(INTERP_KERNEL::Exception)
      {
        DataArrayDouble *ret(self->getData());
        if(ret)
          ret->incrRef();
        return ret;
      }

      DenseMatrix *__add__(const DenseMatrix *other) throw(INTERP_KERNEL::Exception)
      {
        return MEDCoupling::DenseMatrix::Add(self,other);
      }

      DenseMatrix *__sub__(const DenseMatrix *other) throw(INTERP_KERNEL::Exception)
      {
        return MEDCoupling::DenseMatrix::Substract(self,other);
      }

      DenseMatrix *__mul__(const DenseMatrix *other) throw(INTERP_KERNEL::Exception)
      {
        return MEDCoupling::DenseMatrix::Multiply(self,other);
      }

      DenseMatrix *__mul__(const DataArrayDouble *other) throw(INTERP_KERNEL::Exception)
      {
        return MEDCoupling::DenseMatrix::Multiply(self,other);
      }

      PyObject *___iadd___(PyObject *trueSelf, const DenseMatrix *other) throw(INTERP_KERNEL::Exception)
      {
        self->addEqual(other);
        Py_XINCREF(trueSelf);
        return trueSelf;
      }

      PyObject *___isub___(PyObject *trueSelf, const DenseMatrix *other) throw(INTERP_KERNEL::Exception)
      {
        self->substractEqual(other);
        Py_XINCREF(trueSelf);
        return trueSelf;
      }
#ifdef WITH_NUMPY
      PyObject *toNumPyMatrix() throw(INTERP_KERNEL::Exception) // not const. It is not a bug !
      {
        PyObject *obj(ToNumPyArrayUnderground<DataArrayDouble,double>(self->getData(),NPY_DOUBLE,"DataArrayDouble",self->getNumberOfRows(),self->getNumberOfCols()));
        return obj;
      }
#endif
    }
  };
}

%pythoncode %{
import os
__filename=os.environ.get('PYTHONSTARTUP')
if __filename and os.path.isfile(__filename):
  exec(open(__filename).read())
  pass
%}
