// Copyright (C) 2013-2016  CEA/DEN, EDF R&D, OPEN CASCADE
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

// This file contains some code used only for
// * generation of documentation for inline methods of DataArray* classes
// * groupping methods into "Basic API", "Advanced" and "Others..." sections


namespace MEDCoupling
{
  //================================================================================
  /*!
   * Checks if \a this and another MEDCouplingMesh are equal without considering
   * textual data like mesh name, names of spatial components etc.
   *  \param [in] other - an instance of MEDCouplingMesh to compare with \a this one.
   *  \param [in] prec - precision value used to compare node coordinates.
   *  \return bool - \c true if the two meshes are equal, \c false else.
   */
  //================================================================================

  bool MEDCouplingMesh::isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const {}

  /*!
   * Checks if \a this and \a other meshes are geometrically equivalent, else an
   * exception is thrown. The meshes are
   * considered equivalent if (1) \a this mesh contains the same nodes as the \a other
   * mesh (with a specified precision) and (2) \a this mesh contains the same cells as
   * the \a other mesh (with use of a specified cell comparison technique). The mapping 
   * from \a other to \a this for nodes and cells is returned via out parameters.
   *  \param [in] other - the mesh to compare with.
   *  \param [in] cellCompPol - id [0-2] of cell comparison method. See meaning of
   *         each method in description of MEDCouplingUMesh::zipConnectivityTraducer().
   *  \param [in] prec - the precision used to compare nodes of the two meshes.
   *  \param [out] cellCor - a cell permutation array in "Old to New" mode. The caller is
   *         to delete this array using decrRef() as it is no more needed.
   *  \param [out] nodeCor - a node permutation array in "Old to New" mode. The caller is
   *         to delete this array using decrRef() as it is no more needed.
   *  \throw If the two meshes do not match.
   */
  void MEDCouplingMesh::checkDeepEquivalWith(const MEDCouplingMesh *other, int cellCompPol, double prec,DataArrayInt *&cellCor, DataArrayInt *&nodeCor) const throw(INTERP_KERNEL::Exception) {}

  /*!
   * Checks if \a this and \a other meshes are geometrically equivalent, else an
   * exception is thrown. The meshes are considered equivalent if (1) they share the same
   * node coordinates array(s) and (2) they contain the same cells (with use of a specified
   * cell comparison technique). The mapping from cells of the \a other to ones of \a this 
   * is returned via an out parameter.
   *  \param [in] other - the mesh to compare with.
   *  \param [in] cellCompPol - id [0-2] of cell comparison method. See the meaning of
   *         each method in description of MEDCouplingUMesh::zipConnectivityTraducer().
   *  \param [in] prec - a not used parameter.
   *  \param [out] cellCor - the permutation array in "Old to New" mode. The caller is
   *         to delete this array using decrRef() as it is no more needed.
   *  \throw If the two meshes do not match.
   */
  void MEDCouplingMesh::checkDeepEquivalOnSameNodesWith(const MEDCouplingMesh *other, int cellCompPol, double prec,DataArrayInt *&cellCor) const throw(INTERP_KERNEL::Exception) {}
}

namespace MEDCoupling
{
//================================================================================
/////////////////////// GROUPPING members of MEDCouplingMesh /////////////////////
//================================================================================
/*! \name Basic API   */
///@{
  MEDCouplingMesh::MergeMeshes(const MEDCouplingMesh *mesh1, const MEDCouplingMesh *mesh2);
  MEDCouplingMesh::MergeMeshes(std::vector<const MEDCouplingMesh *>& meshes);
  MEDCouplingMesh::checkDeepEquivalOnSameNodesWith(const MEDCouplingMesh *other, int cellCompPol, double prec,DataArrayInt *&cellCor) const;
  MEDCouplingMesh::checkDeepEquivalWith(const MEDCouplingMesh *other, int cellCompPol, double prec,DataArrayInt *&cellCor, DataArrayInt *&nodeCor) const;
  MEDCouplingMesh::fillFromAnalytic(TypeOfField t, int nbOfComp, FunctionToEvaluate func) const;
  MEDCouplingMesh::fillFromAnalytic(TypeOfField t, int nbOfComp, const std::string& func) const;
  MEDCouplingMesh::fillFromAnalytic2(TypeOfField t, int nbOfComp, const std::string& func) const;
  MEDCouplingMesh::fillFromAnalytic3(TypeOfField t, int nbOfComp, const std::vector<std::string>& varsOrder, const std::string& func) const;
  MEDCouplingMesh::getCellsContainingPoint(const double *pos, double eps, std::vector<int>& elts) const;
  MEDCouplingMesh::getCellsContainingPoints(const double *pos, int nbOfPoints, double eps, std::vector<int>& elts, std::vector<int>& eltsIndex) const;
  MEDCouplingMesh::isEqual(const MEDCouplingMesh *other, double prec) const;
  MEDCouplingMesh::isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const = 0;
  MEDCouplingMesh::writeVTK(const char *fileName) const;
///@} 

/*! \name   Advanced API   */
///@{
  MEDCouplingMesh::getCellIdsFullyIncludedInNodeIds(const int *partBg, const int *partEnd) const;
///@} 

/*! \name Others... */
///@{
  MEDCouplingMesh::GetDimensionOfGeometricType(INTERP_KERNEL::NormalizedCellType type);
  MEDCouplingMesh::GetReprOfGeometricType(INTERP_KERNEL::NormalizedCellType type);
  MEDCouplingMesh::MEDCouplingMesh();
  MEDCouplingMesh::~MEDCouplingMesh();
  MEDCouplingMesh::MEDCouplingMesh(const MEDCouplingMesh& other);
  MEDCouplingMesh::advancedRepr() const = 0;
  MEDCouplingMesh::areCompatibleForMerge(const MEDCouplingMesh *other) const;
  MEDCouplingMesh::buildOrthogonalField() const = 0;
  MEDCouplingMesh::buildPart(const int *start, const int *end) const = 0;
  MEDCouplingMesh::buildPartAndReduceNodes(const int *start, const int *end, DataArrayInt*& arr) const = 0;
  MEDCouplingMesh::buildUnstructured() const;
  MEDCouplingMesh::checkCoherency() const;
  MEDCouplingMesh::checkCoherency1(double eps=1e-12) const;
  MEDCouplingMesh::checkFastEquivalWith(const MEDCouplingMesh *other, double prec) const;
  MEDCouplingMesh::checkGeoEquivalWith(const MEDCouplingMesh *other, int levOfCheck, double prec,DataArrayInt *&cellCor, DataArrayInt *&nodeCor) const;
  MEDCouplingMesh::checkTypeConsistencyAndContig(const std::vector<int>& code, const std::vector<const DataArrayInt *>& idsPerType) const;
  MEDCouplingMesh::computeIsoBarycenterOfNodesPerCell() const;
  MEDCouplingMesh::computeNbOfNodesPerCell() const;
  MEDCouplingMesh::copyTinyInfoFrom(const MEDCouplingMesh *other);
  MEDCouplingMesh::copyTinyStringsFrom(const MEDCouplingMesh *other);
  MEDCouplingMesh::deepCpy() const = 0;
  MEDCouplingMesh::getAllGeoTypes() const = 0;
  MEDCouplingMesh::getBarycenterAndOwner() const = 0;
  MEDCouplingMesh::getBoundingBox(double *bbox) const = 0;
  MEDCouplingMesh::getCellContainingPoint(const double *pos, double eps) const = 0;
  MEDCouplingMesh::getCoordinatesAndOwner() const = 0;
  MEDCouplingMesh::getCoordinatesOfNode(int nodeId, std::vector<double>& coo) const;
  MEDCouplingMesh::getDescription() const;
  MEDCouplingMesh::getDistributionOfTypes() const;
//  MEDCouplingMesh::getHeapMemorySize() const;
  MEDCouplingMesh::getMeasureField(bool isAbs) const = 0;
  MEDCouplingMesh::getMeasureFieldOnNode(bool isAbs) const = 0;
  MEDCouplingMesh::getMeshDimension() const = 0;
  MEDCouplingMesh::getName() const;
  MEDCouplingMesh::getNodeIdsOfCell(int cellId, std::vector<int>& conn) const = 0;
  MEDCouplingMesh::getNumberOfCells() const = 0;
  MEDCouplingMesh::getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType type) const = 0;
  MEDCouplingMesh::getNumberOfNodes() const = 0;
  MEDCouplingMesh::getSpaceDimension() const = 0;
  MEDCouplingMesh::getTime(int& iteration, int& order) const;
  MEDCouplingMesh::getTimeUnit() const;
  MEDCouplingMesh::getTinySerializationInformation(std::vector<double>& tinyInfoD, std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const = 0;
  MEDCouplingMesh::getType() const = 0;
  MEDCouplingMesh::getTypeOfCell(int cellId) const = 0;
  MEDCouplingMesh::getVTKDataSetType() const;
  MEDCouplingMesh::giveCellsWithType(INTERP_KERNEL::NormalizedCellType type) const;
  MEDCouplingMesh::isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const;
  MEDCouplingMesh::isStructured() const;
  MEDCouplingMesh::mergeMyselfWith(const MEDCouplingMesh *other) const = 0;
  MEDCouplingMesh::renumberCells(const int *old2NewBg, bool check=true);
  MEDCouplingMesh::reprQuickOverview(std::ostream& stream) const;
  MEDCouplingMesh::resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings) const = 0;
  MEDCouplingMesh::rotate(const double *center, const double *vector, double angle) = 0;
  MEDCouplingMesh::scale(const double *point, double factor) = 0;
  MEDCouplingMesh::serialize(DataArrayInt *&a1, DataArrayDouble *&a2) const = 0;
  MEDCouplingMesh::setDescription(const char *descr);
  MEDCouplingMesh::setName(const char *name);
  MEDCouplingMesh::setTime(double val, int iteration, int order);
  MEDCouplingMesh::setTimeUnit(const char *unit);
  MEDCouplingMesh::simpleRepr() const = 0;
  MEDCouplingMesh::simplexize(int policy);
  MEDCouplingMesh::splitProfilePerType(const DataArrayInt *profile, std::vector<int>& code, std::vector<DataArrayInt *>& idsInPflPerType, std::vector<DataArrayInt *>& idsPerType) const;
  MEDCouplingMesh::translate(const double *vector) = 0;
  MEDCouplingMesh::unserialization(const std::vector<double>& tinyInfoD, const std::vector<int>& tinyInfo, const DataArrayInt *a1, DataArrayDouble *a2,const std::vector<std::string>& littleStrings) = 0;
  //MEDCouplingMesh::writeVTKAdvanced(const char *fileName, const std::string& cda, const std::string& pda) const;
  MEDCouplingMesh::writeVTKLL(std::ostream& ofs, const std::string& cellData, const std::string& pointData) const;
  MEDCouplingMesh::_description;
  MEDCouplingMesh::_iteration;
  MEDCouplingMesh::_name;
  MEDCouplingMesh::_order;
  MEDCouplingMesh::_time;
  MEDCouplingMesh::_time_unit;
///@}
}
