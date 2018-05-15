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
// * generation of documentation for inline methods,
// * groupping methods into "Basic API", "Advanced" and "Others..." sections


namespace MEDCoupling
{
  /*!
   * Returns a new MEDCouplingFieldDouble containing sum values of corresponding values of
   * \a this and a given field ( _f_ [ i, j ] = _this_ [ i, j ] + _other_ [ i, j ] ).
   * Number of tuples and components in the two fields must be the same.
   *  \param [in] other - the input field.
   *  \return MEDCouplingFieldDouble * - the new instance of MEDCouplingFieldDouble.
   *          The caller is to delete this result field using decrRef() as it is no more
   *          needed.
   *  \throw If the fields are not strictly compatible (areStrictlyCompatible()), i.e. they
   *         differ not only in values.
   */
  MEDCouplingFieldDouble *MEDCouplingFieldDouble::operator+(const MEDCouplingFieldDouble& other) const throw(INTERP_KERNEL::Exception) {}
  /*!
   * Returns a new MEDCouplingFieldDouble containing subtraction of corresponding values of
   * \a this and a given field ( _f_ [ i, j ] = _this_ [ i, j ] - _other_ [ i, j ] ).
   * Number of tuples and components in the two fields must be the same.
   *  \param [in] other - the field to subtract from \a this one.
   *  \return MEDCouplingFieldDouble * - the new instance of MEDCouplingFieldDouble.
   *          The caller is to delete this result field using decrRef() as it is no more
   *          needed.
   *  \throw If the fields are not strictly compatible (areStrictlyCompatible()), i.e. they
   *         differ not only in values.
   */
  MEDCouplingFieldDouble *MEDCouplingFieldDouble::operator-(const MEDCouplingFieldDouble& other) const throw(INTERP_KERNEL::Exception) {}
  /*!
   * Returns a new MEDCouplingFieldDouble containing product values of \a this and a
   * given field. There are 2 valid cases.
   * 1.  The fields have same number of tuples and components. Then each value of
   *   the result field (_f_) is a product of the corresponding values of _this_ and
   *   _other_, i.e. _f_ [ i, j ] = _this_ [ i, j ] * _other_ [ i, j ].
   * 2.  The fields have same number of tuples and one field, say _other_, has one
   *   component. Then
   *   _f_ [ i, j ] = _this_ [ i, j ] * _other_ [ i, 0 ].
   *
   * The two fields must have same number of tuples and same underlying mesh.
   *  \param [in] other - a factor field.
   *  \return MEDCouplingFieldDouble * - the new instance of MEDCouplingFieldDouble.
   *          The caller is to delete this result field using decrRef() as it is no more
   *          needed.
   *  \throw If the fields are not compatible for production (areCompatibleForMul()),
   *         i.e. they differ not only in values and possibly number of components.
   */
  MEDCouplingFieldDouble *MEDCouplingFieldDouble::operator*(const MEDCouplingFieldDouble& other) const throw(INTERP_KERNEL::Exception) {}
  /*!
   * Returns a new MEDCouplingFieldDouble containing division of \a this and a given
   * field. There are 2 valid cases.
   * 1.  The fields have same number of tuples and components. Then each value of
   *   the result field (_f_) is a division of the corresponding values of \a this and
   *   \a other, i.e. _f_ [ i, j ] = _this_ [ i, j ] / _other_ [ i, j ].
   * 2.  The fields have same number of tuples and _other_ has one component. Then
   *   _f_ [ i, j ] = _this_ [ i, j ] / _other_ [ i, 0 ].
   *
   *  \param [in] other - a denominator field.
   *  \return MEDCouplingFieldDouble * - the new instance of MEDCouplingFieldDouble.
   *          The caller is to delete this result field using decrRef() as it is no more
   *          needed.
   *  \throw If the fields are not compatible for division (areCompatibleForDiv()),
   *         i.e. they differ not only in values and possibly in number of components.
   */
  MEDCouplingFieldDouble *MEDCouplingFieldDouble::operator/(const MEDCouplingFieldDouble& other) const throw(INTERP_KERNEL::Exception) {}
  /*!
   * Returns a new MEDCouplingFieldDouble containing a dot product of \a this and a given field, 
   * so that the i-th tuple of the result field (_f_) is a sum of products of j-th components of
   * i-th tuples of two fields (\f$ f_i = \sum_ {}^n f1_j * f2_j \f$). 
   * Number of tuples and components in the two fields must be the same.
   *  \param [in] other - the input field.
   *  \return MEDCouplingFieldDouble * - the new instance of MEDCouplingFieldDouble.
   *          The caller is to delete this result field using decrRef() as it is no more
   *          needed.
   *  \throw If the fields are not strictly compatible (areStrictlyCompatible()), i.e. they
   *         differ not only in values.
   */
  MEDCouplingFieldDouble *MEDCouplingFieldDouble::dot(const MEDCouplingFieldDouble& other) const throw(INTERP_KERNEL::Exception) {}
  /*!
   * Returns a new MEDCouplingFieldDouble containing a cross product of \a this and 
   * a given field, so that the i-th tuple of the result field is a 3D vector which 
   * is a cross product of two vectors defined by the i-th tuples of the two fields.
   * Number of tuples in the fields must be the same.
   * Number of components in the fields must be 3.
   *  \param [in] other - the input field.
   *  \return MEDCouplingFieldDouble * - the new instance of MEDCouplingFieldDouble.
   *          The caller is to delete this result field using decrRef() as it is no more
   *          needed.
   *  \throw If \a this->getNumberOfComponents() != 3
   *  \throw If \a other->getNumberOfComponents() != 3
   *  \throw If the fields are not strictly compatible (areStrictlyCompatible()), i.e. they
   *         differ not only in values.
   */
  MEDCouplingFieldDouble *MEDCouplingFieldDouble::crossProduct(const MEDCouplingFieldDouble& other) const throw(INTERP_KERNEL::Exception) {}
  /*!
   * Returns a new MEDCouplingFieldDouble containing maximal values of \a this and a
   * given field. Number of tuples and components in the two fields must be the same.
   *  \param [in] other - the field to compare values with \a this one.
   *  \return MEDCouplingFieldDouble * - the new instance of MEDCouplingFieldDouble.
   *          The caller is to delete this result field using decrRef() as it is no more
   *          needed.
   *  \throw If the fields are not strictly compatible (areStrictlyCompatible()), i.e. they
   *         differ not only in values.
   */
  MEDCouplingFieldDouble *MEDCouplingFieldDouble::max(const MEDCouplingFieldDouble& other) const throw(INTERP_KERNEL::Exception) {}
  /*!
   * Returns a new MEDCouplingFieldDouble containing minimal values of \a this and a
   * given field. Number of tuples and components in the two fields must be the same.
   *  \param [in] other - the field to compare values with \a this one.
   *  \return MEDCouplingFieldDouble * - the new instance of MEDCouplingFieldDouble.
   *          The caller is to delete this result field using decrRef() as it is no more
   *          needed.
   *  \throw If the fields are not strictly compatible (areStrictlyCompatible()), i.e. they
   *         differ not only in values.
   */
  MEDCouplingFieldDouble *MEDCouplingFieldDouble::min(const MEDCouplingFieldDouble& other) const throw(INTERP_KERNEL::Exception) {}
  /*!
   * Returns the data array of \a this field.
   *  \return const DataArrayDouble * - a const pointer to the data array of \a this field.
   */
  const DataArrayDouble *MEDCouplingFieldDouble::getArray() const {}
  /*!
   * Returns the data array of \a this field apt for modification.
   *  \return DataArrayDouble * - a non-const pointer to the data array of \a this field.
   */
  DataArrayDouble *MEDCouplingFieldDouble::getArray() {}
  /*!
   * Sets a precision used to compare time values.
   *  \param [in] val - the precision value.
   */
  void MEDCouplingFieldDouble::setTimeTolerance(double val) {}
  /*!
   * Returns a precision used to compare time values.
   *  \return double - the precision value.
   */
  double MEDCouplingFieldDouble::getTimeTolerance() const {}
  /*!
   * Sets the number of iteration where the data array of \a this field has been calculated.
   * For examples of field construction, see \ref MEDCouplingFirstSteps3.
   *  \param [in] it - the iteration number.
   */
  void MEDCouplingFieldDouble::setIteration(int it) throw(INTERP_KERNEL::Exception) {}
  /*!
   * Sets the number of iteration where the second data array of \a this field has been calculated.
   * For examples of field construction, see \ref MEDCouplingFirstSteps3.
   *  \param [in] it - the iteration number.
   */
  void MEDCouplingFieldDouble::setEndIteration(int it) throw(INTERP_KERNEL::Exception) {}
  /*!
   * Sets the order number of iteration where the data array of \a this field has been calculated.
   * For examples of field construction, see \ref MEDCouplingFirstSteps3.
   *  \param [in] order - the order number.
   */
  void MEDCouplingFieldDouble::setOrder(int order) throw(INTERP_KERNEL::Exception) {}
  /*!
   * Sets the order number of iteration where the second data array of \a this field has
   * been calculated. 
   *  \param [in] order - the order number.
   */
  void MEDCouplingFieldDouble::setEndOrder(int order) throw(INTERP_KERNEL::Exception) {}
  /*!
   * Sets the time when the data array of \a this field has been calculated.
   * For examples of field construction, see \ref MEDCouplingFirstSteps3.
   *  \param [in] val - the time value.
   */
  void MEDCouplingFieldDouble::setTimeValue(double val) throw(INTERP_KERNEL::Exception) {}
  /*!
   * Sets the time when the second data array of \a this field has been calculated.
   *  \param [in] val - the time value.
   */
  void MEDCouplingFieldDouble::setEndTimeValue(double val) throw(INTERP_KERNEL::Exception) {}
  /*!
   * Sets time, number of iteration and order number of iteration when the data array
   * of \a this field has been calculated.
   * For examples of field construction, see \ref MEDCouplingFirstSteps3.
   *  \param [in] val - the time value.
   *  \param [in] iteration - the iteration number.
   *  \param [in] order - the order number.
   */
  void MEDCouplingFieldDouble::setTime(double val, int iteration, int order) {}
  /*!
   * Returns time, number of iteration and order number of iteration when the data array
   * of \a this field has been calculated.
   * For examples of field construction, see \ref MEDCouplingFirstSteps3.
   *  \param [out] iteration - the iteration number.
   *  \param [out] order - the order number.
   *  \return double - the time value.
   */
  double MEDCouplingFieldDouble::getTime(int& iteration, int& order) const {}
  /*!
   * Returns a value indexed by a tuple id and a component id.
   *  \param [in] tupleId - the id of the tuple of interest.
   *  \param [in] compoId - the id of the component of interest.
   *  \return double - the field value.
   */
  double MEDCouplingFieldDouble::getIJ(int tupleId, int compoId) const {}
}

namespace MEDCoupling
{
/*! \name Basic API   */
///@{
MEDCouplingFieldDouble::AddFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
MEDCouplingFieldDouble::CrossProductFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
MEDCouplingFieldDouble::DivideFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
MEDCouplingFieldDouble::DotFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
MEDCouplingFieldDouble::MaxFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
MEDCouplingFieldDouble::MeldFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
MEDCouplingFieldDouble::MergeFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
MEDCouplingFieldDouble::MergeFields(const std::vector<const MEDCouplingFieldDouble *>& a);
MEDCouplingFieldDouble::MinFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
MEDCouplingFieldDouble::MultiplyFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
MEDCouplingFieldDouble::New(TypeOfField type, TypeOfTimeDiscretization td=ONE_TIME);
MEDCouplingFieldDouble::New(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td=ONE_TIME);
MEDCouplingFieldDouble::SubstractFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
MEDCouplingFieldDouble::WriteVTK(const char *fileName, const std::vector<const MEDCouplingFieldDouble *>& fs);
MEDCouplingFieldDouble::accumulate(double *res) const;
MEDCouplingFieldDouble::accumulate(int compId) const;
MEDCouplingFieldDouble::advancedRepr() const;
MEDCouplingFieldDouble::applyFunc(const std::string &func);
MEDCouplingFieldDouble::applyFunc(int nbOfComp, FunctionToEvaluate func);
MEDCouplingFieldDouble::applyFunc(int nbOfComp, const std::string &func);
MEDCouplingFieldDouble::applyFunc(int nbOfComp, double val);
MEDCouplingFieldDouble::applyFunc2(int nbOfComp, const std::string &func);
MEDCouplingFieldDouble::applyFunc3(int nbOfComp, const std::vector<std::string>& varsOrder, const std::string &func);
MEDCouplingFieldDouble::applyLin(double a, double b, int compoId);
MEDCouplingFieldDouble::buildNewTimeReprFromThis(TypeOfTimeDiscretization td, bool deepCopy) const;
MEDCouplingFieldDouble::buildSubPart(const DataArrayInt *part) const;
MEDCouplingFieldDouble::buildSubPart(const int *partBg, const int *partEnd) const;
MEDCouplingFieldDouble::changeNbOfComponents(int newNbOfComp, double dftValue=0.);
MEDCouplingFieldDouble::changeUnderlyingMesh(const MEDCouplingMesh *other, int levOfCheck, double precOnMesh, double eps=1e-15);
MEDCouplingFieldDouble::checkCoherency() const;
MEDCouplingFieldDouble::clone(bool recDeepCpy) const;
MEDCouplingFieldDouble::cloneWithMesh(bool recDeepCpy) const;
MEDCouplingFieldDouble::copyTinyAttrFrom(const MEDCouplingFieldDouble *other);
MEDCouplingFieldDouble::copyTinyStringsFrom(const MEDCouplingField *other);
MEDCouplingFieldDouble::crossProduct(const MEDCouplingFieldDouble& other) const;
MEDCouplingFieldDouble::deepCpy() const;
MEDCouplingFieldDouble::determinant() const;
MEDCouplingFieldDouble::deviator() const;
MEDCouplingFieldDouble::dot(const MEDCouplingFieldDouble& other) const;
MEDCouplingFieldDouble::doublyContractedProduct() const;
MEDCouplingFieldDouble::eigenValues() const;
MEDCouplingFieldDouble::eigenVectors() const;
MEDCouplingFieldDouble::fillFromAnalytic(int nbOfComp, FunctionToEvaluate func);
MEDCouplingFieldDouble::fillFromAnalytic(int nbOfComp, const std::string &func);
MEDCouplingFieldDouble::fillFromAnalytic2(int nbOfComp, const std::string &func);
MEDCouplingFieldDouble::fillFromAnalytic3(int nbOfComp, const std::vector<std::string>& varsOrder, const std::string &func);
MEDCouplingFieldDouble::getArray() const;
MEDCouplingFieldDouble::getArray();
MEDCouplingFieldDouble::getAverageValue() const;
MEDCouplingFieldDouble::getIJ(int tupleId, int compoId) const;
MEDCouplingFieldDouble::getIJK(int cellId, int nodeIdInCell, int compoId) const;
MEDCouplingFieldDouble::getIdsInRange(double vmin, double vmax) const;
MEDCouplingFieldDouble::getMaxValue() const;
MEDCouplingFieldDouble::getMaxValue2(DataArrayInt*& tupleIds) const;
MEDCouplingFieldDouble::getMinValue() const;
MEDCouplingFieldDouble::getMinValue2(DataArrayInt*& tupleIds) const;
MEDCouplingFieldDouble::getNumberOfComponents() const;
MEDCouplingFieldDouble::getNumberOfTuples() const;
MEDCouplingFieldDouble::getNumberOfValues() const;
MEDCouplingFieldDouble::getTime(int& iteration, int& order) const;
MEDCouplingFieldDouble::getTimeDiscretization() const;
MEDCouplingFieldDouble::getTimeTolerance() const;
MEDCouplingFieldDouble::getTimeUnit() const;
MEDCouplingFieldDouble::getValueOn(const double *spaceLoc, double *res) const;
MEDCouplingFieldDouble::getValueOn(const double *spaceLoc, double time, double *res) const;
MEDCouplingFieldDouble::getValueOnMulti(const double *spaceLoc, int nbOfPoints) const;
MEDCouplingFieldDouble::getValueOnPos(int i, int j, int k, double *res) const;
MEDCouplingFieldDouble::getWeightedAverageValue(double *res, bool isWAbs=true) const;
MEDCouplingFieldDouble::getWeightedAverageValue(int compId, bool isWAbs=true) const;
MEDCouplingFieldDouble::integral(bool isWAbs, double *res) const;
MEDCouplingFieldDouble::integral(int compId, bool isWAbs) const;
MEDCouplingFieldDouble::inverse() const;
MEDCouplingFieldDouble::isEqualWithoutConsideringStr(const MEDCouplingField *other, double meshPrec, double valsPrec) const;
MEDCouplingFieldDouble::keepSelectedComponents(const std::vector<int>& compoIds) const;
MEDCouplingFieldDouble::magnitude() const;
MEDCouplingFieldDouble::max(const MEDCouplingFieldDouble& other) const;
MEDCouplingFieldDouble::maxPerTuple() const;
MEDCouplingFieldDouble::mergeNodes(double eps, double epsOnVals=1e-15);
MEDCouplingFieldDouble::mergeNodes2(double eps, double epsOnVals=1e-15);
MEDCouplingFieldDouble::min(const MEDCouplingFieldDouble& other) const;
MEDCouplingFieldDouble::norm2() const;
MEDCouplingFieldDouble::normL1(double *res) const;
MEDCouplingFieldDouble::normL1(int compId) const;
MEDCouplingFieldDouble::normL2(double *res) const;
MEDCouplingFieldDouble::normL2(int compId) const;
MEDCouplingFieldDouble::normMax() const;
MEDCouplingFieldDouble::renumberCells(const int *old2NewBg, bool check=true);
MEDCouplingFieldDouble::renumberNodes(const int *old2NewBg, double eps=1e-15);
MEDCouplingFieldDouble::setArray(DataArrayDouble *array);
MEDCouplingFieldDouble::setArrays(const std::vector<DataArrayDouble *>& arrs);
MEDCouplingFieldDouble::setEndArray(DataArrayDouble *array);
MEDCouplingFieldDouble::setEndIteration(int it);
MEDCouplingFieldDouble::setIteration(int it);
MEDCouplingFieldDouble::setNature(NatureOfField nat);
MEDCouplingFieldDouble::setOrder(int order);
MEDCouplingFieldDouble::setSelectedComponents(const MEDCouplingFieldDouble *f, const std::vector<int>& compoIds);
MEDCouplingFieldDouble::setTime(double val, int iteration, int order);
MEDCouplingFieldDouble::setTimeTolerance(double val);
MEDCouplingFieldDouble::setTimeUnit(const char *unit);
MEDCouplingFieldDouble::setTimeValue(double val);
MEDCouplingFieldDouble::simpleRepr() const;
MEDCouplingFieldDouble::simplexize(int policy);
MEDCouplingFieldDouble::sortPerTuple(bool asc);
MEDCouplingFieldDouble::substractInPlaceDM(const MEDCouplingFieldDouble *f, int levOfCheck, double precOnMesh, double eps=1e-15);
MEDCouplingFieldDouble::trace() const;
MEDCouplingFieldDouble::updateTime() const;
MEDCouplingFieldDouble::writeVTK(const char *fileName) const;
MEDCouplingFieldDouble::zipConnectivity(int compType, double epsOnVals=1e-15);
MEDCouplingFieldDouble::zipCoords(double epsOnVals=1e-15);
    MEDCouplingFieldDouble & MEDCouplingFieldDouble::operator=(double value);
    MEDCouplingFieldDouble * MEDCouplingFieldDouble::operator*(const MEDCouplingFieldDouble& other) const;
    MEDCouplingFieldDouble * MEDCouplingFieldDouble::operator+(const MEDCouplingFieldDouble& other) const;
    MEDCouplingFieldDouble * MEDCouplingFieldDouble::operator-(const MEDCouplingFieldDouble& other) const;
    MEDCouplingFieldDouble * MEDCouplingFieldDouble::operator/(const MEDCouplingFieldDouble& other) const;
    const MEDCouplingFieldDouble & MEDCouplingFieldDouble::operator*=(const MEDCouplingFieldDouble& other);
    const MEDCouplingFieldDouble & MEDCouplingFieldDouble::operator+=(const MEDCouplingFieldDouble& other);
    const MEDCouplingFieldDouble & MEDCouplingFieldDouble::operator-=(const MEDCouplingFieldDouble& other);
    const MEDCouplingFieldDouble & MEDCouplingFieldDouble::operator/=(const MEDCouplingFieldDouble& other);
///@} 
/*! \name   Advanced API   */
///@{
MEDCouplingFieldDouble::renumberCellsWithoutMesh(const int *old2NewBg, bool check=true);
MEDCouplingFieldDouble::renumberNodesWithoutMesh(const int *old2NewBg, int newNbOfNodes, double eps=1e-15);
///@} 

/*! \name Others... */
///@{
  MEDCouplingFieldDouble::negate() const;
  MEDCouplingFieldDouble::operator^(const MEDCouplingFieldDouble& other) const;
  MEDCouplingFieldDouble::operator^=(const MEDCouplingFieldDouble& other);
  MEDCouplingFieldDouble::PowFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2);
MEDCouplingFieldDouble::buildSubPartRange(int begin, int end, int step) const;
MEDCouplingFieldDouble::MEDCouplingFieldDouble(NatureOfField n, MEDCouplingTimeDiscretization *td, MEDCouplingFieldDiscretization *type);
MEDCouplingFieldDouble::MEDCouplingFieldDouble(TypeOfField type, TypeOfTimeDiscretization td);
MEDCouplingFieldDouble::MEDCouplingFieldDouble(const MEDCouplingFieldDouble& other, bool deepCopy);
MEDCouplingFieldDouble::MEDCouplingFieldDouble(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td);
MEDCouplingFieldDouble::applyFuncFast32(const std::string &func);
MEDCouplingFieldDouble::applyFuncFast64(const std::string &func);
MEDCouplingFieldDouble::areCompatibleForDiv(const MEDCouplingField *other) const;
MEDCouplingFieldDouble::areCompatibleForMeld(const MEDCouplingFieldDouble *other) const;
MEDCouplingFieldDouble::areCompatibleForMerge(const MEDCouplingField *other) const;
MEDCouplingFieldDouble::areCompatibleForMul(const MEDCouplingField *other) const;
MEDCouplingFieldDouble::areStrictlyCompatible(const MEDCouplingField *other) const;
MEDCouplingFieldDouble::copyAllTinyAttrFrom(const MEDCouplingFieldDouble *other);
MEDCouplingFieldDouble::extractSlice3D(const double *origin, const double *vec, double eps) const;
MEDCouplingFieldDouble::finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD, const std::vector<std::string>& tinyInfoS);
MEDCouplingFieldDouble::getArrays() const;
MEDCouplingFieldDouble::getEndArray() const;
MEDCouplingFieldDouble::getEndArray();
MEDCouplingFieldDouble::getEndTime(int& iteration, int& order) const;
MEDCouplingFieldDouble::getStartTime(int& iteration, int& order) const;
MEDCouplingFieldDouble::getTimeDiscretizationUnderGround() const;
MEDCouplingFieldDouble::getTimeDiscretizationUnderGround();
MEDCouplingFieldDouble::getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const;
MEDCouplingFieldDouble::getTinySerializationIntInformation(std::vector<int>& tinyInfo) const;
MEDCouplingFieldDouble::getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const;
MEDCouplingFieldDouble::isEqualIfNotWhy(const MEDCouplingField *other, double meshPrec, double valsPrec, std::string& reason) const;
MEDCouplingFieldDouble::reprQuickOverview(std::ostream& stream) const;
MEDCouplingFieldDouble::resizeForUnserialization(const std::vector<int>& tinyInfoI, DataArrayInt *&dataInt, std::vector<DataArrayDouble *>& arrays);
MEDCouplingFieldDouble::serialize(DataArrayInt *&dataInt, std::vector<DataArrayDouble *>& arrays) const;
MEDCouplingFieldDouble::setEndOrder(int order);
MEDCouplingFieldDouble::setEndTime(double val, int iteration, int order);
MEDCouplingFieldDouble::setEndTimeValue(double val);
MEDCouplingFieldDouble::setStartTime(double val, int iteration, int order);
MEDCouplingFieldDouble::synchronizeTimeWithMesh();
MEDCouplingFieldDouble::synchronizeTimeWithSupport();
MEDCouplingFieldDouble::~MEDCouplingFieldDouble();
MEDCouplingFieldDouble::_time_discr;
///@} 
}
