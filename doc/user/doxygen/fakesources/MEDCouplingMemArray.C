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
/*!
 * Returns the attribute \a _name of \a this array.
 * See \ref MEDCouplingArrayBasicsName "DataArrays infos" for more information.
 *  \return std::string - array name
 */
std::string DataArray::getName() const {}

/*!
 * Returns number of components of \a this array.
 * See \ref MEDCouplingArrayBasicsTuplesAndCompo "DataArrays infos" for more information.
 *  \return int - number of components
 */
int DataArray::getNumberOfComponents() const {}

/*!
 * Returns number of tuples of \a this array.
 * See \ref MEDCouplingArrayBasicsName "DataArrays infos" for more information.
 *  \return int - number of tuples
 */
int DataArray::getNumberOfTuples() const {}

/*!
 * Returns number of elements of \a this array.
 * See \ref MEDCouplingArrayBasicsName "DataArrays infos" for more information.
 *  \return int - number of elements == <em> this->getNumberOfTuples() * 
 *          this->getNumberOfComponents() </em>
 */
int DataArray::getNbOfElems() const {}

/*!
 * Throws an exception if number of elements in \a this array differs from a given one.
 * See \ref MEDCouplingArrayBasicsName "DataArrays infos" for more information.
 *  \param [in] nbOfElems - expected array size.
 *  \param [in] msg - message to return within the thrown exception.
 *  \throw if <em> this->getNbOfElems() != nbOfElems</em>.
 */
void DataArray::checkNbOfElems(int nbOfElems, const char *msg) const {}

/*!
 * Returns values of a specified tuple.
 *  \param [in] tupleId - index of the tuple of interest.
 *  \param [out] res - C array returning values of the \a tupleId-th tuple. The \a res
 *         must be allocated by the caller and be of size not less than \a
 *         this->getNumberOfComponents().
 */
void DataArrayDouble::getTuple(int tupleId, double *res) const {}

/*!
 * Returns a value of a specified element of \a this array.
 *  \param [in] tupleId - index of the tuple of interest.
 *  \param [in] compoId - index of the component of interest.
 *  \return double - the value of \a compoId-th component of \a tupleId-th tuple.
 */
double DataArrayDouble::getIJ(int tupleId, int compoId) const {}

/*!
 * Returns a pointer to the first element of the raw data of \a this array.
 *  \return double* - the pointer to the value of 0-th tuple and 0-th component.
 */
double * DataArrayDouble::getPointer() {}

/*!
 * Returns a const pointer to the first element of the raw data of \a this array.
 *  \return const double* - the pointer to the value of 0-th tuple and 0-th component.
 */
const double * DataArrayDouble::getConstPointer() const {}

/*!
 * Assigns a given value to a specified element of \a this array.
 *  \param [in] tupleId - index of the tuple to modify.
 *  \param [in] compoId - index of the component to modify.
 *  \param [in] newVal - the value to assign to the value at \a compoId-th component
 *              of \a tupleId-th tuple. 
 */
void DataArrayDouble::setIJ(int tupleId, int compoId, double newVal) {}

/*!
 * Assigns a given value to a specified element of \a this array which is not marked
 * as changed.
 *  \param [in] tupleId - index of the tuple to modify.
 *  \param [in] compoId - index of the component to modify.
 *  \param [in] newVal - the value to assign to the value at \a compoId-th component
 *              of \a tupleId-th tuple.
 */
void DataArrayDouble::setIJSilent(int tupleId, int compoId, double newVal) {}

/*!
 * Copies values from another array starting from a specified element of \a this.
 *  \param [in] id - index of element to assign the value \a element0 to.
 *  \param [in] element0 - value to assign to the \a id-th element of \a this.
 *  \param [in] others - values to assign to elements following the \a id-th
 *         element of \a this.
 *  \param [in] sizeOfOthers - number of values to copy from \a others.
 */
void DataArrayDouble::writeOnPlace(int id, double element0, const double *others, int sizeOfOthers) {}

/*!
 * Does nothing because this class does not aggregate any TimeLabel instance.
 */
void DataArrayDouble::updateTime() const {}

/*!
 * Returns values of a specified tuple.
 *  \param [in] tupleId - index of the tuple of interest.
 *  \param [out] res - C array returning values of the \a tupleId-th tuple. The \a res
 *         must be allocated by the caller and be of size not less than \a
 *         this->getNumberOfComponents().
 *  \if ENABLE_EXAMPLES
 *  \ref py_mcdataarrayint_getTuple "Here is a Python example".
 *  \endif
 */
void DataArrayInt::getTuple(int tupleId, int *res) const {}

/*!
 * Returns a value of a specified element of \a this array.
 *  \param [in] tupleId - index of the tuple of interest.
 *  \param [in] compoId - index of the component of interest.
 *  \return int - the value of \a compoId-th component of \a tupleId-th tuple.
 */
int DataArrayInt::getIJ(int tupleId, int compoId) const {}


/*!
 * Assigns a given value to a specified element of \a this array.
 *  \param [in] tupleId - index of the tuple to modify.
 *  \param [in] compoId - index of the component to modify.
 *  \param [in] newVal - the value to assign to the value at \a compoId-th component
 *              of \a tupleId-th tuple.
 * \warning As this method declares \a this array as modified, it is more optimal to use
 *          setIJSilent() for modification of muliple values of array and to call
 *          declareAsNew() after the modification is done.
 */
void DataArrayInt::setIJ(int tupleId, int compoId, int newVal) {}


/*!
 * Assigns a given value to a specified element of \a this array which is \b not marked
 * as changed.
 *  \param [in] tupleId - index of the tuple to modify.
 *  \param [in] compoId - index of the component to modify.
 *  \param [in] newVal - the value to assign to the value at \a compoId-th component
 *              of \a tupleId-th tuple.
 */
void DataArrayInt::setIJSilent(int tupleId, int compoId, int newVal) {}

/*!
 * Returns a pointer to the first element of the raw data of \a this array.
 *  \return int* - the pointer to the value of 0-th tuple and 0-th component.
 */
int * DataArrayInt::getPointer() {}

/*!
 * Returns a const pointer to the first element of the raw data of \a this array.
 *  \return const int* - the pointer to the value of 0-th tuple and 0-th component.
 */
const int * DataArrayInt::getConstPointer() const {}

/*!
 * Copies values from another array starting from a given element of \a this.
 *  \param [in] id - index of element to assign the value \a element0 to.
 *  \param [in] element0 - value to assign to the \a id-th element of \a this.
 *  \param [in] others - values to assign to elements following the \a id-th
 *         element of \a this.
 *  \param [in] sizeOfOthers - number of values to copy from \a others.
 */
void DataArrayInt::writeOnPlace(int id, int element0, const int *others, int sizeOfOthers) {}

}

namespace MEDCoupling
{
//================================================================================
/////////////////////// DataArray GROUPPING //////////////////////////////////////
//================================================================================

/*! \name Basic API   */
///@{
DataArray::setName(const char *name);
DataArray::copyStringInfoFrom(const DataArray& other);
DataArray::areInfoEquals(const DataArray& other) const;
DataArray::getName() const;
DataArray::setInfoOnComponents(const std::vector<std::string>& info);
DataArray::getVarOnComponent(int i) const;
DataArray::getUnitOnComponent(int i) const;
DataArray::setInfoOnComponent(int i, const char *info);
DataArray::getNumberOfComponents() const;
DataArray::getNumberOfTuples() const;
DataArray::getNbOfElems() const;
DataArray::checkNbOfElems(int nbOfElems, const char *msg) const;
DataArray::GetVarNameFromInfo(const std::string& info);
DataArray::GetUnitFromInfo(const std::string& info);
///@} 

/*! \name Others... */
///@{
DataArray::getHeapMemorySizeWithoutChildren() const;
DataArray::copyPartOfStringInfoFrom(const DataArray& other, const std::vector<int>& compoIds);
DataArray::copyPartOfStringInfoFrom2(const std::vector<int>& compoIds, const DataArray& other);
DataArray::areInfoEqualsIfNotWhy(const DataArray& other, std::string& reason) const;
DataArray::reprWithoutNameStream(std::ostream& stream) const;
DataArray::cppRepr(const char *varName) const;
DataArray::getInfoOnComponents() const;
DataArray::getInfoOnComponents();
DataArray::getVarsOnComponent() const;
DataArray::getUnitsOnComponent() const;
DataArray::getInfoOnComponent(int i) const;
DataArray::checkNbOfTuples(int nbOfTuples, const char *msg) const;
DataArray::checkNbOfComps(int nbOfCompo, const char *msg) const;
//DataArray::checkNbOfTuplesAndComp(const DataArray& other, const char *msg) const throw(INTERP_KERNEL::Exception);
//DataArray::checkNbOfTuplesAndComp(int nbOfTuples, int nbOfCompo, const char *msg) const throw(INTERP_KERNEL::Exception);
DataArray::GetNumberOfItemGivenBES(int begin, int end, int step, const char *msg);
DataArray::GetNumberOfItemGivenBESRelative(int begin, int end, int step, const char *msg);
DataArray::GetPosOfItemGivenBESRelativeNoThrow(int value, int begin, int end, int step);
DataArray::reprCppStream(const char *varName, std::ostream& stream) const;
DataArray::DataArray();
DataArray::CheckValueInRange(int ref, int value, const char *msg);
DataArray::CheckValueInRangeEx(int value, int start, int end, const char *msg);
DataArray::CheckClosingParInRange(int ref, int value, const char *msg);
std::string              DataArray::_name;
std::vector<std::string> DataArray::_info_on_compo;
///@} 

//================================================================================
/////////////////////// DataArrayDouble GROUPPING ////////////////////////////////
//================================================================================

/*! \name Basic API   */
///@{
DataArrayDouble::isAllocated() const;
//DataArrayDouble::setInfoAndChangeNbOfCompo(const std::vector<std::string>& info);
DataArrayDouble::doubleValue() const;
DataArrayDouble::empty() const;
DataArrayDouble::deepCpy() const;
DataArrayDouble::performCpy(bool deepCpy) const;
DataArrayDouble::cpyFrom(const DataArrayDouble& other);
DataArrayDouble::alloc(int nbOfTuple, int nbOfCompo);
DataArrayDouble::allocIfNecessary(int nbOfTuple, int nbOfCompo);
DataArrayDouble::fillWithZero();
DataArrayDouble::fillWithValue(double val);
DataArrayDouble::iota(double init=0.);
DataArrayDouble::isUniform(double val, double eps) const;
DataArrayDouble::sort(bool asc=true);
DataArrayDouble::reverse();
DataArrayDouble::checkMonotonic(bool increasing, double eps) const;
DataArrayDouble::isMonotonic(bool increasing, double eps) const;
DataArrayDouble::repr() const;
DataArrayDouble::isEqual(const DataArrayDouble& other, double prec) const;
DataArrayDouble::isEqualWithoutConsideringStr(const DataArrayDouble& other, double prec) const;
DataArrayDouble::reAlloc(int nbOfTuples);
DataArrayDouble::convertToIntArr() const;
DataArrayDouble::fromNoInterlace() const;
DataArrayDouble::toNoInterlace() const;
DataArrayDouble::renumberInPlace(const int* old2New);
DataArrayDouble::renumberInPlaceR(const int* new2Old);
DataArrayDouble::renumber(const int* old2New) const;
DataArrayDouble::renumberR(const int* new2Old) const;
DataArrayDouble::renumberAndReduce(const int* old2New, int newNbOfTuple) const;
DataArrayDouble::selectByTupleId(const int* new2OldBg, const int* new2OldEnd) const;
DataArrayDouble::selectByTupleIdSafe(const int* new2OldBg, const int* new2OldEnd) const;
DataArrayDouble::selectByTupleId2(int bg, int end2, int step) const;
DataArrayDouble::selectByTupleRanges(const std::vector<std::pair<int,int> >& ranges) const;
DataArrayDouble::subArray(int tupleIdBg, int tupleIdEnd=-1) const;
DataArrayDouble::rearrange(int newNbOfCompo);
DataArrayDouble::transpose();
DataArrayDouble::changeNbOfComponents(int newNbOfComp, double dftValue) const;
DataArrayDouble::keepSelectedComponents(const std::vector<int>& compoIds) const;
DataArrayDouble::meldWith(const DataArrayDouble* other);
DataArrayDouble::findCommonTuples(double prec, int limitTupleId, DataArrayInt *&comm, DataArrayInt *&commIndex) const;
DataArrayDouble::getDifferentValues(double prec, int limitTupleId=-1) const;
DataArrayDouble::setSelectedComponents(const DataArrayDouble* a, const std::vector<int>& compoIds);
DataArrayDouble::setPartOfValues1(const DataArrayDouble* a, int bgTuples, int endTuples, int stepTuples, int bgComp, int endComp, int stepComp, bool strictCompoCompare=true);
DataArrayDouble::setPartOfValuesSimple1(double a, int bgTuples, int endTuples, int stepTuples, int bgComp, int endComp, int stepComp);
DataArrayDouble::setPartOfValues2(const DataArrayDouble* a, const int* bgTuples, const int* endTuples, const int* bgComp, const int* endComp, bool strictCompoCompare=true);
DataArrayDouble::setPartOfValuesSimple2(double a, const int* bgTuples, const int* endTuples, const int* bgComp, const int* endComp);
DataArrayDouble::setPartOfValues3(const DataArrayDouble* a, const int* bgTuples, const int* endTuples, int bgComp, int endComp, int stepComp, bool strictCompoCompare=true);
DataArrayDouble::setPartOfValuesSimple3(double a, const int* bgTuples, const int* endTuples, int bgComp, int endComp, int stepComp);
DataArrayDouble::getTuple(int tupleId, double* res) const;
DataArrayDouble::getIJ(int tupleId, int compoId) const;
DataArrayDouble::back() const;
DataArrayDouble::getIJSafe(int tupleId, int compoId) const;
DataArrayDouble::setIJ(int tupleId, int compoId, double newVal);
DataArrayDouble::setIJSilent(int tupleId, int compoId, double newVal);
DataArrayDouble::writeOnPlace(int id, double element0, const double* others, int sizeOfOthers);
DataArrayDouble::checkNoNullValues() const;
DataArrayDouble::getMinMaxPerComponent(double* bounds) const;
DataArrayDouble::getMaxValue(int& tupleId) const;
DataArrayDouble::getMaxValueInArray() const;
DataArrayDouble::getMinValue(int& tupleId) const;
DataArrayDouble::getMinValueInArray() const;
DataArrayDouble::getMaxValue2(DataArrayInt*& tupleIds) const;
DataArrayDouble::getMinValue2(DataArrayInt*& tupleIds) const;
DataArrayDouble::getAverageValue() const;
DataArrayDouble::norm2() const;
DataArrayDouble::normMax() const;
DataArrayDouble::accumulate(double* res) const;
DataArrayDouble::accumulate(int compId) const;
DataArrayDouble::fromPolarToCart() const;
DataArrayDouble::fromCylToCart() const;
DataArrayDouble::fromSpherToCart() const;
DataArrayDouble::doublyContractedProduct() const;
DataArrayDouble::determinant() const;
DataArrayDouble::eigenValues() const;
DataArrayDouble::eigenVectors() const;
DataArrayDouble::inverse() const;
DataArrayDouble::trace() const;
DataArrayDouble::deviator() const;
DataArrayDouble::magnitude() const;
DataArrayDouble::maxPerTuple() const;
DataArrayDouble::sortPerTuple(bool asc);
DataArrayDouble::abs();
DataArrayDouble::applyLin(double a, double b, int compoId);
DataArrayDouble::applyLin(double a, double b);
DataArrayDouble::applyInv(double numerator);
DataArrayDouble::negate() const;
DataArrayDouble::applyFunc(int nbOfComp, const std::string& func, bool isSafe=true) const;
DataArrayDouble::applyFunc(const std::string& func, bool isSafe=true) const;
DataArrayDouble::applyFunc2(int nbOfComp, const std::string& func, bool isSafe=true) const;
DataArrayDouble::applyFunc3(int nbOfComp, const std::vector<std::string>& varsOrder, const std::string& func, bool isSafe=true) const;
DataArrayDouble::getIdsInRange(double vmin, double vmax) const;
DataArrayDouble::addEqual(const DataArrayDouble* other);
DataArrayDouble::substractEqual(const DataArrayDouble* other);
DataArrayDouble::multiplyEqual(const DataArrayDouble* other);
DataArrayDouble::divideEqual(const DataArrayDouble* other);
DataArrayDouble::updateTime() const;
DataArrayDouble::New();
DataArrayDouble::Aggregate(const DataArrayDouble* a1, const DataArrayDouble* a2);
DataArrayDouble::Aggregate(const std::vector<const DataArrayDouble *>& arr);
DataArrayDouble::Meld(const DataArrayDouble* a1, const DataArrayDouble* a2);
DataArrayDouble::Meld(const std::vector<const DataArrayDouble *>& arr);
DataArrayDouble::Dot(const DataArrayDouble* a1, const DataArrayDouble* a2);
DataArrayDouble::CrossProduct(const DataArrayDouble* a1, const DataArrayDouble* a2);
DataArrayDouble::Max(const DataArrayDouble* a1, const DataArrayDouble* a2);
DataArrayDouble::Min(const DataArrayDouble* a1, const DataArrayDouble* a2);
DataArrayDouble::Add(const DataArrayDouble* a1, const DataArrayDouble* a2);
DataArrayDouble::Substract(const DataArrayDouble* a1, const DataArrayDouble* a2);
DataArrayDouble::Multiply(const DataArrayDouble* a1, const DataArrayDouble* a2);
DataArrayDouble::Divide(const DataArrayDouble* a1, const DataArrayDouble* a2);
///@}

/*! \name   Advanced API   */
///@{
DataArrayDouble::checkAllocated() const;
DataArrayDouble::setPartOfValuesAdv(const DataArrayDouble* a, const DataArrayInt* tuplesSelec);
DataArrayDouble::setContigPartOfSelectedValues(int tupleIdStart, const DataArrayDouble* a, const DataArrayInt* tuplesSelec);
DataArrayDouble::setContigPartOfSelectedValues2(int tupleIdStart, const DataArrayDouble* a, int bg, int end2, int step);
DataArrayDouble::applyFunc(int nbOfComp, FunctionToEvaluate func) const;
///@} 

/*! \name Others... */
///@{
DataArrayDouble::getNumberOfTuples() const;
DataArrayDouble::getNbOfElems() const;
//DataArrayDouble::getHeapMemorySize() const;
DataArrayDouble::reserve(int nbOfElems);
DataArrayDouble::pushBackSilent(double val);
DataArrayDouble::popBackSilent();
DataArrayDouble::pack() const;
DataArrayDouble::getNbOfElemAllocated() const;
DataArrayDouble::reprZip() const;
DataArrayDouble::writeVTK(std::ostream& ofs, int indent, const char* nameInFile) const;
DataArrayDouble::reprStream(std::ostream& stream) const;
DataArrayDouble::reprZipStream(std::ostream& stream) const;
DataArrayDouble::reprWithoutNameStream(std::ostream& stream) const;
DataArrayDouble::reprZipWithoutNameStream(std::ostream& stream) const;
DataArrayDouble::reprCppStream(const char* varName, std::ostream& stream) const;
DataArrayDouble::isEqualIfNotWhy(const DataArrayDouble& other, double prec, std::string& reason) const;
DataArrayDouble::duplicateEachTupleNTimes(int nbTimes) const;
DataArrayDouble::setPartOfValues4(const DataArrayDouble* a, int bgTuples, int endTuples, int stepTuples, const int* bgComp, const int* endComp, bool strictCompoCompare=true);
DataArrayDouble::setPartOfValuesSimple4(double a, int bgTuples, int endTuples, int stepTuples, const int* bgComp, const int* endComp);
DataArrayDouble::getPointer();
DataArrayDouble::SetArrayIn(DataArrayDouble* newArray, DataArrayDouble* &arrayToSet);
DataArrayDouble::iterator();
DataArrayDouble::getConstPointer() const;
DataArrayDouble::begin() const;
DataArrayDouble::end() const;
DataArrayDouble::useArray(const double* array, bool ownership, DeallocType type, int nbOfTuple, int nbOfCompo);
DataArrayDouble::useExternalArrayWithRWAccess(const double* array, int nbOfTuple, int nbOfCompo);
DataArrayDouble::insertAtTheEnd(InputIterator first, InputIterator last);
DataArrayDouble::computeBBoxPerTuple(double epsilon=0.0) const;
DataArrayDouble::computeTupleIdsNearTuples(const DataArrayDouble* other, double eps, DataArrayInt *& c, DataArrayInt *& cI) const;
DataArrayDouble::recenterForMaxPrecision(double eps);
DataArrayDouble::distanceToTuple(const double* tupleBg, const double* tupleEnd, int& tupleId) const;
DataArrayDouble::buildEuclidianDistanceDenseMatrix() const;
DataArrayDouble::buildEuclidianDistanceDenseMatrixWith(const DataArrayDouble* other) const;
DataArrayDouble::applyFuncFast32(const std::string& func);
DataArrayDouble::applyFuncFast64(const std::string& func);
DataArrayDouble::getTinySerializationIntInformation(std::vector<int>& tinyInfo) const;
DataArrayDouble::getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const;
DataArrayDouble::resizeForUnserialization(const std::vector<int>& tinyInfoI);
DataArrayDouble::finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<std::string>& tinyInfoS);

DataArrayDouble::findCommonTuplesAlg(const double* bbox, int nbNodes, int limitNodeId, double prec, DataArrayInt* c, DataArrayInt* cI) const;

DataArrayDouble::FindTupleIdsNearTuplesAlg(const BBTree<SPACEDIM,int>& myTree, const double* pos, int nbOfTuples, double eps, DataArrayInt* c, DataArrayInt* cI);
DataArrayDouble::DataArrayDouble();
MemArray<double> DataArrayDouble::_mem;
///@} 

//================================================================================
/////////////////////// DataArrayInt GROUPPING ///////////////////////////////////
//================================================================================

/*! \name Advanced API  */
///@{
DataArrayInt::checkAllocated() const;
DataArrayInt::getHashCode() const;
DataArrayInt::buildPermutationArr(const DataArrayInt& other) const;
DataArrayInt::transformWithIndArr(const int* indArrBg, const int* indArrEnd);
DataArrayInt::transformWithIndArrR(const int* indArrBg, const int* indArrEnd) const;
DataArrayInt::splitByValueRange(const int* arrBg, const int* arrEnd,DataArrayInt *& castArr, DataArrayInt *& rankInsideCast, DataArrayInt *& castsPresent) const;
DataArrayInt::setPartOfValuesAdv(const DataArrayInt* a, const DataArrayInt* tuplesSelec);
DataArrayInt::setContigPartOfSelectedValues(int tupleIdStart, const DataArrayInt*a, const DataArrayInt* tuplesSelec);
DataArrayInt::setContigPartOfSelectedValues2(int tupleIdStart, const DataArrayInt* a, int bg, int end2, int step);
DataArrayInt::SetArrayIn(DataArrayInt* newArray, DataArrayInt* &arrayToSet);
///@} 

/*! \name Basic API   */
///@{
DataArrayInt::isAllocated() const;
//DataArrayInt::setInfoAndChangeNbOfCompo(const std::vector<std::string>& info);
DataArrayInt::intValue() const;
DataArrayInt::empty() const;
DataArrayInt::deepCpy() const;
DataArrayInt::performCpy(bool deepCpy) const;
DataArrayInt::cpyFrom(const DataArrayInt& other);
DataArrayInt::alloc(int nbOfTuple, int nbOfCompo);
DataArrayInt::allocIfNecessary(int nbOfTuple, int nbOfCompo);
DataArrayInt::isEqual(const DataArrayInt& other) const;
DataArrayInt::isEqualWithoutConsideringStr(const DataArrayInt& other) const;
DataArrayInt::isEqualWithoutConsideringStrAndOrder(const DataArrayInt& other) const;
DataArrayInt::sort(bool asc=true);
DataArrayInt::reverse();
DataArrayInt::fillWithZero();
DataArrayInt::fillWithValue(int val);
DataArrayInt::iota(int init=0);
DataArrayInt::repr() const;
DataArrayInt::invertArrayO2N2N2O(int newNbOfElem) const;
DataArrayInt::invertArrayN2O2O2N(int oldNbOfElem) const;
DataArrayInt::reAlloc(int nbOfTuples);
DataArrayInt::convertToDblArr() const;
DataArrayInt::fromNoInterlace() const;
DataArrayInt::toNoInterlace() const;
DataArrayInt::renumberInPlace(const int* old2New);
DataArrayInt::renumberInPlaceR(const int* new2Old);
DataArrayInt::renumber(const int* old2New) const;
DataArrayInt::renumberR(const int* new2Old) const;
DataArrayInt::renumberAndReduce(const int* old2NewBg, int newNbOfTuple) const;
DataArrayInt::selectByTupleId(const int* new2OldBg, const int* new2OldEnd) const;
DataArrayInt::selectByTupleIdSafe(const int* new2OldBg, const int* new2OldEnd) const;
DataArrayInt::selectByTupleId2(int bg, int end, int step) const;
DataArrayInt::selectByTupleRanges(const std::vector<std::pair<int,int> >& ranges) const;
DataArrayInt::checkAndPreparePermutation() const;
DataArrayInt::changeSurjectiveFormat(int targetNb, DataArrayInt *&arr, DataArrayInt *&arrI) const;
DataArrayInt::buildPermArrPerLevel() const;
DataArrayInt::isIdentity2(int sizeExpected) const;
DataArrayInt::isUniform(int val) const;
DataArrayInt::subArray(int tupleIdBg, int tupleIdEnd=-1) const;
DataArrayInt::rearrange(int newNbOfCompo);
DataArrayInt::transpose();
DataArrayInt::changeNbOfComponents(int newNbOfComp, int dftValue) const;
DataArrayInt::keepSelectedComponents(const std::vector<int>& compoIds) const;
DataArrayInt::meldWith(const DataArrayInt* other);
DataArrayInt::setSelectedComponents(const DataArrayInt* a, const std::vector<int>& compoIds);
DataArrayInt::setPartOfValues1(const DataArrayInt* a, int bgTuples, int endTuples, int stepTuples, int bgComp, int endComp, int stepComp, bool strictCompoCompare=true);
DataArrayInt::setPartOfValuesSimple1(int a, int bgTuples, int endTuples, int stepTuples, int bgComp, int endComp, int stepComp);
DataArrayInt::setPartOfValues2(const DataArrayInt* a, const int* bgTuples, const int* endTuples, const int* bgComp, const int* endComp, bool strictCompoCompare=true);
DataArrayInt::setPartOfValuesSimple2(int a, const int* bgTuples, const int* endTuples, const int* bgComp, const int* endComp);
DataArrayInt::setPartOfValues3(const DataArrayInt* a, const int* bgTuples, const int* endTuples, int bgComp, int endComp, int stepComp, bool strictCompoCompare=true);
DataArrayInt::setPartOfValuesSimple3(int a, const int* bgTuples, const int* endTuples, int bgComp, int endComp, int stepComp);
DataArrayInt::getTuple(int tupleId, int* res) const;
DataArrayInt::getIJ(int tupleId, int compoId) const;
DataArrayInt::getIJSafe(int tupleId, int compoId) const;
DataArrayInt::back() const;
DataArrayInt::setIJ(int tupleId, int compoId, int newVal);
DataArrayInt::setIJSilent(int tupleId, int compoId, int newVal);
DataArrayInt::getPointer();
DataArrayInt::getConstPointer() const;
DataArrayInt::getIdsEqual(int val) const;
DataArrayInt::getIdsNotEqual(int val) const;
DataArrayInt::getIdsEqualList(const int* valsBg, const int* valsEnd) const;
DataArrayInt::getIdsNotEqualList(const int* valsBg, const int* valsEnd) const;
DataArrayInt::changeValue(int oldValue, int newValue);
DataArrayInt::presenceOfValue(int value) const;
DataArrayInt::presenceOfValue(const std::vector<int>& vals) const;
DataArrayInt::getMaxValue(int& tupleId) const;
DataArrayInt::getMaxValueInArray() const;
DataArrayInt::getMinValue(int& tupleId) const;
DataArrayInt::getMinValueInArray() const;
DataArrayInt::abs();
DataArrayInt::applyLin(int a, int b, int compoId);
DataArrayInt::applyLin(int a, int b);
DataArrayInt::applyInv(int numerator);
DataArrayInt::negate() const;
DataArrayInt::applyDivideBy(int val);
DataArrayInt::applyModulus(int val);
DataArrayInt::applyRModulus(int val);
DataArrayInt::buildComplement(int nbOfElement) const;
DataArrayInt::buildSubstraction(const DataArrayInt* other) const;
DataArrayInt::buildUnion(const DataArrayInt* other) const;
DataArrayInt::buildIntersection(const DataArrayInt* other) const;
DataArrayInt::deltaShiftIndex() const;
DataArrayInt::computeOffsets();
DataArrayInt::computeOffsets2();
DataArrayInt::buildExplicitArrByRanges(const DataArrayInt* offsets) const;
DataArrayInt::useArray(const int* array, bool ownership, DeallocType type, int nbOfTuple, int nbOfCompo);
DataArrayInt::writeOnPlace(int id, int element0, const int* others, int sizeOfOthers);
DataArrayInt::addEqual(const DataArrayInt* other);
DataArrayInt::substractEqual(const DataArrayInt* other);
DataArrayInt::multiplyEqual(const DataArrayInt* other);
DataArrayInt::divideEqual(const DataArrayInt* other);
DataArrayInt::modulusEqual(const DataArrayInt* other);
DataArrayInt::updateTime() const;
DataArrayInt::New();
DataArrayInt::BuildOld2NewArrayFromSurjectiveFormat2(int nbOfOldTuples, const int* arr, const int* arrIBg, const int* arrIEnd, int &newNbOfTuples);
DataArrayInt::Aggregate(const DataArrayInt* a1, const DataArrayInt* a2, int offsetA2);
DataArrayInt::Aggregate(const std::vector<const DataArrayInt *>& arr);
DataArrayInt::Meld(const DataArrayInt* a1, const DataArrayInt* a2);
DataArrayInt::Meld(const std::vector<const DataArrayInt *>& arr);
DataArrayInt::MakePartition(const std::vector<const DataArrayInt *>& groups, int newNb, std::vector< std::vector<int> >& fidsOfGroups);
DataArrayInt::BuildUnion(const std::vector<const DataArrayInt *>& arr);
DataArrayInt::BuildIntersection(const std::vector<const DataArrayInt *>& arr);
DataArrayInt::Add(const DataArrayInt* a1, const DataArrayInt* a2);
DataArrayInt::Substract(const DataArrayInt* a1, const DataArrayInt* a2);
DataArrayInt::Multiply(const DataArrayInt* a1, const DataArrayInt* a2);
DataArrayInt::Divide(const DataArrayInt* a1, const DataArrayInt* a2);
DataArrayInt::Modulus(const DataArrayInt* a1, const DataArrayInt* a2);
DataArrayInt::CheckAndPreparePermutation(const int* start, const int* end);
DataArrayInt::Range(int begin, int end, int step);
///@} 

/*! \name Others... */
///@{
DataArrayInt::getNumberOfTuples() const;
DataArrayInt::getNbOfElems() const;
//DataArrayInt::getHeapMemorySize() const;
DataArrayInt::reserve(int nbOfElems);
DataArrayInt::pushBackSilent(int val);
DataArrayInt::popBackSilent();
DataArrayInt::pack() const;
DataArrayInt::getNbOfElemAllocated() const;
DataArrayInt::isEqualIfNotWhy(const DataArrayInt& other, std::string& reason) const;
DataArrayInt::checkMonotonic(bool increasing) const;
DataArrayInt::isMonotonic(bool increasing) const;
DataArrayInt::checkStrictlyMonotonic(bool increasing) const;
DataArrayInt::isStrictlyMonotonic(bool increasing) const;
DataArrayInt::reprZip() const;
DataArrayInt::writeVTK(std::ostream& ofs, int indent, const char* type, const char* nameInFile) const;
DataArrayInt::reprStream(std::ostream& stream) const;
DataArrayInt::reprZipStream(std::ostream& stream) const;
DataArrayInt::reprWithoutNameStream(std::ostream& stream) const;
DataArrayInt::reprZipWithoutNameStream(std::ostream& stream) const;
DataArrayInt::reprCppStream(const char* varName, std::ostream& stream) const;
DataArrayInt::invertArrayO2N2N2OBis(int newNbOfElem) const;
DataArrayInt::setPartOfValues4(const DataArrayInt* a, int bgTuples, int endTuples, int stepTuples, const int* bgComp, const int* endComp, bool strictCompoCompare=true);
DataArrayInt::setPartOfValuesSimple4(int a, int bgTuples, int endTuples, int stepTuples, const int* bgComp, const int* endComp);
DataArrayInt::iterator();
DataArrayInt::begin() const;
DataArrayInt::end() const;
DataArrayInt::locateTuple(const std::vector<int>& tupl) const;
DataArrayInt::locateValue(int value) const;
DataArrayInt::locateValue(const std::vector<int>& vals) const;
DataArrayInt::findIdSequence(const std::vector<int>& vals) const;
DataArrayInt::presenceOfTuple(const std::vector<int>& tupl) const;
DataArrayInt::accumulate(int* res) const;
DataArrayInt::accumulate(int compId) const;
DataArrayInt::getIdsInRange(int vmin, int vmax) const;
DataArrayInt::buildSubstractionOptimized(const DataArrayInt* other) const;
DataArrayInt::buildUnique() const;
DataArrayInt::findRangeIdForEachTuple(const DataArrayInt* ranges) const;
DataArrayInt::findIdInRangeForEachTuple(const DataArrayInt* ranges) const;
DataArrayInt::duplicateEachTupleNTimes(int nbTimes) const;
DataArrayInt::getDifferentValues() const;
> DataArrayInt::partitionByDifferentValues(std::vector<int>& differentIds) const;
DataArrayInt::useExternalArrayWithRWAccess(const int* array, int nbOfTuple, int nbOfCompo);
tor>
DataArrayInt::insertAtTheEnd(InputIterator first, InputIterator last);
DataArrayInt::getTinySerializationIntInformation(std::vector<int>& tinyInfo) const;
DataArrayInt::getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const;
DataArrayInt::resizeForUnserialization(const std::vector<int>& tinyInfoI);
DataArrayInt::finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<std::string>& tinyInfoS);
DataArrayInt::DataArrayInt();
MemArray<int> DataArrayInt::_mem;
///@}

}
