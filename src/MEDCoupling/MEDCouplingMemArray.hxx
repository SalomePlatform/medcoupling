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

#ifndef __MEDCOUPLING_MEDCOUPLINGMEMARRAY_HXX__
#define __MEDCOUPLING_MEDCOUPLINGMEMARRAY_HXX__

#include "MEDCoupling.hxx"
#include "MCAuto.hxx"
#include "MEDCouplingTimeLabel.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "InterpKernelException.hxx"
#include "MEDCouplingTraits.hxx"
#include "BBTreePts.txx"

#include <string>
#include <vector>
#include <iterator>

namespace MEDCoupling
{
  typedef enum
    {
      AX_CART = 3,
      AX_CYL = 4,
      AX_SPHER = 5
    } MEDCouplingAxisType;
  // -- WARNING this enum must be synchronized with MEDCouplingCommon.i file ! --

  class PartDefinition;
  
  template<class T>
  class MEDCouplingPointer
  {
  public:
    MEDCouplingPointer():_internal(0),_external(0) { }
    void null() { _internal=0; _external=0; }
    bool isNull() const { return _internal==0 && _external==0; }
    void setInternal(T *pointer);
    void setExternal(const T *pointer);
    const T *getConstPointer() const { if(_internal) return _internal; else return _external; }
    const T *getConstPointerLoc(std::size_t offset) const { if(_internal) return _internal+offset; else return _external+offset; }
    T *getPointer() { if(_internal) return _internal; if(_external) throw INTERP_KERNEL::Exception("Trying to write on an external pointer."); else return 0; }
  private:
    T *_internal;
    const T *_external;
  };

  template<class T>
  class MemArray
  {
  public:
    typedef void (*Deallocator)(void *,void *);
  public:
    MemArray():_nb_of_elem(0),_nb_of_elem_alloc(0),_ownership(false),_dealloc(0),_param_for_deallocator(0) { }
    MemArray(const MemArray<T>& other);
    bool isNull() const { return _pointer.isNull(); }
    const T *getConstPointerLoc(std::size_t offset) const { return _pointer.getConstPointerLoc(offset); }
    const T *getConstPointer() const { return _pointer.getConstPointer(); }
    std::size_t getNbOfElem() const { return _nb_of_elem; }
    std::size_t getNbOfElemAllocated() const { return _nb_of_elem_alloc; }
    T *getPointer() { return _pointer.getPointer(); }
     MemArray<T> &operator=(const MemArray<T>& other);
    T operator[](std::size_t id) const { return _pointer.getConstPointer()[id]; }
    T& operator[](std::size_t id) { return _pointer.getPointer()[id]; }
    bool isEqual(const MemArray<T>& other, T prec, std::string& reason) const;
    void repr(int sl, std::ostream& stream) const;
    bool reprHeader(int sl, std::ostream& stream) const;
    void reprZip(int sl, std::ostream& stream) const;
    void reprNotTooLong(int sl, std::ostream& stream) const;
    void fillWithValue(const T& val);
    T *fromNoInterlace(int nbOfComp) const;
    T *toNoInterlace(int nbOfComp) const;
    void sort(bool asc);
    void reverse(int nbOfComp);
    void alloc(std::size_t nbOfElements);
    void reserve(std::size_t newNbOfElements);
    void reAlloc(std::size_t newNbOfElements);
    void useArray(const T *array, bool ownership, DeallocType type, std::size_t nbOfElem);
    void useExternalArrayWithRWAccess(const T *array, std::size_t nbOfElem);
    void writeOnPlace(std::size_t id, T element0, const T *others, std::size_t sizeOfOthers);
    template<class InputIterator>
    void insertAtTheEnd(InputIterator first, InputIterator last);
    void pushBack(T elem);
    T popBack();
    void pack() const;
    bool isDeallocatorCalled() const { return _ownership; }
    Deallocator getDeallocator() const { return _dealloc; }
    void setSpecificDeallocator(Deallocator dealloc) { _dealloc=dealloc; }
    void setParameterForDeallocator(void *param) { _param_for_deallocator=param; }
    void *getParameterForDeallocator() const { return _param_for_deallocator; }
    void destroy();
    ~MemArray() { destroy(); }
  public:
    static void CPPDeallocator(void *pt, void *param);
    static void CDeallocator(void *pt, void *param);
  private:
    static void DestroyPointer(T *pt, Deallocator dealloc, void *param);
    static Deallocator BuildFromType(DeallocType type);
  private:
    std::size_t _nb_of_elem;
    std::size_t _nb_of_elem_alloc;
    bool _ownership;
    MEDCouplingPointer<T> _pointer;
    Deallocator _dealloc;
    void *_param_for_deallocator;
  };

  class DataArrayInt;
  class DataArrayByte;

  class MEDCOUPLING_EXPORT DataArray : public RefCountObject, public TimeLabel
  {
  public:
     std::size_t getHeapMemorySizeWithoutChildren() const;
     std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
     void setName(const std::string& name);
     void copyStringInfoFrom(const DataArray& other);
     void copyPartOfStringInfoFrom(const DataArray& other, const std::vector<int>& compoIds);
     void copyPartOfStringInfoFrom2(const std::vector<int>& compoIds, const DataArray& other);
     bool areInfoEqualsIfNotWhy(const DataArray& other, std::string& reason) const;
     bool areInfoEquals(const DataArray& other) const;
     std::string cppRepr(const std::string& varName) const;
     std::string getName() const { return _name; }
     const std::vector<std::string> &getInfoOnComponents() const { return _info_on_compo; }
     std::vector<std::string> &getInfoOnComponents() { return _info_on_compo; }
     void setInfoOnComponents(const std::vector<std::string>& info);
     void setInfoAndChangeNbOfCompo(const std::vector<std::string>& info);
     std::vector<std::string> getVarsOnComponent() const;
     std::vector<std::string> getUnitsOnComponent() const;
     std::string getInfoOnComponent(int i) const;
     std::string getVarOnComponent(int i) const;
     std::string getUnitOnComponent(int i) const;
     void setInfoOnComponent(int i, const std::string& info);
     std::size_t getNumberOfComponents() const { return _info_on_compo.size(); }
     void setPartOfValuesBase3(const DataArray *aBase, const int *bgTuples, const int *endTuples, int bgComp, int endComp, int stepComp, bool strictCompoCompare=true);
     virtual void *getVoidStarPointer() = 0;
     virtual DataArray *deepCopy() const = 0;
     virtual DataArray *buildNewEmptyInstance() const = 0;
     virtual bool isAllocated() const = 0;
     virtual void checkAllocated() const = 0;
     virtual void desallocate() = 0;
     virtual int getNumberOfTuples() const = 0;
     virtual std::size_t getNbOfElems() const = 0;
     virtual std::size_t getNbOfElemAllocated() const = 0;
     virtual void alloc(std::size_t nbOfTuple, std::size_t nbOfCompo=1) = 0;
     virtual void reAlloc(std::size_t newNbOfTuple) = 0;
     virtual void renumberInPlace(const int *old2New) = 0;
     virtual void renumberInPlaceR(const int *new2Old) = 0;
     virtual void setContigPartOfSelectedValues(int tupleIdStart, const DataArray *aBase, const DataArrayInt *tuplesSelec) = 0;
     virtual void setContigPartOfSelectedValuesSlice(int tupleIdStart, const DataArray *aBase, int bg, int end2, int step) = 0;
     virtual DataArray *selectByTupleRanges(const std::vector<std::pair<int,int> >& ranges) const = 0;
     virtual DataArray *keepSelectedComponents(const std::vector<int>& compoIds) const = 0;
     virtual DataArray *selectByTupleId(const int *new2OldBg, const int *new2OldEnd) const = 0;
     virtual DataArray *selectByTupleIdSafe(const int *new2OldBg, const int *new2OldEnd) const = 0;
     virtual DataArray *selectByTupleIdSafeSlice(int bg, int end2, int step) const = 0;
     virtual void rearrange(int newNbOfCompo) = 0;
     virtual void circularPermutation(int nbOfShift=1) = 0;
     virtual void circularPermutationPerTuple(int nbOfShift=1) = 0;
     virtual void reversePerTuple() = 0;
     void checkNbOfTuples(int nbOfTuples, const std::string& msg) const;
     void checkNbOfComps(int nbOfCompo, const std::string& msg) const;
     void checkNbOfTuplesAndComp(const DataArray& other, const std::string& msg) const;
     void checkNbOfTuplesAndComp(int nbOfTuples, int nbOfCompo, const std::string& msg) const;
     void checkNbOfElems(std::size_t nbOfElems, const std::string& msg) const;
     static void GetSlice(int start, int stop, int step, int sliceId, int nbOfSlices, int& startSlice, int& stopSlice);
     static int GetNumberOfItemGivenBES(int begin, int end, int step, const std::string& msg);
     static int GetNumberOfItemGivenBESRelative(int begin, int end, int step, const std::string& msg);
     static int GetPosOfItemGivenBESRelativeNoThrow(int value, int begin, int end, int step);
     static std::string GetVarNameFromInfo(const std::string& info);
     static std::string GetUnitFromInfo(const std::string& info);
     static std::string BuildInfoFromVarAndUnit(const std::string& var, const std::string& unit);
     static std::string GetAxisTypeRepr(MEDCouplingAxisType at);
     static DataArray *Aggregate(const std::vector<const DataArray *>& arrs);
     virtual void reprStream(std::ostream& stream) const = 0;
     virtual void reprZipStream(std::ostream& stream) const = 0;
     virtual void reprWithoutNameStream(std::ostream& stream) const;
     virtual void reprZipWithoutNameStream(std::ostream& stream) const = 0;
     virtual void reprCppStream(const std::string& varName, std::ostream& stream) const = 0;
     virtual void reprQuickOverview(std::ostream& stream) const = 0;
     virtual void reprQuickOverviewData(std::ostream& stream, std::size_t maxNbOfByteInRepr) const = 0;
  protected:
    DataArray() { }
    ~DataArray() { }
  protected:
     static void CheckValueInRange(int ref, int value, const std::string& msg);
     static void CheckValueInRangeEx(int value, int start, int end, const std::string& msg);
     static void CheckClosingParInRange(int ref, int value, const std::string& msg);
     static int EffectiveCircPerm(int nbOfShift, int nbOfTuples);
  protected:
    std::string _name;
    std::vector<std::string> _info_on_compo;
  };
}

namespace MEDCoupling
{
  class DataArrayInt;

  template<class T>
  class DataArrayTemplate : public DataArray
  {
  public:
     static MCAuto< typename Traits<T>::ArrayTypeCh > NewFromStdVector(const typename std::vector<T>& v);
     std::vector< MCAuto< typename Traits<T>::ArrayTypeCh > > explodeComponents() const;
    //
     std::size_t getHeapMemorySizeWithoutChildren() const;
    //
     int getNumberOfTuples() const { return _info_on_compo.empty()?0:_mem.getNbOfElem()/getNumberOfComponents(); }
     std::size_t getNbOfElems() const { return _mem.getNbOfElem(); }
     bool empty() const;
     void *getVoidStarPointer() { return getPointer(); }
     const T *getConstPointer() const { return _mem.getConstPointer(); }
     const T *begin() const { return getConstPointer(); }
     const T *end() const { return getConstPointer()+getNbOfElems(); }
     void alloc(std::size_t nbOfTuple, std::size_t nbOfCompo=1);
     void useArray(const T *array, bool ownership, DeallocType type, int nbOfTuple, int nbOfCompo);
     void useExternalArrayWithRWAccess(const T *array, int nbOfTuple, int nbOfCompo);
     T getIJSafe(int tupleId, int compoId) const;
     T getIJ(int tupleId, int compoId) const { return _mem[tupleId*_info_on_compo.size()+compoId]; }
     void setIJ(int tupleId, int compoId, T newVal) { _mem[tupleId*_info_on_compo.size()+compoId]=newVal; declareAsNew(); }
     void setIJSilent(int tupleId, int compoId, T newVal) { _mem[tupleId*_info_on_compo.size()+compoId]=newVal; }
     T *getPointer() { return _mem.getPointer(); declareAsNew(); }
     void pack() const;
     bool isAllocated() const;
     void checkAllocated() const;
     void desallocate();
     void reserve(std::size_t nbOfElems);
     void rearrange(int newNbOfCompo);
     void transpose();
     void pushBackSilent(T val);
     void pushBackValsSilent(const T *valsBg, const T *valsEnd);
     T popBackSilent();
     T front() const;
     T back() const;
     std::size_t getNbOfElemAllocated() const { return _mem.getNbOfElemAllocated(); }
     void allocIfNecessary(int nbOfTuple, int nbOfCompo);
     void deepCopyFrom(const DataArrayTemplate<T>& other);
     void reverse();
     void fillWithValue(T val);
     void reAlloc(std::size_t newNbOfTuple);
     void renumberInPlace(const int *old2New);
     void renumberInPlaceR(const int *new2Old);
     void sort(bool asc=true);
     typename Traits<T>::ArrayType *renumber(const int *old2New) const;
     typename Traits<T>::ArrayType *renumberR(const int *new2Old) const;
     typename Traits<T>::ArrayType *renumberAndReduce(const int *old2New, int newNbOfTuple) const;
     typename Traits<T>::ArrayType *changeNbOfComponents(int newNbOfComp, T dftValue) const;
     typename Traits<T>::ArrayType *subArray(int tupleIdBg, int tupleIdEnd=-1) const;
     MCAuto<typename Traits<T>::ArrayTypeCh> selectPartDef(const PartDefinition* pd) const;
     void circularPermutation(int nbOfShift=1);
     void circularPermutationPerTuple(int nbOfShift=1);
     void reversePerTuple();
     void setPartOfValues1(const typename Traits<T>::ArrayType *a, int bgTuples, int endTuples, int stepTuples, int bgComp, int endComp, int stepComp, bool strictCompoCompare=true);
     void setPartOfValuesSimple1(T a, int bgTuples, int endTuples, int stepTuples, int bgComp, int endComp, int stepComp);
     void setPartOfValues2(const typename Traits<T>::ArrayType *a, const int *bgTuples, const int *endTuples, const int *bgComp, const int *endComp, bool strictCompoCompare=true);
     void setPartOfValuesSimple2(T a, const int *bgTuples, const int *endTuples, const int *bgComp, const int *endComp);
     void setPartOfValues3(const typename Traits<T>::ArrayType *a, const int *bgTuples, const int *endTuples, int bgComp, int endComp, int stepComp, bool strictCompoCompare=true);
     void setPartOfValuesSimple3(T a, const int *bgTuples, const int *endTuples, int bgComp, int endComp, int stepComp);
     void setPartOfValues4(const typename Traits<T>::ArrayType *a, int bgTuples, int endTuples, int stepTuples, const int *bgComp, const int *endComp, bool strictCompoCompare=true);
     void setPartOfValuesSimple4(T a, int bgTuples, int endTuples, int stepTuples, const int *bgComp, const int *endComp);
     void setPartOfValuesAdv(const typename Traits<T>::ArrayType *a, const DataArrayInt *tuplesSelec);
     void setContigPartOfSelectedValues(int tupleIdStart, const DataArray *aBase, const DataArrayInt *tuplesSelec);
     void setContigPartOfSelectedValuesSlice(int tupleIdStart, const DataArray *aBase, int bg, int end2, int step);
     T getMaxValue(int& tupleId) const;
     T getMaxValueInArray() const;
     T getMinValue(int& tupleId) const;
     T getMinValueInArray() const;
  protected:
    typename Traits<T>::ArrayType *mySelectByTupleId(const int *new2OldBg, const int *new2OldEnd) const;
    typename Traits<T>::ArrayType *mySelectByTupleId(const DataArrayInt& di) const;
    typename Traits<T>::ArrayType *mySelectByTupleIdSafe(const int *new2OldBg, const int *new2OldEnd) const;
    typename Traits<T>::ArrayType *myKeepSelectedComponents(const std::vector<int>& compoIds) const;
    typename Traits<T>::ArrayType *mySelectByTupleIdSafeSlice(int bg, int end2, int step) const;
    typename Traits<T>::ArrayType *mySelectByTupleRanges(const std::vector<std::pair<int,int> >& ranges) const;
  protected:
    MemArray<T> _mem;
  };
}

namespace MEDCoupling
{
  class DataArrayDoubleIterator;
  class MEDCOUPLING_EXPORT DataArrayDouble : public DataArrayTemplate<double>
  {
  public:
     static DataArrayDouble *New();
     double doubleValue() const;
     DataArrayDouble *deepCopy() const;
     DataArrayDouble *buildNewEmptyInstance() const { return DataArrayDouble::New(); }
     DataArrayDouble *performCopyOrIncrRef(bool deepCopy) const;
     void fillWithZero();
     void iota(double init=0.);
     bool isUniform(double val, double eps) const;
     void checkMonotonic(bool increasing, double eps) const;
     bool isMonotonic(bool increasing, double eps) const;
     std::string repr() const;
     std::string reprZip() const;
     std::string reprNotTooLong() const;
     void writeVTK(std::ostream& ofs, int indent, const std::string& nameInFile, DataArrayByte *byteArr) const;
     void reprStream(std::ostream& stream) const;
     void reprZipStream(std::ostream& stream) const;
     void reprNotTooLongStream(std::ostream& stream) const;
     void reprWithoutNameStream(std::ostream& stream) const;
     void reprZipWithoutNameStream(std::ostream& stream) const;
     void reprNotTooLongWithoutNameStream(std::ostream& stream) const;
     void reprCppStream(const std::string& varName, std::ostream& stream) const;
     void reprQuickOverview(std::ostream& stream) const;
     void reprQuickOverviewData(std::ostream& stream, std::size_t maxNbOfByteInRepr) const;
     bool isEqual(const DataArrayDouble& other, double prec) const;
     bool isEqualIfNotWhy(const DataArrayDouble& other, double prec, std::string& reason) const;
     bool isEqualWithoutConsideringStr(const DataArrayDouble& other, double prec) const;
     DataArrayInt *convertToIntArr() const;
     DataArrayDouble *fromNoInterlace() const;
     DataArrayDouble *toNoInterlace() const;
     DataArrayDouble *selectByTupleId(const int *new2OldBg, const int *new2OldEnd) const { return DataArrayTemplate<double>::mySelectByTupleId(new2OldBg,new2OldEnd); }
     DataArrayDouble *selectByTupleId(const DataArrayInt& di) const { return DataArrayTemplate<double>::mySelectByTupleId(di); }
     DataArrayDouble *selectByTupleIdSafe(const int *new2OldBg, const int *new2OldEnd) const { return DataArrayTemplate<double>::mySelectByTupleIdSafe(new2OldBg,new2OldEnd); }
     DataArrayDouble *keepSelectedComponents(const std::vector<int>& compoIds) const { return DataArrayTemplate<double>::myKeepSelectedComponents(compoIds); }
     DataArrayDouble *selectByTupleIdSafeSlice(int bg, int end2, int step) const { return DataArrayTemplate<double>::mySelectByTupleIdSafeSlice(bg,end2,step); }
     DataArrayDouble *selectByTupleRanges(const std::vector<std::pair<int,int> >& ranges) const { return DataArrayTemplate<double>::mySelectByTupleRanges(ranges); }
     void meldWith(const DataArrayDouble *other);
     bool areIncludedInMe(const DataArrayDouble *other, double prec, DataArrayInt *&tupleIds) const;
     void findCommonTuples(double prec, int limitTupleId, DataArrayInt *&comm, DataArrayInt *&commIndex) const;
     double minimalDistanceTo(const DataArrayDouble *other, int& thisTupleId, int& otherTupleId) const;
     DataArrayDouble *duplicateEachTupleNTimes(int nbTimes) const;
     DataArrayDouble *getDifferentValues(double prec, int limitTupleId=-1) const;
     DataArrayInt *findClosestTupleId(const DataArrayDouble *other) const;
     DataArrayInt *computeNbOfInteractionsWith(const DataArrayDouble *otherBBoxFrmt, double eps) const;
     void setSelectedComponents(const DataArrayDouble *a, const std::vector<int>& compoIds);
     void getTuple(int tupleId, double *res) const { std::copy(_mem.getConstPointerLoc(tupleId*_info_on_compo.size()),_mem.getConstPointerLoc((tupleId+1)*_info_on_compo.size()),res); }
     static void SetArrayIn(DataArrayDouble *newArray, DataArrayDouble* &arrayToSet);
     DataArrayDoubleIterator *iterator();
     template<class InputIterator>
     void insertAtTheEnd(InputIterator first, InputIterator last);
     void aggregate(const DataArrayDouble *other);
     void writeOnPlace(std::size_t id, double element0, const double *others, int sizeOfOthers) { _mem.writeOnPlace(id,element0,others,sizeOfOthers); }
     void checkNoNullValues() const;
     void getMinMaxPerComponent(double *bounds) const;
     DataArrayDouble *computeBBoxPerTuple(double epsilon=0.0) const;
     void computeTupleIdsNearTuples(const DataArrayDouble *other, double eps, DataArrayInt *& c, DataArrayInt *& cI) const;
     void recenterForMaxPrecision(double eps);
     double getMaxValue2(DataArrayInt*& tupleIds) const;
     double getMinValue2(DataArrayInt*& tupleIds) const;
     int count(double value, double eps) const;
     double getAverageValue() const;
     double norm2() const;
     double normMax() const;
     double normMin() const;
     void accumulate(double *res) const;
     double accumulate(int compId) const;
     DataArrayDouble *accumulatePerChunck(const int *bgOfIndex, const int *endOfIndex) const;
     MCAuto<DataArrayDouble> cumSum() const;
     double distanceToTuple(const double *tupleBg, const double *tupleEnd, int& tupleId) const;
     DataArrayDouble *fromPolarToCart() const;
     DataArrayDouble *fromCylToCart() const;
     DataArrayDouble *fromSpherToCart() const;
     DataArrayDouble *cartesianize(MEDCouplingAxisType atOfThis) const;
     DataArrayDouble *fromCartToPolar() const;
     DataArrayDouble *fromCartToCyl() const;
     DataArrayDouble *fromCartToSpher() const;
     DataArrayDouble *fromCartToCylGiven(const DataArrayDouble *coords, const double center[3], const double vect[3]) const;
     DataArrayDouble *doublyContractedProduct() const;
     DataArrayDouble *determinant() const;
     DataArrayDouble *eigenValues() const;
     DataArrayDouble *eigenVectors() const;
     DataArrayDouble *inverse() const;
     DataArrayDouble *trace() const;
     DataArrayDouble *deviator() const;
     DataArrayDouble *magnitude() const;
     DataArrayDouble *sumPerTuple() const;
     DataArrayDouble *maxPerTuple() const;
     DataArrayDouble *maxPerTupleWithCompoId(DataArrayInt* &compoIdOfMaxPerTuple) const;
     DataArrayDouble *buildEuclidianDistanceDenseMatrix() const;
     DataArrayDouble *buildEuclidianDistanceDenseMatrixWith(const DataArrayDouble *other) const;
     void sortPerTuple(bool asc);
     void abs();
     DataArrayDouble *computeAbs() const;
     void applyLin(double a, double b, int compoId);
     void applyLin(double a, double b);
     void applyInv(double numerator);
     void applyPow(double val);
     void applyRPow(double val);
     DataArrayDouble *negate() const;
     DataArrayDouble *applyFunc(int nbOfComp, FunctionToEvaluate func) const;
     DataArrayDouble *applyFunc(int nbOfComp, const std::string& func, bool isSafe=true) const;
     DataArrayDouble *applyFunc(const std::string& func, bool isSafe=true) const;
     void applyFuncOnThis(const std::string& func, bool isSafe=true);
     DataArrayDouble *applyFuncCompo(int nbOfComp, const std::string& func, bool isSafe=true) const;
     DataArrayDouble *applyFuncNamedCompo(int nbOfComp, const std::vector<std::string>& varsOrder, const std::string& func, bool isSafe=true) const;
     void applyFuncFast32(const std::string& func);
     void applyFuncFast64(const std::string& func);
     MCAuto<DataArrayDouble> symmetry3DPlane(const double point[3], const double normalVector[3]) const;
     DataArrayInt *findIdsInRange(double vmin, double vmax) const;
     DataArrayInt *findIdsNotInRange(double vmin, double vmax) const;
     static DataArrayDouble *Aggregate(const DataArrayDouble *a1, const DataArrayDouble *a2);
     static DataArrayDouble *Aggregate(const std::vector<const DataArrayDouble *>& arr);
     static DataArrayDouble *Meld(const DataArrayDouble *a1, const DataArrayDouble *a2);
     static DataArrayDouble *Meld(const std::vector<const DataArrayDouble *>& arr);
     static DataArrayDouble *Dot(const DataArrayDouble *a1, const DataArrayDouble *a2);
     static DataArrayDouble *CrossProduct(const DataArrayDouble *a1, const DataArrayDouble *a2);
     static DataArrayDouble *Max(const DataArrayDouble *a1, const DataArrayDouble *a2);
     static DataArrayDouble *Min(const DataArrayDouble *a1, const DataArrayDouble *a2);
     static DataArrayDouble *Add(const DataArrayDouble *a1, const DataArrayDouble *a2);
     void addEqual(const DataArrayDouble *other);
     static DataArrayDouble *Substract(const DataArrayDouble *a1, const DataArrayDouble *a2);
     void substractEqual(const DataArrayDouble *other);
     static DataArrayDouble *Multiply(const DataArrayDouble *a1, const DataArrayDouble *a2);
     void multiplyEqual(const DataArrayDouble *other);
     static DataArrayDouble *Divide(const DataArrayDouble *a1, const DataArrayDouble *a2);
     void divideEqual(const DataArrayDouble *other);
     static DataArrayDouble *Pow(const DataArrayDouble *a1, const DataArrayDouble *a2);
     void powEqual(const DataArrayDouble *other);
     void updateTime() const { }
     MemArray<double>& accessToMemArray() { return _mem; }
     const MemArray<double>& accessToMemArray() const { return _mem; }
     std::vector<bool> toVectorOfBool(double eps) const;
     static void Rotate2DAlg(const double *center, double angle, int nbNodes, const double *coordsIn, double *coordsOut);
     static void Rotate3DAlg(const double *center, const double *vect, double angle, int nbNodes, const double *coordsIn, double *coordsOut);
     static void Symmetry3DPlane(const double point[3], const double normalVector[3], int nbNodes, const double *coordsIn, double *coordsOut);
     static void GiveBaseForPlane(const double normalVector[3], double baseOfPlane[9]);
  public:
     void getTinySerializationIntInformation(std::vector<int>& tinyInfo) const;
     void getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const;
     bool resizeForUnserialization(const std::vector<int>& tinyInfoI);
     void finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<std::string>& tinyInfoS);
  public:
    template<int SPACEDIM>
    void findCommonTuplesAlg(const double *bbox, int nbNodes, int limitNodeId, double prec, DataArrayInt *c, DataArrayInt *cI) const;
    template<int SPACEDIM>
    static void FindClosestTupleIdAlg(const BBTreePts<SPACEDIM,int>& myTree, double dist, const double *pos, int nbOfTuples, const double *thisPt, int thisNbOfTuples, int *res);
    template<int SPACEDIM>
    static void FindTupleIdsNearTuplesAlg(const BBTreePts<SPACEDIM,int>& myTree, const double *pos, int nbOfTuples, double eps,
                                          DataArrayInt *c, DataArrayInt *cI);
  private:
    ~DataArrayDouble() { }
    DataArrayDouble() { }
  };

  class DataArrayDoubleTuple;

  class MEDCOUPLING_EXPORT DataArrayDoubleIterator
  {
  public:
     DataArrayDoubleIterator(DataArrayDouble *da);
     ~DataArrayDoubleIterator();
     DataArrayDoubleTuple *nextt();
  private:
    DataArrayDouble *_da;
    double *_pt;
    int _tuple_id;
    int _nb_comp;
    int _nb_tuple;
  };

  class MEDCOUPLING_EXPORT DataArrayDoubleTuple
  {
  public:
     DataArrayDoubleTuple(double *pt, int nbOfComp);
     std::string repr() const;
     int getNumberOfCompo() const { return _nb_of_compo; }
     const double *getConstPointer() const { return  _pt; }
     double *getPointer() { return _pt; }
     double doubleValue() const;
     DataArrayDouble *buildDADouble(int nbOfTuples, int nbOfCompo) const;
  private:
    double *_pt;
    int _nb_of_compo;
  };

  class DataArrayIntIterator;

  class MEDCOUPLING_EXPORT DataArrayInt : public DataArrayTemplate<int>
  {
  public:
     static DataArrayInt *New();
     int intValue() const;
     int getHashCode() const;
     DataArrayInt *deepCopy() const;
     DataArrayInt *performCopyOrIncrRef(bool deepCopy) const;
     DataArrayInt *buildNewEmptyInstance() const { return DataArrayInt::New(); }
     bool isEqual(const DataArrayInt& other) const;
     bool isEqualIfNotWhy(const DataArrayInt& other, std::string& reason) const;
     bool isEqualWithoutConsideringStr(const DataArrayInt& other) const;
     bool isEqualWithoutConsideringStrAndOrder(const DataArrayInt& other) const;
     bool isFittingWith(const std::vector<bool>& v) const;
     void switchOnTupleEqualTo(int val, std::vector<bool>& vec) const;
     void switchOnTupleNotEqualTo(int val, std::vector<bool>& vec) const;
     DataArrayInt *buildPermutationArr(const DataArrayInt& other) const;
     DataArrayInt *indicesOfSubPart(const DataArrayInt& partOfThis) const;
     DataArrayInt *sumPerTuple() const;
     void checkMonotonic(bool increasing) const;
     bool isMonotonic(bool increasing) const;
     void checkStrictlyMonotonic(bool increasing) const;
     bool isStrictlyMonotonic(bool increasing) const;
     void fillWithZero();
     void iota(int init=0);
     std::string repr() const;
     std::string reprZip() const;
     std::string reprNotTooLong() const;
     void writeVTK(std::ostream& ofs, int indent, const std::string& type, const std::string& nameInFile, DataArrayByte *byteArr) const;
     void reprStream(std::ostream& stream) const;
     void reprZipStream(std::ostream& stream) const;
     void reprNotTooLongStream(std::ostream& stream) const;
     void reprWithoutNameStream(std::ostream& stream) const;
     void reprZipWithoutNameStream(std::ostream& stream) const;
     void reprNotTooLongWithoutNameStream(std::ostream& stream) const;
     void reprCppStream(const std::string& varName, std::ostream& stream) const;
     void reprQuickOverview(std::ostream& stream) const;
     void reprQuickOverviewData(std::ostream& stream, std::size_t maxNbOfByteInRepr) const;
     void transformWithIndArr(const int *indArrBg, const int *indArrEnd);
     DataArrayInt *transformWithIndArrR(const int *indArrBg, const int *indArrEnd) const;
     void splitByValueRange(const int *arrBg, const int *arrEnd,
                                              DataArrayInt *& castArr, DataArrayInt *& rankInsideCast, DataArrayInt *& castsPresent) const;
     bool isRange(int& strt, int& sttoopp, int& stteepp) const;
     DataArrayInt *invertArrayO2N2N2O(int newNbOfElem) const;
     DataArrayInt *invertArrayN2O2O2N(int oldNbOfElem) const;
     DataArrayInt *invertArrayO2N2N2OBis(int newNbOfElem) const;
     DataArrayDouble *convertToDblArr() const;
     DataArrayInt *fromNoInterlace() const;
     DataArrayInt *toNoInterlace() const;
     DataArrayInt *selectByTupleId(const int *new2OldBg, const int *new2OldEnd) const { return DataArrayTemplate<int>::mySelectByTupleId(new2OldBg,new2OldEnd); }
     DataArrayInt *selectByTupleId(const DataArrayInt& di) const { return DataArrayTemplate<int>::mySelectByTupleId(di); }
     DataArrayInt *selectByTupleIdSafe(const int *new2OldBg, const int *new2OldEnd) const { return DataArrayTemplate<int>::mySelectByTupleIdSafe(new2OldBg,new2OldEnd); }
     DataArrayInt *keepSelectedComponents(const std::vector<int>& compoIds) const { return DataArrayTemplate<int>::myKeepSelectedComponents(compoIds); }
     DataArrayInt *selectByTupleIdSafeSlice(int bg, int end2, int step) const { return DataArrayTemplate<int>::mySelectByTupleIdSafeSlice(bg,end2,step); }
     DataArrayInt *selectByTupleRanges(const std::vector<std::pair<int,int> >& ranges) const { return DataArrayTemplate<int>::mySelectByTupleRanges(ranges); }
     DataArrayInt *checkAndPreparePermutation() const;
     static DataArrayInt *FindPermutationFromFirstToSecond(const DataArrayInt *ids1, const DataArrayInt *ids2);
     void changeSurjectiveFormat(int targetNb, DataArrayInt *&arr, DataArrayInt *&arrI) const;
     static DataArrayInt *ConvertIndexArrayToO2N(int nbOfOldTuples, const int *arr, const int *arrIBg, const int *arrIEnd, int &newNbOfTuples);
     DataArrayInt *buildPermArrPerLevel() const;
     bool isIota(int sizeExpected) const;
     bool isUniform(int val) const;
     bool hasUniqueValues() const;
     void meldWith(const DataArrayInt *other);
     void setSelectedComponents(const DataArrayInt *a, const std::vector<int>& compoIds);
     void getTuple(int tupleId, int *res) const { std::copy(_mem.getConstPointerLoc(tupleId*_info_on_compo.size()),_mem.getConstPointerLoc((tupleId+1)*_info_on_compo.size()),res); }
     static void SetArrayIn(DataArrayInt *newArray, DataArrayInt* &arrayToSet);
     DataArrayIntIterator *iterator();
     DataArrayInt *findIdsEqual(int val) const;
     DataArrayInt *findIdsNotEqual(int val) const;
     DataArrayInt *findIdsEqualList(const int *valsBg, const int *valsEnd) const;
     DataArrayInt *findIdsNotEqualList(const int *valsBg, const int *valsEnd) const;
     DataArrayInt *findIdsEqualTuple(const int *tupleBg, const int *tupleEnd) const;
     int changeValue(int oldValue, int newValue);
     int findIdFirstEqualTuple(const std::vector<int>& tupl) const;
     int findIdFirstEqual(int value) const;
     int findIdFirstEqual(const std::vector<int>& vals) const;
     int findIdSequence(const std::vector<int>& vals) const;
     bool presenceOfTuple(const std::vector<int>& tupl) const;
     bool presenceOfValue(int value) const;
     bool presenceOfValue(const std::vector<int>& vals) const;
     int count(int value) const;
     void accumulate(int *res) const;
     int accumulate(int compId) const;
     DataArrayInt *accumulatePerChunck(const int *bgOfIndex, const int *endOfIndex) const;
     void getMinMaxValues(int& minValue, int& maxValue) const;
     void abs();
     DataArrayInt *computeAbs() const;
     void applyLin(int a, int b, int compoId);
     void applyLin(int a, int b);
     void applyInv(int numerator);
     DataArrayInt *negate() const;
     void applyDivideBy(int val);
     void applyModulus(int val);
     void applyRModulus(int val);
     void applyPow(int val);
     void applyRPow(int val);
     DataArrayInt *findIdsInRange(int vmin, int vmax) const;
     DataArrayInt *findIdsNotInRange(int vmin, int vmax) const;
     DataArrayInt *findIdsStricltyNegative() const;
     bool checkAllIdsInRange(int vmin, int vmax) const;
     static DataArrayInt *Aggregate(const DataArrayInt *a1, const DataArrayInt *a2, int offsetA2);
     static DataArrayInt *Aggregate(const std::vector<const DataArrayInt *>& arr);
     static DataArrayInt *AggregateIndexes(const std::vector<const DataArrayInt *>& arrs);
     static DataArrayInt *Meld(const DataArrayInt *a1, const DataArrayInt *a2);
     static DataArrayInt *Meld(const std::vector<const DataArrayInt *>& arr);
     static DataArrayInt *MakePartition(const std::vector<const DataArrayInt *>& groups, int newNb, std::vector< std::vector<int> >& fidsOfGroups);
     static DataArrayInt *BuildUnion(const std::vector<const DataArrayInt *>& arr);
     static DataArrayInt *BuildIntersection(const std::vector<const DataArrayInt *>& arr);
     static DataArrayInt *BuildListOfSwitchedOn(const std::vector<bool>& v);
     static DataArrayInt *BuildListOfSwitchedOff(const std::vector<bool>& v);
     static void PutIntoToSkylineFrmt(const std::vector< std::vector<int> >& v, DataArrayInt *& data, DataArrayInt *& dataIndex);
     DataArrayInt *buildComplement(int nbOfElement) const;
     DataArrayInt *buildSubstraction(const DataArrayInt *other) const;
     DataArrayInt *buildSubstractionOptimized(const DataArrayInt *other) const;
     DataArrayInt *buildUnion(const DataArrayInt *other) const;
     DataArrayInt *buildIntersection(const DataArrayInt *other) const;
     DataArrayInt *buildUnique() const;
     DataArrayInt *buildUniqueNotSorted() const;
     DataArrayInt *deltaShiftIndex() const;
     void computeOffsets();
     void computeOffsetsFull();
     void findIdsRangesInListOfIds(const DataArrayInt *listOfIds, DataArrayInt *& rangeIdsFetched, DataArrayInt *& idsInInputListThatFetch) const;
     DataArrayInt *buildExplicitArrByRanges(const DataArrayInt *offsets) const;
     DataArrayInt *buildExplicitArrOfSliceOnScaledArr(int begin, int stop, int step) const;
     DataArrayInt *findRangeIdForEachTuple(const DataArrayInt *ranges) const;
     DataArrayInt *findIdInRangeForEachTuple(const DataArrayInt *ranges) const;
     void sortEachPairToMakeALinkedList();
     MCAuto<DataArrayInt> fromLinkedListOfPairToList() const;
     DataArrayInt *duplicateEachTupleNTimes(int nbTimes) const;
     DataArrayInt *getDifferentValues() const;
     std::vector<DataArrayInt *> partitionByDifferentValues(std::vector<int>& differentIds) const;
     std::vector< std::pair<int,int> > splitInBalancedSlices(int nbOfSlices) const;
     template<class InputIterator>
     void insertAtTheEnd(InputIterator first, InputIterator last);
     void aggregate(const DataArrayInt *other);
     void writeOnPlace(std::size_t id, int element0, const int *others, int sizeOfOthers) { _mem.writeOnPlace(id,element0,others,sizeOfOthers); }
     static DataArrayInt *Add(const DataArrayInt *a1, const DataArrayInt *a2);
     void addEqual(const DataArrayInt *other);
     static DataArrayInt *Substract(const DataArrayInt *a1, const DataArrayInt *a2);
     void substractEqual(const DataArrayInt *other);
     static DataArrayInt *Multiply(const DataArrayInt *a1, const DataArrayInt *a2);
     void multiplyEqual(const DataArrayInt *other);
     static DataArrayInt *Divide(const DataArrayInt *a1, const DataArrayInt *a2);
     void divideEqual(const DataArrayInt *other);
     static DataArrayInt *Modulus(const DataArrayInt *a1, const DataArrayInt *a2);
     void modulusEqual(const DataArrayInt *other);
     static DataArrayInt *Pow(const DataArrayInt *a1, const DataArrayInt *a2);
     void powEqual(const DataArrayInt *other);
     void updateTime() const { }
     MemArray<int>& accessToMemArray() { return _mem; }
     const MemArray<int>& accessToMemArray() const { return _mem; }
  public:
     static int *CheckAndPreparePermutation(const int *start, const int *end);
     static DataArrayInt *Range(int begin, int end, int step);
  public:
     void getTinySerializationIntInformation(std::vector<int>& tinyInfo) const;
     void getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const;
     bool resizeForUnserialization(const std::vector<int>& tinyInfoI);
     void finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<std::string>& tinyInfoS);
  private:
    ~DataArrayInt() { }
    DataArrayInt() { }
  };

  class DataArrayIntTuple;

  class MEDCOUPLING_EXPORT DataArrayIntIterator
  {
  public:
     DataArrayIntIterator(DataArrayInt *da);
     ~DataArrayIntIterator();
     DataArrayIntTuple *nextt();
  private:
    DataArrayInt *_da;
    int *_pt;
    int _tuple_id;
    int _nb_comp;
    int _nb_tuple;
  };

  class MEDCOUPLING_EXPORT DataArrayIntTuple
  {
  public:
     DataArrayIntTuple(int *pt, int nbOfComp);
     std::string repr() const;
     int getNumberOfCompo() const { return _nb_of_compo; }
     const int *getConstPointer() const { return  _pt; }
     int *getPointer() { return _pt; }
     int intValue() const;
     DataArrayInt *buildDAInt(int nbOfTuples, int nbOfCompo) const;
  private:
    int *_pt;
    int _nb_of_compo;
  };

  class MEDCOUPLING_EXPORT DataArrayChar : public DataArrayTemplate<char>
  {
  public:
     virtual DataArrayChar *buildEmptySpecializedDAChar() const = 0;
     int getHashCode() const;
     bool isEqual(const DataArrayChar& other) const;
     virtual bool isEqualIfNotWhy(const DataArrayChar& other, std::string& reason) const;
     bool isEqualWithoutConsideringStr(const DataArrayChar& other) const;
     void fillWithZero();
     std::string repr() const;
     std::string reprZip() const;
     DataArrayInt *convertToIntArr() const;
     DataArrayChar *selectByTupleId(const int *new2OldBg, const int *new2OldEnd) const { return DataArrayTemplate<char>::mySelectByTupleId(new2OldBg,new2OldEnd); }
     DataArrayChar *selectByTupleId(const DataArrayInt& di) const { return DataArrayTemplate<char>::mySelectByTupleId(di); }
     DataArrayChar *selectByTupleIdSafe(const int *new2OldBg, const int *new2OldEnd) const { return DataArrayTemplate<char>::mySelectByTupleIdSafe(new2OldBg,new2OldEnd); }
     DataArrayChar *keepSelectedComponents(const std::vector<int>& compoIds) const { return DataArrayTemplate<char>::myKeepSelectedComponents(compoIds); }
     DataArrayChar *selectByTupleIdSafeSlice(int bg, int end2, int step) const { return DataArrayTemplate<char>::mySelectByTupleIdSafeSlice(bg,end2,step); }
     bool isUniform(char val) const;
     void meldWith(const DataArrayChar *other);
     DataArray *selectByTupleRanges(const std::vector<std::pair<int,int> >& ranges) const { return DataArrayTemplate<char>::mySelectByTupleRanges(ranges); }
     void getTuple(int tupleId, char *res) const { std::copy(_mem.getConstPointerLoc(tupleId*_info_on_compo.size()),_mem.getConstPointerLoc((tupleId+1)*_info_on_compo.size()),res); }
     DataArrayInt *findIdsEqual(char val) const;
     DataArrayInt *findIdsNotEqual(char val) const;
     int findIdSequence(const std::vector<char>& vals) const;
     int findIdFirstEqualTuple(const std::vector<char>& tupl) const;
     int findIdFirstEqual(char value) const;
     int findIdFirstEqual(const std::vector<char>& vals) const;
     bool presenceOfTuple(const std::vector<char>& tupl) const;
     bool presenceOfValue(char value) const;
     bool presenceOfValue(const std::vector<char>& vals) const;
     DataArrayInt *findIdsInRange(char vmin, char vmax) const;
     static DataArrayChar *Aggregate(const DataArrayChar *a1, const DataArrayChar *a2);
     static DataArrayChar *Aggregate(const std::vector<const DataArrayChar *>& arr);
     static DataArrayChar *Meld(const DataArrayChar *a1, const DataArrayChar *a2);
     static DataArrayChar *Meld(const std::vector<const DataArrayChar *>& arr);
     template<class InputIterator>
     void insertAtTheEnd(InputIterator first, InputIterator last);
     void updateTime() const { }
     MemArray<char>& accessToMemArray() { return _mem; }
     const MemArray<char>& accessToMemArray() const { return _mem; }
  public:
    // void getTinySerializationIntInformation(std::vector<int>& tinyInfo) const;
    // void getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const;
    // bool resizeForUnserialization(const std::vector<int>& tinyInfoI);
    // void finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<std::string>& tinyInfoS);
  protected:
    DataArrayChar() { }
  };

  class DataArrayByteIterator;

  class MEDCOUPLING_EXPORT DataArrayByte : public DataArrayChar
  {
  public:
     static DataArrayByte *New();
     DataArrayChar *buildEmptySpecializedDAChar() const;
     DataArrayByteIterator *iterator();
     DataArrayByte *deepCopy() const;
     DataArrayByte *performCopyOrIncrRef(bool deepCopy) const;
     DataArrayByte *buildNewEmptyInstance() const { return DataArrayByte::New(); }
     char byteValue() const;
     void reprStream(std::ostream& stream) const;
     void reprZipStream(std::ostream& stream) const;
     void reprWithoutNameStream(std::ostream& stream) const;
     void reprZipWithoutNameStream(std::ostream& stream) const;
     void reprCppStream(const std::string& varName, std::ostream& stream) const;
     void reprQuickOverview(std::ostream& stream) const;
     void reprQuickOverviewData(std::ostream& stream, std::size_t maxNbOfByteInRepr) const;
     bool isEqualIfNotWhy(const DataArrayChar& other, std::string& reason) const;
     std::vector<bool> toVectorOfBool() const;
  private:
    ~DataArrayByte() { }
    DataArrayByte() { }
  };

  class DataArrayByteTuple;

  class MEDCOUPLING_EXPORT DataArrayByteIterator
  {
  public:
     DataArrayByteIterator(DataArrayByte *da);
     ~DataArrayByteIterator();
     DataArrayByteTuple *nextt();
  private:
    DataArrayByte *_da;
    char *_pt;
    int _tuple_id;
    int _nb_comp;
    int _nb_tuple;
  };

  class MEDCOUPLING_EXPORT DataArrayByteTuple
  {
  public:
     DataArrayByteTuple(char *pt, int nbOfComp);
     std::string repr() const;
     int getNumberOfCompo() const { return _nb_of_compo; }
     const char *getConstPointer() const { return  _pt; }
     char *getPointer() { return _pt; }
     char byteValue() const;
     DataArrayByte *buildDAByte(int nbOfTuples, int nbOfCompo) const;
  private:
    char *_pt;
    int _nb_of_compo;
  };

  class DataArrayAsciiCharIterator;

  class MEDCOUPLING_EXPORT DataArrayAsciiChar : public DataArrayChar
  {
  public:
     static DataArrayAsciiChar *New();
     static DataArrayAsciiChar *New(const std::string& st);
     static DataArrayAsciiChar *New(const std::vector<std::string>& vst, char defaultChar);
     DataArrayChar *buildEmptySpecializedDAChar() const;
     DataArrayAsciiCharIterator *iterator();
     DataArrayAsciiChar *deepCopy() const;
     DataArrayAsciiChar *performCopyOrIncrRef(bool deepCopy) const;
     DataArrayAsciiChar *buildNewEmptyInstance() const { return DataArrayAsciiChar::New(); }
     char asciiCharValue() const;
     void reprStream(std::ostream& stream) const;
     void reprZipStream(std::ostream& stream) const;
     void reprWithoutNameStream(std::ostream& stream) const;
     void reprZipWithoutNameStream(std::ostream& stream) const;
     void reprCppStream(const std::string& varName, std::ostream& stream) const;
     void reprQuickOverview(std::ostream& stream) const;
     void reprQuickOverviewData(std::ostream& stream, std::size_t maxNbOfByteInRepr) const;
     bool isEqualIfNotWhy(const DataArrayChar& other, std::string& reason) const;
  private:
    ~DataArrayAsciiChar() { }
    DataArrayAsciiChar() { }
    DataArrayAsciiChar(const std::string& st);
    DataArrayAsciiChar(const std::vector<std::string>& vst, char defaultChar);
  };

  class DataArrayAsciiCharTuple;

  class MEDCOUPLING_EXPORT DataArrayAsciiCharIterator
  {
  public:
     DataArrayAsciiCharIterator(DataArrayAsciiChar *da);
     ~DataArrayAsciiCharIterator();
     DataArrayAsciiCharTuple *nextt();
  private:
    DataArrayAsciiChar *_da;
    char *_pt;
    int _tuple_id;
    int _nb_comp;
    int _nb_tuple;
  };

  class MEDCOUPLING_EXPORT DataArrayAsciiCharTuple
  {
  public:
     DataArrayAsciiCharTuple(char *pt, int nbOfComp);
     std::string repr() const;
     int getNumberOfCompo() const { return _nb_of_compo; }
     const char *getConstPointer() const { return  _pt; }
     char *getPointer() { return _pt; }
     char asciiCharValue() const;
     DataArrayAsciiChar *buildDAAsciiChar(int nbOfTuples, int nbOfCompo) const;
  private:
    char *_pt;
    int _nb_of_compo;
  };
}

namespace MEDCoupling
{
  template<class T>
  template<class InputIterator>
  void MemArray<T>::insertAtTheEnd(InputIterator first, InputIterator last)
  {
    T *pointer=_pointer.getPointer();
    while(first!=last)
      {
        if(_nb_of_elem>=_nb_of_elem_alloc)
          {
            reserve(_nb_of_elem_alloc>0?2*_nb_of_elem_alloc:1);
            pointer=_pointer.getPointer();
          }
        pointer[_nb_of_elem++]=*first++;
      }
  }
  
  template<class InputIterator>
  void DataArrayDouble::insertAtTheEnd(InputIterator first, InputIterator last)
  {
    int nbCompo(getNumberOfComponents());
    if(nbCompo==1)
      _mem.insertAtTheEnd(first,last);
    else if(nbCompo==0)
      {
        _info_on_compo.resize(1);
        _mem.insertAtTheEnd(first,last);
      }
    else
      throw INTERP_KERNEL::Exception("DataArrayDouble::insertAtTheEnd : not available for DataArrayDouble with number of components different than 1 !");
  }

  template<class InputIterator>
  void DataArrayInt::insertAtTheEnd(InputIterator first, InputIterator last)
  {
    int nbCompo(getNumberOfComponents());
    if(nbCompo==1)
      _mem.insertAtTheEnd(first,last);
    else if(nbCompo==0)
      {
        _info_on_compo.resize(1);
        _mem.insertAtTheEnd(first,last);
      }
    else
      throw INTERP_KERNEL::Exception("DataArrayInt::insertAtTheEnd : not available for DataArrayInt with number of components different than 1 !");
  }

  template<class InputIterator>
  void DataArrayChar::insertAtTheEnd(InputIterator first, InputIterator last)
  {
    int nbCompo(getNumberOfComponents());
    if(nbCompo==1)
      _mem.insertAtTheEnd(first,last);
    else if(nbCompo==0)
      {
        _info_on_compo.resize(1);
        _mem.insertAtTheEnd(first,last);
      }
    else
      throw INTERP_KERNEL::Exception("DataArrayChar::insertAtTheEnd : not available for DataArrayChar with number of components different than 1 !");
  }
}

#endif
