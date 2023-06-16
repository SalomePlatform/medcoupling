// Copyright (C) 2007-2023  CEA/DEN, EDF R&D
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

#pragma once

#include "MEDCoupling.hxx"
#include "MCType.hxx"
#include "MCAuto.hxx"
#include "MEDCouplingTimeLabel.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "InterpKernelException.hxx"
#include "MEDCouplingTraits.hxx"
#include "MEDCouplingMap.hxx"
#include "BBTreePts.txx"

#include <string>
#include <vector>
#include <iterator>
#include <functional>

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
    void repr(mcIdType sl, std::ostream& stream) const;
    bool reprHeader(mcIdType sl, std::ostream& stream) const;
    void reprZip(mcIdType sl, std::ostream& stream) const;
    void reprNotTooLong(mcIdType sl, std::ostream& stream) const;
    void fillWithValue(const T& val);
    T *fromNoInterlace(std::size_t nbOfComp) const;
    T *toNoInterlace(std::size_t nbOfComp) const;
    void sort(bool asc);
    void reverse(std::size_t nbOfComp);
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
    static void COffsetDeallocator(void *pt, void *param);
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

  template <class T> class DataArrayTools
  {
  public:
    static void GetSlice(T start, T stop, T step, mcIdType sliceId, mcIdType nbOfSlices, T& startSlice, T& stopSlice);
    static mcIdType GetNumberOfItemGivenBES(T begin, T end, T step, const std::string& msg);
    static mcIdType GetNumberOfItemGivenBESRelative(T begin, T end, T step, const std::string& msg);
    static mcIdType GetPosOfItemGivenBESRelativeNoThrow(T value, T begin, T end, T step);
  };

  class DataArray;
  class DataArrayByte;

  MEDCOUPLING_EXPORT void DACheckNbOfTuplesAndComp(const DataArray *da, mcIdType nbOfTuples, std::size_t nbOfCompo, const std::string& msg);

  class MEDCOUPLING_EXPORT DataArray : public RefCountObject, public TimeLabel
  {
  public:
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    void setName(const std::string& name);
    void copyStringInfoFrom(const DataArray& other);
    void copyPartOfStringInfoFrom(const DataArray& other, const std::vector<std::size_t>& compoIds);
    void copyPartOfStringInfoFrom2(const std::vector<std::size_t>& compoIds, const DataArray& other);
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
    std::string getInfoOnComponent(std::size_t i) const;
    std::string getVarOnComponent(std::size_t i) const;
    std::string getUnitOnComponent(std::size_t i) const;
    void setInfoOnComponent(std::size_t i, const std::string& info);
    std::size_t getNumberOfComponents() const { return _info_on_compo.size(); }
    void setPartOfValuesBase3(const DataArray *aBase, const mcIdType *bgTuples, const mcIdType *endTuples, mcIdType bgComp, mcIdType endComp, mcIdType stepComp, bool strictCompoCompare=true);
    virtual void *getVoidStarPointer() = 0;
    virtual DataArray *deepCopy() const = 0;
    virtual DataArray *copySorted(bool asc=true) const = 0;
    virtual DataArray *buildNewEmptyInstance() const = 0;
    virtual bool isAllocated() const = 0;
    virtual void checkAllocated() const = 0;
    virtual void desallocate() = 0;
    virtual mcIdType getNumberOfTuples() const = 0;
    virtual mcIdType getNbOfElems() const = 0;
    virtual std::size_t getNbOfElemAllocated() const = 0;
    virtual void alloc(std::size_t nbOfTuple, std::size_t nbOfCompo=1) = 0;
    virtual void reAlloc(std::size_t newNbOfTuple) = 0;
    virtual void renumberInPlace(const mcIdType *old2New) = 0;
    virtual void renumberInPlaceR(const mcIdType *new2Old) = 0;
    virtual void setContigPartOfSelectedValues(mcIdType tupleIdStart, const DataArray *aBase, const DataArrayIdType *tuplesSelec) = 0;
    virtual void setContigPartOfSelectedValuesSlice(mcIdType tupleIdStart, const DataArray *aBase, mcIdType bg, mcIdType end2, mcIdType step) = 0;
    virtual DataArray *selectByTupleRanges(const std::vector<std::pair<mcIdType,mcIdType> >& ranges) const = 0;
    virtual DataArray *keepSelectedComponents(const std::vector<std::size_t>& compoIds) const = 0;
    virtual DataArray *selectByTupleId(const mcIdType *new2OldBg, const mcIdType *new2OldEnd) const = 0;
    virtual DataArray *selectByTupleIdSafe(const mcIdType *new2OldBg, const mcIdType *new2OldEnd) const = 0;
    virtual DataArray *selectByTupleIdSafeSlice(mcIdType bg, mcIdType end2, mcIdType step) const = 0;
    virtual void rearrange(std::size_t newNbOfCompo) = 0;
    virtual void circularPermutation(mcIdType nbOfShift=1) = 0;
    virtual void circularPermutationPerTuple(mcIdType nbOfShift=1) = 0;
    virtual void reversePerTuple() = 0;
    void checkNbOfTuples(mcIdType nbOfTuples, const std::string& msg) const;
    void checkNbOfComps(std::size_t nbOfCompo, const std::string& msg) const;
    void checkNbOfTuplesAndComp(const DataArray& other, const std::string& msg) const;
    void checkNbOfTuplesAndComp(mcIdType nbOfTuples, std::size_t nbOfCompo, const std::string& msg) const;
    void checkNbOfElems(mcIdType nbOfElems, const std::string& msg) const;
    static void GetSlice(mcIdType start, mcIdType stop, mcIdType step, mcIdType sliceId, mcIdType nbOfSlices, mcIdType& startSlice, mcIdType& stopSlice);
    static mcIdType GetNumberOfItemGivenBES(mcIdType begin, mcIdType end, mcIdType step, const std::string& msg);
    static mcIdType GetNumberOfItemGivenBESRelative(mcIdType begin, mcIdType end, mcIdType step, const std::string& msg);
    static mcIdType GetPosOfItemGivenBESRelativeNoThrow(mcIdType value, mcIdType begin, mcIdType end, mcIdType step);
    static std::vector<std::string> SplitStringInChuncks(const std::string st, std::size_t sz);
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
    static void CheckValueInRange(mcIdType ref, mcIdType value, const std::string& msg);
    static void CheckValueInRangeEx(mcIdType value, mcIdType start, mcIdType end, const std::string& msg);
    static void CheckClosingParInRange(mcIdType ref, mcIdType value, const std::string& msg);
    static mcIdType EffectiveCircPerm(mcIdType nbOfShift, mcIdType nbOfTuples);
  protected:
    std::string _name;
    std::vector<std::string> _info_on_compo;
  };
}

namespace MEDCoupling
{
  template<class T>
  class DataArrayTemplate : public DataArray
  {
  public:
    typedef T Type;
  public:
    static MCAuto< typename Traits<T>::ArrayTypeCh > NewFromStdVector(const typename std::vector<T>& v);
    static MCAuto< typename Traits<T>::ArrayTypeCh > NewFromArray(const T *arrBegin, const T *arrEnd);
    std::vector< MCAuto< typename Traits<T>::ArrayTypeCh > > explodeComponents() const;
    //
    void printForDebug(std::ostream& oss) const
    {
      this->checkAllocated();
      char comma[3] = {'\0',' ','\0'};
      std::for_each(this->begin(),this->end(),[&comma,&oss](const T& elt) { oss << comma << elt; comma[0]=','; } );
      oss << std::endl;
    }
    std::size_t getHeapMemorySizeWithoutChildren() const;
    void updateTime() const { }
    //
    mcIdType getNumberOfTuples() const { return ToIdType(_info_on_compo.empty()?0:_mem.getNbOfElem()/getNumberOfComponents()); }
    mcIdType getNbOfElems() const { return ToIdType(_mem.getNbOfElem()); }
    bool empty() const;
    void *getVoidStarPointer() { return getPointer(); }
    const T *getConstPointer() const { return _mem.getConstPointer(); }
    const T *begin() const { return getConstPointer(); }
    const T *end() const { return getConstPointer()+getNbOfElems(); }
    T *rwBegin() { return getPointer(); }
    T *rwEnd() { return getPointer()+getNbOfElems(); }
    void alloc(std::size_t nbOfTuple, std::size_t nbOfCompo=1);
    void useArray(const T *array, bool ownership, DeallocType type, std::size_t nbOfTuple, std::size_t nbOfCompo);
    void useExternalArrayWithRWAccess(const T *array, std::size_t nbOfTuple, std::size_t nbOfCompo);
    T getIJSafe(std::size_t tupleId, std::size_t compoId) const;
    T getIJ(std::size_t tupleId, std::size_t compoId) const { return _mem[tupleId*_info_on_compo.size()+compoId]; }
    void setIJ(std::size_t tupleId, std::size_t compoId, T newVal) { _mem[tupleId*_info_on_compo.size()+compoId]=newVal; declareAsNew(); }
    void setIJSilent(std::size_t tupleId, std::size_t compoId, T newVal) { _mem[tupleId*_info_on_compo.size()+compoId]=newVal; }
    T *getPointer() { declareAsNew(); return getPointerSilent(); }
    T *getPointerSilent() { return _mem.getPointer(); }
    void pack() const;
    bool isAllocated() const override;
    void checkAllocated() const;
    void desallocate();
    void reserve(std::size_t nbOfElems);
    void rearrange(std::size_t newNbOfCompo);
    void transpose();
    void pushBackSilent(T val);
    void pushBackValsSilent(const T *valsBg, const T *valsEnd);
    T popBackSilent();
    T front() const;
    T back() const;
    std::size_t getNbOfElemAllocated() const { return _mem.getNbOfElemAllocated(); }
    void allocIfNecessary(std::size_t nbOfTuple, std::size_t nbOfCompo);
    void deepCopyFrom(const DataArrayTemplate<T>& other);
    void reverse();
    void fillWithValue(T val);
    void reAlloc(std::size_t newNbOfTuple);
    void renumberInPlace(const mcIdType *old2New);
    void renumberInPlaceR(const mcIdType *new2Old);
    void sort(bool asc=true);
    typename Traits<T>::ArrayType *renumber(const mcIdType *old2New) const;
    typename Traits<T>::ArrayType *renumberR(const mcIdType *new2Old) const;
    typename Traits<T>::ArrayType *renumberAndReduce(const mcIdType *old2New, mcIdType newNbOfTuple) const;
    typename Traits<T>::ArrayType *changeNbOfComponents(std::size_t newNbOfComp, T dftValue) const;
    typename Traits<T>::ArrayType *subArray(mcIdType tupleIdBg, mcIdType tupleIdEnd=-1) const;
    MCAuto<typename Traits<T>::ArrayTypeCh> selectPartDef(const PartDefinition* pd) const;
    void circularPermutation(mcIdType nbOfShift=1);
    void circularPermutationPerTuple(mcIdType nbOfShift=1);
    void reversePerTuple();
    void setPartOfValues1(const typename Traits<T>::ArrayType *a, mcIdType bgTuples, mcIdType endTuples, mcIdType stepTuples, mcIdType bgComp, mcIdType endComp, mcIdType stepComp, bool strictCompoCompare=true);
    void setPartOfValuesSimple1(T a, mcIdType bgTuples, mcIdType endTuples, mcIdType stepTuples, mcIdType bgComp, mcIdType endComp, mcIdType stepComp);
    void setPartOfValues2(const typename Traits<T>::ArrayType *a, const mcIdType *bgTuples, const mcIdType *endTuples, const mcIdType *bgComp, const mcIdType *endComp, bool strictCompoCompare=true);
    void setPartOfValuesSimple2(T a, const mcIdType *bgTuples, const mcIdType *endTuples, const mcIdType *bgComp, const mcIdType *endComp);
    void setPartOfValues3(const typename Traits<T>::ArrayType *a, const mcIdType *bgTuples, const mcIdType *endTuples, mcIdType bgComp, mcIdType endComp, mcIdType stepComp, bool strictCompoCompare=true);
    void setPartOfValuesSimple3(T a, const mcIdType *bgTuples, const mcIdType *endTuples, mcIdType bgComp, mcIdType endComp, mcIdType stepComp);
    void setPartOfValues4(const typename Traits<T>::ArrayType *a, mcIdType bgTuples, mcIdType endTuples, mcIdType stepTuples, const mcIdType *bgComp, const mcIdType *endComp, bool strictCompoCompare=true);
    void setPartOfValuesSimple4(T a, mcIdType bgTuples, mcIdType endTuples, mcIdType stepTuples, const mcIdType *bgComp, const mcIdType *endComp);
    void setPartOfValuesAdv(const typename Traits<T>::ArrayType *a, const DataArrayIdType *tuplesSelec);
    void setContigPartOfSelectedValues(mcIdType tupleIdStart, const DataArray *aBase, const DataArrayIdType *tuplesSelec);
    void setContigPartOfSelectedValuesSlice(mcIdType tupleIdStart, const DataArray *aBase, mcIdType bg, mcIdType end2, mcIdType step);
    T getMaxValue(mcIdType& tupleId) const;
    T getMaxValueInArray() const;
    T getMaxAbsValue(std::size_t& tupleId) const;
    T getMaxAbsValueInArray() const;
    T getMinValue(mcIdType& tupleId) const;
    T getMinValueInArray() const;
    void getTuple(mcIdType tupleId, T *res) const { std::copy(_mem.getConstPointerLoc(tupleId*_info_on_compo.size()),_mem.getConstPointerLoc((tupleId+1)*_info_on_compo.size()),res); }
    template<class InputIterator>
    void insertAtTheEnd(InputIterator first, InputIterator last);
    static void SetArrayIn(typename Traits<T>::ArrayType *newArray, typename Traits<T>::ArrayType* &arrayToSet);
    void writeOnPlace(std::size_t id, T element0, const T *others, mcIdType sizeOfOthers) { _mem.writeOnPlace(id,element0,others,sizeOfOthers); }
    void fillWithZero();
  public:
    MemArray<T>& accessToMemArray() { return _mem; }
    const MemArray<T>& accessToMemArray() const { return _mem; }
  protected:
    typename Traits<T>::ArrayTypeCh *copySortedImpl(bool asc) const;
    typename Traits<T>::ArrayType *mySelectByTupleId(const mcIdType *new2OldBg, const mcIdType *new2OldEnd) const;
    typename Traits<T>::ArrayType *mySelectByTupleId(const DataArrayIdType& di) const;
    typename Traits<T>::ArrayType *mySelectByTupleIdSafe(const mcIdType *new2OldBg, const mcIdType *new2OldEnd) const;
    typename Traits<T>::ArrayType *myKeepSelectedComponents(const std::vector<std::size_t>& compoIds) const;
    typename Traits<T>::ArrayType *mySelectByTupleIdSafeSlice(mcIdType bg, mcIdType end2, mcIdType step) const;
    typename Traits<T>::ArrayType *mySelectByTupleRanges(const std::vector<std::pair<mcIdType,mcIdType> >& ranges) const;
  protected:
    MemArray<T> _mem;
  };

  template<class T>
  class DataArrayTemplateClassic : public DataArrayTemplate<T>
  {
  public:
    MCAuto<DataArrayDouble> convertToDblArr() const;
    MCAuto<DataArrayInt> convertToIntArr() const;
    MCAuto<DataArrayInt64> convertToInt64Arr() const;
    MCAuto<DataArrayFloat> convertToFloatArr() const;
    void applyLin(T a, T b, std::size_t compoId);
    void applyLin(T a, T b);
    typename Traits<T>::ArrayType *negate() const;
    void addEqual(const typename Traits<T>::ArrayType *other);
    void substractEqual(const typename Traits<T>::ArrayType *other);
    void multiplyEqual(const typename Traits<T>::ArrayType *other);
    void divideEqual(const typename Traits<T>::ArrayType *other);
    static typename Traits<T>::ArrayType *Substract(const typename Traits<T>::ArrayType *a1, const typename Traits<T>::ArrayType *a2);
    static typename Traits<T>::ArrayType *Divide(const typename Traits<T>::ArrayType *a1, const typename Traits<T>::ArrayType *a2);
    static typename Traits<T>::ArrayType *Add(const typename Traits<T>::ArrayType *a1, const typename Traits<T>::ArrayType *a2);
    static typename Traits<T>::ArrayType *Multiply(const typename Traits<T>::ArrayType *a1, const typename Traits<T>::ArrayType *a2);
    static typename Traits<T>::ArrayType *Meld(const typename Traits<T>::ArrayType *a1, const typename Traits<T>::ArrayType *a2);
    static typename Traits<T>::ArrayType *Meld(const std::vector<const typename Traits<T>::ArrayType *>& arr);
    MCAuto<DataArrayIdType> findIdsGreaterOrEqualTo(T val) const;
    MCAuto<DataArrayIdType> findIdsGreaterThan(T val) const;
    MCAuto<DataArrayIdType> findIdsLowerOrEqualTo(T val) const;
    MCAuto<DataArrayIdType> findIdsLowerThan(T val) const;
    DataArrayIdType *findIdsStrictlyNegative() const;
    typename Traits<T>::ArrayType *fromNoInterlace() const;
    typename Traits<T>::ArrayType *toNoInterlace() const;
    void meldWith(const typename Traits<T>::ArrayType *other);
    typename Traits<T>::ArrayType *duplicateEachTupleNTimes(mcIdType nbTimes) const;
    void aggregate(const typename Traits<T>::ArrayType *other);
    void abs();
    typename Traits<T>::ArrayType *computeAbs() const;
    typename Traits<T>::ArrayType *performCopyOrIncrRef(bool dCpy) const;
    typename Traits<T>::ArrayType *sumPerTuple() const;
    void iota(T init=(T)0);
    void reprStream(std::ostream& stream) const;
    void reprZipStream(std::ostream& stream) const;
    void reprNotTooLongStream(std::ostream& stream) const;
    void reprWithoutNameStream(std::ostream& stream) const;
    void reprZipWithoutNameStream(std::ostream& stream) const;
    void reprNotTooLongWithoutNameStream(std::ostream& stream) const;
    std::string repr() const;
    std::string reprZip() const;
    std::string reprNotTooLong() const;
    template<class U>
    MCAuto< typename Traits<U>::ArrayType > convertToOtherTypeOfArr() const;
  protected:
    static typename Traits<T>::ArrayType *PerformCopyOrIncrRef(bool dCpy, const typename Traits<T>::ArrayType& self);
    template<class OP>
    MCAuto<DataArrayIdType> findIdsAdv(const OP& op) const;
  private:
    template<class FCT>
    void somethingEqual(const typename Traits<T>::ArrayType *other);
  };

  template<class T>
  class DataArrayTemplateFP : public DataArrayTemplateClassic<T>
  {
  public:
    bool isUniform(T val, T eps) const;
  };
}

namespace MEDCoupling
{
  class DataArrayFloatIterator;
  class MEDCOUPLING_EXPORT DataArrayFloat : public DataArrayTemplateFP<float>
  {
  public:
    static DataArrayFloat *New();
  public:// abstract method overload
    DataArrayFloat *deepCopy() const;
    DataArrayFloat *copySorted(bool asc=true) const override { return this->copySortedImpl(asc); }
    std::string getClassName() const override { return std::string("DataArrayFloat"); }
    DataArrayFloat *buildNewEmptyInstance() const { return DataArrayFloat::New(); }
    DataArrayFloat *selectByTupleRanges(const std::vector<std::pair<mcIdType,mcIdType> >& ranges) const { return DataArrayTemplateFP<float>::mySelectByTupleRanges(ranges); }
    DataArrayFloat *keepSelectedComponents(const std::vector<std::size_t>& compoIds) const { return DataArrayTemplateFP<float>::myKeepSelectedComponents(compoIds); }
    DataArrayFloat *selectByTupleId(const mcIdType *new2OldBg, const mcIdType *new2OldEnd) const { return this->mySelectByTupleId(new2OldBg,new2OldEnd); }
    DataArrayFloat *selectByTupleIdSafe(const mcIdType *new2OldBg, const mcIdType *new2OldEnd) const { return DataArrayTemplateFP<float>::mySelectByTupleIdSafe(new2OldBg,new2OldEnd); }
    DataArrayFloat *selectByTupleIdSafeSlice(mcIdType bg, mcIdType end2, mcIdType step) const { return DataArrayTemplateFP<float>::mySelectByTupleIdSafeSlice(bg,end2,step); }
    void reprCppStream(const std::string& varName, std::ostream& stream) const;
    void reprQuickOverview(std::ostream& stream) const;
    void reprQuickOverviewData(std::ostream& stream, std::size_t maxNbOfByteInRepr) const;
  public:// non abstract but essential
    bool isEqual(const DataArrayFloat& other, float prec) const;
    bool isEqualIfNotWhy(const DataArrayFloat& other, float prec, std::string& reason) const;
    bool isEqualWithoutConsideringStr(const DataArrayFloat& other, float prec) const;
  public:
    DataArrayFloatIterator *iterator();
  private:
    ~DataArrayFloat() { }
    DataArrayFloat() { }
  };
}

namespace MEDCoupling
{
  class DataArrayDoubleIterator;
  class MEDCOUPLING_EXPORT DataArrayDouble : public DataArrayTemplateFP<double>
  {
  public:
    static DataArrayDouble *New();
    double doubleValue() const;
    DataArrayDouble *deepCopy() const;
    DataArrayDouble *copySorted(bool asc=true) const override { return this->copySortedImpl(asc); }
    std::string getClassName() const override { return std::string("DataArrayDouble"); }
    DataArrayDouble *buildNewEmptyInstance() const { return DataArrayDouble::New(); }
    void checkMonotonic(bool increasing, double eps) const;
    bool isMonotonic(bool increasing, double eps) const;
    void writeVTK(std::ostream& ofs, mcIdType indent, const std::string& nameInFile, DataArrayByte *byteArr) const;
    void reprCppStream(const std::string& varName, std::ostream& stream) const;
    void reprQuickOverview(std::ostream& stream) const;
    void reprQuickOverviewData(std::ostream& stream, std::size_t maxNbOfByteInRepr) const;
    bool isEqual(const DataArrayDouble& other, double prec) const;
    bool isEqualIfNotWhy(const DataArrayDouble& other, double prec, std::string& reason) const;
    bool isEqualWithoutConsideringStr(const DataArrayDouble& other, double prec) const;
    DataArrayDouble *selectByTupleId(const mcIdType *new2OldBg, const mcIdType *new2OldEnd) const { return this->mySelectByTupleId(new2OldBg,new2OldEnd); }
    DataArrayDouble *selectByTupleId(const DataArrayIdType& di) const { return this->mySelectByTupleId(di); }
    DataArrayDouble *selectByTupleIdSafe(const mcIdType *new2OldBg, const mcIdType *new2OldEnd) const { return DataArrayTemplateFP<double>::mySelectByTupleIdSafe(new2OldBg,new2OldEnd); }
    DataArrayDouble *keepSelectedComponents(const std::vector<std::size_t>& compoIds) const { return DataArrayTemplateFP<double>::myKeepSelectedComponents(compoIds); }
    DataArrayDouble *selectByTupleIdSafeSlice(mcIdType bg, mcIdType end2, mcIdType step) const { return DataArrayTemplateFP<double>::mySelectByTupleIdSafeSlice(bg,end2,step); }
    DataArrayDouble *selectByTupleRanges(const std::vector<std::pair<mcIdType,mcIdType> >& ranges) const { return DataArrayTemplateFP<double>::mySelectByTupleRanges(ranges); }
    bool areIncludedInMe(const DataArrayDouble *other, double prec, DataArrayIdType *&tupleIds) const;
    void findCommonTuples(double prec, mcIdType limitTupleId, DataArrayIdType *&comm, DataArrayIdType *&commIndex) const;
    double minimalDistanceTo(const DataArrayDouble *other, mcIdType& thisTupleId, mcIdType& otherTupleId) const;
    DataArrayDouble *getDifferentValues(double prec, mcIdType limitTupleId=-1) const;
    DataArrayIdType *findClosestTupleId(const DataArrayDouble *other) const;
    DataArrayIdType *computeNbOfInteractionsWith(const DataArrayDouble *otherBBoxFrmt, double eps) const;
    void setSelectedComponents(const DataArrayDouble *a, const std::vector<std::size_t>& compoIds);
    DataArrayDoubleIterator *iterator();
    void checkNoNullValues() const;
    void getMinMaxPerComponent(double *bounds) const;
    DataArrayDouble *computeBBoxPerTuple(double epsilon=0.0) const;
    void computeTupleIdsNearTuples(const DataArrayDouble *other, double eps, DataArrayIdType *& c, DataArrayIdType *& cI) const;
    void recenterForMaxPrecision(double eps);
    double getMaxValue2(DataArrayIdType*& tupleIds) const;
    double getMinValue2(DataArrayIdType*& tupleIds) const;
    mcIdType count(double value, double eps) const;
    double getAverageValue() const;
    double norm2() const;
    double normMax() const;
    void normMaxPerComponent(double * res) const;
    double normMin() const;
    void accumulate(double *res) const;
    double accumulate(std::size_t compId) const;
    DataArrayDouble *accumulatePerChunck(const mcIdType *bgOfIndex, const mcIdType *endOfIndex) const;
    MCAuto<DataArrayDouble> cumSum() const;
    double distanceToTuple(const double *tupleBg, const double *tupleEnd, mcIdType& tupleId) const;
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
    DataArrayDouble *minPerTuple() const;
    DataArrayDouble *maxPerTuple() const;
    DataArrayDouble *maxPerTupleWithCompoId(DataArrayIdType* &compoIdOfMaxPerTuple) const;
    DataArrayDouble *buildEuclidianDistanceDenseMatrix() const;
    DataArrayDouble *buildEuclidianDistanceDenseMatrixWith(const DataArrayDouble *other) const;
    void asArcOfCircle(double center[2], double& radius, double& ang) const;
    void sortPerTuple(bool asc);
    void applyInv(double numerator);
    void applyPow(double val);
    void applyRPow(double val);
    DataArrayDouble *applyFunc(std::size_t nbOfComp, FunctionToEvaluate func) const;
    DataArrayDouble *applyFunc(std::size_t nbOfComp, const std::string& func, bool isSafe=true) const;
    DataArrayDouble *applyFunc(const std::string& func, bool isSafe=true) const;
    void applyFuncOnThis(const std::string& func, bool isSafe=true);
    DataArrayDouble *applyFuncCompo(std::size_t nbOfComp, const std::string& func, bool isSafe=true) const;
    DataArrayDouble *applyFuncNamedCompo(std::size_t nbOfComp, const std::vector<std::string>& varsOrder, const std::string& func, bool isSafe=true) const;
    void applyFuncFast32(const std::string& func);
    void applyFuncFast64(const std::string& func);
    MCAuto<DataArrayDouble> symmetry3DPlane(const double point[3], const double normalVector[3]) const;
    DataArrayIdType *findIdsInRange(double vmin, double vmax) const;
    DataArrayIdType *findIdsNotInRange(double vmin, double vmax) const;
    static DataArrayDouble *Aggregate(const DataArrayDouble *a1, const DataArrayDouble *a2);
    static DataArrayDouble *Aggregate(const std::vector<const DataArrayDouble *>& arr);
    static DataArrayDouble *Dot(const DataArrayDouble *a1, const DataArrayDouble *a2);
    static DataArrayDouble *CrossProduct(const DataArrayDouble *a1, const DataArrayDouble *a2);
    static DataArrayDouble *Max(const DataArrayDouble *a1, const DataArrayDouble *a2);
    static DataArrayDouble *Min(const DataArrayDouble *a1, const DataArrayDouble *a2);
    static DataArrayDouble *Pow(const DataArrayDouble *a1, const DataArrayDouble *a2);
    void powEqual(const DataArrayDouble *other);
    std::vector<bool> toVectorOfBool(double eps) const;
    static void Rotate2DAlg(const double *center, double angle, mcIdType nbNodes, const double *coordsIn, double *coordsOut);
    static void Rotate3DAlg(const double *center, const double *vect, double angle, mcIdType nbNodes, const double *coordsIn, double *coordsOut);
    static void Symmetry3DPlane(const double point[3], const double normalVector[3], mcIdType nbNodes, const double *coordsIn, double *coordsOut);
    static void GiveBaseForPlane(const double normalVector[3], double baseOfPlane[9]);
    static void ComputeIntegralOfSeg2IntoTri3(const double seg2[4], const double tri3[6], double coeffs[3], double& length);
  public:
    void getTinySerializationIntInformation(std::vector<mcIdType>& tinyInfo) const;
    void getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const;
    bool resizeForUnserialization(const std::vector<mcIdType>& tinyInfoI);
    void finishUnserialization(const std::vector<mcIdType>& tinyInfoI, const std::vector<std::string>& tinyInfoS);
  public:
    template<int SPACEDIM>
    void findCommonTuplesAlg(const double *bbox, mcIdType nbNodes, mcIdType limitNodeId, double prec, DataArrayIdType *c, DataArrayIdType *cI) const;
    template<int SPACEDIM>
    static void FindClosestTupleIdAlg(const BBTreePts<SPACEDIM,mcIdType>& myTree, double dist, const double *pos, mcIdType nbOfTuples, const double *thisPt, mcIdType thisNbOfTuples, mcIdType *res);
    template<int SPACEDIM>
    static void FindTupleIdsNearTuplesAlg(const BBTreePts<SPACEDIM,mcIdType>& myTree, const double *pos, mcIdType nbOfTuples, double eps,
                                          DataArrayIdType *c, DataArrayIdType *cI);
  private:
    DataArrayDouble *operatePerTuple(std::function<double(const double *bg, const double *endd)> func) const;
  private:
    ~DataArrayDouble() { }
    DataArrayDouble() { }
  };
}

namespace MEDCoupling
{
  template<class T>
  class DataArrayDiscrete : public DataArrayTemplateClassic<T>
  {
  public:
    using DataArrayType = typename Traits<T>::ArrayType;
  public:
    static DataArrayType *New();
    T intValue() const;
    bool isEqual(const DataArrayDiscrete<T>& other) const;
    bool isEqualIfNotWhy(const DataArrayDiscrete<T>& other, std::string& reason) const;
    bool isEqualWithoutConsideringStr(const DataArrayDiscrete<T>& other) const;
    bool isEqualWithoutConsideringStrAndOrder(const typename Traits<T>::ArrayType& other) const;
    void switchOnTupleEqualTo(T val, std::vector<bool>& vec) const;
    void switchOnTupleNotEqualTo(T val, std::vector<bool>& vec) const;
    DataArrayIdType *occurenceRankInThis() const;
    DataArrayIdType *buildPermutationArr(const DataArrayDiscrete<T>& other) const;
    DataArrayIdType *indicesOfSubPart(const DataArrayDiscrete<T>& partOfThis) const;
    void checkMonotonic(bool increasing) const;
    bool isMonotonic(bool increasing) const;
    void checkStrictlyMonotonic(bool increasing) const;
    bool isStrictlyMonotonic(bool increasing) const;
    mcIdType getHashCode() const;
    void reprCppStream(const std::string& varName, std::ostream& stream) const;
    void reprQuickOverview(std::ostream& stream) const;
    void reprQuickOverviewData(std::ostream& stream, std::size_t maxNbOfByteInRepr) const;
    void writeVTK(std::ostream& ofs, mcIdType indent, const std::string& type, const std::string& nameInFile, DataArrayByte *byteArr) const;
    void transformWithIndArr(const T *indArrBg, const T *indArrEnd);
    void transformWithIndArr(const MapKeyVal<T, T>& m);
    DataArrayIdType *findIdsEqual(T val) const;
    DataArrayIdType *transformWithIndArrR(const T *indArr2Bg, const T *indArrEnd) const;
    void splitByValueRange(const T *arrBg, const T *arrEnd,
                           DataArrayType *& castArr, DataArrayType *& rankInsideCast, DataArrayType *& castsPresent) const;
    bool isRange(T& strt, T& sttoopp, T& stteepp) const;
    DataArrayIdType *invertArrayO2N2N2O(mcIdType newNbOfElem) const;
    DataArrayIdType *invertArrayN2O2O2N(mcIdType oldNbOfElem) const;
    DataArrayIdType *invertArrayO2N2N2OBis(mcIdType newNbOfElem) const;
    MCAuto< MapKeyVal<T, mcIdType> > invertArrayN2O2O2NOptimized() const;
    MCAuto< MapKeyVal<mcIdType, T> > giveN2OOptimized() const;
    MCAuto<DataArrayIdType> findIdForEach(const T *valsBg, const T *valsEnd) const;
    DataArrayIdType *checkAndPreparePermutation() const;
    void changeSurjectiveFormat(T targetNb, DataArrayIdType *&arr, DataArrayIdType *&arrI) const;
    DataArrayIdType *buildPermArrPerLevel() const;
    bool isIota(mcIdType sizeExpected) const;
    bool isUniform(T val) const;
    T checkUniformAndGuess() const;
    bool hasUniqueValues() const;
    void setSelectedComponents(const DataArrayType *a, const std::vector<std::size_t>& compoIds);
    DataArrayIdType *locateComponentId(const DataArrayType *valToSearchIntoTuples, const DataArrayIdType *tupleIdHint) const;
    DataArrayIdType *findIdsNotEqual(T val) const;
    DataArrayIdType *findIdsEqualTuple(const T *tupleBg, const T *tupleEnd) const;
    DataArrayIdType *findIdsEqualList(const T *valsBg, const T *valsEnd) const;
    DataArrayIdType *findIdsNotEqualList(const T *valsBg, const T *valsEnd) const;
    mcIdType findIdFirstEqual(T value) const;
    mcIdType findIdFirstEqual(const std::vector<T>& vals) const;
    mcIdType findIdFirstEqualTuple(const std::vector<T>& tupl) const;
    mcIdType findIdSequence(const std::vector<T>& vals) const;
    mcIdType changeValue(T oldValue, T newValue);
    mcIdType count(T value) const;
    bool presenceOfTuple(const std::vector<T>& tupl) const;
    bool presenceOfValue(T value) const;
    bool presenceOfValue(const std::vector<T>& vals) const;
    void accumulate(T *res) const;
    T accumulate(std::size_t compId) const;
    DataArrayType *accumulatePerChunck(const mcIdType *bgOfIndex, const mcIdType *endOfIndex) const;
    void getMinMaxValues(T& minValue, T& maxValue) const;
    void applyInv(T numerator);
    void applyDivideBy(T val);
    void applyModulus(T val);
    void applyRModulus(T val);
    void applyPow(T val);
    void applyRPow(T val);
    DataArrayIdType *findIdsInRange(T vmin, T vmax) const;
    DataArrayIdType *findIdsNotInRange(T vmin, T vmax) const;
    bool checkAllIdsInRange(T vmin, T vmax) const;
    static DataArrayType *Aggregate(const DataArrayType *a1, const DataArrayType *a2, T offsetA2);
    static DataArrayType *Aggregate(const std::vector<const DataArrayType *>& arr);
    static DataArrayType *AggregateIndexes(const std::vector<const DataArrayType *>& arrs);
    static DataArrayType *BuildUnion(const std::vector<const DataArrayType *>& arr);
    static DataArrayType *BuildIntersection(const std::vector<const DataArrayType *>& arr);
    static void PutIntoToSkylineFrmt(const std::vector< std::vector<T> >& v, DataArrayType *& data, DataArrayIdType *& dataIndex);
    DataArrayIdType *buildComplement(mcIdType nbOfElement) const;
    DataArrayType *buildSubstraction(const DataArrayType *other) const;
    DataArrayType *buildSubstractionOptimized(const DataArrayType *other) const;
    DataArrayType *buildUnion(const DataArrayType *other) const;
    DataArrayType *buildIntersection(const DataArrayType *other) const;
    DataArrayIdType *indexOfSameConsecutiveValueGroups() const;
    DataArrayType *buildUnique() const;
    DataArrayType *buildUniqueNotSorted() const;
    DataArrayType *deltaShiftIndex() const;
    void computeOffsets();
    void computeOffsetsFull();
    void findIdsRangesInListOfIds(const DataArrayType *listOfIds, DataArrayIdType *& rangeIdsFetched, DataArrayType *& idsInInputListThatFetch) const;
    DataArrayType *buildExplicitArrByRanges(const DataArrayType *offsets) const;
    DataArrayType *buildExplicitArrOfSliceOnScaledArr(T begin, T stop, T step) const;
    DataArrayIdType *findRangeIdForEachTuple(const DataArrayType *ranges) const;
    DataArrayType *findIdInRangeForEachTuple(const DataArrayType *ranges) const;
    void sortEachPairToMakeALinkedList();
    void sortToHaveConsecutivePairs();
    MCAuto<DataArrayType> fromLinkedListOfPairToList() const;
    DataArrayType *getDifferentValues() const;
    std::vector<DataArrayIdType *> partitionByDifferentValues(std::vector<T>& differentIds) const;
    std::vector< std::pair<mcIdType,mcIdType> > splitInBalancedSlices(mcIdType nbOfSlices) const;
    static DataArrayType *Modulus(const DataArrayType *a1, const DataArrayType *a2);
    void modulusEqual(const DataArrayType *other);
    static DataArrayType *Pow(const DataArrayType *a1, const DataArrayType *a2);
    void powEqual(const DataArrayType *other);
    //MemArray<T>& accessToMemArray() { return _mem; }
    //const MemArray<T>& accessToMemArray() const { return _mem; }
  public:
    static DataArrayIdType *FindPermutationFromFirstToSecond(const DataArrayType *ids1, const DataArrayType *ids2);
    static DataArrayIdType *FindPermutationFromFirstToSecondDuplicate(const DataArrayType *ids1, const DataArrayType *ids2);
    static mcIdType *CheckAndPreparePermutation(const T *start, const T *end);
    static DataArrayType *BuildListOfSwitchedOn(const std::vector<bool>& v);
    static DataArrayType *BuildListOfSwitchedOff(const std::vector<bool>& v);
    static DataArrayIdType *ConvertIndexArrayToO2N(mcIdType nbOfOldTuples, const mcIdType *arr, const mcIdType *arrIBg, const mcIdType *arrIEnd, mcIdType &newNbOfTuples);
    static DataArrayIdType *MakePartition(const std::vector<const DataArrayType *>& groups, mcIdType newNb, std::vector< std::vector<mcIdType> >& fidsOfGroups);
  public:
    static void ExtractFromIndexedArrays(const mcIdType *idsOfSelectBg, const mcIdType *idsOfSelectEnd,
                                                            const DataArrayType *arrIn, const DataArrayIdType *arrIndxIn,
                                                            DataArrayType* &arrOut, DataArrayIdType* &arrIndexOut);
    static void ExtractFromIndexedArraysSlice(mcIdType idsOfSelectStart, mcIdType idsOfSelectStop, mcIdType idsOfSelectStep,
                                                                 const DataArrayType *arrIn, const DataArrayIdType *arrIndxIn,
                                                                 DataArrayType* &arrOut, DataArrayIdType* &arrIndexOut);
    static void SetPartOfIndexedArrays(const mcIdType *idsOfSelectBg, const mcIdType *idsOfSelectEnd,
                                                          const DataArrayType *arrIn, const DataArrayIdType *arrIndxIn,
                                                          const DataArrayType *srcArr, const DataArrayIdType *srcArrIndex,
                                                          DataArrayType* &arrOut, DataArrayIdType* &arrIndexOut);
    static void SetPartOfIndexedArraysSlice(mcIdType start, mcIdType end, mcIdType step,
                                                               const DataArrayType *arrIn, const DataArrayIdType *arrIndxIn,
                                                               const DataArrayType *srcArr, const DataArrayIdType *srcArrIndex,
                                                               DataArrayType* &arrOut, DataArrayIdType* &arrIndexOut);
    static void SetPartOfIndexedArraysSameIdx(const mcIdType *idsOfSelectBg, const mcIdType *idsOfSelectEnd,
                                                                 DataArrayType *arrInOut, const DataArrayIdType *arrIndxIn,
                                                                 const DataArrayType *srcArr, const DataArrayIdType *srcArrIndex);
    static void SetPartOfIndexedArraysSameIdxSlice(mcIdType start, mcIdType end, mcIdType step,
                                                                      DataArrayType *arrInOut, const DataArrayIdType *arrIndxIn,
                                                                      const DataArrayType *srcArr, const DataArrayIdType *srcArrIndex);
    static bool RemoveIdsFromIndexedArrays(const T *idsToRemoveBg, const T *idsToRemoveEnd,
                                                              DataArrayType *arr, DataArrayIdType *arrIndx, mcIdType offsetForRemoval=0);
    static DataArrayType *Range(T begin, T end, T step);
  public:
    void getTinySerializationIntInformation(std::vector<mcIdType>& tinyInfo) const;
    void getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const;
    bool resizeForUnserialization(const std::vector<mcIdType>& tinyInfoI);
    void finishUnserialization(const std::vector<mcIdType>& tinyInfoI, const std::vector<std::string>& tinyInfoS);
  protected:
    template<class ALG>
    void switchOnTupleAlg(T val, std::vector<bool>& vec, ALG algo) const;
  protected:
    ~DataArrayDiscrete() { }
  };

  template<class T>
  class DataArrayDiscreteSigned : public DataArrayDiscrete<T>
  {
  public:
    bool isFittingWith(const std::vector<bool>& v) const;
  protected:
    ~DataArrayDiscreteSigned() { }
  };

  class DataArrayInt32Iterator;

  class MEDCOUPLING_EXPORT DataArrayInt32 : public DataArrayDiscreteSigned<Int32>
  {
    friend class DataArrayDiscrete<Int32>;
  public:
    DataArrayInt32 *deepCopy() const;
    DataArrayInt32 *copySorted(bool asc=true) const override { return this->copySortedImpl(asc); }
    DataArrayInt32 *buildNewEmptyInstance() const { return DataArrayInt32::New(); }
    MCAuto<DataArrayInt64> convertToInt64Arr() const;
  public:
    DataArrayInt32 *selectByTupleId(const mcIdType *new2OldBg, const mcIdType *new2OldEnd) const { return this->mySelectByTupleId(new2OldBg,new2OldEnd); }
    DataArrayInt32 *selectByTupleId(const DataArrayIdType& di) const { return this->mySelectByTupleId(di); }
    DataArrayInt32 *selectByTupleIdSafe(const mcIdType *new2OldBg, const mcIdType *new2OldEnd) const { return this->mySelectByTupleIdSafe(new2OldBg,new2OldEnd); }
    DataArrayInt32 *keepSelectedComponents(const std::vector<std::size_t>& compoIds) const { return this->myKeepSelectedComponents(compoIds); }
    DataArrayInt32 *selectByTupleIdSafeSlice(mcIdType bg, mcIdType end2, mcIdType step) const { return this->mySelectByTupleIdSafeSlice(bg,end2,step); }
    DataArrayInt32 *selectByTupleRanges(const std::vector<std::pair<mcIdType,mcIdType> >& ranges) const { return this->mySelectByTupleRanges(ranges); }
    std::string getClassName() const override { return std::string("DataArrayInt32"); }
  public:
    DataArrayInt32Iterator *iterator();
  private:
    ~DataArrayInt32() { }
    DataArrayInt32() { }
  };

  class MEDCOUPLING_EXPORT DataArrayInt64 : public DataArrayDiscreteSigned<Int64>
  {
    friend class DataArrayDiscrete<Int64>;
  public:
    DataArrayInt64 *deepCopy() const;
    DataArrayInt64 *copySorted(bool asc=true) const override { return this->copySortedImpl(asc); }
    DataArrayInt64 *buildNewEmptyInstance() const { return DataArrayInt64::New(); }//ok
    MCAuto<DataArrayInt32> convertToInt32Arr() const;
  public:
    DataArrayInt64 *selectByTupleId(const mcIdType *new2OldBg, const mcIdType *new2OldEnd) const { return this->mySelectByTupleId(new2OldBg,new2OldEnd); }
    DataArrayInt64 *selectByTupleId(const DataArrayIdType& di) const { return this->mySelectByTupleId(di); }
    DataArrayInt64 *selectByTupleIdSafe(const mcIdType *new2OldBg, const mcIdType *new2OldEnd) const { return DataArrayTemplate<Int64>::mySelectByTupleIdSafe(new2OldBg,new2OldEnd); }
    DataArrayInt64 *keepSelectedComponents(const std::vector<std::size_t>& compoIds) const { return DataArrayTemplate<Int64>::myKeepSelectedComponents(compoIds); }
    DataArrayInt64 *selectByTupleIdSafeSlice(mcIdType bg, mcIdType end2, mcIdType step) const { return DataArrayTemplate<Int64>::mySelectByTupleIdSafeSlice(bg,end2,step); }
    DataArrayInt64 *selectByTupleRanges(const std::vector<std::pair<mcIdType,mcIdType> >& ranges) const { return DataArrayTemplate<Int64>::mySelectByTupleRanges(ranges); }
    std::string getClassName() const override { return std::string("DataArrayInt64"); }
  public:
    DataArrayInt64Iterator *iterator();
  private:
    ~DataArrayInt64() { }
    DataArrayInt64() { }
  };
}

namespace MEDCoupling
{

  template<class T>
  template<class OP>
  MCAuto<DataArrayIdType> DataArrayTemplateClassic<T>::findIdsAdv(const OP& op) const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::findIdsAdv : this must have exactly one component !");
    const T *cptr(this->begin());
    MCAuto<DataArrayIdType> ret(DataArrayIdType::New()); ret->alloc(0,1);
    mcIdType nbOfTuples=this->getNumberOfTuples();
    for(mcIdType i=0;i<nbOfTuples;i++,cptr++)
      if(op(*cptr))
        ret->pushBackSilent(i);
    return ret;
  }

  class MEDCOUPLING_EXPORT DataArrayChar : public DataArrayTemplate<char>
  {
  public:
    virtual DataArrayChar *buildEmptySpecializedDAChar() const = 0;
    mcIdType getHashCode() const;
    bool isEqual(const DataArrayChar& other) const;
    virtual bool isEqualIfNotWhy(const DataArrayChar& other, std::string& reason) const;
    bool isEqualWithoutConsideringStr(const DataArrayChar& other) const;
    std::string repr() const;
    std::string reprZip() const;
    DataArrayInt *convertToIntArr() const;
    DataArrayChar *selectByTupleId(const mcIdType *new2OldBg, const mcIdType *new2OldEnd) const { return this->mySelectByTupleId(new2OldBg,new2OldEnd); }
    DataArrayChar *selectByTupleId(const DataArrayIdType& di) const { return this->mySelectByTupleId(di); }
    DataArrayChar *selectByTupleIdSafe(const mcIdType *new2OldBg, const mcIdType *new2OldEnd) const { return DataArrayTemplate<char>::mySelectByTupleIdSafe(new2OldBg,new2OldEnd); }
    DataArrayChar *keepSelectedComponents(const std::vector<std::size_t>& compoIds) const { return DataArrayTemplate<char>::myKeepSelectedComponents(compoIds); }
    DataArrayChar *selectByTupleIdSafeSlice(mcIdType bg, mcIdType end2, mcIdType step) const { return DataArrayTemplate<char>::mySelectByTupleIdSafeSlice(bg,end2,step); }
    bool isUniform(char val) const;
    void meldWith(const DataArrayChar *other);
    DataArray *selectByTupleRanges(const std::vector<std::pair<mcIdType,mcIdType> >& ranges) const { return DataArrayTemplate<char>::mySelectByTupleRanges(ranges); }
    DataArrayIdType *findIdsEqual(char val) const;
    DataArrayIdType *findIdsNotEqual(char val) const;
    mcIdType findIdSequence(const std::vector<char>& vals) const;
    mcIdType findIdFirstEqualTuple(const std::vector<char>& tupl) const;
    mcIdType findIdFirstEqual(char value) const;
    mcIdType findIdFirstEqual(const std::vector<char>& vals) const;
    bool presenceOfTuple(const std::vector<char>& tupl) const;
    bool presenceOfValue(char value) const;
    bool presenceOfValue(const std::vector<char>& vals) const;
    DataArrayIdType *findIdsInRange(char vmin, char vmax) const;
    static DataArrayChar *Aggregate(const DataArrayChar *a1, const DataArrayChar *a2);
    static DataArrayChar *Aggregate(const std::vector<const DataArrayChar *>& arr);
    static DataArrayChar *Meld(const DataArrayChar *a1, const DataArrayChar *a2);
    static DataArrayChar *Meld(const std::vector<const DataArrayChar *>& arr);
    MemArray<char>& accessToMemArray() { return _mem; }
    const MemArray<char>& accessToMemArray() const { return _mem; }
  public:
    //void getTinySerializationIntInformation(std::vector<mcIdType>& tinyInfo) const;
    //void getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const;
    //bool resizeForUnserialization(const std::vector<mcIdType>& tinyInfoI);
    //void finishUnserialization(const std::vector<mcIdType>& tinyInfoI, const std::vector<std::string>& tinyInfoS);
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
    DataArrayByte *copySorted(bool asc=true) const override { return this->copySortedImpl(asc); }
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
    std::string getClassName() const override { return std::string("DataArrayByte"); }
  private:
    ~DataArrayByte() { }
    DataArrayByte() { }
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
    DataArrayAsciiChar *copySorted(bool asc=true) const override { (void)asc;throw INTERP_KERNEL::Exception("DataArrayAsciiChar::copySorted : not implemented for DataArrayByte"); }
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
    std::string getClassName() const override { return std::string("DataArrayAsciiChar"); }
  private:
    ~DataArrayAsciiChar() { }
    DataArrayAsciiChar() { }
    DataArrayAsciiChar(const std::string& st);
    DataArrayAsciiChar(const std::vector<std::string>& vst, char defaultChar);
  };

  template<class T>
  class DataArrayIterator
  {
  public:
    DataArrayIterator(typename Traits<T>::ArrayType *da);
    ~DataArrayIterator();
    typename Traits<T>::ArrayTuple *nextt();
  private:
    typename Traits<T>::ArrayType *_da;
    T *_pt;
    mcIdType _tuple_id;
    std::size_t _nb_comp;
    mcIdType _nb_tuple;
  };

  template<class T>
  class DataArrayTuple
  {
  public:
    DataArrayTuple(T *pt, std::size_t nbOfComp);
    std::string repr() const;
    std::size_t getNumberOfCompo() const { return _nb_of_compo; }
    const T *getConstPointer() const { return  _pt; }
    T *getPointer() { return _pt; }
    typename Traits<T>::ArrayType *buildDA(std::size_t nbOfTuples, std::size_t nbOfCompo) const;
  protected:
    T zeValue() const;
  protected:
    T *_pt;
    std::size_t _nb_of_compo;
  };

  class DataArrayDoubleTuple;

  class MEDCOUPLING_EXPORT DataArrayDoubleIterator : public DataArrayIterator<double>
  {
  public:
    DataArrayDoubleIterator(DataArrayDouble *da);
    ~DataArrayDoubleIterator() { }
  };

  class MEDCOUPLING_EXPORT DataArrayDoubleTuple : public DataArrayTuple<double>
  {
  public:
    DataArrayDoubleTuple(double *pt, std::size_t nbOfComp);
    std::string repr() const;
    double doubleValue() const;
    DataArrayDouble *buildDADouble(std::size_t nbOfTuples, std::size_t nbOfCompo) const;
  };

  class DataArrayFloatTuple;

  class MEDCOUPLING_EXPORT DataArrayFloatIterator : public DataArrayIterator<float>
  {
  public:
    DataArrayFloatIterator(DataArrayFloat *da);
    ~DataArrayFloatIterator() { }
  };

  class MEDCOUPLING_EXPORT DataArrayFloatTuple : public DataArrayTuple<float>
  {
  public:
    DataArrayFloatTuple(float *pt, std::size_t nbOfComp);
    std::string repr() const;
    float floatValue() const;
    DataArrayFloat *buildDAFloat(std::size_t nbOfTuples, std::size_t nbOfCompo) const;
  };

  class MEDCOUPLING_EXPORT DataArrayInt32Iterator : public DataArrayIterator<Int32>
  {
  public:
    DataArrayInt32Iterator(DataArrayInt32 *da);
    ~DataArrayInt32Iterator() { }
  };

  class MEDCOUPLING_EXPORT DataArrayInt64Iterator : public DataArrayIterator<Int64>
  {
  public:
    DataArrayInt64Iterator(DataArrayInt64 *da);
    ~DataArrayInt64Iterator() { }
  };

  class MEDCOUPLING_EXPORT  DataArrayInt32Tuple : public DataArrayTuple<Int32>
  {
  public:
     DataArrayInt32Tuple(Int32 *pt, std::size_t nbOfComp);
    std::string repr() const;
    Int32 intValue() const;
    DataArrayInt32 *buildDAInt(std::size_t nbOfTuples, std::size_t nbOfCompo) const;
  };

  class MEDCOUPLING_EXPORT DataArrayInt64Tuple : public DataArrayTuple<Int64>
  {
  public:
     DataArrayInt64Tuple(Int64 *pt, std::size_t nbOfComp);
    std::string repr() const;
    Int64 intValue() const;
    DataArrayInt64 *buildDAInt(std::size_t nbOfTuples, std::size_t nbOfCompo) const;
  };

  typedef DataArrayInt32Tuple DataArrayIntTuple;

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
    mcIdType _tuple_id;
    std::size_t _nb_comp;
    mcIdType _nb_tuple;
  };

  class MEDCOUPLING_EXPORT DataArrayAsciiCharTuple
  {
  public:
    DataArrayAsciiCharTuple(char *pt, std::size_t nbOfComp);
    std::string repr() const;
    std::size_t getNumberOfCompo() const { return _nb_of_compo; }
    const char *getConstPointer() const { return  _pt; }
    char *getPointer() { return _pt; }
    char asciiCharValue() const;
    DataArrayAsciiChar *buildDAAsciiChar(std::size_t nbOfTuples, std::size_t nbOfCompo) const;
  private:
    char *_pt;
    std::size_t _nb_of_compo;
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
    mcIdType _tuple_id;
    std::size_t _nb_comp;
    mcIdType _nb_tuple;
  };

  class MEDCOUPLING_EXPORT DataArrayByteTuple
  {
  public:
    DataArrayByteTuple(char *pt, std::size_t nbOfComp);
    std::string repr() const;
    std::size_t getNumberOfCompo() const { return _nb_of_compo; }
    const char *getConstPointer() const { return  _pt; }
    char *getPointer() { return _pt; }
    char byteValue() const;
    DataArrayByte *buildDAByte(std::size_t nbOfTuples, std::size_t nbOfCompo) const;
  private:
    char *_pt;
    std::size_t _nb_of_compo;
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

  template<class T>
  template<class InputIterator>
  void DataArrayTemplate<T>::insertAtTheEnd(InputIterator first, InputIterator last)
  {
    std::size_t nbCompo(this->getNumberOfComponents());
    if(nbCompo==1)
      this->_mem.insertAtTheEnd(first,last);
    else if(nbCompo==0)
      {
        this->_info_on_compo.resize(1);
        this->_mem.insertAtTheEnd(first,last);
      }
    else
      throw INTERP_KERNEL::Exception("DataArrayDouble::insertAtTheEnd : not available for DataArrayDouble with number of components different than 1 !");
  }
}
