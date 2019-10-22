// Copyright (C) 2007-2019  CEA/DEN, EDF R&D
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

#ifndef __MEDCOUPLING_MEDCOUPLINGMEMARRAY_HXX__
#define __MEDCOUPLING_MEDCOUPLINGMEMARRAY_HXX__

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
    MEDCOUPLING_EXPORT MEDCouplingPointer():_internal(0),_external(0) { }
    MEDCOUPLING_EXPORT void null() { _internal=0; _external=0; }
    MEDCOUPLING_EXPORT bool isNull() const { return _internal==0 && _external==0; }
    MEDCOUPLING_EXPORT void setInternal(T *pointer);
    MEDCOUPLING_EXPORT void setExternal(const T *pointer);
    MEDCOUPLING_EXPORT const T *getConstPointer() const { if(_internal) return _internal; else return _external; }
    MEDCOUPLING_EXPORT const T *getConstPointerLoc(std::size_t offset) const { if(_internal) return _internal+offset; else return _external+offset; }
    MEDCOUPLING_EXPORT T *getPointer() { if(_internal) return _internal; if(_external) throw INTERP_KERNEL::Exception("Trying to write on an external pointer."); else return 0; }
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
    MEDCOUPLING_EXPORT MemArray():_nb_of_elem(0),_nb_of_elem_alloc(0),_ownership(false),_dealloc(0),_param_for_deallocator(0) { }
    MEDCOUPLING_EXPORT MemArray(const MemArray<T>& other);
    MEDCOUPLING_EXPORT bool isNull() const { return _pointer.isNull(); }
    MEDCOUPLING_EXPORT const T *getConstPointerLoc(std::size_t offset) const { return _pointer.getConstPointerLoc(offset); }
    MEDCOUPLING_EXPORT const T *getConstPointer() const { return _pointer.getConstPointer(); }
    MEDCOUPLING_EXPORT std::size_t getNbOfElem() const { return _nb_of_elem; }
    MEDCOUPLING_EXPORT std::size_t getNbOfElemAllocated() const { return _nb_of_elem_alloc; }
    MEDCOUPLING_EXPORT T *getPointer() { return _pointer.getPointer(); }
    MEDCOUPLING_EXPORT MemArray<T> &operator=(const MemArray<T>& other);
    MEDCOUPLING_EXPORT T operator[](std::size_t id) const { return _pointer.getConstPointer()[id]; }
    MEDCOUPLING_EXPORT T& operator[](std::size_t id) { return _pointer.getPointer()[id]; }
    MEDCOUPLING_EXPORT bool isEqual(const MemArray<T>& other, T prec, std::string& reason) const;
    MEDCOUPLING_EXPORT void repr(mcIdType sl, std::ostream& stream) const;
    MEDCOUPLING_EXPORT bool reprHeader(mcIdType sl, std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprZip(mcIdType sl, std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprNotTooLong(mcIdType sl, std::ostream& stream) const;
    MEDCOUPLING_EXPORT void fillWithValue(const T& val);
    MEDCOUPLING_EXPORT T *fromNoInterlace(std::size_t nbOfComp) const;
    MEDCOUPLING_EXPORT T *toNoInterlace(std::size_t nbOfComp) const;
    MEDCOUPLING_EXPORT void sort(bool asc);
    MEDCOUPLING_EXPORT void reverse(std::size_t nbOfComp);
    MEDCOUPLING_EXPORT void alloc(std::size_t nbOfElements);
    MEDCOUPLING_EXPORT void reserve(std::size_t newNbOfElements);
    MEDCOUPLING_EXPORT void reAlloc(std::size_t newNbOfElements);
    MEDCOUPLING_EXPORT void useArray(const T *array, bool ownership, DeallocType type, std::size_t nbOfElem);
    MEDCOUPLING_EXPORT void useExternalArrayWithRWAccess(const T *array, std::size_t nbOfElem);
    MEDCOUPLING_EXPORT void writeOnPlace(std::size_t id, T element0, const T *others, std::size_t sizeOfOthers);
    template<class InputIterator>
    void insertAtTheEnd(InputIterator first, InputIterator last);
    MEDCOUPLING_EXPORT void pushBack(T elem);
    MEDCOUPLING_EXPORT T popBack();
    MEDCOUPLING_EXPORT void pack() const;
    MEDCOUPLING_EXPORT bool isDeallocatorCalled() const { return _ownership; }
    MEDCOUPLING_EXPORT Deallocator getDeallocator() const { return _dealloc; }
    MEDCOUPLING_EXPORT void setSpecificDeallocator(Deallocator dealloc) { _dealloc=dealloc; }
    MEDCOUPLING_EXPORT void setParameterForDeallocator(void *param) { _param_for_deallocator=param; }
    MEDCOUPLING_EXPORT void *getParameterForDeallocator() const { return _param_for_deallocator; }
    MEDCOUPLING_EXPORT void destroy();
    MEDCOUPLING_EXPORT ~MemArray() { destroy(); }
  public:
    MEDCOUPLING_EXPORT static void CPPDeallocator(void *pt, void *param);
    MEDCOUPLING_EXPORT static void CDeallocator(void *pt, void *param);
    MEDCOUPLING_EXPORT static void COffsetDeallocator(void *pt, void *param);
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
    MEDCOUPLING_EXPORT static void GetSlice(T start, T stop, T step, mcIdType sliceId, mcIdType nbOfSlices, T& startSlice, T& stopSlice);
    MEDCOUPLING_EXPORT static mcIdType GetNumberOfItemGivenBES(T begin, T end, T step, const std::string& msg);
    MEDCOUPLING_EXPORT static mcIdType GetNumberOfItemGivenBESRelative(T begin, T end, T step, const std::string& msg);
    MEDCOUPLING_EXPORT static mcIdType GetPosOfItemGivenBESRelativeNoThrow(T value, T begin, T end, T step);
  };

  class DataArrayByte;

  class DataArray : public RefCountObject, public TimeLabel
  {
  public:
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDCOUPLING_EXPORT void setName(const std::string& name);
    MEDCOUPLING_EXPORT void copyStringInfoFrom(const DataArray& other);
    MEDCOUPLING_EXPORT void copyPartOfStringInfoFrom(const DataArray& other, const std::vector<std::size_t>& compoIds);
    MEDCOUPLING_EXPORT void copyPartOfStringInfoFrom2(const std::vector<std::size_t>& compoIds, const DataArray& other);
    MEDCOUPLING_EXPORT bool areInfoEqualsIfNotWhy(const DataArray& other, std::string& reason) const;
    MEDCOUPLING_EXPORT bool areInfoEquals(const DataArray& other) const;
    MEDCOUPLING_EXPORT std::string cppRepr(const std::string& varName) const;
    MEDCOUPLING_EXPORT std::string getName() const { return _name; }
    MEDCOUPLING_EXPORT const std::vector<std::string> &getInfoOnComponents() const { return _info_on_compo; }
    MEDCOUPLING_EXPORT std::vector<std::string> &getInfoOnComponents() { return _info_on_compo; }
    MEDCOUPLING_EXPORT void setInfoOnComponents(const std::vector<std::string>& info);
    MEDCOUPLING_EXPORT void setInfoAndChangeNbOfCompo(const std::vector<std::string>& info);
    MEDCOUPLING_EXPORT std::vector<std::string> getVarsOnComponent() const;
    MEDCOUPLING_EXPORT std::vector<std::string> getUnitsOnComponent() const;
    MEDCOUPLING_EXPORT std::string getInfoOnComponent(std::size_t i) const;
    MEDCOUPLING_EXPORT std::string getVarOnComponent(std::size_t i) const;
    MEDCOUPLING_EXPORT std::string getUnitOnComponent(std::size_t i) const;
    MEDCOUPLING_EXPORT void setInfoOnComponent(std::size_t i, const std::string& info);
    MEDCOUPLING_EXPORT std::size_t getNumberOfComponents() const { return _info_on_compo.size(); }
    MEDCOUPLING_EXPORT void setPartOfValuesBase3(const DataArray *aBase, const mcIdType *bgTuples, const mcIdType *endTuples, mcIdType bgComp, mcIdType endComp, mcIdType stepComp, bool strictCompoCompare=true);
    MEDCOUPLING_EXPORT virtual void *getVoidStarPointer() = 0;
    MEDCOUPLING_EXPORT virtual DataArray *deepCopy() const = 0;
    MEDCOUPLING_EXPORT virtual DataArray *buildNewEmptyInstance() const = 0;
    MEDCOUPLING_EXPORT virtual bool isAllocated() const = 0;
    MEDCOUPLING_EXPORT virtual void checkAllocated() const = 0;
    MEDCOUPLING_EXPORT virtual void desallocate() = 0;
    MEDCOUPLING_EXPORT virtual mcIdType getNumberOfTuples() const = 0;
    MEDCOUPLING_EXPORT virtual mcIdType getNbOfElems() const = 0;
    MEDCOUPLING_EXPORT virtual std::size_t getNbOfElemAllocated() const = 0;
    MEDCOUPLING_EXPORT virtual void alloc(std::size_t nbOfTuple, std::size_t nbOfCompo=1) = 0;
    MEDCOUPLING_EXPORT virtual void reAlloc(std::size_t newNbOfTuple) = 0;
    MEDCOUPLING_EXPORT virtual void renumberInPlace(const mcIdType *old2New) = 0;
    MEDCOUPLING_EXPORT virtual void renumberInPlaceR(const mcIdType *new2Old) = 0;
    MEDCOUPLING_EXPORT virtual void setContigPartOfSelectedValues(mcIdType tupleIdStart, const DataArray *aBase, const DataArrayIdType *tuplesSelec) = 0;
    MEDCOUPLING_EXPORT virtual void setContigPartOfSelectedValuesSlice(mcIdType tupleIdStart, const DataArray *aBase, mcIdType bg, mcIdType end2, mcIdType step) = 0;
    MEDCOUPLING_EXPORT virtual DataArray *selectByTupleRanges(const std::vector<std::pair<mcIdType,mcIdType> >& ranges) const = 0;
    MEDCOUPLING_EXPORT virtual DataArray *keepSelectedComponents(const std::vector<std::size_t>& compoIds) const = 0;
    MEDCOUPLING_EXPORT virtual DataArray *selectByTupleId(const mcIdType *new2OldBg, const mcIdType *new2OldEnd) const = 0;
    MEDCOUPLING_EXPORT virtual DataArray *selectByTupleIdSafe(const mcIdType *new2OldBg, const mcIdType *new2OldEnd) const = 0;
    MEDCOUPLING_EXPORT virtual DataArray *selectByTupleIdSafeSlice(mcIdType bg, mcIdType end2, mcIdType step) const = 0;
    MEDCOUPLING_EXPORT virtual void rearrange(std::size_t newNbOfCompo) = 0;
    MEDCOUPLING_EXPORT virtual void circularPermutation(mcIdType nbOfShift=1) = 0;
    MEDCOUPLING_EXPORT virtual void circularPermutationPerTuple(mcIdType nbOfShift=1) = 0;
    MEDCOUPLING_EXPORT virtual void reversePerTuple() = 0;
    MEDCOUPLING_EXPORT void checkNbOfTuples(mcIdType nbOfTuples, const std::string& msg) const;
    MEDCOUPLING_EXPORT void checkNbOfComps(std::size_t nbOfCompo, const std::string& msg) const;
    MEDCOUPLING_EXPORT void checkNbOfTuplesAndComp(const DataArray& other, const std::string& msg) const;
    MEDCOUPLING_EXPORT void checkNbOfTuplesAndComp(mcIdType nbOfTuples, std::size_t nbOfCompo, const std::string& msg) const;
    MEDCOUPLING_EXPORT void checkNbOfElems(mcIdType nbOfElems, const std::string& msg) const;
    MEDCOUPLING_EXPORT static void GetSlice(mcIdType start, mcIdType stop, mcIdType step, mcIdType sliceId, mcIdType nbOfSlices, mcIdType& startSlice, mcIdType& stopSlice);
    MEDCOUPLING_EXPORT static mcIdType GetNumberOfItemGivenBES(mcIdType begin, mcIdType end, mcIdType step, const std::string& msg);
    MEDCOUPLING_EXPORT static mcIdType GetNumberOfItemGivenBESRelative(mcIdType begin, mcIdType end, mcIdType step, const std::string& msg);
    MEDCOUPLING_EXPORT static mcIdType GetPosOfItemGivenBESRelativeNoThrow(mcIdType value, mcIdType begin, mcIdType end, mcIdType step);
    MEDCOUPLING_EXPORT static std::string GetVarNameFromInfo(const std::string& info);
    MEDCOUPLING_EXPORT static std::string GetUnitFromInfo(const std::string& info);
    MEDCOUPLING_EXPORT static std::string BuildInfoFromVarAndUnit(const std::string& var, const std::string& unit);
    MEDCOUPLING_EXPORT static std::string GetAxisTypeRepr(MEDCouplingAxisType at);
    MEDCOUPLING_EXPORT static DataArray *Aggregate(const std::vector<const DataArray *>& arrs);
    MEDCOUPLING_EXPORT virtual void reprStream(std::ostream& stream) const = 0;
    MEDCOUPLING_EXPORT virtual void reprZipStream(std::ostream& stream) const = 0;
    MEDCOUPLING_EXPORT virtual void reprWithoutNameStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT virtual void reprZipWithoutNameStream(std::ostream& stream) const = 0;
    MEDCOUPLING_EXPORT virtual void reprCppStream(const std::string& varName, std::ostream& stream) const = 0;
    MEDCOUPLING_EXPORT virtual void reprQuickOverview(std::ostream& stream) const = 0;
    MEDCOUPLING_EXPORT virtual void reprQuickOverviewData(std::ostream& stream, std::size_t maxNbOfByteInRepr) const = 0;
  protected:
    DataArray() { }
    ~DataArray() { }
  protected:
    MEDCOUPLING_EXPORT static void CheckValueInRange(mcIdType ref, mcIdType value, const std::string& msg);
    MEDCOUPLING_EXPORT static void CheckValueInRangeEx(mcIdType value, mcIdType start, mcIdType end, const std::string& msg);
    MEDCOUPLING_EXPORT static void CheckClosingParInRange(mcIdType ref, mcIdType value, const std::string& msg);
    MEDCOUPLING_EXPORT static mcIdType EffectiveCircPerm(mcIdType nbOfShift, mcIdType nbOfTuples);
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
    MEDCOUPLING_EXPORT static MCAuto< typename Traits<T>::ArrayTypeCh > NewFromStdVector(const typename std::vector<T>& v);
    MEDCOUPLING_EXPORT std::vector< MCAuto< typename Traits<T>::ArrayTypeCh > > explodeComponents() const;
    //
    std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT void updateTime() const { }
    //
    MEDCOUPLING_EXPORT mcIdType getNumberOfTuples() const { return ToIdType(_info_on_compo.empty()?0:_mem.getNbOfElem()/getNumberOfComponents()); }
    MEDCOUPLING_EXPORT mcIdType getNbOfElems() const { return ToIdType(_mem.getNbOfElem()); }
    MEDCOUPLING_EXPORT bool empty() const;
    MEDCOUPLING_EXPORT void *getVoidStarPointer() { return getPointer(); }
    MEDCOUPLING_EXPORT const T *getConstPointer() const { return _mem.getConstPointer(); }
    MEDCOUPLING_EXPORT const T *begin() const { return getConstPointer(); }
    MEDCOUPLING_EXPORT const T *end() const { return getConstPointer()+getNbOfElems(); }
    MEDCOUPLING_EXPORT T *rwBegin() { return getPointer(); }
    MEDCOUPLING_EXPORT T *rwEnd() { return getPointer()+getNbOfElems(); }
    MEDCOUPLING_EXPORT void alloc(std::size_t nbOfTuple, std::size_t nbOfCompo=1);
    MEDCOUPLING_EXPORT void useArray(const T *array, bool ownership, DeallocType type, std::size_t nbOfTuple, std::size_t nbOfCompo);
    MEDCOUPLING_EXPORT void useExternalArrayWithRWAccess(const T *array, std::size_t nbOfTuple, std::size_t nbOfCompo);
    MEDCOUPLING_EXPORT T getIJSafe(std::size_t tupleId, std::size_t compoId) const;
    MEDCOUPLING_EXPORT T getIJ(std::size_t tupleId, std::size_t compoId) const { return _mem[tupleId*_info_on_compo.size()+compoId]; }
    MEDCOUPLING_EXPORT void setIJ(std::size_t tupleId, std::size_t compoId, T newVal) { _mem[tupleId*_info_on_compo.size()+compoId]=newVal; declareAsNew(); }
    MEDCOUPLING_EXPORT void setIJSilent(std::size_t tupleId, std::size_t compoId, T newVal) { _mem[tupleId*_info_on_compo.size()+compoId]=newVal; }
    MEDCOUPLING_EXPORT T *getPointer() { return _mem.getPointer(); declareAsNew(); }
    MEDCOUPLING_EXPORT void pack() const;
    MEDCOUPLING_EXPORT bool isAllocated() const override;
    MEDCOUPLING_EXPORT void checkAllocated() const;
    MEDCOUPLING_EXPORT void desallocate();
    MEDCOUPLING_EXPORT void reserve(std::size_t nbOfElems);
    MEDCOUPLING_EXPORT void rearrange(std::size_t newNbOfCompo);
    MEDCOUPLING_EXPORT void transpose();
    MEDCOUPLING_EXPORT void pushBackSilent(T val);
    MEDCOUPLING_EXPORT void pushBackValsSilent(const T *valsBg, const T *valsEnd);
    MEDCOUPLING_EXPORT T popBackSilent();
    MEDCOUPLING_EXPORT T front() const;
    MEDCOUPLING_EXPORT T back() const;
    MEDCOUPLING_EXPORT std::size_t getNbOfElemAllocated() const { return _mem.getNbOfElemAllocated(); }
    MEDCOUPLING_EXPORT void allocIfNecessary(std::size_t nbOfTuple, std::size_t nbOfCompo);
    MEDCOUPLING_EXPORT void deepCopyFrom(const DataArrayTemplate<T>& other);
    MEDCOUPLING_EXPORT void reverse();
    MEDCOUPLING_EXPORT void fillWithValue(T val);
    MEDCOUPLING_EXPORT void reAlloc(std::size_t newNbOfTuple);
    MEDCOUPLING_EXPORT void renumberInPlace(const mcIdType *old2New);
    MEDCOUPLING_EXPORT void renumberInPlaceR(const mcIdType *new2Old);
    MEDCOUPLING_EXPORT void sort(bool asc=true);
    MEDCOUPLING_EXPORT typename Traits<T>::ArrayType *renumber(const mcIdType *old2New) const;
    MEDCOUPLING_EXPORT typename Traits<T>::ArrayType *renumberR(const mcIdType *new2Old) const;
    MEDCOUPLING_EXPORT typename Traits<T>::ArrayType *renumberAndReduce(const mcIdType *old2New, mcIdType newNbOfTuple) const;
    MEDCOUPLING_EXPORT typename Traits<T>::ArrayType *changeNbOfComponents(std::size_t newNbOfComp, T dftValue) const;
    MEDCOUPLING_EXPORT typename Traits<T>::ArrayType *subArray(mcIdType tupleIdBg, mcIdType tupleIdEnd=-1) const;
    MEDCOUPLING_EXPORT MCAuto<typename Traits<T>::ArrayTypeCh> selectPartDef(const PartDefinition* pd) const;
    MEDCOUPLING_EXPORT void circularPermutation(mcIdType nbOfShift=1);
    MEDCOUPLING_EXPORT void circularPermutationPerTuple(mcIdType nbOfShift=1);
    MEDCOUPLING_EXPORT void reversePerTuple();
    MEDCOUPLING_EXPORT void setPartOfValues1(const typename Traits<T>::ArrayType *a, mcIdType bgTuples, mcIdType endTuples, mcIdType stepTuples, mcIdType bgComp, mcIdType endComp, mcIdType stepComp, bool strictCompoCompare=true);
    MEDCOUPLING_EXPORT void setPartOfValuesSimple1(T a, mcIdType bgTuples, mcIdType endTuples, mcIdType stepTuples, mcIdType bgComp, mcIdType endComp, mcIdType stepComp);
    MEDCOUPLING_EXPORT void setPartOfValues2(const typename Traits<T>::ArrayType *a, const mcIdType *bgTuples, const mcIdType *endTuples, const mcIdType *bgComp, const mcIdType *endComp, bool strictCompoCompare=true);
    MEDCOUPLING_EXPORT void setPartOfValuesSimple2(T a, const mcIdType *bgTuples, const mcIdType *endTuples, const mcIdType *bgComp, const mcIdType *endComp);
    MEDCOUPLING_EXPORT void setPartOfValues3(const typename Traits<T>::ArrayType *a, const mcIdType *bgTuples, const mcIdType *endTuples, mcIdType bgComp, mcIdType endComp, mcIdType stepComp, bool strictCompoCompare=true);
    MEDCOUPLING_EXPORT void setPartOfValuesSimple3(T a, const mcIdType *bgTuples, const mcIdType *endTuples, mcIdType bgComp, mcIdType endComp, mcIdType stepComp);
    MEDCOUPLING_EXPORT void setPartOfValues4(const typename Traits<T>::ArrayType *a, mcIdType bgTuples, mcIdType endTuples, mcIdType stepTuples, const mcIdType *bgComp, const mcIdType *endComp, bool strictCompoCompare=true);
    MEDCOUPLING_EXPORT void setPartOfValuesSimple4(T a, mcIdType bgTuples, mcIdType endTuples, mcIdType stepTuples, const mcIdType *bgComp, const mcIdType *endComp);
    MEDCOUPLING_EXPORT void setPartOfValuesAdv(const typename Traits<T>::ArrayType *a, const DataArrayIdType *tuplesSelec);
    MEDCOUPLING_EXPORT void setContigPartOfSelectedValues(mcIdType tupleIdStart, const DataArray *aBase, const DataArrayIdType *tuplesSelec);
    MEDCOUPLING_EXPORT void setContigPartOfSelectedValuesSlice(mcIdType tupleIdStart, const DataArray *aBase, mcIdType bg, mcIdType end2, mcIdType step);
    MEDCOUPLING_EXPORT T getMaxValue(mcIdType& tupleId) const;
    MEDCOUPLING_EXPORT T getMaxValueInArray() const;
    MEDCOUPLING_EXPORT T getMaxAbsValue(std::size_t& tupleId) const;
    MEDCOUPLING_EXPORT T getMaxAbsValueInArray() const;
    MEDCOUPLING_EXPORT T getMinValue(mcIdType& tupleId) const;
    MEDCOUPLING_EXPORT T getMinValueInArray() const;
    MEDCOUPLING_EXPORT void getTuple(mcIdType tupleId, T *res) const { std::copy(_mem.getConstPointerLoc(tupleId*_info_on_compo.size()),_mem.getConstPointerLoc((tupleId+1)*_info_on_compo.size()),res); }
    template<class InputIterator>
    void insertAtTheEnd(InputIterator first, InputIterator last);
    MEDCOUPLING_EXPORT static void SetArrayIn(typename Traits<T>::ArrayType *newArray, typename Traits<T>::ArrayType* &arrayToSet);
    MEDCOUPLING_EXPORT void writeOnPlace(std::size_t id, T element0, const T *others, mcIdType sizeOfOthers) { _mem.writeOnPlace(id,element0,others,sizeOfOthers); }
    MEDCOUPLING_EXPORT void fillWithZero();
  public:
    MEDCOUPLING_EXPORT MemArray<T>& accessToMemArray() { return _mem; }
    MEDCOUPLING_EXPORT const MemArray<T>& accessToMemArray() const { return _mem; }
  protected:
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
    MEDCOUPLING_EXPORT MCAuto<DataArrayDouble> convertToDblArr() const;
    MEDCOUPLING_EXPORT MCAuto<DataArrayInt> convertToIntArr() const;
    MEDCOUPLING_EXPORT MCAuto<DataArrayFloat> convertToFloatArr() const;
    MEDCOUPLING_EXPORT void applyLin(T a, T b, std::size_t compoId);
    MEDCOUPLING_EXPORT void applyLin(T a, T b);
    MEDCOUPLING_EXPORT typename Traits<T>::ArrayType *negate() const;
    MEDCOUPLING_EXPORT void addEqual(const typename Traits<T>::ArrayType *other);
    MEDCOUPLING_EXPORT void substractEqual(const typename Traits<T>::ArrayType *other);
    MEDCOUPLING_EXPORT void multiplyEqual(const typename Traits<T>::ArrayType *other);
    MEDCOUPLING_EXPORT void divideEqual(const typename Traits<T>::ArrayType *other);
    MEDCOUPLING_EXPORT static typename Traits<T>::ArrayType *Substract(const typename Traits<T>::ArrayType *a1, const typename Traits<T>::ArrayType *a2);
    MEDCOUPLING_EXPORT static typename Traits<T>::ArrayType *Divide(const typename Traits<T>::ArrayType *a1, const typename Traits<T>::ArrayType *a2);
    MEDCOUPLING_EXPORT static typename Traits<T>::ArrayType *Add(const typename Traits<T>::ArrayType *a1, const typename Traits<T>::ArrayType *a2);
    MEDCOUPLING_EXPORT static typename Traits<T>::ArrayType *Multiply(const typename Traits<T>::ArrayType *a1, const typename Traits<T>::ArrayType *a2);
    MEDCOUPLING_EXPORT static typename Traits<T>::ArrayType *Meld(const typename Traits<T>::ArrayType *a1, const typename Traits<T>::ArrayType *a2);
    MEDCOUPLING_EXPORT static typename Traits<T>::ArrayType *Meld(const std::vector<const typename Traits<T>::ArrayType *>& arr);
    MEDCOUPLING_EXPORT MCAuto<DataArrayIdType> findIdsGreaterOrEqualTo(T val) const;
    MEDCOUPLING_EXPORT MCAuto<DataArrayIdType> findIdsGreaterThan(T val) const;
    MEDCOUPLING_EXPORT MCAuto<DataArrayIdType> findIdsLowerOrEqualTo(T val) const;
    MEDCOUPLING_EXPORT MCAuto<DataArrayIdType> findIdsLowerThan(T val) const;
    MEDCOUPLING_EXPORT DataArrayIdType *findIdsStrictlyNegative() const;
    MEDCOUPLING_EXPORT typename Traits<T>::ArrayType *fromNoInterlace() const;
    MEDCOUPLING_EXPORT typename Traits<T>::ArrayType *toNoInterlace() const;
    MEDCOUPLING_EXPORT void meldWith(const typename Traits<T>::ArrayType *other);
    MEDCOUPLING_EXPORT typename Traits<T>::ArrayType *duplicateEachTupleNTimes(mcIdType nbTimes) const;
    MEDCOUPLING_EXPORT void aggregate(const typename Traits<T>::ArrayType *other);
    MEDCOUPLING_EXPORT void abs();
    MEDCOUPLING_EXPORT typename Traits<T>::ArrayType *computeAbs() const;
    MEDCOUPLING_EXPORT typename Traits<T>::ArrayType *performCopyOrIncrRef(bool dCpy) const;
    MEDCOUPLING_EXPORT typename Traits<T>::ArrayType *sumPerTuple() const;
    MEDCOUPLING_EXPORT void iota(T init=(T)0);
    MEDCOUPLING_EXPORT void reprStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprZipStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprNotTooLongStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprWithoutNameStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprZipWithoutNameStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprNotTooLongWithoutNameStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT std::string repr() const;
    MEDCOUPLING_EXPORT std::string reprZip() const;
    MEDCOUPLING_EXPORT std::string reprNotTooLong() const;
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
    MEDCOUPLING_EXPORT bool isUniform(T val, T eps) const;
  };
}

namespace MEDCoupling
{
  class DataArrayFloatIterator;
  class DataArrayFloat : public DataArrayTemplateFP<float>
  {
  public:
    MEDCOUPLING_EXPORT static DataArrayFloat *New();
  public:// abstract method overload
    MEDCOUPLING_EXPORT DataArrayFloat *deepCopy() const;
    MEDCOUPLING_EXPORT DataArrayFloat *buildNewEmptyInstance() const { return DataArrayFloat::New(); }
    MEDCOUPLING_EXPORT DataArrayFloat *selectByTupleRanges(const std::vector<std::pair<mcIdType,mcIdType> >& ranges) const { return DataArrayTemplateFP<float>::mySelectByTupleRanges(ranges); }
    MEDCOUPLING_EXPORT DataArrayFloat *keepSelectedComponents(const std::vector<std::size_t>& compoIds) const { return DataArrayTemplateFP<float>::myKeepSelectedComponents(compoIds); }
    MEDCOUPLING_EXPORT DataArrayFloat *selectByTupleId(const mcIdType *new2OldBg, const mcIdType *new2OldEnd) const { return this->mySelectByTupleId(new2OldBg,new2OldEnd); }
    MEDCOUPLING_EXPORT DataArrayFloat *selectByTupleIdSafe(const mcIdType *new2OldBg, const mcIdType *new2OldEnd) const { return DataArrayTemplateFP<float>::mySelectByTupleIdSafe(new2OldBg,new2OldEnd); }
    MEDCOUPLING_EXPORT DataArrayFloat *selectByTupleIdSafeSlice(mcIdType bg, mcIdType end2, mcIdType step) const { return DataArrayTemplateFP<float>::mySelectByTupleIdSafeSlice(bg,end2,step); }
    MEDCOUPLING_EXPORT void reprCppStream(const std::string& varName, std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprQuickOverview(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprQuickOverviewData(std::ostream& stream, std::size_t maxNbOfByteInRepr) const;
  public:// non abstract but essential
    MEDCOUPLING_EXPORT bool isEqual(const DataArrayFloat& other, float prec) const;
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const DataArrayFloat& other, float prec, std::string& reason) const;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const DataArrayFloat& other, float prec) const;
  public:
    MEDCOUPLING_EXPORT DataArrayFloatIterator *iterator();
  private:
    ~DataArrayFloat() { }
    DataArrayFloat() { }
  };
}

namespace MEDCoupling
{
  class DataArrayDoubleIterator;
  class DataArrayDouble : public DataArrayTemplateFP<double>
  {
  public:
    MEDCOUPLING_EXPORT static DataArrayDouble *New();
    MEDCOUPLING_EXPORT double doubleValue() const;
    MEDCOUPLING_EXPORT DataArrayDouble *deepCopy() const;
    MEDCOUPLING_EXPORT DataArrayDouble *buildNewEmptyInstance() const { return DataArrayDouble::New(); }
    MEDCOUPLING_EXPORT void checkMonotonic(bool increasing, double eps) const;
    MEDCOUPLING_EXPORT bool isMonotonic(bool increasing, double eps) const;
    MEDCOUPLING_EXPORT void writeVTK(std::ostream& ofs, mcIdType indent, const std::string& nameInFile, DataArrayByte *byteArr) const;
    MEDCOUPLING_EXPORT void reprCppStream(const std::string& varName, std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprQuickOverview(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprQuickOverviewData(std::ostream& stream, std::size_t maxNbOfByteInRepr) const;
    MEDCOUPLING_EXPORT bool isEqual(const DataArrayDouble& other, double prec) const;
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const DataArrayDouble& other, double prec, std::string& reason) const;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const DataArrayDouble& other, double prec) const;
    MEDCOUPLING_EXPORT DataArrayDouble *selectByTupleId(const mcIdType *new2OldBg, const mcIdType *new2OldEnd) const { return this->mySelectByTupleId(new2OldBg,new2OldEnd); }
    MEDCOUPLING_EXPORT DataArrayDouble *selectByTupleId(const DataArrayIdType& di) const { return this->mySelectByTupleId(di); }
    MEDCOUPLING_EXPORT DataArrayDouble *selectByTupleIdSafe(const mcIdType *new2OldBg, const mcIdType *new2OldEnd) const { return DataArrayTemplateFP<double>::mySelectByTupleIdSafe(new2OldBg,new2OldEnd); }
    MEDCOUPLING_EXPORT DataArrayDouble *keepSelectedComponents(const std::vector<std::size_t>& compoIds) const { return DataArrayTemplateFP<double>::myKeepSelectedComponents(compoIds); }
    MEDCOUPLING_EXPORT DataArrayDouble *selectByTupleIdSafeSlice(mcIdType bg, mcIdType end2, mcIdType step) const { return DataArrayTemplateFP<double>::mySelectByTupleIdSafeSlice(bg,end2,step); }
    MEDCOUPLING_EXPORT DataArrayDouble *selectByTupleRanges(const std::vector<std::pair<mcIdType,mcIdType> >& ranges) const { return DataArrayTemplateFP<double>::mySelectByTupleRanges(ranges); }
    MEDCOUPLING_EXPORT bool areIncludedInMe(const DataArrayDouble *other, double prec, DataArrayIdType *&tupleIds) const;
    MEDCOUPLING_EXPORT void findCommonTuples(double prec, mcIdType limitTupleId, DataArrayIdType *&comm, DataArrayIdType *&commIndex) const;
    MEDCOUPLING_EXPORT double minimalDistanceTo(const DataArrayDouble *other, mcIdType& thisTupleId, mcIdType& otherTupleId) const;
    MEDCOUPLING_EXPORT DataArrayDouble *getDifferentValues(double prec, mcIdType limitTupleId=-1) const;
    MEDCOUPLING_EXPORT DataArrayIdType *findClosestTupleId(const DataArrayDouble *other) const;
    MEDCOUPLING_EXPORT DataArrayIdType *computeNbOfInteractionsWith(const DataArrayDouble *otherBBoxFrmt, double eps) const;
    MEDCOUPLING_EXPORT void setSelectedComponents(const DataArrayDouble *a, const std::vector<std::size_t>& compoIds);
    MEDCOUPLING_EXPORT DataArrayDoubleIterator *iterator();
    MEDCOUPLING_EXPORT void checkNoNullValues() const;
    MEDCOUPLING_EXPORT void getMinMaxPerComponent(double *bounds) const;
    MEDCOUPLING_EXPORT DataArrayDouble *computeBBoxPerTuple(double epsilon=0.0) const;
    MEDCOUPLING_EXPORT void computeTupleIdsNearTuples(const DataArrayDouble *other, double eps, DataArrayIdType *& c, DataArrayIdType *& cI) const;
    MEDCOUPLING_EXPORT void recenterForMaxPrecision(double eps);
    MEDCOUPLING_EXPORT double getMaxValue2(DataArrayIdType*& tupleIds) const;
    MEDCOUPLING_EXPORT double getMinValue2(DataArrayIdType*& tupleIds) const;
    MEDCOUPLING_EXPORT mcIdType count(double value, double eps) const;
    MEDCOUPLING_EXPORT double getAverageValue() const;
    MEDCOUPLING_EXPORT double norm2() const;
    MEDCOUPLING_EXPORT double normMax() const;
    MEDCOUPLING_EXPORT void normMaxPerComponent(double * res) const;
    MEDCOUPLING_EXPORT double normMin() const;
    MEDCOUPLING_EXPORT void accumulate(double *res) const;
    MEDCOUPLING_EXPORT double accumulate(std::size_t compId) const;
    MEDCOUPLING_EXPORT DataArrayDouble *accumulatePerChunck(const mcIdType *bgOfIndex, const mcIdType *endOfIndex) const;
    MEDCOUPLING_EXPORT MCAuto<DataArrayDouble> cumSum() const;
    MEDCOUPLING_EXPORT double distanceToTuple(const double *tupleBg, const double *tupleEnd, mcIdType& tupleId) const;
    MEDCOUPLING_EXPORT DataArrayDouble *fromPolarToCart() const;
    MEDCOUPLING_EXPORT DataArrayDouble *fromCylToCart() const;
    MEDCOUPLING_EXPORT DataArrayDouble *fromSpherToCart() const;
    MEDCOUPLING_EXPORT DataArrayDouble *cartesianize(MEDCouplingAxisType atOfThis) const;
    MEDCOUPLING_EXPORT DataArrayDouble *fromCartToPolar() const;
    MEDCOUPLING_EXPORT DataArrayDouble *fromCartToCyl() const;
    MEDCOUPLING_EXPORT DataArrayDouble *fromCartToSpher() const;
    MEDCOUPLING_EXPORT DataArrayDouble *fromCartToCylGiven(const DataArrayDouble *coords, const double center[3], const double vect[3]) const;
    MEDCOUPLING_EXPORT DataArrayDouble *doublyContractedProduct() const;
    MEDCOUPLING_EXPORT DataArrayDouble *determinant() const;
    MEDCOUPLING_EXPORT DataArrayDouble *eigenValues() const;
    MEDCOUPLING_EXPORT DataArrayDouble *eigenVectors() const;
    MEDCOUPLING_EXPORT DataArrayDouble *inverse() const;
    MEDCOUPLING_EXPORT DataArrayDouble *trace() const;
    MEDCOUPLING_EXPORT DataArrayDouble *deviator() const;
    MEDCOUPLING_EXPORT DataArrayDouble *magnitude() const;
    MEDCOUPLING_EXPORT DataArrayDouble *maxPerTuple() const;
    MEDCOUPLING_EXPORT DataArrayDouble *maxPerTupleWithCompoId(DataArrayIdType* &compoIdOfMaxPerTuple) const;
    MEDCOUPLING_EXPORT DataArrayDouble *buildEuclidianDistanceDenseMatrix() const;
    MEDCOUPLING_EXPORT DataArrayDouble *buildEuclidianDistanceDenseMatrixWith(const DataArrayDouble *other) const;
    MEDCOUPLING_EXPORT void asArcOfCircle(double center[2], double& radius, double& ang) const;
    MEDCOUPLING_EXPORT void sortPerTuple(bool asc);
    MEDCOUPLING_EXPORT void applyInv(double numerator);
    MEDCOUPLING_EXPORT void applyPow(double val);
    MEDCOUPLING_EXPORT void applyRPow(double val);
    MEDCOUPLING_EXPORT DataArrayDouble *applyFunc(std::size_t nbOfComp, FunctionToEvaluate func) const;
    MEDCOUPLING_EXPORT DataArrayDouble *applyFunc(std::size_t nbOfComp, const std::string& func, bool isSafe=true) const;
    MEDCOUPLING_EXPORT DataArrayDouble *applyFunc(const std::string& func, bool isSafe=true) const;
    MEDCOUPLING_EXPORT void applyFuncOnThis(const std::string& func, bool isSafe=true);
    MEDCOUPLING_EXPORT DataArrayDouble *applyFuncCompo(std::size_t nbOfComp, const std::string& func, bool isSafe=true) const;
    MEDCOUPLING_EXPORT DataArrayDouble *applyFuncNamedCompo(std::size_t nbOfComp, const std::vector<std::string>& varsOrder, const std::string& func, bool isSafe=true) const;
    MEDCOUPLING_EXPORT void applyFuncFast32(const std::string& func);
    MEDCOUPLING_EXPORT void applyFuncFast64(const std::string& func);
    MEDCOUPLING_EXPORT MCAuto<DataArrayDouble> symmetry3DPlane(const double point[3], const double normalVector[3]) const;
    MEDCOUPLING_EXPORT DataArrayIdType *findIdsInRange(double vmin, double vmax) const;
    MEDCOUPLING_EXPORT DataArrayIdType *findIdsNotInRange(double vmin, double vmax) const;
    MEDCOUPLING_EXPORT static DataArrayDouble *Aggregate(const DataArrayDouble *a1, const DataArrayDouble *a2);
    MEDCOUPLING_EXPORT static DataArrayDouble *Aggregate(const std::vector<const DataArrayDouble *>& arr);
    MEDCOUPLING_EXPORT static DataArrayDouble *Dot(const DataArrayDouble *a1, const DataArrayDouble *a2);
    MEDCOUPLING_EXPORT static DataArrayDouble *CrossProduct(const DataArrayDouble *a1, const DataArrayDouble *a2);
    MEDCOUPLING_EXPORT static DataArrayDouble *Max(const DataArrayDouble *a1, const DataArrayDouble *a2);
    MEDCOUPLING_EXPORT static DataArrayDouble *Min(const DataArrayDouble *a1, const DataArrayDouble *a2);
    MEDCOUPLING_EXPORT static DataArrayDouble *Pow(const DataArrayDouble *a1, const DataArrayDouble *a2);
    MEDCOUPLING_EXPORT void powEqual(const DataArrayDouble *other);
    MEDCOUPLING_EXPORT std::vector<bool> toVectorOfBool(double eps) const;
    MEDCOUPLING_EXPORT static void Rotate2DAlg(const double *center, double angle, mcIdType nbNodes, const double *coordsIn, double *coordsOut);
    MEDCOUPLING_EXPORT static void Rotate3DAlg(const double *center, const double *vect, double angle, mcIdType nbNodes, const double *coordsIn, double *coordsOut);
    MEDCOUPLING_EXPORT static void Symmetry3DPlane(const double point[3], const double normalVector[3], mcIdType nbNodes, const double *coordsIn, double *coordsOut);
    MEDCOUPLING_EXPORT static void GiveBaseForPlane(const double normalVector[3], double baseOfPlane[9]);
    MEDCOUPLING_EXPORT static void ComputeIntegralOfSeg2IntoTri3(const double seg2[4], const double tri3[6], double coeffs[3], double& length);
  public:
    MEDCOUPLING_EXPORT void getTinySerializationIntInformation(std::vector<mcIdType>& tinyInfo) const;
    MEDCOUPLING_EXPORT void getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const;
    MEDCOUPLING_EXPORT bool resizeForUnserialization(const std::vector<mcIdType>& tinyInfoI);
    MEDCOUPLING_EXPORT void finishUnserialization(const std::vector<mcIdType>& tinyInfoI, const std::vector<std::string>& tinyInfoS);
  public:
    template<mcIdType SPACEDIM>
    void findCommonTuplesAlg(const double *bbox, mcIdType nbNodes, mcIdType limitNodeId, double prec, DataArrayIdType *c, DataArrayIdType *cI) const;
    template<mcIdType SPACEDIM>
    static void FindClosestTupleIdAlg(const BBTreePts<SPACEDIM,mcIdType>& myTree, double dist, const double *pos, mcIdType nbOfTuples, const double *thisPt, mcIdType thisNbOfTuples, mcIdType *res);
    template<mcIdType SPACEDIM>
    static void FindTupleIdsNearTuplesAlg(const BBTreePts<SPACEDIM,mcIdType>& myTree, const double *pos, mcIdType nbOfTuples, double eps,
                                          DataArrayIdType *c, DataArrayIdType *cI);
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
    typedef typename Traits<T>::ArrayType DataArrayType;
  public:
    MEDCOUPLING_EXPORT static DataArrayType *New();
    MEDCOUPLING_EXPORT T intValue() const;
    MEDCOUPLING_EXPORT bool isEqual(const DataArrayDiscrete<T>& other) const;
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const DataArrayDiscrete<T>& other, std::string& reason) const;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const DataArrayDiscrete<T>& other) const;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStrAndOrder(const typename Traits<T>::ArrayType& other) const;
    MEDCOUPLING_EXPORT void switchOnTupleEqualTo(T val, std::vector<bool>& vec) const;
    MEDCOUPLING_EXPORT void switchOnTupleNotEqualTo(T val, std::vector<bool>& vec) const;
    MEDCOUPLING_EXPORT DataArrayIdType *buildPermutationArr(const DataArrayDiscrete<T>& other) const;
    MEDCOUPLING_EXPORT DataArrayIdType *indicesOfSubPart(const DataArrayDiscrete<T>& partOfThis) const;
    MEDCOUPLING_EXPORT void checkMonotonic(bool increasing) const;
    MEDCOUPLING_EXPORT bool isMonotonic(bool increasing) const;
    MEDCOUPLING_EXPORT void checkStrictlyMonotonic(bool increasing) const;
    MEDCOUPLING_EXPORT bool isStrictlyMonotonic(bool increasing) const;
    MEDCOUPLING_EXPORT mcIdType getHashCode() const;
    MEDCOUPLING_EXPORT void reprCppStream(const std::string& varName, std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprQuickOverview(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprQuickOverviewData(std::ostream& stream, std::size_t maxNbOfByteInRepr) const;
    MEDCOUPLING_EXPORT void writeVTK(std::ostream& ofs, mcIdType indent, const std::string& type, const std::string& nameInFile, DataArrayByte *byteArr) const;
    MEDCOUPLING_EXPORT void transformWithIndArr(const T *indArrBg, const T *indArrEnd);
    MEDCOUPLING_EXPORT void transformWithIndArr(const MapKeyVal<T, T>& m);
    MEDCOUPLING_EXPORT DataArrayIdType *findIdsEqual(T val) const;
    MEDCOUPLING_EXPORT DataArrayIdType *transformWithIndArrR(const T *indArr2Bg, const T *indArrEnd) const;
    MEDCOUPLING_EXPORT void splitByValueRange(const T *arrBg, const T *arrEnd,
                                              DataArrayType *& castArr, DataArrayType *& rankInsideCast, DataArrayType *& castsPresent) const;
    MEDCOUPLING_EXPORT bool isRange(T& strt, T& sttoopp, T& stteepp) const;
    MEDCOUPLING_EXPORT DataArrayIdType *invertArrayO2N2N2O(mcIdType newNbOfElem) const;
    MEDCOUPLING_EXPORT DataArrayIdType *invertArrayN2O2O2N(mcIdType oldNbOfElem) const;
    MEDCOUPLING_EXPORT DataArrayIdType *invertArrayO2N2N2OBis(mcIdType newNbOfElem) const;
    MEDCOUPLING_EXPORT MCAuto< MapKeyVal<T, mcIdType> > invertArrayN2O2O2NOptimized() const;
    MEDCOUPLING_EXPORT MCAuto< MapKeyVal<mcIdType, T> > giveN2OOptimized() const;
    MEDCOUPLING_EXPORT MCAuto<DataArrayIdType> findIdForEach(const T *valsBg, const T *valsEnd) const;
    MEDCOUPLING_EXPORT DataArrayIdType *checkAndPreparePermutation() const;
    MEDCOUPLING_EXPORT void changeSurjectiveFormat(T targetNb, DataArrayIdType *&arr, DataArrayIdType *&arrI) const;
    MEDCOUPLING_EXPORT DataArrayIdType *buildPermArrPerLevel() const;
    MEDCOUPLING_EXPORT bool isIota(mcIdType sizeExpected) const;
    MEDCOUPLING_EXPORT bool isUniform(T val) const;
    MEDCOUPLING_EXPORT T checkUniformAndGuess() const;
    MEDCOUPLING_EXPORT bool hasUniqueValues() const;
    MEDCOUPLING_EXPORT void setSelectedComponents(const DataArrayType *a, const std::vector<std::size_t>& compoIds);
    MEDCOUPLING_EXPORT DataArrayIdType *findIdsNotEqual(T val) const;
    MEDCOUPLING_EXPORT DataArrayIdType *findIdsEqualTuple(const T *tupleBg, const T *tupleEnd) const;
    MEDCOUPLING_EXPORT DataArrayIdType *findIdsEqualList(const T *valsBg, const T *valsEnd) const;
    MEDCOUPLING_EXPORT DataArrayIdType *findIdsNotEqualList(const T *valsBg, const T *valsEnd) const;
    MEDCOUPLING_EXPORT mcIdType findIdFirstEqual(T value) const;
    MEDCOUPLING_EXPORT mcIdType findIdFirstEqual(const std::vector<T>& vals) const;
    MEDCOUPLING_EXPORT mcIdType findIdFirstEqualTuple(const std::vector<T>& tupl) const;
    MEDCOUPLING_EXPORT mcIdType findIdSequence(const std::vector<T>& vals) const;
    MEDCOUPLING_EXPORT mcIdType changeValue(T oldValue, T newValue);
    MEDCOUPLING_EXPORT mcIdType count(T value) const;
    MEDCOUPLING_EXPORT bool presenceOfTuple(const std::vector<T>& tupl) const;
    MEDCOUPLING_EXPORT bool presenceOfValue(T value) const;
    MEDCOUPLING_EXPORT bool presenceOfValue(const std::vector<T>& vals) const;
    MEDCOUPLING_EXPORT void accumulate(T *res) const;
    MEDCOUPLING_EXPORT T accumulate(std::size_t compId) const;
    MEDCOUPLING_EXPORT DataArrayType *accumulatePerChunck(const mcIdType *bgOfIndex, const mcIdType *endOfIndex) const;
    MEDCOUPLING_EXPORT void getMinMaxValues(T& minValue, T& maxValue) const;
    MEDCOUPLING_EXPORT void applyInv(T numerator);
    MEDCOUPLING_EXPORT void applyDivideBy(T val);
    MEDCOUPLING_EXPORT void applyModulus(T val);
    MEDCOUPLING_EXPORT void applyRModulus(T val);
    MEDCOUPLING_EXPORT void applyPow(T val);
    MEDCOUPLING_EXPORT void applyRPow(T val);
    MEDCOUPLING_EXPORT DataArrayIdType *findIdsInRange(T vmin, T vmax) const;
    MEDCOUPLING_EXPORT DataArrayIdType *findIdsNotInRange(T vmin, T vmax) const;
    MEDCOUPLING_EXPORT bool checkAllIdsInRange(T vmin, T vmax) const;
    MEDCOUPLING_EXPORT static DataArrayType *Aggregate(const DataArrayType *a1, const DataArrayType *a2, T offsetA2);
    MEDCOUPLING_EXPORT static DataArrayType *Aggregate(const std::vector<const DataArrayType *>& arr);
    MEDCOUPLING_EXPORT static DataArrayType *AggregateIndexes(const std::vector<const DataArrayType *>& arrs);
    MEDCOUPLING_EXPORT static DataArrayType *BuildUnion(const std::vector<const DataArrayType *>& arr);
    MEDCOUPLING_EXPORT static DataArrayType *BuildIntersection(const std::vector<const DataArrayType *>& arr);
    MEDCOUPLING_EXPORT static void PutIntoToSkylineFrmt(const std::vector< std::vector<T> >& v, DataArrayType *& data, DataArrayIdType *& dataIndex);
    MEDCOUPLING_EXPORT DataArrayIdType *buildComplement(mcIdType nbOfElement) const;
    MEDCOUPLING_EXPORT DataArrayType *buildSubstraction(const DataArrayType *other) const;
    MEDCOUPLING_EXPORT DataArrayType *buildSubstractionOptimized(const DataArrayType *other) const;
    MEDCOUPLING_EXPORT DataArrayType *buildUnion(const DataArrayType *other) const;
    MEDCOUPLING_EXPORT DataArrayType *buildIntersection(const DataArrayType *other) const;
    MEDCOUPLING_EXPORT DataArrayType *buildUnique() const;
    MEDCOUPLING_EXPORT DataArrayType *buildUniqueNotSorted() const;
    MEDCOUPLING_EXPORT DataArrayType *deltaShiftIndex() const;
    MEDCOUPLING_EXPORT void computeOffsets();
    MEDCOUPLING_EXPORT void computeOffsetsFull();
    MEDCOUPLING_EXPORT void findIdsRangesInListOfIds(const DataArrayType *listOfIds, DataArrayIdType *& rangeIdsFetched, DataArrayType *& idsInInputListThatFetch) const;
    MEDCOUPLING_EXPORT DataArrayType *buildExplicitArrByRanges(const DataArrayType *offsets) const;
    MEDCOUPLING_EXPORT DataArrayType *buildExplicitArrOfSliceOnScaledArr(T begin, T stop, T step) const;
    MEDCOUPLING_EXPORT DataArrayIdType *findRangeIdForEachTuple(const DataArrayType *ranges) const;
    MEDCOUPLING_EXPORT DataArrayType *findIdInRangeForEachTuple(const DataArrayType *ranges) const;
    MEDCOUPLING_EXPORT void sortEachPairToMakeALinkedList();
    MEDCOUPLING_EXPORT MCAuto<DataArrayType> fromLinkedListOfPairToList() const;
    MEDCOUPLING_EXPORT DataArrayType *getDifferentValues() const;
    MEDCOUPLING_EXPORT std::vector<DataArrayIdType *> partitionByDifferentValues(std::vector<T>& differentIds) const;
    MEDCOUPLING_EXPORT std::vector< std::pair<mcIdType,mcIdType> > splitInBalancedSlices(mcIdType nbOfSlices) const;
    MEDCOUPLING_EXPORT static DataArrayType *Modulus(const DataArrayType *a1, const DataArrayType *a2);
    MEDCOUPLING_EXPORT void modulusEqual(const DataArrayType *other);
    MEDCOUPLING_EXPORT static DataArrayType *Pow(const DataArrayType *a1, const DataArrayType *a2);
    MEDCOUPLING_EXPORT void powEqual(const DataArrayType *other);
    //MEDCOUPLING_EXPORT MemArray<T>& accessToMemArray() { return _mem; }
    //MEDCOUPLING_EXPORT const MemArray<T>& accessToMemArray() const { return _mem; }
  public:
    MEDCOUPLING_EXPORT static DataArrayIdType *FindPermutationFromFirstToSecond(const DataArrayType *ids1, const DataArrayType *ids2);
    MEDCOUPLING_EXPORT static mcIdType *CheckAndPreparePermutation(const T *start, const T *end);
    MEDCOUPLING_EXPORT static DataArrayType *BuildListOfSwitchedOn(const std::vector<bool>& v);
    MEDCOUPLING_EXPORT static DataArrayType *BuildListOfSwitchedOff(const std::vector<bool>& v);
    MEDCOUPLING_EXPORT static DataArrayIdType *ConvertIndexArrayToO2N(mcIdType nbOfOldTuples, const mcIdType *arr, const mcIdType *arrIBg, const mcIdType *arrIEnd, mcIdType &newNbOfTuples);
    MEDCOUPLING_EXPORT static DataArrayIdType *MakePartition(const std::vector<const DataArrayType *>& groups, mcIdType newNb, std::vector< std::vector<mcIdType> >& fidsOfGroups);
  public:
    MEDCOUPLING_EXPORT static void ExtractFromIndexedArrays(const mcIdType *idsOfSelectBg, const mcIdType *idsOfSelectEnd,
                                                            const DataArrayType *arrIn, const DataArrayIdType *arrIndxIn,
                                                            DataArrayType* &arrOut, DataArrayIdType* &arrIndexOut);
    MEDCOUPLING_EXPORT static void ExtractFromIndexedArraysSlice(mcIdType idsOfSelectStart, mcIdType idsOfSelectStop, mcIdType idsOfSelectStep,
                                                                 const DataArrayType *arrIn, const DataArrayIdType *arrIndxIn,
                                                                 DataArrayType* &arrOut, DataArrayIdType* &arrIndexOut);
    MEDCOUPLING_EXPORT static void SetPartOfIndexedArrays(const mcIdType *idsOfSelectBg, const mcIdType *idsOfSelectEnd,
                                                          const DataArrayType *arrIn, const DataArrayIdType *arrIndxIn,
                                                          const DataArrayType *srcArr, const DataArrayIdType *srcArrIndex,
                                                          DataArrayType* &arrOut, DataArrayIdType* &arrIndexOut);
    MEDCOUPLING_EXPORT static void SetPartOfIndexedArraysSlice(mcIdType start, mcIdType end, mcIdType step,
                                                               const DataArrayType *arrIn, const DataArrayIdType *arrIndxIn,
                                                               const DataArrayType *srcArr, const DataArrayIdType *srcArrIndex,
                                                               DataArrayType* &arrOut, DataArrayIdType* &arrIndexOut);
    MEDCOUPLING_EXPORT static void SetPartOfIndexedArraysSameIdx(const mcIdType *idsOfSelectBg, const mcIdType *idsOfSelectEnd,
                                                                 DataArrayType *arrInOut, const DataArrayIdType *arrIndxIn,
                                                                 const DataArrayType *srcArr, const DataArrayIdType *srcArrIndex);
    MEDCOUPLING_EXPORT static void SetPartOfIndexedArraysSameIdxSlice(mcIdType start, mcIdType end, mcIdType step,
                                                                      DataArrayType *arrInOut, const DataArrayIdType *arrIndxIn,
                                                                      const DataArrayType *srcArr, const DataArrayIdType *srcArrIndex);
    MEDCOUPLING_EXPORT static bool RemoveIdsFromIndexedArrays(const T *idsToRemoveBg, const T *idsToRemoveEnd,
                                                              DataArrayType *arr, DataArrayIdType *arrIndx, mcIdType offsetForRemoval=0);
    MEDCOUPLING_EXPORT static DataArrayType *Range(T begin, T end, T step);
  public:
    MEDCOUPLING_EXPORT void getTinySerializationIntInformation(std::vector<mcIdType>& tinyInfo) const;
    MEDCOUPLING_EXPORT void getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const;
    MEDCOUPLING_EXPORT bool resizeForUnserialization(const std::vector<mcIdType>& tinyInfoI);
    MEDCOUPLING_EXPORT void finishUnserialization(const std::vector<mcIdType>& tinyInfoI, const std::vector<std::string>& tinyInfoS);
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
    MEDCOUPLING_EXPORT bool isFittingWith(const std::vector<bool>& v) const;
  protected:
    ~DataArrayDiscreteSigned() { }
  };

  class DataArrayInt32Iterator;

  class DataArrayInt32 : public DataArrayDiscreteSigned<Int32>
  {
    friend class DataArrayDiscrete<Int32>;
  public:
    MEDCOUPLING_EXPORT DataArrayInt32 *deepCopy() const;//ok
    MEDCOUPLING_EXPORT DataArrayInt32 *buildNewEmptyInstance() const { return DataArrayInt32::New(); }//ok
  public:
    MEDCOUPLING_EXPORT DataArrayInt32 *selectByTupleId(const mcIdType *new2OldBg, const mcIdType *new2OldEnd) const { return this->mySelectByTupleId(new2OldBg,new2OldEnd); }
    MEDCOUPLING_EXPORT DataArrayInt32 *selectByTupleId(const DataArrayIdType& di) const { return this->mySelectByTupleId(di); }
    MEDCOUPLING_EXPORT DataArrayInt32 *selectByTupleIdSafe(const mcIdType *new2OldBg, const mcIdType *new2OldEnd) const { return this->mySelectByTupleIdSafe(new2OldBg,new2OldEnd); }
    MEDCOUPLING_EXPORT DataArrayInt32 *keepSelectedComponents(const std::vector<std::size_t>& compoIds) const { return this->myKeepSelectedComponents(compoIds); }
    MEDCOUPLING_EXPORT DataArrayInt32 *selectByTupleIdSafeSlice(mcIdType bg, mcIdType end2, mcIdType step) const { return this->mySelectByTupleIdSafeSlice(bg,end2,step); }
    MEDCOUPLING_EXPORT DataArrayInt32 *selectByTupleRanges(const std::vector<std::pair<mcIdType,mcIdType> >& ranges) const { return this->mySelectByTupleRanges(ranges); }
  public:
    MEDCOUPLING_EXPORT DataArrayInt32Iterator *iterator();
  private:
    ~DataArrayInt32() { }
    DataArrayInt32() { }
  };

  class DataArrayInt64 : public DataArrayDiscreteSigned<Int64>
  {
    friend class DataArrayDiscrete<Int64>;
  public:
    MEDCOUPLING_EXPORT DataArrayInt64 *deepCopy() const;
    MEDCOUPLING_EXPORT DataArrayInt64 *buildNewEmptyInstance() const { return DataArrayInt64::New(); }//ok
  public:
    MEDCOUPLING_EXPORT DataArrayInt64 *selectByTupleId(const mcIdType *new2OldBg, const mcIdType *new2OldEnd) const { return this->mySelectByTupleId(new2OldBg,new2OldEnd); }
    MEDCOUPLING_EXPORT DataArrayInt64 *selectByTupleId(const DataArrayIdType& di) const { return this->mySelectByTupleId(di); }
    MEDCOUPLING_EXPORT DataArrayInt64 *selectByTupleIdSafe(const mcIdType *new2OldBg, const mcIdType *new2OldEnd) const { return DataArrayTemplate<Int64>::mySelectByTupleIdSafe(new2OldBg,new2OldEnd); }
    MEDCOUPLING_EXPORT DataArrayInt64 *keepSelectedComponents(const std::vector<std::size_t>& compoIds) const { return DataArrayTemplate<Int64>::myKeepSelectedComponents(compoIds); }
    MEDCOUPLING_EXPORT DataArrayInt64 *selectByTupleIdSafeSlice(mcIdType bg, mcIdType end2, mcIdType step) const { return DataArrayTemplate<Int64>::mySelectByTupleIdSafeSlice(bg,end2,step); }
    MEDCOUPLING_EXPORT DataArrayInt64 *selectByTupleRanges(const std::vector<std::pair<mcIdType,mcIdType> >& ranges) const { return DataArrayTemplate<Int64>::mySelectByTupleRanges(ranges); }
  public:
    MEDCOUPLING_EXPORT DataArrayInt64Iterator *iterator();
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

  class DataArrayChar : public DataArrayTemplate<char>
  {
  public:
    MEDCOUPLING_EXPORT virtual DataArrayChar *buildEmptySpecializedDAChar() const = 0;
    MEDCOUPLING_EXPORT mcIdType getHashCode() const;
    MEDCOUPLING_EXPORT bool isEqual(const DataArrayChar& other) const;
    MEDCOUPLING_EXPORT virtual bool isEqualIfNotWhy(const DataArrayChar& other, std::string& reason) const;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const DataArrayChar& other) const;
    MEDCOUPLING_EXPORT std::string repr() const;
    MEDCOUPLING_EXPORT std::string reprZip() const;
    MEDCOUPLING_EXPORT DataArrayInt *convertToIntArr() const;
    MEDCOUPLING_EXPORT DataArrayChar *selectByTupleId(const mcIdType *new2OldBg, const mcIdType *new2OldEnd) const { return this->mySelectByTupleId(new2OldBg,new2OldEnd); }
    MEDCOUPLING_EXPORT DataArrayChar *selectByTupleId(const DataArrayIdType& di) const { return this->mySelectByTupleId(di); }
    MEDCOUPLING_EXPORT DataArrayChar *selectByTupleIdSafe(const mcIdType *new2OldBg, const mcIdType *new2OldEnd) const { return DataArrayTemplate<char>::mySelectByTupleIdSafe(new2OldBg,new2OldEnd); }
    MEDCOUPLING_EXPORT DataArrayChar *keepSelectedComponents(const std::vector<std::size_t>& compoIds) const { return DataArrayTemplate<char>::myKeepSelectedComponents(compoIds); }
    MEDCOUPLING_EXPORT DataArrayChar *selectByTupleIdSafeSlice(mcIdType bg, mcIdType end2, mcIdType step) const { return DataArrayTemplate<char>::mySelectByTupleIdSafeSlice(bg,end2,step); }
    MEDCOUPLING_EXPORT bool isUniform(char val) const;
    MEDCOUPLING_EXPORT void meldWith(const DataArrayChar *other);
    MEDCOUPLING_EXPORT DataArray *selectByTupleRanges(const std::vector<std::pair<mcIdType,mcIdType> >& ranges) const { return DataArrayTemplate<char>::mySelectByTupleRanges(ranges); }
    MEDCOUPLING_EXPORT DataArrayIdType *findIdsEqual(char val) const;
    MEDCOUPLING_EXPORT DataArrayIdType *findIdsNotEqual(char val) const;
    MEDCOUPLING_EXPORT mcIdType findIdSequence(const std::vector<char>& vals) const;
    MEDCOUPLING_EXPORT mcIdType findIdFirstEqualTuple(const std::vector<char>& tupl) const;
    MEDCOUPLING_EXPORT mcIdType findIdFirstEqual(char value) const;
    MEDCOUPLING_EXPORT mcIdType findIdFirstEqual(const std::vector<char>& vals) const;
    MEDCOUPLING_EXPORT bool presenceOfTuple(const std::vector<char>& tupl) const;
    MEDCOUPLING_EXPORT bool presenceOfValue(char value) const;
    MEDCOUPLING_EXPORT bool presenceOfValue(const std::vector<char>& vals) const;
    MEDCOUPLING_EXPORT DataArrayIdType *findIdsInRange(char vmin, char vmax) const;
    MEDCOUPLING_EXPORT static DataArrayChar *Aggregate(const DataArrayChar *a1, const DataArrayChar *a2);
    MEDCOUPLING_EXPORT static DataArrayChar *Aggregate(const std::vector<const DataArrayChar *>& arr);
    MEDCOUPLING_EXPORT static DataArrayChar *Meld(const DataArrayChar *a1, const DataArrayChar *a2);
    MEDCOUPLING_EXPORT static DataArrayChar *Meld(const std::vector<const DataArrayChar *>& arr);
    MEDCOUPLING_EXPORT MemArray<char>& accessToMemArray() { return _mem; }
    MEDCOUPLING_EXPORT const MemArray<char>& accessToMemArray() const { return _mem; }
  public:
    //MEDCOUPLING_EXPORT void getTinySerializationIntInformation(std::vector<mcIdType>& tinyInfo) const;
    //MEDCOUPLING_EXPORT void getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const;
    //MEDCOUPLING_EXPORT bool resizeForUnserialization(const std::vector<mcIdType>& tinyInfoI);
    //MEDCOUPLING_EXPORT void finishUnserialization(const std::vector<mcIdType>& tinyInfoI, const std::vector<std::string>& tinyInfoS);
  protected:
    DataArrayChar() { }
  };

  class DataArrayByteIterator;

  class DataArrayByte : public DataArrayChar
  {
  public:
    MEDCOUPLING_EXPORT static DataArrayByte *New();
    MEDCOUPLING_EXPORT DataArrayChar *buildEmptySpecializedDAChar() const;
    MEDCOUPLING_EXPORT DataArrayByteIterator *iterator();
    MEDCOUPLING_EXPORT DataArrayByte *deepCopy() const;
    MEDCOUPLING_EXPORT DataArrayByte *performCopyOrIncrRef(bool deepCopy) const;
    MEDCOUPLING_EXPORT DataArrayByte *buildNewEmptyInstance() const { return DataArrayByte::New(); }
    MEDCOUPLING_EXPORT char byteValue() const;
    MEDCOUPLING_EXPORT void reprStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprZipStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprWithoutNameStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprZipWithoutNameStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprCppStream(const std::string& varName, std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprQuickOverview(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprQuickOverviewData(std::ostream& stream, std::size_t maxNbOfByteInRepr) const;
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const DataArrayChar& other, std::string& reason) const;
    MEDCOUPLING_EXPORT std::vector<bool> toVectorOfBool() const;
  private:
    ~DataArrayByte() { }
    DataArrayByte() { }
  };

  class DataArrayAsciiCharIterator;

  class DataArrayAsciiChar : public DataArrayChar
  {
  public:
    MEDCOUPLING_EXPORT static DataArrayAsciiChar *New();
    MEDCOUPLING_EXPORT static DataArrayAsciiChar *New(const std::string& st);
    MEDCOUPLING_EXPORT static DataArrayAsciiChar *New(const std::vector<std::string>& vst, char defaultChar);
    MEDCOUPLING_EXPORT DataArrayChar *buildEmptySpecializedDAChar() const;
    MEDCOUPLING_EXPORT DataArrayAsciiCharIterator *iterator();
    MEDCOUPLING_EXPORT DataArrayAsciiChar *deepCopy() const;
    MEDCOUPLING_EXPORT DataArrayAsciiChar *performCopyOrIncrRef(bool deepCopy) const;
    MEDCOUPLING_EXPORT DataArrayAsciiChar *buildNewEmptyInstance() const { return DataArrayAsciiChar::New(); }
    MEDCOUPLING_EXPORT char asciiCharValue() const;
    MEDCOUPLING_EXPORT void reprStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprZipStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprWithoutNameStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprZipWithoutNameStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprCppStream(const std::string& varName, std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprQuickOverview(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprQuickOverviewData(std::ostream& stream, std::size_t maxNbOfByteInRepr) const;
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const DataArrayChar& other, std::string& reason) const;
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
    MEDCOUPLING_EXPORT typename Traits<T>::ArrayTuple *nextt();
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
    MEDCOUPLING_EXPORT DataArrayTuple(T *pt, std::size_t nbOfComp);
    //MEDCOUPLING_EXPORT std::string repr() const;
    MEDCOUPLING_EXPORT std::size_t getNumberOfCompo() const { return _nb_of_compo; }
    MEDCOUPLING_EXPORT const T *getConstPointer() const { return  _pt; }
    MEDCOUPLING_EXPORT T *getPointer() { return _pt; }
    MEDCOUPLING_EXPORT typename Traits<T>::ArrayType *buildDA(std::size_t nbOfTuples, std::size_t nbOfCompo) const;
  protected:
    T zeValue() const;
  protected:
    T *_pt;
    std::size_t _nb_of_compo;
  };

  class DataArrayDoubleTuple;

  class DataArrayDoubleIterator : public DataArrayIterator<double>
  {
  public:
    MEDCOUPLING_EXPORT DataArrayDoubleIterator(DataArrayDouble *da);
    MEDCOUPLING_EXPORT ~DataArrayDoubleIterator() { }
  };

  class DataArrayDoubleTuple : public DataArrayTuple<double>
  {
  public:
    MEDCOUPLING_EXPORT DataArrayDoubleTuple(double *pt, std::size_t nbOfComp);
    MEDCOUPLING_EXPORT std::string repr() const;
    MEDCOUPLING_EXPORT double doubleValue() const;
    MEDCOUPLING_EXPORT DataArrayDouble *buildDADouble(std::size_t nbOfTuples, std::size_t nbOfCompo) const;
  };

  class DataArrayFloatTuple;

  class DataArrayFloatIterator : public DataArrayIterator<float>
  {
  public:
    MEDCOUPLING_EXPORT DataArrayFloatIterator(DataArrayFloat *da);
    MEDCOUPLING_EXPORT ~DataArrayFloatIterator() { }
  };

  class DataArrayFloatTuple : public DataArrayTuple<float>
  {
  public:
    MEDCOUPLING_EXPORT DataArrayFloatTuple(float *pt, std::size_t nbOfComp);
    MEDCOUPLING_EXPORT std::string repr() const;
    MEDCOUPLING_EXPORT float floatValue() const;
    MEDCOUPLING_EXPORT DataArrayFloat *buildDAFloat(std::size_t nbOfTuples, std::size_t nbOfCompo) const;
  };
  
  class DataArrayInt32Iterator : public DataArrayIterator<Int32>
  {
  public:
    MEDCOUPLING_EXPORT DataArrayInt32Iterator(DataArrayInt32 *da);
    MEDCOUPLING_EXPORT ~DataArrayInt32Iterator() { }
  };

  class DataArrayInt64Iterator : public DataArrayIterator<Int64>
  {
  public:
    MEDCOUPLING_EXPORT DataArrayInt64Iterator(DataArrayInt64 *da);
    MEDCOUPLING_EXPORT ~DataArrayInt64Iterator() { }
  };

  class DataArrayInt32Tuple : public DataArrayTuple<Int32>
  {
  public:
    MEDCOUPLING_EXPORT DataArrayInt32Tuple(Int32 *pt, std::size_t nbOfComp);
    MEDCOUPLING_EXPORT std::string repr() const;
    MEDCOUPLING_EXPORT Int32 intValue() const;
    MEDCOUPLING_EXPORT DataArrayInt32 *buildDAInt(std::size_t nbOfTuples, std::size_t nbOfCompo) const;
  };

  class DataArrayInt64Tuple : public DataArrayTuple<Int64>
  {
  public:
    MEDCOUPLING_EXPORT DataArrayInt64Tuple(Int64 *pt, std::size_t nbOfComp);
    MEDCOUPLING_EXPORT std::string repr() const;
    MEDCOUPLING_EXPORT Int64 intValue() const;
    MEDCOUPLING_EXPORT DataArrayInt64 *buildDAInt(std::size_t nbOfTuples, std::size_t nbOfCompo) const;
  };

  typedef DataArrayInt32Tuple DataArrayIntTuple;

  class DataArrayAsciiCharTuple;

  class DataArrayAsciiCharIterator
  {
  public:
    MEDCOUPLING_EXPORT DataArrayAsciiCharIterator(DataArrayAsciiChar *da);
    MEDCOUPLING_EXPORT ~DataArrayAsciiCharIterator();
    MEDCOUPLING_EXPORT DataArrayAsciiCharTuple *nextt();
  private:
    DataArrayAsciiChar *_da;
    char *_pt;
    mcIdType _tuple_id;
    std::size_t _nb_comp;
    mcIdType _nb_tuple;
  };

  class DataArrayAsciiCharTuple
  {
  public:
    MEDCOUPLING_EXPORT DataArrayAsciiCharTuple(char *pt, std::size_t nbOfComp);
    MEDCOUPLING_EXPORT std::string repr() const;
    MEDCOUPLING_EXPORT std::size_t getNumberOfCompo() const { return _nb_of_compo; }
    MEDCOUPLING_EXPORT const char *getConstPointer() const { return  _pt; }
    MEDCOUPLING_EXPORT char *getPointer() { return _pt; }
    MEDCOUPLING_EXPORT char asciiCharValue() const;
    MEDCOUPLING_EXPORT DataArrayAsciiChar *buildDAAsciiChar(std::size_t nbOfTuples, std::size_t nbOfCompo) const;
  private:
    char *_pt;
    std::size_t _nb_of_compo;
  };
  
  class DataArrayByteTuple;

  class DataArrayByteIterator
  {
  public:
    MEDCOUPLING_EXPORT DataArrayByteIterator(DataArrayByte *da);
    MEDCOUPLING_EXPORT ~DataArrayByteIterator();
    MEDCOUPLING_EXPORT DataArrayByteTuple *nextt();
  private:
    DataArrayByte *_da;
    char *_pt;
    mcIdType _tuple_id;
    std::size_t _nb_comp;
    mcIdType _nb_tuple;
  };

  class DataArrayByteTuple
  {
  public:
    MEDCOUPLING_EXPORT DataArrayByteTuple(char *pt, std::size_t nbOfComp);
    MEDCOUPLING_EXPORT std::string repr() const;
    MEDCOUPLING_EXPORT std::size_t getNumberOfCompo() const { return _nb_of_compo; }
    MEDCOUPLING_EXPORT const char *getConstPointer() const { return  _pt; }
    MEDCOUPLING_EXPORT char *getPointer() { return _pt; }
    MEDCOUPLING_EXPORT char byteValue() const;
    MEDCOUPLING_EXPORT DataArrayByte *buildDAByte(std::size_t nbOfTuples, std::size_t nbOfCompo) const;
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

#endif
