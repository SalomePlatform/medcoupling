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

  class DataArrayInt32;
  class DataArrayByte;

  class DataArray : public RefCountObject, public TimeLabel
  {
  public:
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDCOUPLING_EXPORT void setName(const std::string& name);
    MEDCOUPLING_EXPORT void copyStringInfoFrom(const DataArray& other);
    MEDCOUPLING_EXPORT void copyPartOfStringInfoFrom(const DataArray& other, const std::vector<int>& compoIds);
    MEDCOUPLING_EXPORT void copyPartOfStringInfoFrom2(const std::vector<int>& compoIds, const DataArray& other);
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
    MEDCOUPLING_EXPORT std::string getInfoOnComponent(int i) const;
    MEDCOUPLING_EXPORT std::string getVarOnComponent(int i) const;
    MEDCOUPLING_EXPORT std::string getUnitOnComponent(int i) const;
    MEDCOUPLING_EXPORT void setInfoOnComponent(int i, const std::string& info);
    MEDCOUPLING_EXPORT std::size_t getNumberOfComponents() const { return _info_on_compo.size(); }
    MEDCOUPLING_EXPORT void setPartOfValuesBase3(const DataArray *aBase, const int *bgTuples, const int *endTuples, int bgComp, int endComp, int stepComp, bool strictCompoCompare=true);
    MEDCOUPLING_EXPORT virtual void *getVoidStarPointer() = 0;
    MEDCOUPLING_EXPORT virtual DataArray *deepCopy() const = 0;
    MEDCOUPLING_EXPORT virtual DataArray *buildNewEmptyInstance() const = 0;
    MEDCOUPLING_EXPORT virtual bool isAllocated() const = 0;
    MEDCOUPLING_EXPORT virtual void checkAllocated() const = 0;
    MEDCOUPLING_EXPORT virtual void desallocate() = 0;
    MEDCOUPLING_EXPORT virtual std::size_t getNumberOfTuples() const = 0;
    MEDCOUPLING_EXPORT virtual std::size_t getNbOfElems() const = 0;
    MEDCOUPLING_EXPORT virtual std::size_t getNbOfElemAllocated() const = 0;
    MEDCOUPLING_EXPORT virtual void alloc(std::size_t nbOfTuple, std::size_t nbOfCompo=1) = 0;
    MEDCOUPLING_EXPORT virtual void reAlloc(std::size_t newNbOfTuple) = 0;
    MEDCOUPLING_EXPORT virtual void renumberInPlace(const int *old2New) = 0;
    MEDCOUPLING_EXPORT virtual void renumberInPlaceR(const int *new2Old) = 0;
    MEDCOUPLING_EXPORT virtual void setContigPartOfSelectedValues(int tupleIdStart, const DataArray *aBase, const DataArrayInt32 *tuplesSelec) = 0;
    MEDCOUPLING_EXPORT virtual void setContigPartOfSelectedValuesSlice(int tupleIdStart, const DataArray *aBase, int bg, int end2, int step) = 0;
    MEDCOUPLING_EXPORT virtual DataArray *selectByTupleRanges(const std::vector<std::pair<int,int> >& ranges) const = 0;
    MEDCOUPLING_EXPORT virtual DataArray *keepSelectedComponents(const std::vector<int>& compoIds) const = 0;
    MEDCOUPLING_EXPORT virtual DataArray *selectByTupleId(const int *new2OldBg, const int *new2OldEnd) const = 0;
    MEDCOUPLING_EXPORT virtual DataArray *selectByTupleIdSafe(const int *new2OldBg, const int *new2OldEnd) const = 0;
    MEDCOUPLING_EXPORT virtual DataArray *selectByTupleIdSafeSlice(int bg, int end2, int step) const = 0;
    MEDCOUPLING_EXPORT virtual void rearrange(int newNbOfCompo) = 0;
    MEDCOUPLING_EXPORT virtual void circularPermutation(int nbOfShift=1) = 0;
    MEDCOUPLING_EXPORT virtual void circularPermutationPerTuple(int nbOfShift=1) = 0;
    MEDCOUPLING_EXPORT virtual void reversePerTuple() = 0;
    MEDCOUPLING_EXPORT void checkNbOfTuples(int nbOfTuples, const std::string& msg) const;
    MEDCOUPLING_EXPORT void checkNbOfComps(int nbOfCompo, const std::string& msg) const;
    MEDCOUPLING_EXPORT void checkNbOfTuplesAndComp(const DataArray& other, const std::string& msg) const;
    MEDCOUPLING_EXPORT void checkNbOfTuplesAndComp(int nbOfTuples, int nbOfCompo, const std::string& msg) const;
    MEDCOUPLING_EXPORT void checkNbOfElems(std::size_t nbOfElems, const std::string& msg) const;
    MEDCOUPLING_EXPORT static void GetSlice(int start, int stop, int step, int sliceId, int nbOfSlices, int& startSlice, int& stopSlice);
    MEDCOUPLING_EXPORT static int GetNumberOfItemGivenBES(int begin, int end, int step, const std::string& msg);
    MEDCOUPLING_EXPORT static int GetNumberOfItemGivenBESRelative(int begin, int end, int step, const std::string& msg);
    MEDCOUPLING_EXPORT static int GetPosOfItemGivenBESRelativeNoThrow(int value, int begin, int end, int step);
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
    MEDCOUPLING_EXPORT static void CheckValueInRange(int ref, int value, const std::string& msg);
    MEDCOUPLING_EXPORT static void CheckValueInRangeEx(int value, int start, int end, const std::string& msg);
    MEDCOUPLING_EXPORT static void CheckClosingParInRange(int ref, int value, const std::string& msg);
    MEDCOUPLING_EXPORT static int EffectiveCircPerm(int nbOfShift, int nbOfTuples);
  protected:
    std::string _name;
    std::vector<std::string> _info_on_compo;
  };
}

namespace MEDCoupling
{
  class DataArrayInt32;

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
    MEDCOUPLING_EXPORT std::size_t getNumberOfTuples() const { return _info_on_compo.empty()?0:_mem.getNbOfElem()/getNumberOfComponents(); }
    MEDCOUPLING_EXPORT std::size_t getNbOfElems() const { return _mem.getNbOfElem(); }
    bool empty() const;
    MEDCOUPLING_EXPORT void *getVoidStarPointer() { return getPointer(); }
    MEDCOUPLING_EXPORT const T *getConstPointer() const { return _mem.getConstPointer(); }
    MEDCOUPLING_EXPORT const T *begin() const { return getConstPointer(); }
    MEDCOUPLING_EXPORT const T *end() const { return getConstPointer()+getNbOfElems(); }
    void alloc(std::size_t nbOfTuple, std::size_t nbOfCompo=1);
    void useArray(const T *array, bool ownership, DeallocType type, int nbOfTuple, int nbOfCompo);
    void useExternalArrayWithRWAccess(const T *array, int nbOfTuple, int nbOfCompo);
    T getIJSafe(int tupleId, int compoId) const;
    MEDCOUPLING_EXPORT T getIJ(int tupleId, int compoId) const { return _mem[tupleId*_info_on_compo.size()+compoId]; }
    MEDCOUPLING_EXPORT void setIJ(int tupleId, int compoId, T newVal) { _mem[tupleId*_info_on_compo.size()+compoId]=newVal; declareAsNew(); }
    MEDCOUPLING_EXPORT void setIJSilent(int tupleId, int compoId, T newVal) { _mem[tupleId*_info_on_compo.size()+compoId]=newVal; }
    MEDCOUPLING_EXPORT T *getPointer() { return _mem.getPointer(); declareAsNew(); }
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
    void setPartOfValuesAdv(const typename Traits<T>::ArrayType *a, const DataArrayInt32 *tuplesSelec);
    void setContigPartOfSelectedValues(int tupleIdStart, const DataArray *aBase, const DataArrayInt32 *tuplesSelec);
    void setContigPartOfSelectedValuesSlice(int tupleIdStart, const DataArray *aBase, int bg, int end2, int step);
    T getMaxValue(int& tupleId) const;
    T getMaxValueInArray() const;
    T getMinValue(int& tupleId) const;
    T getMinValueInArray() const;
    MEDCOUPLING_EXPORT void getTuple(int tupleId, T *res) const { std::copy(_mem.getConstPointerLoc(tupleId*_info_on_compo.size()),_mem.getConstPointerLoc((tupleId+1)*_info_on_compo.size()),res); }
    template<class InputIterator>
    void insertAtTheEnd(InputIterator first, InputIterator last);
    MEDCOUPLING_EXPORT static void SetArrayIn(typename Traits<T>::ArrayType *newArray, typename Traits<T>::ArrayType* &arrayToSet);
    MEDCOUPLING_EXPORT void writeOnPlace(std::size_t id, T element0, const T *others, int sizeOfOthers) { _mem.writeOnPlace(id,element0,others,sizeOfOthers); }
    MEDCOUPLING_EXPORT void fillWithZero();
  public:
    MEDCOUPLING_EXPORT MemArray<T>& accessToMemArray() { return _mem; }
    MEDCOUPLING_EXPORT const MemArray<T>& accessToMemArray() const { return _mem; }
  protected:
    typename Traits<T>::ArrayType *mySelectByTupleId(const int *new2OldBg, const int *new2OldEnd) const;
    typename Traits<T>::ArrayType *mySelectByTupleId(const DataArrayInt32& di) const;
    typename Traits<T>::ArrayType *mySelectByTupleIdSafe(const int *new2OldBg, const int *new2OldEnd) const;
    typename Traits<T>::ArrayType *myKeepSelectedComponents(const std::vector<int>& compoIds) const;
    typename Traits<T>::ArrayType *mySelectByTupleIdSafeSlice(int bg, int end2, int step) const;
    typename Traits<T>::ArrayType *mySelectByTupleRanges(const std::vector<std::pair<int,int> >& ranges) const;
  protected:
    MemArray<T> _mem;
  };

  template<class T>
  class DataArrayTemplateClassic : public DataArrayTemplate<T>
  {
  public:
    MEDCOUPLING_EXPORT MCAuto<DataArrayDouble> convertToDblArr() const;
    MEDCOUPLING_EXPORT MCAuto<DataArrayInt32> convertToIntArr() const;
    MEDCOUPLING_EXPORT MCAuto<DataArrayFloat> convertToFloatArr() const;
    MEDCOUPLING_EXPORT void applyLin(T a, T b, int compoId);
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
    MEDCOUPLING_EXPORT MCAuto<DataArrayInt32> findIdsGreaterOrEqualTo(T val) const;
    MEDCOUPLING_EXPORT MCAuto<DataArrayInt32> findIdsGreaterThan(T val) const;
    MEDCOUPLING_EXPORT MCAuto<DataArrayInt32> findIdsLowerOrEqualTo(T val) const;
    MEDCOUPLING_EXPORT MCAuto<DataArrayInt32> findIdsLowerThan(T val) const;
    MEDCOUPLING_EXPORT DataArrayInt32 *findIdsStrictlyNegative() const;
    MEDCOUPLING_EXPORT typename Traits<T>::ArrayType *fromNoInterlace() const;
    MEDCOUPLING_EXPORT typename Traits<T>::ArrayType *toNoInterlace() const;
    MEDCOUPLING_EXPORT void meldWith(const typename Traits<T>::ArrayType *other);
    MEDCOUPLING_EXPORT typename Traits<T>::ArrayType *duplicateEachTupleNTimes(int nbTimes) const;
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
    MEDCOUPLING_EXPORT std::string reprNotTooLong() const;
    template<class U>
    MCAuto< typename Traits<U>::ArrayType > convertToOtherTypeOfArr() const;
  protected:
    static typename Traits<T>::ArrayType *PerformCopyOrIncrRef(bool dCpy, const typename Traits<T>::ArrayType& self);
    template<class OP>
    MCAuto<DataArrayInt32> findIdsAdv(const OP& op) const;
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
    MEDCOUPLING_EXPORT DataArrayFloat *selectByTupleRanges(const std::vector<std::pair<int,int> >& ranges) const { return DataArrayTemplateFP<float>::mySelectByTupleRanges(ranges); }
    MEDCOUPLING_EXPORT DataArrayFloat *keepSelectedComponents(const std::vector<int>& compoIds) const { return DataArrayTemplateFP<float>::myKeepSelectedComponents(compoIds); }
    MEDCOUPLING_EXPORT DataArrayFloat *selectByTupleId(const int *new2OldBg, const int *new2OldEnd) const { return DataArrayTemplateFP<float>::mySelectByTupleId(new2OldBg,new2OldEnd); }
    MEDCOUPLING_EXPORT DataArrayFloat *selectByTupleIdSafe(const int *new2OldBg, const int *new2OldEnd) const { return DataArrayTemplateFP<float>::mySelectByTupleIdSafe(new2OldBg,new2OldEnd); }
    MEDCOUPLING_EXPORT DataArrayFloat *selectByTupleIdSafeSlice(int bg, int end2, int step) const { return DataArrayTemplateFP<float>::mySelectByTupleIdSafeSlice(bg,end2,step); }
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
    MEDCOUPLING_EXPORT std::string repr() const;
    MEDCOUPLING_EXPORT std::string reprZip() const;
    MEDCOUPLING_EXPORT void writeVTK(std::ostream& ofs, int indent, const std::string& nameInFile, DataArrayByte *byteArr) const;
    MEDCOUPLING_EXPORT void reprCppStream(const std::string& varName, std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprQuickOverview(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprQuickOverviewData(std::ostream& stream, std::size_t maxNbOfByteInRepr) const;
    MEDCOUPLING_EXPORT bool isEqual(const DataArrayDouble& other, double prec) const;
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const DataArrayDouble& other, double prec, std::string& reason) const;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const DataArrayDouble& other, double prec) const;
    MEDCOUPLING_EXPORT DataArrayDouble *selectByTupleId(const int *new2OldBg, const int *new2OldEnd) const { return DataArrayTemplateFP<double>::mySelectByTupleId(new2OldBg,new2OldEnd); }
    MEDCOUPLING_EXPORT DataArrayDouble *selectByTupleId(const DataArrayInt32& di) const { return DataArrayTemplateFP<double>::mySelectByTupleId(di); }
    MEDCOUPLING_EXPORT DataArrayDouble *selectByTupleIdSafe(const int *new2OldBg, const int *new2OldEnd) const { return DataArrayTemplateFP<double>::mySelectByTupleIdSafe(new2OldBg,new2OldEnd); }
    MEDCOUPLING_EXPORT DataArrayDouble *keepSelectedComponents(const std::vector<int>& compoIds) const { return DataArrayTemplateFP<double>::myKeepSelectedComponents(compoIds); }
    MEDCOUPLING_EXPORT DataArrayDouble *selectByTupleIdSafeSlice(int bg, int end2, int step) const { return DataArrayTemplateFP<double>::mySelectByTupleIdSafeSlice(bg,end2,step); }
    MEDCOUPLING_EXPORT DataArrayDouble *selectByTupleRanges(const std::vector<std::pair<int,int> >& ranges) const { return DataArrayTemplateFP<double>::mySelectByTupleRanges(ranges); }
    MEDCOUPLING_EXPORT bool areIncludedInMe(const DataArrayDouble *other, double prec, DataArrayInt32 *&tupleIds) const;
    MEDCOUPLING_EXPORT void findCommonTuples(double prec, int limitTupleId, DataArrayInt32 *&comm, DataArrayInt32 *&commIndex) const;
    MEDCOUPLING_EXPORT double minimalDistanceTo(const DataArrayDouble *other, int& thisTupleId, int& otherTupleId) const;
    MEDCOUPLING_EXPORT DataArrayDouble *getDifferentValues(double prec, int limitTupleId=-1) const;
    MEDCOUPLING_EXPORT DataArrayInt32 *findClosestTupleId(const DataArrayDouble *other) const;
    MEDCOUPLING_EXPORT DataArrayInt32 *computeNbOfInteractionsWith(const DataArrayDouble *otherBBoxFrmt, double eps) const;
    MEDCOUPLING_EXPORT void setSelectedComponents(const DataArrayDouble *a, const std::vector<int>& compoIds);
    MEDCOUPLING_EXPORT DataArrayDoubleIterator *iterator();
    MEDCOUPLING_EXPORT void checkNoNullValues() const;
    MEDCOUPLING_EXPORT void getMinMaxPerComponent(double *bounds) const;
    MEDCOUPLING_EXPORT DataArrayDouble *computeBBoxPerTuple(double epsilon=0.0) const;
    MEDCOUPLING_EXPORT void computeTupleIdsNearTuples(const DataArrayDouble *other, double eps, DataArrayInt32 *& c, DataArrayInt32 *& cI) const;
    MEDCOUPLING_EXPORT void recenterForMaxPrecision(double eps);
    MEDCOUPLING_EXPORT double getMaxValue2(DataArrayInt32*& tupleIds) const;
    MEDCOUPLING_EXPORT double getMinValue2(DataArrayInt32*& tupleIds) const;
    MEDCOUPLING_EXPORT int count(double value, double eps) const;
    MEDCOUPLING_EXPORT double getAverageValue() const;
    MEDCOUPLING_EXPORT double norm2() const;
    MEDCOUPLING_EXPORT double normMax() const;
    MEDCOUPLING_EXPORT double normMin() const;
    MEDCOUPLING_EXPORT void accumulate(double *res) const;
    MEDCOUPLING_EXPORT double accumulate(int compId) const;
    MEDCOUPLING_EXPORT DataArrayDouble *accumulatePerChunck(const int *bgOfIndex, const int *endOfIndex) const;
    MEDCOUPLING_EXPORT MCAuto<DataArrayDouble> cumSum() const;
    MEDCOUPLING_EXPORT double distanceToTuple(const double *tupleBg, const double *tupleEnd, int& tupleId) const;
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
    MEDCOUPLING_EXPORT DataArrayDouble *maxPerTupleWithCompoId(DataArrayInt32* &compoIdOfMaxPerTuple) const;
    MEDCOUPLING_EXPORT DataArrayDouble *buildEuclidianDistanceDenseMatrix() const;
    MEDCOUPLING_EXPORT DataArrayDouble *buildEuclidianDistanceDenseMatrixWith(const DataArrayDouble *other) const;
    MEDCOUPLING_EXPORT void sortPerTuple(bool asc);
    MEDCOUPLING_EXPORT void applyInv(double numerator);
    MEDCOUPLING_EXPORT void applyPow(double val);
    MEDCOUPLING_EXPORT void applyRPow(double val);
    MEDCOUPLING_EXPORT DataArrayDouble *applyFunc(int nbOfComp, FunctionToEvaluate func) const;
    MEDCOUPLING_EXPORT DataArrayDouble *applyFunc(int nbOfComp, const std::string& func, bool isSafe=true) const;
    MEDCOUPLING_EXPORT DataArrayDouble *applyFunc(const std::string& func, bool isSafe=true) const;
    MEDCOUPLING_EXPORT void applyFuncOnThis(const std::string& func, bool isSafe=true);
    MEDCOUPLING_EXPORT DataArrayDouble *applyFuncCompo(int nbOfComp, const std::string& func, bool isSafe=true) const;
    MEDCOUPLING_EXPORT DataArrayDouble *applyFuncNamedCompo(int nbOfComp, const std::vector<std::string>& varsOrder, const std::string& func, bool isSafe=true) const;
    MEDCOUPLING_EXPORT void applyFuncFast32(const std::string& func);
    MEDCOUPLING_EXPORT void applyFuncFast64(const std::string& func);
    MEDCOUPLING_EXPORT MCAuto<DataArrayDouble> symmetry3DPlane(const double point[3], const double normalVector[3]) const;
    MEDCOUPLING_EXPORT DataArrayInt32 *findIdsInRange(double vmin, double vmax) const;
    MEDCOUPLING_EXPORT DataArrayInt32 *findIdsNotInRange(double vmin, double vmax) const;
    MEDCOUPLING_EXPORT static DataArrayDouble *Aggregate(const DataArrayDouble *a1, const DataArrayDouble *a2);
    MEDCOUPLING_EXPORT static DataArrayDouble *Aggregate(const std::vector<const DataArrayDouble *>& arr);
    MEDCOUPLING_EXPORT static DataArrayDouble *Dot(const DataArrayDouble *a1, const DataArrayDouble *a2);
    MEDCOUPLING_EXPORT static DataArrayDouble *CrossProduct(const DataArrayDouble *a1, const DataArrayDouble *a2);
    MEDCOUPLING_EXPORT static DataArrayDouble *Max(const DataArrayDouble *a1, const DataArrayDouble *a2);
    MEDCOUPLING_EXPORT static DataArrayDouble *Min(const DataArrayDouble *a1, const DataArrayDouble *a2);
    MEDCOUPLING_EXPORT static DataArrayDouble *Pow(const DataArrayDouble *a1, const DataArrayDouble *a2);
    MEDCOUPLING_EXPORT void powEqual(const DataArrayDouble *other);
    MEDCOUPLING_EXPORT std::vector<bool> toVectorOfBool(double eps) const;
    MEDCOUPLING_EXPORT static void Rotate2DAlg(const double *center, double angle, int nbNodes, const double *coordsIn, double *coordsOut);
    MEDCOUPLING_EXPORT static void Rotate3DAlg(const double *center, const double *vect, double angle, int nbNodes, const double *coordsIn, double *coordsOut);
    MEDCOUPLING_EXPORT static void Symmetry3DPlane(const double point[3], const double normalVector[3], int nbNodes, const double *coordsIn, double *coordsOut);
    MEDCOUPLING_EXPORT static void GiveBaseForPlane(const double normalVector[3], double baseOfPlane[9]);
  public:
    MEDCOUPLING_EXPORT void getTinySerializationIntInformation(std::vector<int>& tinyInfo) const;
    MEDCOUPLING_EXPORT void getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const;
    MEDCOUPLING_EXPORT bool resizeForUnserialization(const std::vector<int>& tinyInfoI);
    MEDCOUPLING_EXPORT void finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<std::string>& tinyInfoS);
  public:
    template<int SPACEDIM>
    void findCommonTuplesAlg(const double *bbox, int nbNodes, int limitNodeId, double prec, DataArrayInt32 *c, DataArrayInt32 *cI) const;
    template<int SPACEDIM>
    static void FindClosestTupleIdAlg(const BBTreePts<SPACEDIM,int>& myTree, double dist, const double *pos, int nbOfTuples, const double *thisPt, int thisNbOfTuples, int *res);
    template<int SPACEDIM>
    static void FindTupleIdsNearTuplesAlg(const BBTreePts<SPACEDIM,int>& myTree, const double *pos, int nbOfTuples, double eps,
                                          DataArrayInt32 *c, DataArrayInt32 *cI);
  private:
    ~DataArrayDouble() { }
    DataArrayDouble() { }
  };

  template<class T>
  class DataArrayDiscrete : public DataArrayTemplateClassic<T>
  {
  public:
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
  public:
    MEDCOUPLING_EXPORT static DataArrayInt32 *New();
    MEDCOUPLING_EXPORT int intValue() const;
    MEDCOUPLING_EXPORT int getHashCode() const;
    MEDCOUPLING_EXPORT DataArrayInt32 *deepCopy() const;//ok
    MEDCOUPLING_EXPORT DataArrayInt32 *buildNewEmptyInstance() const { return DataArrayInt32::New(); }//ok
    MEDCOUPLING_EXPORT std::string repr() const;
    MEDCOUPLING_EXPORT std::string reprZip() const;
    MEDCOUPLING_EXPORT void writeVTK(std::ostream& ofs, int indent, const std::string& type, const std::string& nameInFile, DataArrayByte *byteArr) const;
    MEDCOUPLING_EXPORT void reprCppStream(const std::string& varName, std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprQuickOverview(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprQuickOverviewData(std::ostream& stream, std::size_t maxNbOfByteInRepr) const;
    MEDCOUPLING_EXPORT void transformWithIndArr(const int *indArrBg, const int *indArrEnd);
    MEDCOUPLING_EXPORT void transformWithIndArr(const MapKeyVal<int>& m);
    MEDCOUPLING_EXPORT DataArrayInt32 *transformWithIndArrR(const int *indArrBg, const int *indArrEnd) const;
    MEDCOUPLING_EXPORT void splitByValueRange(const int *arrBg, const int *arrEnd,
                                              DataArrayInt32 *& castArr, DataArrayInt32 *& rankInsideCast, DataArrayInt32 *& castsPresent) const;
    MEDCOUPLING_EXPORT bool isRange(int& strt, int& sttoopp, int& stteepp) const;
    MEDCOUPLING_EXPORT DataArrayInt32 *invertArrayO2N2N2O(int newNbOfElem) const;
    MEDCOUPLING_EXPORT DataArrayInt32 *invertArrayN2O2O2N(int oldNbOfElem) const;
    MEDCOUPLING_EXPORT MCAuto< MapKeyVal<int> > invertArrayN2O2O2NOptimized() const;
    MEDCOUPLING_EXPORT DataArrayInt32 *invertArrayO2N2N2OBis(int newNbOfElem) const;
    MEDCOUPLING_EXPORT DataArrayInt32 *selectByTupleId(const int *new2OldBg, const int *new2OldEnd) const { return DataArrayTemplate<int>::mySelectByTupleId(new2OldBg,new2OldEnd); }
    MEDCOUPLING_EXPORT DataArrayInt32 *selectByTupleId(const DataArrayInt32& di) const { return DataArrayTemplate<int>::mySelectByTupleId(di); }
    MEDCOUPLING_EXPORT DataArrayInt32 *selectByTupleIdSafe(const int *new2OldBg, const int *new2OldEnd) const { return DataArrayTemplate<int>::mySelectByTupleIdSafe(new2OldBg,new2OldEnd); }
    MEDCOUPLING_EXPORT DataArrayInt32 *keepSelectedComponents(const std::vector<int>& compoIds) const { return DataArrayTemplate<int>::myKeepSelectedComponents(compoIds); }
    MEDCOUPLING_EXPORT DataArrayInt32 *selectByTupleIdSafeSlice(int bg, int end2, int step) const { return DataArrayTemplate<int>::mySelectByTupleIdSafeSlice(bg,end2,step); }
    MEDCOUPLING_EXPORT DataArrayInt32 *selectByTupleRanges(const std::vector<std::pair<int,int> >& ranges) const { return DataArrayTemplate<int>::mySelectByTupleRanges(ranges); }
    MEDCOUPLING_EXPORT DataArrayInt32 *checkAndPreparePermutation() const;
    MEDCOUPLING_EXPORT static DataArrayInt32 *FindPermutationFromFirstToSecond(const DataArrayInt32 *ids1, const DataArrayInt32 *ids2);
    MEDCOUPLING_EXPORT void changeSurjectiveFormat(int targetNb, DataArrayInt32 *&arr, DataArrayInt32 *&arrI) const;
    MEDCOUPLING_EXPORT static DataArrayInt32 *ConvertIndexArrayToO2N(int nbOfOldTuples, const int *arr, const int *arrIBg, const int *arrIEnd, int &newNbOfTuples);
    MEDCOUPLING_EXPORT DataArrayInt32 *buildPermArrPerLevel() const;
    MEDCOUPLING_EXPORT bool isIota(int sizeExpected) const;
    MEDCOUPLING_EXPORT bool isUniform(int val) const;
    MEDCOUPLING_EXPORT int checkUniformAndGuess() const;
    MEDCOUPLING_EXPORT bool hasUniqueValues() const;
    MEDCOUPLING_EXPORT void setSelectedComponents(const DataArrayInt32 *a, const std::vector<int>& compoIds);
    MEDCOUPLING_EXPORT DataArrayInt32Iterator *iterator();
    MEDCOUPLING_EXPORT DataArrayInt32 *findIdsEqual(int val) const;
    MEDCOUPLING_EXPORT DataArrayInt32 *findIdsNotEqual(int val) const;
    MEDCOUPLING_EXPORT DataArrayInt32 *findIdsEqualList(const int *valsBg, const int *valsEnd) const;
    MEDCOUPLING_EXPORT DataArrayInt32 *findIdsNotEqualList(const int *valsBg, const int *valsEnd) const;
    MEDCOUPLING_EXPORT DataArrayInt32 *findIdsEqualTuple(const int *tupleBg, const int *tupleEnd) const;
    MEDCOUPLING_EXPORT int changeValue(int oldValue, int newValue);
    MEDCOUPLING_EXPORT int findIdFirstEqualTuple(const std::vector<int>& tupl) const;
    MEDCOUPLING_EXPORT int findIdFirstEqual(int value) const;
    MEDCOUPLING_EXPORT int findIdFirstEqual(const std::vector<int>& vals) const;
    MEDCOUPLING_EXPORT int findIdSequence(const std::vector<int>& vals) const;
    MEDCOUPLING_EXPORT bool presenceOfTuple(const std::vector<int>& tupl) const;
    MEDCOUPLING_EXPORT bool presenceOfValue(int value) const;
    MEDCOUPLING_EXPORT bool presenceOfValue(const std::vector<int>& vals) const;
    MEDCOUPLING_EXPORT int count(int value) const;
    MEDCOUPLING_EXPORT void accumulate(int *res) const;
    MEDCOUPLING_EXPORT int accumulate(int compId) const;
    MEDCOUPLING_EXPORT DataArrayInt32 *accumulatePerChunck(const int *bgOfIndex, const int *endOfIndex) const;
    MEDCOUPLING_EXPORT void getMinMaxValues(int& minValue, int& maxValue) const;
    MEDCOUPLING_EXPORT void applyInv(int numerator);
    MEDCOUPLING_EXPORT void applyDivideBy(int val);
    MEDCOUPLING_EXPORT void applyModulus(int val);
    MEDCOUPLING_EXPORT void applyRModulus(int val);
    MEDCOUPLING_EXPORT void applyPow(int val);
    MEDCOUPLING_EXPORT void applyRPow(int val);
    MEDCOUPLING_EXPORT DataArrayInt32 *findIdsInRange(int vmin, int vmax) const;
    MEDCOUPLING_EXPORT DataArrayInt32 *findIdsNotInRange(int vmin, int vmax) const;
    MEDCOUPLING_EXPORT bool checkAllIdsInRange(int vmin, int vmax) const;
    MEDCOUPLING_EXPORT static DataArrayInt32 *Aggregate(const DataArrayInt32 *a1, const DataArrayInt32 *a2, int offsetA2);
    MEDCOUPLING_EXPORT static DataArrayInt32 *Aggregate(const std::vector<const DataArrayInt32 *>& arr);
    MEDCOUPLING_EXPORT static DataArrayInt32 *AggregateIndexes(const std::vector<const DataArrayInt32 *>& arrs);
    MEDCOUPLING_EXPORT static DataArrayInt32 *MakePartition(const std::vector<const DataArrayInt32 *>& groups, int newNb, std::vector< std::vector<int> >& fidsOfGroups);
    MEDCOUPLING_EXPORT static DataArrayInt32 *BuildUnion(const std::vector<const DataArrayInt32 *>& arr);
    MEDCOUPLING_EXPORT static DataArrayInt32 *BuildIntersection(const std::vector<const DataArrayInt32 *>& arr);
    MEDCOUPLING_EXPORT static DataArrayInt32 *BuildListOfSwitchedOn(const std::vector<bool>& v);
    MEDCOUPLING_EXPORT static DataArrayInt32 *BuildListOfSwitchedOff(const std::vector<bool>& v);
    MEDCOUPLING_EXPORT static void PutIntoToSkylineFrmt(const std::vector< std::vector<int> >& v, DataArrayInt32 *& data, DataArrayInt32 *& dataIndex);
    MEDCOUPLING_EXPORT DataArrayInt32 *buildComplement(int nbOfElement) const;
    MEDCOUPLING_EXPORT DataArrayInt32 *buildSubstraction(const DataArrayInt32 *other) const;
    MEDCOUPLING_EXPORT DataArrayInt32 *buildSubstractionOptimized(const DataArrayInt32 *other) const;
    MEDCOUPLING_EXPORT DataArrayInt32 *buildUnion(const DataArrayInt32 *other) const;
    MEDCOUPLING_EXPORT DataArrayInt32 *buildIntersection(const DataArrayInt32 *other) const;
    MEDCOUPLING_EXPORT DataArrayInt32 *buildUnique() const;
    MEDCOUPLING_EXPORT DataArrayInt32 *buildUniqueNotSorted() const;
    MEDCOUPLING_EXPORT DataArrayInt32 *deltaShiftIndex() const;
    MEDCOUPLING_EXPORT void computeOffsets();
    MEDCOUPLING_EXPORT void computeOffsetsFull();
    MEDCOUPLING_EXPORT void findIdsRangesInListOfIds(const DataArrayInt32 *listOfIds, DataArrayInt32 *& rangeIdsFetched, DataArrayInt32 *& idsInInputListThatFetch) const;
    MEDCOUPLING_EXPORT DataArrayInt32 *buildExplicitArrByRanges(const DataArrayInt32 *offsets) const;
    MEDCOUPLING_EXPORT DataArrayInt32 *buildExplicitArrOfSliceOnScaledArr(int begin, int stop, int step) const;
    MEDCOUPLING_EXPORT DataArrayInt32 *findRangeIdForEachTuple(const DataArrayInt32 *ranges) const;
    MEDCOUPLING_EXPORT DataArrayInt32 *findIdInRangeForEachTuple(const DataArrayInt32 *ranges) const;
    MEDCOUPLING_EXPORT void sortEachPairToMakeALinkedList();
    MEDCOUPLING_EXPORT MCAuto<DataArrayInt32> fromLinkedListOfPairToList() const;
    MEDCOUPLING_EXPORT DataArrayInt32 *getDifferentValues() const;
    MEDCOUPLING_EXPORT std::vector<DataArrayInt32 *> partitionByDifferentValues(std::vector<int>& differentIds) const;
    MEDCOUPLING_EXPORT std::vector< std::pair<int,int> > splitInBalancedSlices(int nbOfSlices) const;
    MEDCOUPLING_EXPORT static DataArrayInt32 *Modulus(const DataArrayInt32 *a1, const DataArrayInt32 *a2);
    MEDCOUPLING_EXPORT void modulusEqual(const DataArrayInt32 *other);
    MEDCOUPLING_EXPORT static DataArrayInt32 *Pow(const DataArrayInt32 *a1, const DataArrayInt32 *a2);
    MEDCOUPLING_EXPORT void powEqual(const DataArrayInt32 *other);
    MEDCOUPLING_EXPORT MemArray<int>& accessToMemArray() { return _mem; }
    MEDCOUPLING_EXPORT const MemArray<int>& accessToMemArray() const { return _mem; }
  public:
    MEDCOUPLING_EXPORT static int *CheckAndPreparePermutation(const int *start, const int *end);
    MEDCOUPLING_EXPORT static DataArrayInt32 *Range(int begin, int end, int step);
  public:
    MEDCOUPLING_EXPORT void getTinySerializationIntInformation(std::vector<int>& tinyInfo) const;
    MEDCOUPLING_EXPORT void getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const;
    MEDCOUPLING_EXPORT bool resizeForUnserialization(const std::vector<int>& tinyInfoI);
    MEDCOUPLING_EXPORT void finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<std::string>& tinyInfoS);
  private:
    ~DataArrayInt32() { }
    DataArrayInt32() { }
  };

  class DataArrayInt64 : public DataArrayDiscrete<Int64>
  {
  };
  
  template<class T>
  template<class OP>
  MCAuto<DataArrayInt> DataArrayTemplateClassic<T>::findIdsAdv(const OP& op) const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::findIdsAdv : this must have exactly one component !");
    const T *cptr(this->begin());
    MCAuto<DataArrayInt> ret(DataArrayInt::New()); ret->alloc(0,1);
    int nbOfTuples(this->getNumberOfTuples());
    for(int i=0;i<nbOfTuples;i++,cptr++)
      if(op(*cptr))
        ret->pushBackSilent(i);
    return ret;
  }

  class DataArrayChar : public DataArrayTemplate<char>
  {
  public:
    MEDCOUPLING_EXPORT virtual DataArrayChar *buildEmptySpecializedDAChar() const = 0;
    MEDCOUPLING_EXPORT int getHashCode() const;
    MEDCOUPLING_EXPORT bool isEqual(const DataArrayChar& other) const;
    MEDCOUPLING_EXPORT virtual bool isEqualIfNotWhy(const DataArrayChar& other, std::string& reason) const;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const DataArrayChar& other) const;
    MEDCOUPLING_EXPORT std::string repr() const;
    MEDCOUPLING_EXPORT std::string reprZip() const;
    MEDCOUPLING_EXPORT DataArrayInt *convertToIntArr() const;
    MEDCOUPLING_EXPORT DataArrayChar *selectByTupleId(const int *new2OldBg, const int *new2OldEnd) const { return DataArrayTemplate<char>::mySelectByTupleId(new2OldBg,new2OldEnd); }
    MEDCOUPLING_EXPORT DataArrayChar *selectByTupleId(const DataArrayInt& di) const { return DataArrayTemplate<char>::mySelectByTupleId(di); }
    MEDCOUPLING_EXPORT DataArrayChar *selectByTupleIdSafe(const int *new2OldBg, const int *new2OldEnd) const { return DataArrayTemplate<char>::mySelectByTupleIdSafe(new2OldBg,new2OldEnd); }
    MEDCOUPLING_EXPORT DataArrayChar *keepSelectedComponents(const std::vector<int>& compoIds) const { return DataArrayTemplate<char>::myKeepSelectedComponents(compoIds); }
    MEDCOUPLING_EXPORT DataArrayChar *selectByTupleIdSafeSlice(int bg, int end2, int step) const { return DataArrayTemplate<char>::mySelectByTupleIdSafeSlice(bg,end2,step); }
    MEDCOUPLING_EXPORT bool isUniform(char val) const;
    MEDCOUPLING_EXPORT void meldWith(const DataArrayChar *other);
    MEDCOUPLING_EXPORT DataArray *selectByTupleRanges(const std::vector<std::pair<int,int> >& ranges) const { return DataArrayTemplate<char>::mySelectByTupleRanges(ranges); }
    MEDCOUPLING_EXPORT DataArrayInt *findIdsEqual(char val) const;
    MEDCOUPLING_EXPORT DataArrayInt *findIdsNotEqual(char val) const;
    MEDCOUPLING_EXPORT int findIdSequence(const std::vector<char>& vals) const;
    MEDCOUPLING_EXPORT int findIdFirstEqualTuple(const std::vector<char>& tupl) const;
    MEDCOUPLING_EXPORT int findIdFirstEqual(char value) const;
    MEDCOUPLING_EXPORT int findIdFirstEqual(const std::vector<char>& vals) const;
    MEDCOUPLING_EXPORT bool presenceOfTuple(const std::vector<char>& tupl) const;
    MEDCOUPLING_EXPORT bool presenceOfValue(char value) const;
    MEDCOUPLING_EXPORT bool presenceOfValue(const std::vector<char>& vals) const;
    MEDCOUPLING_EXPORT DataArrayInt *findIdsInRange(char vmin, char vmax) const;
    MEDCOUPLING_EXPORT static DataArrayChar *Aggregate(const DataArrayChar *a1, const DataArrayChar *a2);
    MEDCOUPLING_EXPORT static DataArrayChar *Aggregate(const std::vector<const DataArrayChar *>& arr);
    MEDCOUPLING_EXPORT static DataArrayChar *Meld(const DataArrayChar *a1, const DataArrayChar *a2);
    MEDCOUPLING_EXPORT static DataArrayChar *Meld(const std::vector<const DataArrayChar *>& arr);
    MEDCOUPLING_EXPORT MemArray<char>& accessToMemArray() { return _mem; }
    MEDCOUPLING_EXPORT const MemArray<char>& accessToMemArray() const { return _mem; }
  public:
    //MEDCOUPLING_EXPORT void getTinySerializationIntInformation(std::vector<int>& tinyInfo) const;
    //MEDCOUPLING_EXPORT void getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const;
    //MEDCOUPLING_EXPORT bool resizeForUnserialization(const std::vector<int>& tinyInfoI);
    //MEDCOUPLING_EXPORT void finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<std::string>& tinyInfoS);
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
    int _tuple_id;
    int _nb_comp;
    int _nb_tuple;
  };

  template<class T>
  class DataArrayTuple
  {
  public:
    MEDCOUPLING_EXPORT DataArrayTuple(T *pt, int nbOfComp);
    MEDCOUPLING_EXPORT std::string repr() const;
    MEDCOUPLING_EXPORT int getNumberOfCompo() const { return _nb_of_compo; }
    MEDCOUPLING_EXPORT const T *getConstPointer() const { return  _pt; }
    MEDCOUPLING_EXPORT T *getPointer() { return _pt; }
    MEDCOUPLING_EXPORT typename Traits<T>::ArrayType *buildDA(int nbOfTuples, int nbOfCompo) const;
  protected:
    T zeValue() const;
  protected:
    T *_pt;
    int _nb_of_compo;
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
    MEDCOUPLING_EXPORT DataArrayDoubleTuple(double *pt, int nbOfComp);
    MEDCOUPLING_EXPORT std::string repr() const;
    MEDCOUPLING_EXPORT double doubleValue() const;
    MEDCOUPLING_EXPORT DataArrayDouble *buildDADouble(int nbOfTuples, int nbOfCompo) const;
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
    MEDCOUPLING_EXPORT DataArrayFloatTuple(float *pt, int nbOfComp);
    MEDCOUPLING_EXPORT std::string repr() const;
    MEDCOUPLING_EXPORT float floatValue() const;
    MEDCOUPLING_EXPORT DataArrayFloat *buildDAFloat(int nbOfTuples, int nbOfCompo) const;
  };
  
  class DataArrayIntIterator : public DataArrayIterator<int>
  {
  public:
    MEDCOUPLING_EXPORT DataArrayIntIterator(DataArrayInt *da);
    MEDCOUPLING_EXPORT ~DataArrayIntIterator() { }
  };

  class DataArrayInt32Tuple : public DataArrayTuple<Int32>
  {
  public:
    MEDCOUPLING_EXPORT DataArrayInt32Tuple(int *pt, int nbOfComp);
    MEDCOUPLING_EXPORT std::string repr() const;
    MEDCOUPLING_EXPORT int intValue() const;
    MEDCOUPLING_EXPORT DataArrayInt32 *buildDAInt(int nbOfTuples, int nbOfCompo) const;
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
    int _tuple_id;
    int _nb_comp;
    int _nb_tuple;
  };

  class DataArrayAsciiCharTuple
  {
  public:
    MEDCOUPLING_EXPORT DataArrayAsciiCharTuple(char *pt, int nbOfComp);
    MEDCOUPLING_EXPORT std::string repr() const;
    MEDCOUPLING_EXPORT int getNumberOfCompo() const { return _nb_of_compo; }
    MEDCOUPLING_EXPORT const char *getConstPointer() const { return  _pt; }
    MEDCOUPLING_EXPORT char *getPointer() { return _pt; }
    MEDCOUPLING_EXPORT char asciiCharValue() const;
    MEDCOUPLING_EXPORT DataArrayAsciiChar *buildDAAsciiChar(int nbOfTuples, int nbOfCompo) const;
  private:
    char *_pt;
    int _nb_of_compo;
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
    int _tuple_id;
    int _nb_comp;
    int _nb_tuple;
  };

  class DataArrayByteTuple
  {
  public:
    MEDCOUPLING_EXPORT DataArrayByteTuple(char *pt, int nbOfComp);
    MEDCOUPLING_EXPORT std::string repr() const;
    MEDCOUPLING_EXPORT int getNumberOfCompo() const { return _nb_of_compo; }
    MEDCOUPLING_EXPORT const char *getConstPointer() const { return  _pt; }
    MEDCOUPLING_EXPORT char *getPointer() { return _pt; }
    MEDCOUPLING_EXPORT char byteValue() const;
    MEDCOUPLING_EXPORT DataArrayByte *buildDAByte(int nbOfTuples, int nbOfCompo) const;
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

  template<class T>
  template<class InputIterator>
  void DataArrayTemplate<T>::insertAtTheEnd(InputIterator first, InputIterator last)
  {
    int nbCompo(this->getNumberOfComponents());
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
