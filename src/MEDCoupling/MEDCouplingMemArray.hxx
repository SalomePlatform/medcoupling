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
// Author : Anthony Geay (CEA/DEN)

#ifndef __PARAMEDMEM_MEDCOUPLINGMEMARRAY_HXX__
#define __PARAMEDMEM_MEDCOUPLINGMEMARRAY_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingTimeLabel.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "InterpKernelException.hxx"
#include "BBTreePts.txx"

#include <string>
#include <vector>
#include <iterator>

namespace ParaMEDMEM
{
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

  class DataArray : public RefCountObject, public TimeLabel
  {
  public:
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildren() const;
    MEDCOUPLING_EXPORT void setName(const char *name);
    MEDCOUPLING_EXPORT void copyStringInfoFrom(const DataArray& other);
    MEDCOUPLING_EXPORT void copyPartOfStringInfoFrom(const DataArray& other, const std::vector<int>& compoIds);
    MEDCOUPLING_EXPORT void copyPartOfStringInfoFrom2(const std::vector<int>& compoIds, const DataArray& other);
    MEDCOUPLING_EXPORT bool areInfoEqualsIfNotWhy(const DataArray& other, std::string& reason) const;
    MEDCOUPLING_EXPORT bool areInfoEquals(const DataArray& other) const;
    MEDCOUPLING_EXPORT std::string cppRepr(const char *varName) const;
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
    MEDCOUPLING_EXPORT void setInfoOnComponent(int i, const char *info);
    MEDCOUPLING_EXPORT int getNumberOfComponents() const { return (int)_info_on_compo.size(); }
    MEDCOUPLING_EXPORT void setPartOfValuesBase3(const DataArray *aBase, const int *bgTuples, const int *endTuples, int bgComp, int endComp, int stepComp, bool strictCompoCompare=true);
    MEDCOUPLING_EXPORT virtual DataArray *deepCpy() const = 0;
    MEDCOUPLING_EXPORT virtual bool isAllocated() const = 0;
    MEDCOUPLING_EXPORT virtual void checkAllocated() const = 0;
    MEDCOUPLING_EXPORT virtual void desallocate() = 0;
    MEDCOUPLING_EXPORT virtual int getNumberOfTuples() const = 0;
    MEDCOUPLING_EXPORT virtual std::size_t getNbOfElems() const = 0;
    MEDCOUPLING_EXPORT virtual std::size_t getNbOfElemAllocated() const = 0;
    MEDCOUPLING_EXPORT virtual void alloc(int nbOfTuple, int nbOfCompo=1) = 0;
    MEDCOUPLING_EXPORT virtual void reAlloc(int newNbOfTuple) = 0;
    MEDCOUPLING_EXPORT virtual void renumberInPlace(const int *old2New) = 0;
    MEDCOUPLING_EXPORT virtual void renumberInPlaceR(const int *new2Old) = 0;
    MEDCOUPLING_EXPORT virtual void setContigPartOfSelectedValues(int tupleIdStart, const DataArray *aBase, const DataArrayInt *tuplesSelec) = 0;
    MEDCOUPLING_EXPORT virtual void setContigPartOfSelectedValues2(int tupleIdStart, const DataArray *aBase, int bg, int end2, int step) = 0;
    MEDCOUPLING_EXPORT virtual DataArray *selectByTupleRanges(const std::vector<std::pair<int,int> >& ranges) const = 0;
    MEDCOUPLING_EXPORT virtual DataArray *keepSelectedComponents(const std::vector<int>& compoIds) const = 0;
    MEDCOUPLING_EXPORT virtual DataArray *selectByTupleId(const int *new2OldBg, const int *new2OldEnd) const = 0;
    MEDCOUPLING_EXPORT virtual DataArray *selectByTupleIdSafe(const int *new2OldBg, const int *new2OldEnd) const = 0;
    MEDCOUPLING_EXPORT virtual DataArray *selectByTupleId2(int bg, int end2, int step) const = 0;
    MEDCOUPLING_EXPORT virtual void rearrange(int newNbOfCompo) = 0;
    MEDCOUPLING_EXPORT void checkNbOfTuples(int nbOfTuples, const char *msg) const;
    MEDCOUPLING_EXPORT void checkNbOfComps(int nbOfCompo, const char *msg) const;
    MEDCOUPLING_EXPORT void checkNbOfTuplesAndComp(const DataArray& other, const char *msg) const;
    MEDCOUPLING_EXPORT void checkNbOfTuplesAndComp(int nbOfTuples, int nbOfCompo, const char *msg) const;
    MEDCOUPLING_EXPORT void checkNbOfElems(std::size_t nbOfElems, const char *msg) const;
    MEDCOUPLING_EXPORT static void GetSlice(int start, int stop, int step, int sliceId, int nbOfSlices, int& startSlice, int& stopSlice);
    MEDCOUPLING_EXPORT static int GetNumberOfItemGivenBES(int begin, int end, int step, const char *msg);
    MEDCOUPLING_EXPORT static int GetNumberOfItemGivenBESRelative(int begin, int end, int step, const char *msg);
    MEDCOUPLING_EXPORT static int GetPosOfItemGivenBESRelativeNoThrow(int value, int begin, int end, int step);
    MEDCOUPLING_EXPORT static std::string GetVarNameFromInfo(const std::string& info);
    MEDCOUPLING_EXPORT static std::string GetUnitFromInfo(const std::string& info);
    MEDCOUPLING_EXPORT static DataArray *Aggregate(const std::vector<const DataArray *>& arrs);
    MEDCOUPLING_EXPORT virtual void reprStream(std::ostream& stream) const = 0;
    MEDCOUPLING_EXPORT virtual void reprZipStream(std::ostream& stream) const = 0;
    MEDCOUPLING_EXPORT virtual void reprWithoutNameStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT virtual void reprZipWithoutNameStream(std::ostream& stream) const = 0;
    MEDCOUPLING_EXPORT virtual void reprCppStream(const char *varName, std::ostream& stream) const = 0;
    MEDCOUPLING_EXPORT virtual void reprQuickOverview(std::ostream& stream) const = 0;
    MEDCOUPLING_EXPORT virtual void reprQuickOverviewData(std::ostream& stream, std::size_t maxNbOfByteInRepr) const = 0;
  protected:
    DataArray() { }
    ~DataArray() { }
  protected:
    static void CheckValueInRange(int ref, int value, const char *msg);
    static void CheckValueInRangeEx(int value, int start, int end, const char *msg);
    static void CheckClosingParInRange(int ref, int value, const char *msg);
  protected:
    std::string _name;
    std::vector<std::string> _info_on_compo;
  };
}

#include "MEDCouplingMemArray.txx"

namespace ParaMEDMEM
{
  class DataArrayInt;
  class DataArrayDoubleIterator;
  class DataArrayDouble : public DataArray
  {
  public:
    MEDCOUPLING_EXPORT static DataArrayDouble *New();
    MEDCOUPLING_EXPORT bool isAllocated() const;
    MEDCOUPLING_EXPORT void checkAllocated() const;
    MEDCOUPLING_EXPORT void desallocate();
    MEDCOUPLING_EXPORT int getNumberOfTuples() const { return _info_on_compo.empty()?0:_mem.getNbOfElem()/getNumberOfComponents(); }
    MEDCOUPLING_EXPORT std::size_t getNbOfElems() const { return _mem.getNbOfElem(); }
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT double doubleValue() const;
    MEDCOUPLING_EXPORT bool empty() const;
    MEDCOUPLING_EXPORT DataArrayDouble *deepCpy() const;
    MEDCOUPLING_EXPORT DataArrayDouble *performCpy(bool deepCpy) const;
    MEDCOUPLING_EXPORT void cpyFrom(const DataArrayDouble& other);
    MEDCOUPLING_EXPORT void reserve(std::size_t nbOfElems);
    MEDCOUPLING_EXPORT void pushBackSilent(double val);
    MEDCOUPLING_EXPORT void pushBackValsSilent(const double *valsBg, const double *valsEnd);
    MEDCOUPLING_EXPORT double popBackSilent();
    MEDCOUPLING_EXPORT void pack() const;
    MEDCOUPLING_EXPORT std::size_t getNbOfElemAllocated() const { return _mem.getNbOfElemAllocated(); }
    MEDCOUPLING_EXPORT void alloc(int nbOfTuple, int nbOfCompo=1);
    MEDCOUPLING_EXPORT void allocIfNecessary(int nbOfTuple, int nbOfCompo);
    MEDCOUPLING_EXPORT void fillWithZero();
    MEDCOUPLING_EXPORT void fillWithValue(double val);
    MEDCOUPLING_EXPORT void iota(double init=0.);
    MEDCOUPLING_EXPORT bool isUniform(double val, double eps) const;
    MEDCOUPLING_EXPORT void sort(bool asc=true);
    MEDCOUPLING_EXPORT void reverse();
    MEDCOUPLING_EXPORT void checkMonotonic(bool increasing, double eps) const;
    MEDCOUPLING_EXPORT bool isMonotonic(bool increasing, double eps) const;
    MEDCOUPLING_EXPORT std::string repr() const;
    MEDCOUPLING_EXPORT std::string reprZip() const;
    MEDCOUPLING_EXPORT void writeVTK(std::ostream& ofs, int indent, const char *nameInFile, DataArrayByte *byteArr) const;
    MEDCOUPLING_EXPORT void reprStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprZipStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprWithoutNameStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprZipWithoutNameStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprCppStream(const char *varName, std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprQuickOverview(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprQuickOverviewData(std::ostream& stream, std::size_t maxNbOfByteInRepr) const;
    MEDCOUPLING_EXPORT bool isEqual(const DataArrayDouble& other, double prec) const;
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const DataArrayDouble& other, double prec, std::string& reason) const;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const DataArrayDouble& other, double prec) const;
    MEDCOUPLING_EXPORT void reAlloc(int nbOfTuples);
    MEDCOUPLING_EXPORT DataArrayInt *convertToIntArr() const;
    MEDCOUPLING_EXPORT DataArrayDouble *fromNoInterlace() const;
    MEDCOUPLING_EXPORT DataArrayDouble *toNoInterlace() const;
    MEDCOUPLING_EXPORT void renumberInPlace(const int *old2New);
    MEDCOUPLING_EXPORT void renumberInPlaceR(const int *new2Old);
    MEDCOUPLING_EXPORT DataArrayDouble *renumber(const int *old2New) const;
    MEDCOUPLING_EXPORT DataArrayDouble *renumberR(const int *new2Old) const;
    MEDCOUPLING_EXPORT DataArrayDouble *renumberAndReduce(const int *old2New, int newNbOfTuple) const;
    MEDCOUPLING_EXPORT DataArrayDouble *selectByTupleId(const int *new2OldBg, const int *new2OldEnd) const;
    MEDCOUPLING_EXPORT DataArrayDouble *selectByTupleIdSafe(const int *new2OldBg, const int *new2OldEnd) const;
    MEDCOUPLING_EXPORT DataArrayDouble *selectByTupleId2(int bg, int end2, int step) const;
    MEDCOUPLING_EXPORT DataArray *selectByTupleRanges(const std::vector<std::pair<int,int> >& ranges) const;
    MEDCOUPLING_EXPORT DataArrayDouble *substr(int tupleIdBg, int tupleIdEnd=-1) const;
    MEDCOUPLING_EXPORT void rearrange(int newNbOfCompo);
    MEDCOUPLING_EXPORT void transpose();
    MEDCOUPLING_EXPORT DataArrayDouble *changeNbOfComponents(int newNbOfComp, double dftValue) const;
    MEDCOUPLING_EXPORT DataArray *keepSelectedComponents(const std::vector<int>& compoIds) const;
    MEDCOUPLING_EXPORT void meldWith(const DataArrayDouble *other);
    MEDCOUPLING_EXPORT bool areIncludedInMe(const DataArrayDouble *other, double prec, DataArrayInt *&tupleIds) const;
    MEDCOUPLING_EXPORT void findCommonTuples(double prec, int limitTupleId, DataArrayInt *&comm, DataArrayInt *&commIndex) const;
    MEDCOUPLING_EXPORT double minimalDistanceTo(const DataArrayDouble *other, int& thisTupleId, int& otherTupleId) const;
    MEDCOUPLING_EXPORT DataArrayDouble *duplicateEachTupleNTimes(int nbTimes) const;
    MEDCOUPLING_EXPORT DataArrayDouble *getDifferentValues(double prec, int limitTupleId=-1) const;
    MEDCOUPLING_EXPORT DataArrayInt *findClosestTupleId(const DataArrayDouble *other) const;
    MEDCOUPLING_EXPORT DataArrayInt *computeNbOfInteractionsWith(const DataArrayDouble *otherBBoxFrmt, double eps) const;
    MEDCOUPLING_EXPORT void setSelectedComponents(const DataArrayDouble *a, const std::vector<int>& compoIds);
    MEDCOUPLING_EXPORT void setPartOfValues1(const DataArrayDouble *a, int bgTuples, int endTuples, int stepTuples, int bgComp, int endComp, int stepComp, bool strictCompoCompare=true);
    MEDCOUPLING_EXPORT void setPartOfValuesSimple1(double a, int bgTuples, int endTuples, int stepTuples, int bgComp, int endComp, int stepComp);
    MEDCOUPLING_EXPORT void setPartOfValues2(const DataArrayDouble *a, const int *bgTuples, const int *endTuples, const int *bgComp, const int *endComp, bool strictCompoCompare=true);
    MEDCOUPLING_EXPORT void setPartOfValuesSimple2(double a, const int *bgTuples, const int *endTuples, const int *bgComp, const int *endComp);
    MEDCOUPLING_EXPORT void setPartOfValues3(const DataArrayDouble *a, const int *bgTuples, const int *endTuples, int bgComp, int endComp, int stepComp, bool strictCompoCompare=true);
    MEDCOUPLING_EXPORT void setPartOfValuesSimple3(double a, const int *bgTuples, const int *endTuples, int bgComp, int endComp, int stepComp);
    MEDCOUPLING_EXPORT void setPartOfValues4(const DataArrayDouble *a, int bgTuples, int endTuples, int stepTuples, const int *bgComp, const int *endComp, bool strictCompoCompare=true);
    MEDCOUPLING_EXPORT void setPartOfValuesSimple4(double a, int bgTuples, int endTuples, int stepTuples, const int *bgComp, const int *endComp);
    MEDCOUPLING_EXPORT void setPartOfValuesAdv(const DataArrayDouble *a, const DataArrayInt *tuplesSelec);
    MEDCOUPLING_EXPORT void setContigPartOfSelectedValues(int tupleIdStart, const DataArray *aBase, const DataArrayInt *tuplesSelec);
    MEDCOUPLING_EXPORT void setContigPartOfSelectedValues2(int tupleIdStart, const DataArray *aBase, int bg, int end2, int step);
    MEDCOUPLING_EXPORT void getTuple(int tupleId, double *res) const { std::copy(_mem.getConstPointerLoc(tupleId*_info_on_compo.size()),_mem.getConstPointerLoc((tupleId+1)*_info_on_compo.size()),res); }
    MEDCOUPLING_EXPORT double getIJ(int tupleId, int compoId) const { return _mem[tupleId*_info_on_compo.size()+compoId]; }
    MEDCOUPLING_EXPORT double front() const;
    MEDCOUPLING_EXPORT double back() const;
    MEDCOUPLING_EXPORT double getIJSafe(int tupleId, int compoId) const;
    MEDCOUPLING_EXPORT void setIJ(int tupleId, int compoId, double newVal) { _mem[tupleId*_info_on_compo.size()+compoId]=newVal; declareAsNew(); }
    MEDCOUPLING_EXPORT void setIJSilent(int tupleId, int compoId, double newVal) { _mem[tupleId*_info_on_compo.size()+compoId]=newVal; }
    MEDCOUPLING_EXPORT double *getPointer() { return _mem.getPointer(); }
    MEDCOUPLING_EXPORT static void SetArrayIn(DataArrayDouble *newArray, DataArrayDouble* &arrayToSet);
    MEDCOUPLING_EXPORT const double *getConstPointer() const { return _mem.getConstPointer(); }
    MEDCOUPLING_EXPORT DataArrayDoubleIterator *iterator();
    MEDCOUPLING_EXPORT const double *begin() const { return getConstPointer(); }
    MEDCOUPLING_EXPORT const double *end() const { return getConstPointer()+getNbOfElems(); }
    MEDCOUPLING_EXPORT void useArray(const double *array, bool ownership, DeallocType type, int nbOfTuple, int nbOfCompo);
    MEDCOUPLING_EXPORT void useExternalArrayWithRWAccess(const double *array, int nbOfTuple, int nbOfCompo);
    template<class InputIterator>
    void insertAtTheEnd(InputIterator first, InputIterator last);
    MEDCOUPLING_EXPORT void writeOnPlace(std::size_t id, double element0, const double *others, int sizeOfOthers) { _mem.writeOnPlace(id,element0,others,sizeOfOthers); }
    MEDCOUPLING_EXPORT void checkNoNullValues() const;
    MEDCOUPLING_EXPORT void getMinMaxPerComponent(double *bounds) const;
    MEDCOUPLING_EXPORT DataArrayDouble *computeBBoxPerTuple(double epsilon=0.0) const;
    MEDCOUPLING_EXPORT void computeTupleIdsNearTuples(const DataArrayDouble *other, double eps, DataArrayInt *& c, DataArrayInt *& cI) const;
    MEDCOUPLING_EXPORT void recenterForMaxPrecision(double eps);
    MEDCOUPLING_EXPORT double getMaxValue(int& tupleId) const;
    MEDCOUPLING_EXPORT double getMaxValueInArray() const;
    MEDCOUPLING_EXPORT double getMinValue(int& tupleId) const;
    MEDCOUPLING_EXPORT double getMinValueInArray() const;
    MEDCOUPLING_EXPORT double getMaxValue2(DataArrayInt*& tupleIds) const;
    MEDCOUPLING_EXPORT double getMinValue2(DataArrayInt*& tupleIds) const;
    MEDCOUPLING_EXPORT int count(double value, double eps) const;
    MEDCOUPLING_EXPORT double getAverageValue() const;
    MEDCOUPLING_EXPORT double norm2() const;
    MEDCOUPLING_EXPORT double normMax() const;
    MEDCOUPLING_EXPORT double normMin() const;
    MEDCOUPLING_EXPORT void accumulate(double *res) const;
    MEDCOUPLING_EXPORT double accumulate(int compId) const;
    MEDCOUPLING_EXPORT DataArrayDouble *accumulatePerChunck(const int *bgOfIndex, const int *endOfIndex) const;
    MEDCOUPLING_EXPORT double distanceToTuple(const double *tupleBg, const double *tupleEnd, int& tupleId) const;
    MEDCOUPLING_EXPORT DataArrayDouble *fromPolarToCart() const;
    MEDCOUPLING_EXPORT DataArrayDouble *fromCylToCart() const;
    MEDCOUPLING_EXPORT DataArrayDouble *fromSpherToCart() const;
    MEDCOUPLING_EXPORT DataArrayDouble *doublyContractedProduct() const;
    MEDCOUPLING_EXPORT DataArrayDouble *determinant() const;
    MEDCOUPLING_EXPORT DataArrayDouble *eigenValues() const;
    MEDCOUPLING_EXPORT DataArrayDouble *eigenVectors() const;
    MEDCOUPLING_EXPORT DataArrayDouble *inverse() const;
    MEDCOUPLING_EXPORT DataArrayDouble *trace() const;
    MEDCOUPLING_EXPORT DataArrayDouble *deviator() const;
    MEDCOUPLING_EXPORT DataArrayDouble *magnitude() const;
    MEDCOUPLING_EXPORT DataArrayDouble *sumPerTuple() const;
    MEDCOUPLING_EXPORT DataArrayDouble *maxPerTuple() const;
    MEDCOUPLING_EXPORT DataArrayDouble *maxPerTupleWithCompoId(DataArrayInt* &compoIdOfMaxPerTuple) const;
    MEDCOUPLING_EXPORT DataArrayDouble *buildEuclidianDistanceDenseMatrix() const;
    MEDCOUPLING_EXPORT DataArrayDouble *buildEuclidianDistanceDenseMatrixWith(const DataArrayDouble *other) const;
    MEDCOUPLING_EXPORT void sortPerTuple(bool asc);
    MEDCOUPLING_EXPORT void abs();
    MEDCOUPLING_EXPORT DataArrayDouble *computeAbs() const;
    MEDCOUPLING_EXPORT void applyLin(double a, double b, int compoId);
    MEDCOUPLING_EXPORT void applyLin(double a, double b);
    MEDCOUPLING_EXPORT void applyInv(double numerator);
    MEDCOUPLING_EXPORT void applyPow(double val);
    MEDCOUPLING_EXPORT void applyRPow(double val);
    MEDCOUPLING_EXPORT DataArrayDouble *negate() const;
    MEDCOUPLING_EXPORT DataArrayDouble *applyFunc(int nbOfComp, FunctionToEvaluate func) const;
    MEDCOUPLING_EXPORT DataArrayDouble *applyFunc(int nbOfComp, const char *func) const;
    MEDCOUPLING_EXPORT DataArrayDouble *applyFunc(const char *func) const;
    MEDCOUPLING_EXPORT DataArrayDouble *applyFunc2(int nbOfComp, const char *func) const;
    MEDCOUPLING_EXPORT DataArrayDouble *applyFunc3(int nbOfComp, const std::vector<std::string>& varsOrder, const char *func) const;
    MEDCOUPLING_EXPORT void applyFuncFast32(const char *func);
    MEDCOUPLING_EXPORT void applyFuncFast64(const char *func);
    MEDCOUPLING_EXPORT DataArrayInt *getIdsInRange(double vmin, double vmax) const;
    MEDCOUPLING_EXPORT DataArrayInt *getIdsNotInRange(double vmin, double vmax) const;
    MEDCOUPLING_EXPORT static DataArrayDouble *Aggregate(const DataArrayDouble *a1, const DataArrayDouble *a2);
    MEDCOUPLING_EXPORT static DataArrayDouble *Aggregate(const std::vector<const DataArrayDouble *>& arr);
    MEDCOUPLING_EXPORT static DataArrayDouble *Meld(const DataArrayDouble *a1, const DataArrayDouble *a2);
    MEDCOUPLING_EXPORT static DataArrayDouble *Meld(const std::vector<const DataArrayDouble *>& arr);
    MEDCOUPLING_EXPORT static DataArrayDouble *Dot(const DataArrayDouble *a1, const DataArrayDouble *a2);
    MEDCOUPLING_EXPORT static DataArrayDouble *CrossProduct(const DataArrayDouble *a1, const DataArrayDouble *a2);
    MEDCOUPLING_EXPORT static DataArrayDouble *Max(const DataArrayDouble *a1, const DataArrayDouble *a2);
    MEDCOUPLING_EXPORT static DataArrayDouble *Min(const DataArrayDouble *a1, const DataArrayDouble *a2);
    MEDCOUPLING_EXPORT static DataArrayDouble *Add(const DataArrayDouble *a1, const DataArrayDouble *a2);
    MEDCOUPLING_EXPORT void addEqual(const DataArrayDouble *other);
    MEDCOUPLING_EXPORT static DataArrayDouble *Substract(const DataArrayDouble *a1, const DataArrayDouble *a2);
    MEDCOUPLING_EXPORT void substractEqual(const DataArrayDouble *other);
    MEDCOUPLING_EXPORT static DataArrayDouble *Multiply(const DataArrayDouble *a1, const DataArrayDouble *a2);
    MEDCOUPLING_EXPORT void multiplyEqual(const DataArrayDouble *other);
    MEDCOUPLING_EXPORT static DataArrayDouble *Divide(const DataArrayDouble *a1, const DataArrayDouble *a2);
    MEDCOUPLING_EXPORT void divideEqual(const DataArrayDouble *other);
    MEDCOUPLING_EXPORT static DataArrayDouble *Pow(const DataArrayDouble *a1, const DataArrayDouble *a2);
    MEDCOUPLING_EXPORT void powEqual(const DataArrayDouble *other);
    MEDCOUPLING_EXPORT void updateTime() const { }
    MEDCOUPLING_EXPORT MemArray<double>& accessToMemArray() { return _mem; }
    MEDCOUPLING_EXPORT const MemArray<double>& accessToMemArray() const { return _mem; }
  public:
    MEDCOUPLING_EXPORT void getTinySerializationIntInformation(std::vector<int>& tinyInfo) const;
    MEDCOUPLING_EXPORT void getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const;
    MEDCOUPLING_EXPORT bool resizeForUnserialization(const std::vector<int>& tinyInfoI);
    MEDCOUPLING_EXPORT void finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<std::string>& tinyInfoS);
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
  private:
    MemArray<double> _mem;
  };

  class DataArrayDoubleTuple;

  class DataArrayDoubleIterator
  {
  public:
    MEDCOUPLING_EXPORT DataArrayDoubleIterator(DataArrayDouble *da);
    MEDCOUPLING_EXPORT ~DataArrayDoubleIterator();
    MEDCOUPLING_EXPORT DataArrayDoubleTuple *nextt();
  private:
    DataArrayDouble *_da;
    double *_pt;
    int _tuple_id;
    int _nb_comp;
    int _nb_tuple;
  };

  class DataArrayDoubleTuple
  {
  public:
    MEDCOUPLING_EXPORT DataArrayDoubleTuple(double *pt, int nbOfComp);
    MEDCOUPLING_EXPORT std::string repr() const;
    MEDCOUPLING_EXPORT int getNumberOfCompo() const { return _nb_of_compo; }
    MEDCOUPLING_EXPORT const double *getConstPointer() const { return  _pt; }
    MEDCOUPLING_EXPORT double *getPointer() { return _pt; }
    MEDCOUPLING_EXPORT double doubleValue() const;
    MEDCOUPLING_EXPORT DataArrayDouble *buildDADouble(int nbOfTuples, int nbOfCompo) const;
  private:
    double *_pt;
    int _nb_of_compo;
  };

  class DataArrayIntIterator;

  class DataArrayInt : public DataArray
  {
  public:
    MEDCOUPLING_EXPORT static DataArrayInt *New();
    MEDCOUPLING_EXPORT bool isAllocated() const;
    MEDCOUPLING_EXPORT void checkAllocated() const;
    MEDCOUPLING_EXPORT void desallocate();
    MEDCOUPLING_EXPORT int getNumberOfTuples() const { return _info_on_compo.empty()?0:_mem.getNbOfElem()/getNumberOfComponents(); }
    MEDCOUPLING_EXPORT std::size_t getNbOfElems() const { return _mem.getNbOfElem(); }
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT int intValue() const;
    MEDCOUPLING_EXPORT int getHashCode() const;
    MEDCOUPLING_EXPORT bool empty() const;
    MEDCOUPLING_EXPORT DataArrayInt *deepCpy() const;
    MEDCOUPLING_EXPORT DataArrayInt *performCpy(bool deepCpy) const;
    MEDCOUPLING_EXPORT void cpyFrom(const DataArrayInt& other);
    MEDCOUPLING_EXPORT void reserve(std::size_t nbOfElems);
    MEDCOUPLING_EXPORT void pushBackSilent(int val);
    MEDCOUPLING_EXPORT void pushBackValsSilent(const int *valsBg, const int *valsEnd);
    MEDCOUPLING_EXPORT int popBackSilent();
    MEDCOUPLING_EXPORT void pack() const;
    MEDCOUPLING_EXPORT std::size_t getNbOfElemAllocated() const { return _mem.getNbOfElemAllocated(); }
    MEDCOUPLING_EXPORT void alloc(int nbOfTuple, int nbOfCompo=1);
    MEDCOUPLING_EXPORT void allocIfNecessary(int nbOfTuple, int nbOfCompo);
    MEDCOUPLING_EXPORT bool isEqual(const DataArrayInt& other) const;
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const DataArrayInt& other, std::string& reason) const;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const DataArrayInt& other) const;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStrAndOrder(const DataArrayInt& other) const;
    MEDCOUPLING_EXPORT bool isFittingWith(const std::vector<bool>& v) const;
    MEDCOUPLING_EXPORT DataArrayInt *buildPermutationArr(const DataArrayInt& other) const;
    MEDCOUPLING_EXPORT DataArrayInt *sumPerTuple() const;
    MEDCOUPLING_EXPORT void sort(bool asc=true);
    MEDCOUPLING_EXPORT void reverse();
    MEDCOUPLING_EXPORT void checkMonotonic(bool increasing) const;
    MEDCOUPLING_EXPORT bool isMonotonic(bool increasing) const;
    MEDCOUPLING_EXPORT void checkStrictlyMonotonic(bool increasing) const;
    MEDCOUPLING_EXPORT bool isStrictlyMonotonic(bool increasing) const;
    MEDCOUPLING_EXPORT void fillWithZero();
    MEDCOUPLING_EXPORT void fillWithValue(int val);
    MEDCOUPLING_EXPORT void iota(int init=0);
    MEDCOUPLING_EXPORT std::string repr() const;
    MEDCOUPLING_EXPORT std::string reprZip() const;
    MEDCOUPLING_EXPORT void writeVTK(std::ostream& ofs, int indent, const char *type, const char *nameInFile, DataArrayByte *byteArr) const;
    MEDCOUPLING_EXPORT void reprStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprZipStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprWithoutNameStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprZipWithoutNameStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprCppStream(const char *varName, std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprQuickOverview(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprQuickOverviewData(std::ostream& stream, std::size_t maxNbOfByteInRepr) const;
    MEDCOUPLING_EXPORT void transformWithIndArr(const int *indArrBg, const int *indArrEnd);
    MEDCOUPLING_EXPORT DataArrayInt *transformWithIndArrR(const int *indArrBg, const int *indArrEnd) const;
    MEDCOUPLING_EXPORT void splitByValueRange(const int *arrBg, const int *arrEnd,
                                              DataArrayInt *& castArr, DataArrayInt *& rankInsideCast, DataArrayInt *& castsPresent) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *invertArrayO2N2N2O(int newNbOfElem) const;
    MEDCOUPLING_EXPORT DataArrayInt *invertArrayN2O2O2N(int oldNbOfElem) const;
    MEDCOUPLING_EXPORT DataArrayInt *invertArrayO2N2N2OBis(int newNbOfElem) const;
    MEDCOUPLING_EXPORT void reAlloc(int nbOfTuples);
    MEDCOUPLING_EXPORT DataArrayDouble *convertToDblArr() const;
    MEDCOUPLING_EXPORT DataArrayInt *fromNoInterlace() const;
    MEDCOUPLING_EXPORT DataArrayInt *toNoInterlace() const;
    MEDCOUPLING_EXPORT void renumberInPlace(const int *old2New);
    MEDCOUPLING_EXPORT void renumberInPlaceR(const int *new2Old);
    MEDCOUPLING_EXPORT DataArrayInt *renumber(const int *old2New) const;
    MEDCOUPLING_EXPORT DataArrayInt *renumberR(const int *new2Old) const;
    MEDCOUPLING_EXPORT DataArrayInt *renumberAndReduce(const int *old2NewBg, int newNbOfTuple) const;
    MEDCOUPLING_EXPORT DataArrayInt *selectByTupleId(const int *new2OldBg, const int *new2OldEnd) const;
    MEDCOUPLING_EXPORT DataArrayInt *selectByTupleIdSafe(const int *new2OldBg, const int *new2OldEnd) const;
    MEDCOUPLING_EXPORT DataArrayInt *selectByTupleId2(int bg, int end, int step) const;
    MEDCOUPLING_EXPORT DataArray *selectByTupleRanges(const std::vector<std::pair<int,int> >& ranges) const;
    MEDCOUPLING_EXPORT DataArrayInt *checkAndPreparePermutation() const;
    MEDCOUPLING_EXPORT static DataArrayInt *FindPermutationFromFirstToSecond(const DataArrayInt *ids1, const DataArrayInt *ids2);
    MEDCOUPLING_EXPORT void changeSurjectiveFormat(int targetNb, DataArrayInt *&arr, DataArrayInt *&arrI) const;
    MEDCOUPLING_EXPORT static DataArrayInt *BuildOld2NewArrayFromSurjectiveFormat2(int nbOfOldTuples, const int *arr, const int *arrIBg, const int *arrIEnd, int &newNbOfTuples);
    MEDCOUPLING_EXPORT DataArrayInt *buildPermArrPerLevel() const;
    MEDCOUPLING_EXPORT bool isIdentity() const;
    MEDCOUPLING_EXPORT bool isUniform(int val) const;
    MEDCOUPLING_EXPORT DataArrayInt *substr(int tupleIdBg, int tupleIdEnd=-1) const;
    MEDCOUPLING_EXPORT void rearrange(int newNbOfCompo);
    MEDCOUPLING_EXPORT void transpose();
    MEDCOUPLING_EXPORT DataArrayInt *changeNbOfComponents(int newNbOfComp, int dftValue) const;
    MEDCOUPLING_EXPORT DataArray *keepSelectedComponents(const std::vector<int>& compoIds) const;
    MEDCOUPLING_EXPORT void meldWith(const DataArrayInt *other);
    MEDCOUPLING_EXPORT void setSelectedComponents(const DataArrayInt *a, const std::vector<int>& compoIds);
    MEDCOUPLING_EXPORT void setPartOfValues1(const DataArrayInt *a, int bgTuples, int endTuples, int stepTuples, int bgComp, int endComp, int stepComp, bool strictCompoCompare=true);
    MEDCOUPLING_EXPORT void setPartOfValuesSimple1(int a, int bgTuples, int endTuples, int stepTuples, int bgComp, int endComp, int stepComp);
    MEDCOUPLING_EXPORT void setPartOfValues2(const DataArrayInt *a, const int *bgTuples, const int *endTuples, const int *bgComp, const int *endComp, bool strictCompoCompare=true);
    MEDCOUPLING_EXPORT void setPartOfValuesSimple2(int a, const int *bgTuples, const int *endTuples, const int *bgComp, const int *endComp);
    MEDCOUPLING_EXPORT void setPartOfValues3(const DataArrayInt *a, const int *bgTuples, const int *endTuples, int bgComp, int endComp, int stepComp, bool strictCompoCompare=true);
    MEDCOUPLING_EXPORT void setPartOfValuesSimple3(int a, const int *bgTuples, const int *endTuples, int bgComp, int endComp, int stepComp);
    MEDCOUPLING_EXPORT void setPartOfValues4(const DataArrayInt *a, int bgTuples, int endTuples, int stepTuples, const int *bgComp, const int *endComp, bool strictCompoCompare=true);
    MEDCOUPLING_EXPORT void setPartOfValuesSimple4(int a, int bgTuples, int endTuples, int stepTuples, const int *bgComp, const int *endComp);
    MEDCOUPLING_EXPORT void setPartOfValuesAdv(const DataArrayInt *a, const DataArrayInt *tuplesSelec);
    MEDCOUPLING_EXPORT void setContigPartOfSelectedValues(int tupleIdStart, const DataArray *aBase, const DataArrayInt *tuplesSelec);
    MEDCOUPLING_EXPORT void setContigPartOfSelectedValues2(int tupleIdStart, const DataArray *aBase, int bg, int end2, int step);
    MEDCOUPLING_EXPORT void getTuple(int tupleId, int *res) const { std::copy(_mem.getConstPointerLoc(tupleId*_info_on_compo.size()),_mem.getConstPointerLoc((tupleId+1)*_info_on_compo.size()),res); }
    MEDCOUPLING_EXPORT int getIJ(int tupleId, int compoId) const { return _mem[tupleId*_info_on_compo.size()+compoId]; }
    MEDCOUPLING_EXPORT int getIJSafe(int tupleId, int compoId) const;
    MEDCOUPLING_EXPORT int front() const;
    MEDCOUPLING_EXPORT int back() const;
    MEDCOUPLING_EXPORT void setIJ(int tupleId, int compoId, int newVal) { _mem[tupleId*_info_on_compo.size()+compoId]=newVal; declareAsNew(); }
    MEDCOUPLING_EXPORT void setIJSilent(int tupleId, int compoId, int newVal) { _mem[tupleId*_info_on_compo.size()+compoId]=newVal; }
    MEDCOUPLING_EXPORT int *getPointer() { return _mem.getPointer(); }
    MEDCOUPLING_EXPORT static void SetArrayIn(DataArrayInt *newArray, DataArrayInt* &arrayToSet);
    MEDCOUPLING_EXPORT const int *getConstPointer() const { return _mem.getConstPointer(); }
    MEDCOUPLING_EXPORT DataArrayIntIterator *iterator();
    MEDCOUPLING_EXPORT const int *begin() const { return getConstPointer(); }
    MEDCOUPLING_EXPORT const int *end() const { return getConstPointer()+getNbOfElems(); }
    MEDCOUPLING_EXPORT DataArrayInt *getIdsEqual(int val) const;
    MEDCOUPLING_EXPORT DataArrayInt *getIdsNotEqual(int val) const;
    MEDCOUPLING_EXPORT DataArrayInt *getIdsEqualList(const int *valsBg, const int *valsEnd) const;
    MEDCOUPLING_EXPORT DataArrayInt *getIdsNotEqualList(const int *valsBg, const int *valsEnd) const;
    MEDCOUPLING_EXPORT DataArrayInt *getIdsEqualTuple(const int *tupleBg, const int *tupleEnd) const;
    MEDCOUPLING_EXPORT int changeValue(int oldValue, int newValue);
    MEDCOUPLING_EXPORT int locateTuple(const std::vector<int>& tupl) const;
    MEDCOUPLING_EXPORT int locateValue(int value) const;
    MEDCOUPLING_EXPORT int locateValue(const std::vector<int>& vals) const;
    MEDCOUPLING_EXPORT int search(const std::vector<int>& vals) const;
    MEDCOUPLING_EXPORT bool presenceOfTuple(const std::vector<int>& tupl) const;
    MEDCOUPLING_EXPORT bool presenceOfValue(int value) const;
    MEDCOUPLING_EXPORT bool presenceOfValue(const std::vector<int>& vals) const;
    MEDCOUPLING_EXPORT int count(int value) const;
    MEDCOUPLING_EXPORT void accumulate(int *res) const;
    MEDCOUPLING_EXPORT int accumulate(int compId) const;
    MEDCOUPLING_EXPORT DataArrayInt *accumulatePerChunck(const int *bgOfIndex, const int *endOfIndex) const;
    MEDCOUPLING_EXPORT int getMaxValue(int& tupleId) const;
    MEDCOUPLING_EXPORT int getMaxValueInArray() const;
    MEDCOUPLING_EXPORT int getMinValue(int& tupleId) const;
    MEDCOUPLING_EXPORT int getMinValueInArray() const;
    MEDCOUPLING_EXPORT void abs();
    MEDCOUPLING_EXPORT DataArrayInt *computeAbs() const;
    MEDCOUPLING_EXPORT void applyLin(int a, int b, int compoId);
    MEDCOUPLING_EXPORT void applyLin(int a, int b);
    MEDCOUPLING_EXPORT void applyInv(int numerator);
    MEDCOUPLING_EXPORT DataArrayInt *negate() const;
    MEDCOUPLING_EXPORT void applyDivideBy(int val);
    MEDCOUPLING_EXPORT void applyModulus(int val);
    MEDCOUPLING_EXPORT void applyRModulus(int val);
    MEDCOUPLING_EXPORT void applyPow(int val);
    MEDCOUPLING_EXPORT void applyRPow(int val);
    MEDCOUPLING_EXPORT DataArrayInt *getIdsInRange(int vmin, int vmax) const;
    MEDCOUPLING_EXPORT DataArrayInt *getIdsNotInRange(int vmin, int vmax) const;
    MEDCOUPLING_EXPORT bool checkAllIdsInRange(int vmin, int vmax) const;
    MEDCOUPLING_EXPORT static DataArrayInt *Aggregate(const DataArrayInt *a1, const DataArrayInt *a2, int offsetA2);
    MEDCOUPLING_EXPORT static DataArrayInt *Aggregate(const std::vector<const DataArrayInt *>& arr);
    MEDCOUPLING_EXPORT static DataArrayInt *AggregateIndexes(const std::vector<const DataArrayInt *>& arrs);
    MEDCOUPLING_EXPORT static DataArrayInt *Meld(const DataArrayInt *a1, const DataArrayInt *a2);
    MEDCOUPLING_EXPORT static DataArrayInt *Meld(const std::vector<const DataArrayInt *>& arr);
    MEDCOUPLING_EXPORT static DataArrayInt *MakePartition(const std::vector<const DataArrayInt *>& groups, int newNb, std::vector< std::vector<int> >& fidsOfGroups);
    MEDCOUPLING_EXPORT static DataArrayInt *BuildUnion(const std::vector<const DataArrayInt *>& arr);
    MEDCOUPLING_EXPORT static DataArrayInt *BuildIntersection(const std::vector<const DataArrayInt *>& arr);
    MEDCOUPLING_EXPORT DataArrayInt *buildComplement(int nbOfElement) const;
    MEDCOUPLING_EXPORT DataArrayInt *buildSubstraction(const DataArrayInt *other) const;
    MEDCOUPLING_EXPORT DataArrayInt *buildSubstractionOptimized(const DataArrayInt *other) const;
    MEDCOUPLING_EXPORT DataArrayInt *buildUnion(const DataArrayInt *other) const;
    MEDCOUPLING_EXPORT DataArrayInt *buildIntersection(const DataArrayInt *other) const;
    MEDCOUPLING_EXPORT DataArrayInt *buildUnique() const;
    MEDCOUPLING_EXPORT DataArrayInt *deltaShiftIndex() const;
    MEDCOUPLING_EXPORT void computeOffsets();
    MEDCOUPLING_EXPORT void computeOffsets2();
    MEDCOUPLING_EXPORT void searchRangesInListOfIds(const DataArrayInt *listOfIds, DataArrayInt *& rangeIdsFetched, DataArrayInt *& idsInInputListThatFetch) const;
    MEDCOUPLING_EXPORT DataArrayInt *buildExplicitArrByRanges(const DataArrayInt *offsets) const;
    MEDCOUPLING_EXPORT DataArrayInt *buildExplicitArrOfSliceOnScaledArr(int begin, int stop, int step) const;
    MEDCOUPLING_EXPORT DataArrayInt *findRangeIdForEachTuple(const DataArrayInt *ranges) const;
    MEDCOUPLING_EXPORT DataArrayInt *findIdInRangeForEachTuple(const DataArrayInt *ranges) const;
    MEDCOUPLING_EXPORT DataArrayInt *duplicateEachTupleNTimes(int nbTimes) const;
    MEDCOUPLING_EXPORT DataArrayInt *getDifferentValues() const;
    MEDCOUPLING_EXPORT std::vector<DataArrayInt *> partitionByDifferentValues(std::vector<int>& differentIds) const;
    MEDCOUPLING_EXPORT std::vector< std::pair<int,int> > splitInBalancedSlices(int nbOfSlices) const;
    MEDCOUPLING_EXPORT void useArray(const int *array, bool ownership, DeallocType type, int nbOfTuple, int nbOfCompo);
    MEDCOUPLING_EXPORT void useExternalArrayWithRWAccess(const int *array, int nbOfTuple, int nbOfCompo);
    template<class InputIterator>
    void insertAtTheEnd(InputIterator first, InputIterator last);
    MEDCOUPLING_EXPORT void writeOnPlace(std::size_t id, int element0, const int *others, int sizeOfOthers) { _mem.writeOnPlace(id,element0,others,sizeOfOthers); }
    MEDCOUPLING_EXPORT static DataArrayInt *Add(const DataArrayInt *a1, const DataArrayInt *a2);
    MEDCOUPLING_EXPORT void addEqual(const DataArrayInt *other);
    MEDCOUPLING_EXPORT static DataArrayInt *Substract(const DataArrayInt *a1, const DataArrayInt *a2);
    MEDCOUPLING_EXPORT void substractEqual(const DataArrayInt *other);
    MEDCOUPLING_EXPORT static DataArrayInt *Multiply(const DataArrayInt *a1, const DataArrayInt *a2);
    MEDCOUPLING_EXPORT void multiplyEqual(const DataArrayInt *other);
    MEDCOUPLING_EXPORT static DataArrayInt *Divide(const DataArrayInt *a1, const DataArrayInt *a2);
    MEDCOUPLING_EXPORT void divideEqual(const DataArrayInt *other);
    MEDCOUPLING_EXPORT static DataArrayInt *Modulus(const DataArrayInt *a1, const DataArrayInt *a2);
    MEDCOUPLING_EXPORT void modulusEqual(const DataArrayInt *other);
    MEDCOUPLING_EXPORT static DataArrayInt *Pow(const DataArrayInt *a1, const DataArrayInt *a2);
    MEDCOUPLING_EXPORT void powEqual(const DataArrayInt *other);
    MEDCOUPLING_EXPORT void updateTime() const { }
    MEDCOUPLING_EXPORT MemArray<int>& accessToMemArray() { return _mem; }
    MEDCOUPLING_EXPORT const MemArray<int>& accessToMemArray() const { return _mem; }
  public:
    MEDCOUPLING_EXPORT static int *CheckAndPreparePermutation(const int *start, const int *end);
    MEDCOUPLING_EXPORT static DataArrayInt *Range(int begin, int end, int step);
  public:
    MEDCOUPLING_EXPORT void getTinySerializationIntInformation(std::vector<int>& tinyInfo) const;
    MEDCOUPLING_EXPORT void getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const;
    MEDCOUPLING_EXPORT bool resizeForUnserialization(const std::vector<int>& tinyInfoI);
    MEDCOUPLING_EXPORT void finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<std::string>& tinyInfoS);
  private:
    ~DataArrayInt() { }
    DataArrayInt() { }
  private:
    MemArray<int> _mem;
  };

  class DataArrayIntTuple;

  class DataArrayIntIterator
  {
  public:
    MEDCOUPLING_EXPORT DataArrayIntIterator(DataArrayInt *da);
    MEDCOUPLING_EXPORT ~DataArrayIntIterator();
    MEDCOUPLING_EXPORT DataArrayIntTuple *nextt();
  private:
    DataArrayInt *_da;
    int *_pt;
    int _tuple_id;
    int _nb_comp;
    int _nb_tuple;
  };

  class DataArrayIntTuple
  {
  public:
    MEDCOUPLING_EXPORT DataArrayIntTuple(int *pt, int nbOfComp);
    MEDCOUPLING_EXPORT std::string repr() const;
    MEDCOUPLING_EXPORT int getNumberOfCompo() const { return _nb_of_compo; }
    MEDCOUPLING_EXPORT const int *getConstPointer() const { return  _pt; }
    MEDCOUPLING_EXPORT int *getPointer() { return _pt; }
    MEDCOUPLING_EXPORT int intValue() const;
    MEDCOUPLING_EXPORT DataArrayInt *buildDAInt(int nbOfTuples, int nbOfCompo) const;
  private:
    int *_pt;
    int _nb_of_compo;
  };

  class DataArrayChar : public DataArray
  {
  public:
    MEDCOUPLING_EXPORT virtual DataArrayChar *buildEmptySpecializedDAChar() const = 0;
    MEDCOUPLING_EXPORT bool isAllocated() const;
    MEDCOUPLING_EXPORT void checkAllocated() const;
    MEDCOUPLING_EXPORT void desallocate();
    MEDCOUPLING_EXPORT int getNumberOfTuples() const { return _info_on_compo.empty()?0:_mem.getNbOfElem()/getNumberOfComponents(); }
    MEDCOUPLING_EXPORT std::size_t getNbOfElems() const { return _mem.getNbOfElem(); }
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT int getHashCode() const;
    MEDCOUPLING_EXPORT bool empty() const;
    MEDCOUPLING_EXPORT void cpyFrom(const DataArrayChar& other);
    MEDCOUPLING_EXPORT void reserve(std::size_t nbOfElems);
    MEDCOUPLING_EXPORT void pushBackSilent(char val);
    MEDCOUPLING_EXPORT void pushBackValsSilent(const char *valsBg, const char *valsEnd);
    MEDCOUPLING_EXPORT char popBackSilent();
    MEDCOUPLING_EXPORT void pack() const;
    MEDCOUPLING_EXPORT std::size_t getNbOfElemAllocated() const { return _mem.getNbOfElemAllocated(); }
    MEDCOUPLING_EXPORT void alloc(int nbOfTuple, int nbOfCompo=1);
    MEDCOUPLING_EXPORT void allocIfNecessary(int nbOfTuple, int nbOfCompo);
    MEDCOUPLING_EXPORT bool isEqual(const DataArrayChar& other) const;
    MEDCOUPLING_EXPORT virtual bool isEqualIfNotWhy(const DataArrayChar& other, std::string& reason) const;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const DataArrayChar& other) const;
    MEDCOUPLING_EXPORT void reverse();
    MEDCOUPLING_EXPORT void fillWithZero();
    MEDCOUPLING_EXPORT void fillWithValue(char val);
    MEDCOUPLING_EXPORT std::string repr() const;
    MEDCOUPLING_EXPORT std::string reprZip() const;
    MEDCOUPLING_EXPORT void reAlloc(int nbOfTuples);
    MEDCOUPLING_EXPORT DataArrayInt *convertToIntArr() const;
    MEDCOUPLING_EXPORT void renumberInPlace(const int *old2New);
    MEDCOUPLING_EXPORT void renumberInPlaceR(const int *new2Old);
    MEDCOUPLING_EXPORT DataArrayChar *renumber(const int *old2New) const;
    MEDCOUPLING_EXPORT DataArrayChar *renumberR(const int *new2Old) const;
    MEDCOUPLING_EXPORT DataArrayChar *renumberAndReduce(const int *old2NewBg, int newNbOfTuple) const;
    MEDCOUPLING_EXPORT DataArrayChar *selectByTupleId(const int *new2OldBg, const int *new2OldEnd) const;
    MEDCOUPLING_EXPORT DataArrayChar *selectByTupleIdSafe(const int *new2OldBg, const int *new2OldEnd) const;
    MEDCOUPLING_EXPORT DataArrayChar *selectByTupleId2(int bg, int end, int step) const;
    MEDCOUPLING_EXPORT bool isUniform(char val) const;
    MEDCOUPLING_EXPORT void rearrange(int newNbOfCompo);
    MEDCOUPLING_EXPORT DataArrayChar *substr(int tupleIdBg, int tupleIdEnd=-1) const;
    MEDCOUPLING_EXPORT DataArrayChar *changeNbOfComponents(int newNbOfComp, char dftValue) const;
    MEDCOUPLING_EXPORT DataArray *keepSelectedComponents(const std::vector<int>& compoIds) const;
    MEDCOUPLING_EXPORT void meldWith(const DataArrayChar *other);
    MEDCOUPLING_EXPORT void setPartOfValues1(const DataArrayChar *a, int bgTuples, int endTuples, int stepTuples, int bgComp, int endComp, int stepComp, bool strictCompoCompare=true);
    MEDCOUPLING_EXPORT void setPartOfValuesSimple1(char a, int bgTuples, int endTuples, int stepTuples, int bgComp, int endComp, int stepComp);
    MEDCOUPLING_EXPORT void setPartOfValues2(const DataArrayChar *a, const int *bgTuples, const int *endTuples, const int *bgComp, const int *endComp, bool strictCompoCompare=true);
    MEDCOUPLING_EXPORT void setPartOfValuesSimple2(char a, const int *bgTuples, const int *endTuples, const int *bgComp, const int *endComp);
    MEDCOUPLING_EXPORT void setPartOfValues3(const DataArrayChar *a, const int *bgTuples, const int *endTuples, int bgComp, int endComp, int stepComp, bool strictCompoCompare=true);
    MEDCOUPLING_EXPORT void setPartOfValuesSimple3(char a, const int *bgTuples, const int *endTuples, int bgComp, int endComp, int stepComp);
    MEDCOUPLING_EXPORT void setPartOfValues4(const DataArrayChar *a, int bgTuples, int endTuples, int stepTuples, const int *bgComp, const int *endComp, bool strictCompoCompare=true);
    MEDCOUPLING_EXPORT void setPartOfValuesSimple4(char a, int bgTuples, int endTuples, int stepTuples, const int *bgComp, const int *endComp);
    MEDCOUPLING_EXPORT void setPartOfValuesAdv(const DataArrayChar *a, const DataArrayChar *tuplesSelec);
    MEDCOUPLING_EXPORT void setContigPartOfSelectedValues(int tupleIdStart, const DataArray *aBase, const DataArrayInt *tuplesSelec);
    MEDCOUPLING_EXPORT void setContigPartOfSelectedValues2(int tupleIdStart, const DataArray *aBase, int bg, int end2, int step);
    MEDCOUPLING_EXPORT DataArray *selectByTupleRanges(const std::vector<std::pair<int,int> >& ranges) const;
    MEDCOUPLING_EXPORT void getTuple(int tupleId, char *res) const { std::copy(_mem.getConstPointerLoc(tupleId*_info_on_compo.size()),_mem.getConstPointerLoc((tupleId+1)*_info_on_compo.size()),res); }
    MEDCOUPLING_EXPORT char getIJ(int tupleId, int compoId) const { return _mem[tupleId*_info_on_compo.size()+compoId]; }
    MEDCOUPLING_EXPORT char getIJSafe(int tupleId, int compoId) const;
    MEDCOUPLING_EXPORT char front() const;
    MEDCOUPLING_EXPORT char back() const;
    MEDCOUPLING_EXPORT void setIJ(int tupleId, int compoId, char newVal) { _mem[tupleId*_info_on_compo.size()+compoId]=newVal; declareAsNew(); }
    MEDCOUPLING_EXPORT void setIJSilent(int tupleId, int compoId, char newVal) { _mem[tupleId*_info_on_compo.size()+compoId]=newVal; }
    MEDCOUPLING_EXPORT char *getPointer() { return _mem.getPointer(); }
    MEDCOUPLING_EXPORT const char *getConstPointer() const { return _mem.getConstPointer(); }
    MEDCOUPLING_EXPORT const char *begin() const { return getConstPointer(); }
    MEDCOUPLING_EXPORT const char *end() const { return getConstPointer()+getNbOfElems(); }
    MEDCOUPLING_EXPORT DataArrayInt *getIdsEqual(char val) const;
    MEDCOUPLING_EXPORT DataArrayInt *getIdsNotEqual(char val) const;
    MEDCOUPLING_EXPORT int search(const std::vector<char>& vals) const;
    MEDCOUPLING_EXPORT int locateTuple(const std::vector<char>& tupl) const;
    MEDCOUPLING_EXPORT int locateValue(char value) const;
    MEDCOUPLING_EXPORT int locateValue(const std::vector<char>& vals) const;
    MEDCOUPLING_EXPORT bool presenceOfTuple(const std::vector<char>& tupl) const;
    MEDCOUPLING_EXPORT bool presenceOfValue(char value) const;
    MEDCOUPLING_EXPORT bool presenceOfValue(const std::vector<char>& vals) const;
    MEDCOUPLING_EXPORT char getMaxValue(int& tupleId) const;
    MEDCOUPLING_EXPORT char getMaxValueInArray() const;
    MEDCOUPLING_EXPORT char getMinValue(int& tupleId) const;
    MEDCOUPLING_EXPORT char getMinValueInArray() const;
    MEDCOUPLING_EXPORT DataArrayInt *getIdsInRange(char vmin, char vmax) const;
    MEDCOUPLING_EXPORT static DataArrayChar *Aggregate(const DataArrayChar *a1, const DataArrayChar *a2);
    MEDCOUPLING_EXPORT static DataArrayChar *Aggregate(const std::vector<const DataArrayChar *>& arr);
    MEDCOUPLING_EXPORT static DataArrayChar *Meld(const DataArrayChar *a1, const DataArrayChar *a2);
    MEDCOUPLING_EXPORT static DataArrayChar *Meld(const std::vector<const DataArrayChar *>& arr);
    MEDCOUPLING_EXPORT void useArray(const char *array, bool ownership, DeallocType type, int nbOfTuple, int nbOfCompo);
    template<class InputIterator>
    void insertAtTheEnd(InputIterator first, InputIterator last);
    MEDCOUPLING_EXPORT void useExternalArrayWithRWAccess(const char *array, int nbOfTuple, int nbOfCompo);
    MEDCOUPLING_EXPORT void updateTime() const { }
    MEDCOUPLING_EXPORT MemArray<char>& accessToMemArray() { return _mem; }
    MEDCOUPLING_EXPORT const MemArray<char>& accessToMemArray() const { return _mem; }
  public:
    //MEDCOUPLING_EXPORT void getTinySerializationIntInformation(std::vector<int>& tinyInfo) const;
    //MEDCOUPLING_EXPORT void getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const;
    //MEDCOUPLING_EXPORT bool resizeForUnserialization(const std::vector<int>& tinyInfoI);
    //MEDCOUPLING_EXPORT void finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<std::string>& tinyInfoS);
  protected:
    DataArrayChar() { }
  protected:
    MemArray<char> _mem;
  };
  
  class DataArrayByteIterator;

  class DataArrayByte : public DataArrayChar
  {
  public:
    MEDCOUPLING_EXPORT static DataArrayByte *New();
    MEDCOUPLING_EXPORT DataArrayChar *buildEmptySpecializedDAChar() const;
    MEDCOUPLING_EXPORT DataArrayByteIterator *iterator();
    MEDCOUPLING_EXPORT DataArrayByte *deepCpy() const;
    MEDCOUPLING_EXPORT DataArrayByte *performCpy(bool deepCpy) const;
    MEDCOUPLING_EXPORT char byteValue() const;
    MEDCOUPLING_EXPORT void reprStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprZipStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprWithoutNameStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprZipWithoutNameStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprCppStream(const char *varName, std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprQuickOverview(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprQuickOverviewData(std::ostream& stream, std::size_t maxNbOfByteInRepr) const;
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const DataArrayChar& other, std::string& reason) const;
  private:
    ~DataArrayByte() { }
    DataArrayByte() { }
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
  
  class DataArrayAsciiCharIterator;
  
  class DataArrayAsciiChar : public DataArrayChar
  {
  public:
    MEDCOUPLING_EXPORT static DataArrayAsciiChar *New();
    MEDCOUPLING_EXPORT static DataArrayAsciiChar *New(const std::string& st);
    MEDCOUPLING_EXPORT static DataArrayAsciiChar *New(const std::vector<std::string>& vst, char defaultChar);
    MEDCOUPLING_EXPORT DataArrayChar *buildEmptySpecializedDAChar() const;
    MEDCOUPLING_EXPORT DataArrayAsciiCharIterator *iterator();
    MEDCOUPLING_EXPORT DataArrayAsciiChar *deepCpy() const;
    MEDCOUPLING_EXPORT DataArrayAsciiChar *performCpy(bool deepCpy) const;
    MEDCOUPLING_EXPORT char asciiCharValue() const;
    MEDCOUPLING_EXPORT void reprStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprZipStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprWithoutNameStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprZipWithoutNameStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprCppStream(const char *varName, std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprQuickOverview(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprQuickOverviewData(std::ostream& stream, std::size_t maxNbOfByteInRepr) const;
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const DataArrayChar& other, std::string& reason) const;
  private:
    ~DataArrayAsciiChar() { }
    DataArrayAsciiChar() { }
    DataArrayAsciiChar(const std::string& st);
    DataArrayAsciiChar(const std::vector<std::string>& vst, char defaultChar);
  };

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
