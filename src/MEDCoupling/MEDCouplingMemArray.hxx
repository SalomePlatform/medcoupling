// Copyright (C) 2007-2012  CEA/DEN, EDF R&D
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

#ifndef __PARAMEDMEM_MEDCOUPLINGMEMARRAY_HXX__
#define __PARAMEDMEM_MEDCOUPLINGMEMARRAY_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingTimeLabel.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "InterpKernelException.hxx"

#include <set>
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
    const T *getConstPointerLoc(int offset) const { if(_internal) return _internal+offset; else return _external+offset; }
    T *getPointer() { if(_internal) return _internal; if(_external) throw INTERP_KERNEL::Exception("Trying to write on an external pointer."); else return 0; }
  private:
    T *_internal;
    const T *_external;
  };

  template<class T>
  class MemArray
  {
  public:
    MemArray():_nb_of_elem(-1),_ownership(false),_dealloc(CPP_DEALLOC) { }
    MemArray(const MemArray<T>& other);
    bool isNull() const { return _pointer.isNull(); }
    const T *getConstPointerLoc(int offset) const { return _pointer.getConstPointerLoc(offset); }
    const T *getConstPointer() const { return _pointer.getConstPointer(); }
    T *getPointer() { return _pointer.getPointer(); }
    MemArray<T> &operator=(const MemArray<T>& other);
    T operator[](int id) const { return _pointer.getConstPointer()[id]; }
    T& operator[](int id) { return _pointer.getPointer()[id]; }
    bool isEqual(const MemArray<T>& other, T prec) const;
    void repr(int sl, std::ostream& stream) const;
    void reprZip(int sl, std::ostream& stream) const;
    void fillWithValue(const T& val);
    T *fromNoInterlace(int nbOfComp) const;
    T *toNoInterlace(int nbOfComp) const;
    void sort(bool asc);
    void reverse();
    void alloc(int nbOfElements) throw(INTERP_KERNEL::Exception);
    void reAlloc(int newNbOfElements) throw(INTERP_KERNEL::Exception);
    void useArray(const T *array, bool ownership, DeallocType type, int nbOfElem);
    void writeOnPlace(int id, T element0, const T *others, int sizeOfOthers);
    ~MemArray() { destroy(); }
  private:
    void destroy();
    static void destroyPointer(T *pt, DeallocType type);
  private:
    int _nb_of_elem;
    bool _ownership;
    MEDCouplingPointer<T> _pointer;
    DeallocType _dealloc;
  };

  class DataArray : public RefCountObject, public TimeLabel
  {
  public:
    MEDCOUPLING_EXPORT void setName(const char *name);
    MEDCOUPLING_EXPORT void copyStringInfoFrom(const DataArray& other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void copyPartOfStringInfoFrom(const DataArray& other, const std::vector<int>& compoIds) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void copyPartOfStringInfoFrom2(const std::vector<int>& compoIds, const DataArray& other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool areInfoEquals(const DataArray& other) const;
    MEDCOUPLING_EXPORT void reprWithoutNameStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT std::string getName() const { return _name; }
    MEDCOUPLING_EXPORT const std::vector<std::string> &getInfoOnComponents() const { return _info_on_compo; }
    MEDCOUPLING_EXPORT std::vector<std::string> &getInfoOnComponents() { return _info_on_compo; }
    MEDCOUPLING_EXPORT void setInfoOnComponents(const std::vector<std::string>& info) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT std::vector<std::string> getVarsOnComponent() const;
    MEDCOUPLING_EXPORT std::vector<std::string> getUnitsOnComponent() const;
    MEDCOUPLING_EXPORT std::string getInfoOnComponent(int i) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT std::string getVarOnComponent(int i) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT std::string getUnitOnComponent(int i) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setInfoOnComponent(int i, const char *info) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT int getNumberOfComponents() const { return (int)_info_on_compo.size(); }
    MEDCOUPLING_EXPORT int getNumberOfTuples() const { return _nb_of_tuples; }
    MEDCOUPLING_EXPORT int getNbOfElems() const { return ((int)_info_on_compo.size())*_nb_of_tuples; }
    MEDCOUPLING_EXPORT void checkNbOfTuplesAndComp(const DataArray& other, const char *msg) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void checkNbOfTuplesAndComp(int nbOfTuples, int nbOfCompo, const char *msg) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void checkNbOfElems(int nbOfElems, const char *msg) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static int GetNumberOfItemGivenBES(int begin, int end, int step, const char *msg) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static int GetNumberOfItemGivenBESRelative(int begin, int end, int step, const char *msg) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static std::string GetVarNameFromInfo(const std::string& info) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static std::string GetUnitFromInfo(const std::string& info) throw(INTERP_KERNEL::Exception);
  protected:
    DataArray():_nb_of_tuples(-1) { }
  protected:
    static void CheckValueInRange(int ref, int value, const char *msg) throw(INTERP_KERNEL::Exception);
    static void CheckClosingParInRange(int ref, int value, const char *msg) throw(INTERP_KERNEL::Exception);
  protected:
    int _nb_of_tuples;
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
    MEDCOUPLING_EXPORT void checkAllocated() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setInfoAndChangeNbOfCompo(const std::vector<std::string>& info) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT double doubleValue() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool empty() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *deepCpy() const;
    MEDCOUPLING_EXPORT DataArrayDouble *performCpy(bool deepCpy) const;
    MEDCOUPLING_EXPORT void cpyFrom(const DataArrayDouble& other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void alloc(int nbOfTuple, int nbOfCompo) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void allocIfNecessary(int nbOfTuple, int nbOfCompo);
    MEDCOUPLING_EXPORT void fillWithZero() throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void fillWithValue(double val) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void iota(double init=0.) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool isUniform(double val, double eps) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void sort(bool asc=true) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void reverse() throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void checkMonotonic(bool increasing, double eps) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool isMonotonic(bool increasing, double eps) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT std::string repr() const;
    MEDCOUPLING_EXPORT std::string reprZip() const;
    MEDCOUPLING_EXPORT void writeVTK(std::ostream& ofs, int indent, const char *nameInFile) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void reprStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprZipStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprWithoutNameStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprZipWithoutNameStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT bool isEqual(const DataArrayDouble& other, double prec) const;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const DataArrayDouble& other, double prec) const;
    //!alloc or useArray should have been called before.
    MEDCOUPLING_EXPORT void reAlloc(int nbOfTuples) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *convertToIntArr() const;
    MEDCOUPLING_EXPORT DataArrayDouble *fromNoInterlace() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *toNoInterlace() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void renumberInPlace(const int *old2New);
    MEDCOUPLING_EXPORT void renumberInPlaceR(const int *new2Old);
    MEDCOUPLING_EXPORT DataArrayDouble *renumber(const int *old2New) const;
    MEDCOUPLING_EXPORT DataArrayDouble *renumberR(const int *new2Old) const;
    MEDCOUPLING_EXPORT DataArrayDouble *renumberAndReduce(const int *old2New, int newNbOfTuple) const;
    MEDCOUPLING_EXPORT DataArrayDouble *selectByTupleId(const int *new2OldBg, const int *new2OldEnd) const;
    MEDCOUPLING_EXPORT DataArrayDouble *selectByTupleIdSafe(const int *new2OldBg, const int *new2OldEnd) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *selectByTupleId2(int bg, int end2, int step) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *selectByTupleRanges(const std::vector<std::pair<int,int> >& ranges) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *substr(int tupleIdBg, int tupleIdEnd=-1) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void rearrange(int newNbOfCompo) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void transpose() throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *changeNbOfComponents(int newNbOfComp, double dftValue) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *keepSelectedComponents(const std::vector<int>& compoIds) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void meldWith(const DataArrayDouble *other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void findCommonTuples(double prec, int limitTupleId, DataArrayInt *&comm, DataArrayInt *&commIndex) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *getDifferentValues(double prec, int limitTupleId=-1) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setSelectedComponents(const DataArrayDouble *a, const std::vector<int>& compoIds) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setPartOfValues1(const DataArrayDouble *a, int bgTuples, int endTuples, int stepTuples, int bgComp, int endComp, int stepComp, bool strictCompoCompare=true) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setPartOfValuesSimple1(double a, int bgTuples, int endTuples, int stepTuples, int bgComp, int endComp, int stepComp) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setPartOfValues2(const DataArrayDouble *a, const int *bgTuples, const int *endTuples, const int *bgComp, const int *endComp, bool strictCompoCompare=true) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setPartOfValuesSimple2(double a, const int *bgTuples, const int *endTuples, const int *bgComp, const int *endComp) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setPartOfValues3(const DataArrayDouble *a, const int *bgTuples, const int *endTuples, int bgComp, int endComp, int stepComp, bool strictCompoCompare=true) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setPartOfValuesSimple3(double a, const int *bgTuples, const int *endTuples, int bgComp, int endComp, int stepComp) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setPartOfValuesAdv(const DataArrayDouble *a, const DataArrayInt *tuplesSelec) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setContigPartOfSelectedValues(int tupleIdStart, const DataArrayDouble *a, const DataArrayInt *tuplesSelec) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setContigPartOfSelectedValues2(int tupleIdStart, const DataArrayDouble *a, int bg, int end2, int step) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getTuple(int tupleId, double *res) const { std::copy(_mem.getConstPointerLoc(tupleId*((int)_info_on_compo.size())),_mem.getConstPointerLoc((tupleId+1)*((int)_info_on_compo.size())),res); }
    MEDCOUPLING_EXPORT double getIJ(int tupleId, int compoId) const { return _mem[tupleId*((int)_info_on_compo.size())+compoId]; }
    MEDCOUPLING_EXPORT double back() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT double getIJSafe(int tupleId, int compoId) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setIJ(int tupleId, int compoId, double newVal) { _mem[tupleId*((int)_info_on_compo.size())+compoId]=newVal; declareAsNew(); }
    MEDCOUPLING_EXPORT void setIJSilent(int tupleId, int compoId, double newVal) { _mem[tupleId*((int)_info_on_compo.size())+compoId]=newVal; }
    MEDCOUPLING_EXPORT double *getPointer() { return _mem.getPointer(); }
    MEDCOUPLING_EXPORT static void SetArrayIn(DataArrayDouble *newArray, DataArrayDouble* &arrayToSet);
    MEDCOUPLING_EXPORT const double *getConstPointer() const { return _mem.getConstPointer(); }
    MEDCOUPLING_EXPORT DataArrayDoubleIterator *iterator();
    MEDCOUPLING_EXPORT const double *begin() const { return getConstPointer(); }
    MEDCOUPLING_EXPORT const double *end() const { return getConstPointer()+getNbOfElems(); }
    MEDCOUPLING_EXPORT void useArray(const double *array, bool ownership, DeallocType type, int nbOfTuple, int nbOfCompo);
    MEDCOUPLING_EXPORT void writeOnPlace(int id, double element0, const double *others, int sizeOfOthers) { _mem.writeOnPlace(id,element0,others,sizeOfOthers); }
    MEDCOUPLING_EXPORT void checkNoNullValues() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getMinMaxPerComponent(double *bounds) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void recenterForMaxPrecision(double eps) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT double getMaxValue(int& tupleId) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT double getMaxValueInArray() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT double getMinValue(int& tupleId) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT double getMinValueInArray() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT double getMaxValue2(DataArrayInt*& tupleIds) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT double getMinValue2(DataArrayInt*& tupleIds) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT double getAverageValue() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT double norm2() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT double normMax() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void accumulate(double *res) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT double accumulate(int compId) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *fromPolarToCart() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *fromCylToCart() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *fromSpherToCart() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *doublyContractedProduct() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *determinant() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *eigenValues() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *eigenVectors() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *inverse() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *trace() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *deviator() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *magnitude() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *maxPerTuple() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void sortPerTuple(bool asc) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void abs() throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void applyLin(double a, double b, int compoId) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void applyLin(double a, double b) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void applyInv(double numerator) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *negate() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *applyFunc(int nbOfComp, FunctionToEvaluate func) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *applyFunc(int nbOfComp, const char *func) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *applyFunc(const char *func) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *applyFunc2(int nbOfComp, const char *func) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *applyFunc3(int nbOfComp, const std::vector<std::string>& varsOrder, const char *func) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void applyFuncFast32(const char *func) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void applyFuncFast64(const char *func) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *getIdsInRange(double vmin, double vmax) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayDouble *Aggregate(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayDouble *Aggregate(const std::vector<const DataArrayDouble *>& arr) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayDouble *Meld(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayDouble *Meld(const std::vector<const DataArrayDouble *>& arr) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayDouble *Dot(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayDouble *CrossProduct(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayDouble *Max(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayDouble *Min(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayDouble *Add(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void addEqual(const DataArrayDouble *other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayDouble *Substract(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void substractEqual(const DataArrayDouble *other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayDouble *Multiply(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void multiplyEqual(const DataArrayDouble *other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayDouble *Divide(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void divideEqual(const DataArrayDouble *other) throw(INTERP_KERNEL::Exception);
    //! nothing to do here because this class does not aggregate any TimeLabel instance.
    MEDCOUPLING_EXPORT void updateTime() const { }
  public:
    MEDCOUPLING_EXPORT void getTinySerializationIntInformation(std::vector<int>& tinyInfo) const;
    MEDCOUPLING_EXPORT void getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const;
    MEDCOUPLING_EXPORT bool resizeForUnserialization(const std::vector<int>& tinyInfoI);
    MEDCOUPLING_EXPORT void finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<std::string>& tinyInfoS);
  protected:
    template<int SPACEDIM>
    void findCommonTuplesAlg(std::vector<double>& bbox,
                             int nbNodes, int limitNodeId, double prec, std::vector<int>& c, std::vector<int>& cI) const;
  private:
    DataArrayDouble() { }
  private:
    MemArray<double> _mem;
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
    double doubleValue() const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *buildDADouble(int nbOfTuples, int nbOfCompo) const throw(INTERP_KERNEL::Exception);
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
    MEDCOUPLING_EXPORT void checkAllocated() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setInfoAndChangeNbOfCompo(const std::vector<std::string>& info) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT int intValue() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT int getHashCode() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool empty() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *deepCpy() const;
    MEDCOUPLING_EXPORT DataArrayInt *performCpy(bool deepCpy) const;
    MEDCOUPLING_EXPORT void cpyFrom(const DataArrayInt& other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void alloc(int nbOfTuple, int nbOfCompo) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void allocIfNecessary(int nbOfTuple, int nbOfCompo);
    MEDCOUPLING_EXPORT bool isEqual(const DataArrayInt& other) const;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const DataArrayInt& other) const;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStrAndOrder(const DataArrayInt& other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *buildPermutationArr(const DataArrayInt& other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void sort(bool asc=true) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void reverse() throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void fillWithZero() throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void fillWithValue(int val) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void iota(int init=0) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT std::string repr() const;
    MEDCOUPLING_EXPORT std::string reprZip() const;
    MEDCOUPLING_EXPORT void writeVTK(std::ostream& ofs, int indent, const char *type, const char *nameInFile) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void reprStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprZipStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprWithoutNameStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprZipWithoutNameStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void transformWithIndArr(const int *indArrBg, const int *indArrEnd) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *transformWithIndArrR(const int *indArrBg, const int *indArrEnd) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void splitByValueRange(const int *arrBg, const int *arrEnd,
                                              DataArrayInt *& castArr, DataArrayInt *& rankInsideCast, DataArrayInt *& castsPresent) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *invertArrayO2N2N2O(int newNbOfElem) const;
    MEDCOUPLING_EXPORT DataArrayInt *invertArrayN2O2O2N(int oldNbOfElem) const;
    //!alloc or useArray should have been called before.
    MEDCOUPLING_EXPORT void reAlloc(int nbOfTuples) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *convertToDblArr() const;
    MEDCOUPLING_EXPORT DataArrayInt *fromNoInterlace() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *toNoInterlace() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void renumberInPlace(const int *old2New);
    MEDCOUPLING_EXPORT void renumberInPlaceR(const int *new2Old);
    MEDCOUPLING_EXPORT DataArrayInt *renumber(const int *old2New) const;
    MEDCOUPLING_EXPORT DataArrayInt *renumberR(const int *new2Old) const;
    MEDCOUPLING_EXPORT DataArrayInt *renumberAndReduce(const int *old2NewBg, int newNbOfTuple) const;
    MEDCOUPLING_EXPORT DataArrayInt *selectByTupleId(const int *new2OldBg, const int *new2OldEnd) const;
    MEDCOUPLING_EXPORT DataArrayInt *selectByTupleIdSafe(const int *new2OldBg, const int *new2OldEnd) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *selectByTupleId2(int bg, int end, int step) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *selectByTupleRanges(const std::vector<std::pair<int,int> >& ranges) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *checkAndPreparePermutation() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void changeSurjectiveFormat(int targetNb, DataArrayInt *&arr, DataArrayInt *&arrI) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayInt *BuildOld2NewArrayFromSurjectiveFormat2(int nbOfOldTuples, const DataArrayInt *arr, const DataArrayInt *arrI, int &newNbOfTuples) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *buildPermArrPerLevel() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool isIdentity() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool isUniform(int val) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *substr(int tupleIdBg, int tupleIdEnd=-1) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void rearrange(int newNbOfCompo) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void transpose() throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *changeNbOfComponents(int newNbOfComp, int dftValue) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *keepSelectedComponents(const std::vector<int>& compoIds) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void meldWith(const DataArrayInt *other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setSelectedComponents(const DataArrayInt *a, const std::vector<int>& compoIds) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setPartOfValues1(const DataArrayInt *a, int bgTuples, int endTuples, int stepTuples, int bgComp, int endComp, int stepComp, bool strictCompoCompare=true) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setPartOfValuesSimple1(int a, int bgTuples, int endTuples, int stepTuples, int bgComp, int endComp, int stepComp) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setPartOfValues2(const DataArrayInt *a, const int *bgTuples, const int *endTuples, const int *bgComp, const int *endComp, bool strictCompoCompare=true) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setPartOfValuesSimple2(int a, const int *bgTuples, const int *endTuples, const int *bgComp, const int *endComp) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setPartOfValues3(const DataArrayInt *a, const int *bgTuples, const int *endTuples, int bgComp, int endComp, int stepComp, bool strictCompoCompare=true) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setPartOfValuesSimple3(int a, const int *bgTuples, const int *endTuples, int bgComp, int endComp, int stepComp) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setPartOfValuesAdv(const DataArrayInt *a, const DataArrayInt *tuplesSelec) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setContigPartOfSelectedValues(int tupleIdStart, const DataArrayInt*a, const DataArrayInt *tuplesSelec) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setContigPartOfSelectedValues2(int tupleIdStart, const DataArrayInt *a, int bg, int end2, int step) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getTuple(int tupleId, int *res) const { std::copy(_mem.getConstPointerLoc(tupleId*((int)_info_on_compo.size())),_mem.getConstPointerLoc((tupleId+1)*((int)_info_on_compo.size())),res); }
    MEDCOUPLING_EXPORT int getIJ(int tupleId, int compoId) const { return _mem[tupleId*((int)_info_on_compo.size())+compoId]; }
    MEDCOUPLING_EXPORT int getIJSafe(int tupleId, int compoId) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT int back() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setIJ(int tupleId, int compoId, int newVal) { _mem[tupleId*((int)_info_on_compo.size())+compoId]=newVal; declareAsNew(); }
    MEDCOUPLING_EXPORT void setIJSilent(int tupleId, int compoId, int newVal) { _mem[tupleId*((int)_info_on_compo.size())+compoId]=newVal; }
    MEDCOUPLING_EXPORT int *getPointer() { return _mem.getPointer(); }
    MEDCOUPLING_EXPORT static void SetArrayIn(DataArrayInt *newArray, DataArrayInt* &arrayToSet);
    MEDCOUPLING_EXPORT const int *getConstPointer() const { return _mem.getConstPointer(); }
    MEDCOUPLING_EXPORT DataArrayIntIterator *iterator();
    MEDCOUPLING_EXPORT const int *begin() const { return getConstPointer(); }
    MEDCOUPLING_EXPORT const int *end() const { return getConstPointer()+getNbOfElems(); }
    MEDCOUPLING_EXPORT DataArrayInt *getIdsEqual(int val) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *getIdsNotEqual(int val) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *getIdsEqualList(const int *valsBg, const int *valsEnd) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *getIdsNotEqualList(const int *valsBg, const int *valsEnd) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT int changeValue(int oldValue, int newValue) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT int locateTuple(const std::vector<int>& tupl) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT int locateValue(int value) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT int locateValue(const std::vector<int>& vals) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT int search(const std::vector<int>& vals) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool presenceOfTuple(const std::vector<int>& tupl) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool presenceOfValue(int value) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool presenceOfValue(const std::vector<int>& vals) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT int getMaxValue(int& tupleId) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT int getMaxValueInArray() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT int getMinValue(int& tupleId) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT int getMinValueInArray() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void abs() throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void applyLin(int a, int b, int compoId) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void applyLin(int a, int b) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void applyInv(int numerator) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *negate() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void applyDivideBy(int val) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void applyModulus(int val) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void applyRModulus(int val) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayInt *Aggregate(const DataArrayInt *a1, const DataArrayInt *a2, int offsetA2);
    MEDCOUPLING_EXPORT static DataArrayInt *Aggregate(const std::vector<const DataArrayInt *>& arr) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayInt *Meld(const DataArrayInt *a1, const DataArrayInt *a2) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayInt *Meld(const std::vector<const DataArrayInt *>& arr) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayInt *MakePartition(const std::vector<const DataArrayInt *>& groups, int newNb, std::vector< std::vector<int> >& fidsOfGroups);
    MEDCOUPLING_EXPORT static DataArrayInt *BuildUnion(const std::vector<const DataArrayInt *>& arr) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayInt *BuildIntersection(const std::vector<const DataArrayInt *>& arr) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *buildComplement(int nbOfElement) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *buildSubstraction(const DataArrayInt *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *buildUnion(const DataArrayInt *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *buildIntersection(const DataArrayInt *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *deltaShiftIndex() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void computeOffsets() throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void computeOffsets2() throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *buildExplicitArrByRanges(const DataArrayInt *offsets) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT std::set<int> getDifferentValues() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void useArray(const int *array, bool ownership, DeallocType type, int nbOfTuple, int nbOfCompo);
    MEDCOUPLING_EXPORT void writeOnPlace(int id, int element0, const int *others, int sizeOfOthers) { _mem.writeOnPlace(id,element0,others,sizeOfOthers); }
    MEDCOUPLING_EXPORT static DataArrayInt *Add(const DataArrayInt *a1, const DataArrayInt *a2) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void addEqual(const DataArrayInt *other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayInt *Substract(const DataArrayInt *a1, const DataArrayInt *a2) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void substractEqual(const DataArrayInt *other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayInt *Multiply(const DataArrayInt *a1, const DataArrayInt *a2) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void multiplyEqual(const DataArrayInt *other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayInt *Divide(const DataArrayInt *a1, const DataArrayInt *a2) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void divideEqual(const DataArrayInt *other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayInt *Modulus(const DataArrayInt *a1, const DataArrayInt *a2) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void modulusEqual(const DataArrayInt *other) throw(INTERP_KERNEL::Exception);
    //! nothing to do here because this class does not aggregate any TimeLabel instance.
    MEDCOUPLING_EXPORT void updateTime() const { }
  public:
    MEDCOUPLING_EXPORT static int *CheckAndPreparePermutation(const int *start, const int *end);
    MEDCOUPLING_EXPORT static DataArrayInt *Range(int begin, int end, int step) throw(INTERP_KERNEL::Exception);
  public:
    MEDCOUPLING_EXPORT void getTinySerializationIntInformation(std::vector<int>& tinyInfo) const;
    MEDCOUPLING_EXPORT void getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const;
    MEDCOUPLING_EXPORT bool resizeForUnserialization(const std::vector<int>& tinyInfoI);
    MEDCOUPLING_EXPORT void finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<std::string>& tinyInfoS);
  private:
    DataArrayInt() { }
  private:
    MemArray<int> _mem;
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
    int intValue() const throw(INTERP_KERNEL::Exception);
    DataArrayInt *buildDAInt(int nbOfTuples, int nbOfCompo) const throw(INTERP_KERNEL::Exception);
  private:
    int *_pt;
    int _nb_of_compo;
  };
}

#endif
