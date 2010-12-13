//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

#ifndef __PARAMEDMEM_MEDCOUPLINGMEMARRAY_HXX__
#define __PARAMEDMEM_MEDCOUPLINGMEMARRAY_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingTimeLabel.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "InterpKernelException.hxx"

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
    T *getPointer() const { if(_internal) return _internal; if(_external) throw INTERP_KERNEL::Exception("Trying to write on an external pointer."); else return 0; }
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
    T *getPointer() const { return _pointer.getPointer(); }
    MemArray<T> &operator=(const MemArray<T>& other);
    T operator[](int id) const { return _pointer.getConstPointer()[id]; }
    T& operator[](int id) { return _pointer.getPointer()[id]; }
    bool isEqual(const MemArray<T>& other, T prec) const;
    void repr(int sl, std::ostream& stream) const;
    void reprZip(int sl, std::ostream& stream) const;
    void fillWithValue(const T& val);
    T *fromNoInterlace(int nbOfComp) const;
    T *toNoInterlace(int nbOfComp) const;
    void alloc(int nbOfElements);
    void reAlloc(int newNbOfElements);
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
    MEDCOUPLING_EXPORT const std::vector<std::string> &getInfoOnComponent() const { return _info_on_compo; }
    MEDCOUPLING_EXPORT std::string getInfoOnComponent(int i) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setInfoOnComponent(int i, const char *info) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT int getNumberOfComponents() const { return _info_on_compo.size(); }
    MEDCOUPLING_EXPORT int getNumberOfTuples() const { return _nb_of_tuples; }
    MEDCOUPLING_EXPORT int getNbOfElems() const { return _info_on_compo.size()*_nb_of_tuples; }
  protected:
    DataArray():_nb_of_tuples(-1) { }
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
  class DataArrayDouble : public DataArray
  {
  public:
    MEDCOUPLING_EXPORT static DataArrayDouble *New();
    MEDCOUPLING_EXPORT bool isAllocated() const;
    MEDCOUPLING_EXPORT void checkAllocated() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *deepCpy() const;
    MEDCOUPLING_EXPORT DataArrayDouble *performCpy(bool deepCpy) const;
    MEDCOUPLING_EXPORT void alloc(int nbOfTuple, int nbOfCompo);
    MEDCOUPLING_EXPORT void fillWithZero() throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void fillWithValue(double val) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void iota(double init=0.) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool isUniform(double val, double eps) const;
    MEDCOUPLING_EXPORT std::string repr() const;
    MEDCOUPLING_EXPORT std::string reprZip() const;
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
    MEDCOUPLING_EXPORT DataArrayDouble *substr(int tupleIdBg, int tupleIdEnd=-1) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *changeNbOfComponents(int newNbOfComp, double dftValue) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *keepSelectedComponents(const std::vector<int>& compoIds) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void meldWith(const DataArrayDouble *other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setSelectedComponents(const DataArrayDouble *a, const std::vector<int>& compoIds) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getTuple(int tupleId, double *res) const { std::copy(_mem.getConstPointerLoc(tupleId*_info_on_compo.size()),_mem.getConstPointerLoc((tupleId+1)*_info_on_compo.size()),res); }
    MEDCOUPLING_EXPORT double getIJ(int tupleId, int compoId) const { return _mem[tupleId*_info_on_compo.size()+compoId]; }
    MEDCOUPLING_EXPORT void setIJ(int tupleId, int compoId, double newVal) { _mem[tupleId*_info_on_compo.size()+compoId]=newVal; declareAsNew(); }
    MEDCOUPLING_EXPORT double *getPointer() const { return _mem.getPointer(); }
    MEDCOUPLING_EXPORT static void setArrayIn(DataArrayDouble *newArray, DataArrayDouble* &arrayToSet);
    MEDCOUPLING_EXPORT const double *getConstPointer() const { return _mem.getConstPointer(); }
    MEDCOUPLING_EXPORT void useArray(const double *array, bool ownership, DeallocType type, int nbOfTuple, int nbOfCompo);
    MEDCOUPLING_EXPORT void writeOnPlace(int id, double element0, const double *others, int sizeOfOthers) { _mem.writeOnPlace(id,element0,others,sizeOfOthers); }
    MEDCOUPLING_EXPORT void checkNoNullValues() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT double getMaxValue(int& tupleId) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT double getMinValue(int& tupleId) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT double getMaxValue2(DataArrayInt*& tupleIds) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT double getMinValue2(DataArrayInt*& tupleIds) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT double getAverageValue() const throw(INTERP_KERNEL::Exception);
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
    MEDCOUPLING_EXPORT void applyLin(double a, double b, int compoId) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *applyFunc(int nbOfComp, FunctionToEvaluate func) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *applyFunc(int nbOfComp, const char *func) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayDouble *applyFunc(const char *func) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void applyFuncFast32(const char *func) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void applyFuncFast64(const char *func) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *getIdsInRange(double vmin, double vmax) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayDouble *aggregate(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayDouble *aggregate(const std::vector<const DataArrayDouble *>& a) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayDouble *meld(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayDouble *meld(const std::vector<const DataArrayDouble *>& a) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayDouble *dot(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayDouble *crossProduct(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayDouble *max(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayDouble *min(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayDouble *add(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void addEqual(const DataArrayDouble *other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayDouble *substract(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void substractEqual(const DataArrayDouble *other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayDouble *multiply(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void multiplyEqual(const DataArrayDouble *other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayDouble *divide(const DataArrayDouble *a1, const DataArrayDouble *a2) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void divideEqual(const DataArrayDouble *other) throw(INTERP_KERNEL::Exception);
    //! nothing to do here because this class does not aggregate any TimeLabel instance.
    MEDCOUPLING_EXPORT void updateTime() { }
  private:
    DataArrayDouble() { }
  private:
    MemArray<double> _mem;
  };

  class DataArrayInt : public DataArray
  {
  public:
    MEDCOUPLING_EXPORT static DataArrayInt *New();
    MEDCOUPLING_EXPORT bool isAllocated() const;
    MEDCOUPLING_EXPORT void checkAllocated() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *deepCpy() const;
    MEDCOUPLING_EXPORT DataArrayInt *performCpy(bool deepCpy) const;
    MEDCOUPLING_EXPORT void alloc(int nbOfTuple, int nbOfCompo);
    MEDCOUPLING_EXPORT bool isEqual(const DataArrayInt& other) const;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const DataArrayInt& other) const;
    MEDCOUPLING_EXPORT void fillWithZero();
    MEDCOUPLING_EXPORT void fillWithValue(int val);
    MEDCOUPLING_EXPORT void iota(int init=0) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT std::string repr() const;
    MEDCOUPLING_EXPORT std::string reprZip() const;
    MEDCOUPLING_EXPORT void reprStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprZipStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprWithoutNameStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void reprZipWithoutNameStream(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void transformWithIndArr(const int *indArr);
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
    MEDCOUPLING_EXPORT bool isIdentity() const;
    MEDCOUPLING_EXPORT bool isUniform(int val) const;
    MEDCOUPLING_EXPORT DataArrayInt *substr(int tupleIdBg, int tupleIdEnd=-1) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *changeNbOfComponents(int newNbOfComp, int dftValue) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *keepSelectedComponents(const std::vector<int>& compoIds) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void meldWith(const DataArrayInt *other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void setSelectedComponents(const DataArrayInt *a, const std::vector<int>& compoIds) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void getTuple(int tupleId, int *res) const { std::copy(_mem.getConstPointerLoc(tupleId*_info_on_compo.size()),_mem.getConstPointerLoc((tupleId+1)*_info_on_compo.size()),res); }
    MEDCOUPLING_EXPORT int getIJ(int tupleId, int compoId) const { return _mem[tupleId*_info_on_compo.size()+compoId]; }
    MEDCOUPLING_EXPORT void setIJ(int tupleId, int compoId, int newVal) { _mem[tupleId*_info_on_compo.size()+compoId]=newVal; }
    MEDCOUPLING_EXPORT int *getPointer() const { return _mem.getPointer(); }
    MEDCOUPLING_EXPORT static void setArrayIn(DataArrayInt *newArray, DataArrayInt* &arrayToSet);
    MEDCOUPLING_EXPORT const int *getConstPointer() const { return _mem.getConstPointer(); }
    MEDCOUPLING_EXPORT DataArrayInt *getIdsEqual(int val) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *getIdsEqualList(const std::vector<int>& vals) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT int getMaxValue(int& tupleId) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT int getMinValue(int& tupleId) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayInt *aggregate(const DataArrayInt *a1, const DataArrayInt *a2, int offsetA2);
    MEDCOUPLING_EXPORT static DataArrayInt *meld(const DataArrayInt *a1, const DataArrayInt *a2) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayInt *meld(const std::vector<const DataArrayInt *>& a) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayInt *makePartition(const std::vector<DataArrayInt *>& groups, int newNb, std::vector< std::vector<int> >& fidsOfGroups);
    MEDCOUPLING_EXPORT DataArrayInt *buildComplement(int nbOfElement) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *buildSubstraction(const DataArrayInt *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *buildUnion(const DataArrayInt *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *buildIntersection(const DataArrayInt *other) const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT DataArrayInt *deltaShiftIndex() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT void useArray(const int *array, bool ownership, DeallocType type, int nbOfTuple, int nbOfCompo);
    MEDCOUPLING_EXPORT void writeOnPlace(int id, int element0, const int *others, int sizeOfOthers) { _mem.writeOnPlace(id,element0,others,sizeOfOthers); }
    //! nothing to do here because this class does not aggregate any TimeLabel instance.
    MEDCOUPLING_EXPORT void updateTime() { }
  public:
    MEDCOUPLING_EXPORT static int *checkAndPreparePermutation(const int *start, const int *end);
  private:
    DataArrayInt() { }
  private:
    MemArray<int> _mem;
  };
}

#endif
