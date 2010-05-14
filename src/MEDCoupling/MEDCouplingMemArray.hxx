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
    const T *getConstPointerLoc(int offset) const { return _pointer.getConstPointerLoc(offset); }
    const T *getConstPointer() const { return _pointer.getConstPointer(); }
    T *getPointer() const { return _pointer.getPointer(); }
    MemArray<T> &operator=(const MemArray<T>& other);
    T operator[](int id) const { return _pointer.getConstPointer()[id]; }
    T& operator[](int id) { return _pointer.getPointer()[id]; }
    bool isEqual(const MemArray<T>& other, T prec) const;
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
    //T *_pointer;
    DeallocType _dealloc;
  };

  class DataArray : public RefCountObject, public TimeLabel
  {
  public:
    MEDCOUPLING_EXPORT void setName(const char *name);
    MEDCOUPLING_EXPORT void copyStringInfoFrom(const DataArray& other) throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT bool areInfoEquals(const DataArray& other) const;
    MEDCOUPLING_EXPORT std::string getName() const { return _name; }
    MEDCOUPLING_EXPORT std::string getInfoOnComponent(int i) const { return _info_on_compo[i]; }
    MEDCOUPLING_EXPORT void setInfoOnComponent(int i, const char *info) { _info_on_compo[i]=info; }
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
  class DataArrayDouble : public DataArray
  {
  public:
    MEDCOUPLING_EXPORT static DataArrayDouble *New();
    MEDCOUPLING_EXPORT DataArrayDouble *deepCopy() const;
    MEDCOUPLING_EXPORT DataArrayDouble *performCpy(bool deepCpy) const;
    MEDCOUPLING_EXPORT void alloc(int nbOfTuple, int nbOfCompo);
    MEDCOUPLING_EXPORT bool isEqual(const DataArrayDouble& other, double prec) const;
    //!alloc or useArray should have been called before.
    MEDCOUPLING_EXPORT void reAlloc(int nbOfTuples);
    MEDCOUPLING_EXPORT void getTuple(int tupleId, double *res) const { std::copy(_mem.getConstPointerLoc(tupleId*_info_on_compo.size()),_mem.getConstPointerLoc((tupleId+1)*_info_on_compo.size()),res); }
    MEDCOUPLING_EXPORT double getIJ(int tupleId, int compoId) const { return _mem[tupleId*_info_on_compo.size()+compoId]; }
    MEDCOUPLING_EXPORT void setIJ(int tupleId, int compoId, double newVal) { _mem[tupleId*_info_on_compo.size()+compoId]=newVal; }
    MEDCOUPLING_EXPORT double *getPointer() const { return _mem.getPointer(); }
    MEDCOUPLING_EXPORT static void setArrayIn(DataArrayDouble *newArray, DataArrayDouble* &arrayToSet);
    MEDCOUPLING_EXPORT const double *getConstPointer() const { return _mem.getConstPointer(); }
    MEDCOUPLING_EXPORT void useArray(const double *array, bool ownership, DeallocType type, int nbOfTuple, int nbOfCompo);
    MEDCOUPLING_EXPORT void writeOnPlace(int id, double element0, const double *others, int sizeOfOthers) { _mem.writeOnPlace(id,element0,others,sizeOfOthers); }
    MEDCOUPLING_EXPORT void checkNoNullValues() const throw(INTERP_KERNEL::Exception);
    MEDCOUPLING_EXPORT static DataArrayDouble *aggregate(const DataArrayDouble *a1, const DataArrayDouble *a2);
    MEDCOUPLING_EXPORT static DataArrayDouble *add(const DataArrayDouble *a1, const DataArrayDouble *a2);
    MEDCOUPLING_EXPORT void addEqual(const DataArrayDouble *other);
    MEDCOUPLING_EXPORT static DataArrayDouble *substract(const DataArrayDouble *a1, const DataArrayDouble *a2);
    MEDCOUPLING_EXPORT void substractEqual(const DataArrayDouble *other);
    MEDCOUPLING_EXPORT static DataArrayDouble *multiply(const DataArrayDouble *a1, const DataArrayDouble *a2);
    MEDCOUPLING_EXPORT void multiplyEqual(const DataArrayDouble *other);
    MEDCOUPLING_EXPORT static DataArrayDouble *divide(const DataArrayDouble *a1, const DataArrayDouble *a2);
    MEDCOUPLING_EXPORT void divideEqual(const DataArrayDouble *other);
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
    MEDCOUPLING_EXPORT DataArrayInt *deepCopy() const;
    MEDCOUPLING_EXPORT DataArrayInt *performCpy(bool deepCpy) const;
    MEDCOUPLING_EXPORT void alloc(int nbOfTuple, int nbOfCompo);
    MEDCOUPLING_EXPORT bool isEqual(const DataArrayInt& other) const;
    //!alloc or useArray should have been called before.
    MEDCOUPLING_EXPORT void reAlloc(int nbOfTuples);
    MEDCOUPLING_EXPORT void getTuple(int tupleId, int *res) const { std::copy(_mem.getConstPointerLoc(tupleId*_info_on_compo.size()),_mem.getConstPointerLoc((tupleId+1)*_info_on_compo.size()),res); }
    MEDCOUPLING_EXPORT int getIJ(int tupleId, int compoId) const { return _mem[tupleId*_info_on_compo.size()+compoId]; }
    MEDCOUPLING_EXPORT void setIJ(int tupleId, int compoId, int newVal) { _mem[tupleId*_info_on_compo.size()+compoId]=newVal; }
    MEDCOUPLING_EXPORT int *getPointer() const { return _mem.getPointer(); }
    MEDCOUPLING_EXPORT static void setArrayIn(DataArrayInt *newArray, DataArrayInt* &arrayToSet);
    MEDCOUPLING_EXPORT const int *getConstPointer() const { return _mem.getConstPointer(); }
    MEDCOUPLING_EXPORT static DataArrayInt *aggregate(const DataArrayInt *a1, const DataArrayInt *a2, int offsetA2);
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
