//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
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

  class MEDCOUPLING_EXPORT DataArray : public RefCountObject, public TimeLabel
  {
  public:
    void setName(const char *name);
    void copyStringInfoFrom(const DataArray& other) throw(INTERP_KERNEL::Exception);
    bool areInfoEquals(const DataArray& other) const;
    std::string getName() const { return _name; }
    std::string getInfoOnComponent(int i) const { return _info_on_compo[i]; }
    void setInfoOnComponent(int i, const char *info) { _info_on_compo[i]=info; }
    int getNumberOfComponents() const { return _info_on_compo.size(); }
    int getNumberOfTuples() const { return _nb_of_tuples; }
    int getNbOfElems() const { return _info_on_compo.size()*_nb_of_tuples; }
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
  class MEDCOUPLING_EXPORT DataArrayDouble : public DataArray
  {
  public:
    static DataArrayDouble *New();
    DataArrayDouble *deepCopy() const;
    DataArrayDouble *performCpy(bool deepCpy) const;
    void alloc(int nbOfTuple, int nbOfCompo);
    bool isEqual(const DataArrayDouble& other, double prec) const;
    //!alloc or useArray should have been called before.
    void reAlloc(int nbOfTuples);
    void getTuple(int tupleId, double *res) const { std::copy(_mem.getConstPointerLoc(tupleId*_info_on_compo.size()),_mem.getConstPointerLoc((tupleId+1)*_info_on_compo.size()),res); }
    double getIJ(int tupleId, int compoId) const { return _mem[tupleId*_info_on_compo.size()+compoId]; }
    void setIJ(int tupleId, int compoId, double newVal) { _mem[tupleId*_info_on_compo.size()+compoId]=newVal; }
    double *getPointer() const { return _mem.getPointer(); }
    static void setArrayIn(DataArrayDouble *newArray, DataArrayDouble* &arrayToSet);
    const double *getConstPointer() const { return _mem.getConstPointer(); }
    void useArray(const double *array, bool ownership, DeallocType type, int nbOfTuple, int nbOfCompo);
    void writeOnPlace(int id, double element0, const double *others, int sizeOfOthers) { _mem.writeOnPlace(id,element0,others,sizeOfOthers); }
    void checkNoNullValues() const throw(INTERP_KERNEL::Exception);
    static DataArrayDouble *aggregate(const DataArrayDouble *a1, const DataArrayDouble *a2);
    static DataArrayDouble *add(const DataArrayDouble *a1, const DataArrayDouble *a2);
    static DataArrayDouble *substract(const DataArrayDouble *a1, const DataArrayDouble *a2);
    static DataArrayDouble *multiply(const DataArrayDouble *a1, const DataArrayDouble *a2);
    static DataArrayDouble *divide(const DataArrayDouble *a1, const DataArrayDouble *a2);
    //! nothing to do here because this class does not aggregate any TimeLabel instance.
    void updateTime() { }
  private:
    DataArrayDouble() { }
  private:
    MemArray<double> _mem;
  };

  class MEDCOUPLING_EXPORT DataArrayInt : public DataArray
  {
  public:
    static DataArrayInt *New();
    DataArrayInt *deepCopy() const;
    DataArrayInt *performCpy(bool deepCpy) const;
    void alloc(int nbOfTuple, int nbOfCompo);
    bool isEqual(const DataArrayInt& other) const;
    //!alloc or useArray should have been called before.
    void reAlloc(int nbOfTuples);
    void getTuple(int tupleId, int *res) const { std::copy(_mem.getConstPointerLoc(tupleId*_info_on_compo.size()),_mem.getConstPointerLoc((tupleId+1)*_info_on_compo.size()),res); }
    int getIJ(int tupleId, int compoId) const { return _mem[tupleId*_info_on_compo.size()+compoId]; }
    void setIJ(int tupleId, int compoId, int newVal) { _mem[tupleId*_info_on_compo.size()+compoId]=newVal; }
    int *getPointer() const { return _mem.getPointer(); }
    static void setArrayIn(DataArrayInt *newArray, DataArrayInt* &arrayToSet);
    const int *getConstPointer() const { return _mem.getConstPointer(); }
    static DataArrayInt *aggregate(const DataArrayInt *a1, const DataArrayInt *a2, int offsetA2);
    void useArray(const int *array, bool ownership, DeallocType type, int nbOfTuple, int nbOfCompo);
    void writeOnPlace(int id, int element0, const int *others, int sizeOfOthers) { _mem.writeOnPlace(id,element0,others,sizeOfOthers); }
    //! nothing to do here because this class does not aggregate any TimeLabel instance.
    void updateTime() { }
  private:
    DataArrayInt() { }
  private:
    MemArray<int> _mem;
  };
}

#endif
