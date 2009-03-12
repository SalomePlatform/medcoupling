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
#ifndef __PARAMEDMEM_MEMARRAY_HXX__
#define __PARAMEDMEM_MEMARRAY_HXX__

#include "MEDCoupling.hxx"
#include "RefCountObject.hxx"

#include <string>
#include <vector>

namespace ParaMEDMEM
{
  template<class T>
  class MemArray
  {
  public:
    MemArray():_nb_of_elem(-1),_ownership(false),_pointer(0),_dealloc(CPP_DEALLOC) { }
    MemArray(const MemArray<T>& other);
    T *getPointer() const { return _pointer; }
    MemArray<T> &operator=(const MemArray<T>& other);
    T operator[](int id) const { return _pointer[id]; }
    T& operator[](int id) { return _pointer[id]; }
    void alloc(int nbOfElements);
    void reAlloc(int newNbOfElements);
    void useArray(void *array, bool ownership, DeallocType type, int nbOfElem);
    void writeOnPlace(int id, T element0, const T *others, int sizeOfOthers);
    ~MemArray() { destroy(); }
  private:
    void destroy();
    static void destroyPointer(T *pt, DeallocType type);
  private:
    int _nb_of_elem;
    bool _ownership;
    T *_pointer;
    DeallocType _dealloc;
  };

  class MEDCOUPLING_EXPORT DataArray : public RefCountObject
  {
  public:
    void setName(const char *name);
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

#include "MemArray.txx"

namespace ParaMEDMEM
{
  class MEDCOUPLING_EXPORT DataArrayDouble : public DataArray
  {
  public:
    static DataArrayDouble *New();
    DataArrayDouble *deepCopy() const;
    DataArrayDouble *performCpy(bool deepCpy) const;
    void alloc(int nbOfTuple, int nbOfCompo);
    bool isEqual(DataArrayDouble *other, double prec) const;
    //!alloc or useArray should have been called before.
    void reAlloc(int nbOfTuples);
    double getIJ(int tupleId, int compoId) const { return _mem[tupleId*_info_on_compo.size()+compoId]; }
    void setIJ(int tupleId, int compoId, double newVal) { _mem[tupleId*_info_on_compo.size()+compoId]=newVal; }
    double *getPointer() const { return _mem.getPointer(); }
    void useArray(double *array, bool ownership, DeallocType type, int nbOfTuple, int nbOfCompo);
    void writeOnPlace(int id, double element0, const double *others, int sizeOfOthers) { _mem.writeOnPlace(id,element0,others,sizeOfOthers); }
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
    //!alloc or useArray should have been called before.
    void reAlloc(int nbOfTuples);
    int getIJ(int tupleId, int compoId) const { return _mem[tupleId*_info_on_compo.size()+compoId]; }
    void setIJ(int tupleId, int compoId, int newVal) { _mem[tupleId*_info_on_compo.size()+compoId]=newVal; }
    int *getPointer() const { return _mem.getPointer(); }
    void useArray(int *array, bool ownership, DeallocType type, int nbOfTuple, int nbOfCompo);
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
