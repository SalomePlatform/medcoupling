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
// Author : Anthony Geay (CEA/DEN)

#ifndef __MEDFILEBASIS_HXX__
#define __MEDFILEBASIS_HXX__

#include "InterpKernelException.hxx"
#include "MEDCouplingMemArray.hxx"

#include <string>
#include <vector>

#include <med.h>

namespace MEDCoupling
{
  class MEDFileString
  {
  public:
    MEDFileString(int maxLgth);
    ~MEDFileString();
    void clear();
    void set(const char *s);
    char *getPointer() { return _content; }
    const char *getReprForWrite() const { return _content; }
    std::string getRepr() const;
  private:
    int _max_lgth;
    char *_content;
  };


  class MEDFileMultiString
  {
  public:
    MEDFileMultiString(int nbOfCompo, int maxLgthPerCompo);
    ~MEDFileMultiString();
    void set(int compoId, const char *s);
    const char *getReprForWrite() const;
    std::vector<std::string> getRepr() const;
    std::string getReprPerComp(int compId) const;
  private:
    int _nb_of_comp;
    int _max_lgth_per_comp;
    char *_content;
  };
}

namespace MEDCoupling
{

  class DataArrayMedInt : public DataArrayDiscreteSigned< med_int >
  {
    friend class DataArrayDiscrete<med_int>;
  public:
    template<class INTARRAY>
    static DataArrayMedInt *Copy( const INTARRAY* array );
    static DataArrayMedInt *New() { return new DataArrayMedInt(); }
    std::string getClassName() const override { return std::string("DataArrayMedInt"); }
    DataArrayMedInt *deepCopy() const { return new DataArrayMedInt(*this); }
    //DataArrayMedInt *buildNewEmptyInstance() const { return new DataArrayMedInt(); }//ko
    DataArray *buildNewEmptyInstance() const { if ( sizeof(med_int)==sizeof(int)) return DataArrayInt32::New(); return DataArrayInt64::New(); }
  public:
    DataArray *selectByTupleId(const mcIdType *new2OldBg, const mcIdType *new2OldEnd) const { return this->mySelectByTupleId(new2OldBg,new2OldEnd); }
    DataArray *selectByTupleId(const DataArrayIdType& di) const { return this->mySelectByTupleId(di); }
    DataArray *selectByTupleIdSafe(const mcIdType *new2OldBg, const mcIdType *new2OldEnd) const { return this->mySelectByTupleIdSafe(new2OldBg,new2OldEnd); }
    DataArray *keepSelectedComponents(const std::vector<std::size_t>& compoIds) const { return this->myKeepSelectedComponents(compoIds); }
    DataArray *selectByTupleIdSafeSlice(mcIdType bg, mcIdType end2, mcIdType step) const { return this->mySelectByTupleIdSafeSlice(bg,end2,step); }
    DataArray *selectByTupleRanges(const std::vector<std::pair<mcIdType,mcIdType> >& ranges) const { return this->mySelectByTupleRanges(ranges); }
  private:
    ~DataArrayMedInt() { }
    DataArrayMedInt() { }
  };

  template< class T1, class T2 >
  MCAuto<T1> StaticCast( const MCAuto< T2 >& array )
  {
    DataArray *src = const_cast< T2* >((const T2*) array );
    T1*        tgt = static_cast<T1*>( src );
    if ( tgt )
      tgt->incrRef();
    return tgt;
  }

  template< class INTARRAY >
  MCAuto<DataArrayMedInt> ToMedIntArray(const MCAuto<INTARRAY>& intArray )
  {
    if ( sizeof( med_int ) == sizeof( typename INTARRAY::Type ))
      return StaticCast< DataArrayMedInt >( intArray );
    return DataArrayMedInt::Copy((const INTARRAY*) intArray );
  }

  template< class INT >
  MCAuto<DataArrayMedInt> ToMedIntArray(const typename MEDCoupling::Traits<INT>::ArrayType* intArray )
  {
    if ( sizeof( med_int ) == sizeof( INT ))
    {
      typedef typename MEDCoupling::Traits<INT>::ArrayType INTARRAY;
      MCAuto< INTARRAY > ia = const_cast< INTARRAY* >( intArray );
      ia->incrRef();
      return StaticCast< DataArrayMedInt >( ia );
    }
    return DataArrayMedInt::Copy( intArray );
  }

  template< class INT >
  MCAuto< typename MEDCoupling::Traits<INT>::ArrayType> FromMedIntArray(MCAuto<DataArrayMedInt>& medIntArray )
  {
    typedef typename MEDCoupling::Traits<INT>::ArrayType INTARRAY;
    if ( sizeof( med_int ) == sizeof( INT ))
      return StaticCast< INTARRAY >( medIntArray );

    INTARRAY* intArray = INTARRAY::New();
    intArray->alloc( medIntArray->getNumberOfTuples(), medIntArray->getNumberOfComponents() );
    intArray->copyStringInfoFrom( * medIntArray.operator->() );
    std::copy( medIntArray->begin(), medIntArray->end(), intArray->getPointer() );
    return intArray;
  }

  template< class INT >
  MCAuto<DataArrayMedInt> ToMedIntArray(const std::vector<INT>& intVec )
  {
    DataArrayMedInt* medIntArray = DataArrayMedInt::New();
    if ( sizeof( med_int ) == sizeof( INT ))
      {
        medIntArray->useArray( reinterpret_cast<const med_int*>(intVec.data()), /*owner=*/false, DeallocType::CPP_DEALLOC, intVec.size(), /*nbComp=*/1 );
      }
    else
      {
        medIntArray->alloc( intVec.size(), 1 );
        std::copy( intVec.begin(), intVec.end(), medIntArray->getPointer() );
      }
    return medIntArray;
  }

  template< class INT >
  med_int ToMedInt( INT i )
  {
    return static_cast< med_int >( i );
  }

  template< class INT >
  INT FromMedInt( med_int mi )
  {
    return static_cast< INT >( mi );
  }

  template<class INTARRAY>
  DataArrayMedInt * DataArrayMedInt::Copy( const INTARRAY* intArray )
  {
    DataArrayMedInt* medIntArray = DataArrayMedInt::New();
    medIntArray->alloc( intArray->getNumberOfTuples(), intArray->getNumberOfComponents() );
    medIntArray->copyStringInfoFrom( *intArray );
    std::copy( intArray->begin(), intArray->end(), medIntArray->getPointer() );
    return medIntArray;
  }
}

#endif
