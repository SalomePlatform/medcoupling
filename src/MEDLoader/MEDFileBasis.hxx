// Copyright (C) 2007-2024  CEA, EDF
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
#include "MEDCouplingMemArray.txx"

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
#ifdef MED_INT_IS_LONG
  using DataArrayMedInt = DataArrayInt64;
#else
  using DataArrayMedInt = DataArrayInt32;
#endif

  template<class INTARRAY>
  DataArrayMedInt * DataArrayMedInt_Copy( const INTARRAY* intArray )
  {
    DataArrayMedInt* medIntArray = DataArrayMedInt::New();
    medIntArray->alloc( intArray->getNumberOfTuples(), intArray->getNumberOfComponents() );
    medIntArray->copyStringInfoFrom( *intArray );
    std::copy( intArray->begin(), intArray->end(), medIntArray->getPointer() );
    return medIntArray;
  }

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
    return DataArrayMedInt_Copy((const INTARRAY*) intArray );
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
    return DataArrayMedInt_Copy( intArray );
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

}

#endif
