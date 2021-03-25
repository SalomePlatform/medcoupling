// Copyright (C) 2016-2021  CEA/DEN, EDF R&D
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

#pragma once

#include "MEDCoupling.hxx"
#include "MCType.hxx"
#include <vector>

namespace MEDCoupling
{
  template<class T>
  struct MEDCOUPLING_EXPORT Traits
  {
    using EltType = T;
  };

  class DataArrayInt32;
  class DataArrayFloat;
  class DataArrayDouble;
  class DataArrayChar;
  class DataArrayByte;
  class DataArrayInt64;
  class MEDCouplingFieldDouble;
  class MEDCouplingFieldFloat;
  class MEDCouplingFieldInt32;
  class MEDCouplingFieldInt64;
  class DataArrayInt32Tuple;
  class DataArrayInt64Tuple;
  class DataArrayFloatTuple;
  class DataArrayDoubleTuple;
  class DataArrayByteTuple;
  class DataArrayInt32Iterator;
  class DataArrayInt64Iterator;
  class DataArrayByteIterator;
  
  template<>
  struct MEDCOUPLING_EXPORT Traits<double>
  {
    static const char ArrayTypeName[];
    static const char FieldTypeName[];
    static const char NPYStr[];
    static const char ReprStr[];
    using ArrayType = DataArrayDouble;
    using ArrayTypeCh = DataArrayDouble;
    using FieldType = MEDCouplingFieldDouble;
    using ArrayTuple = DataArrayDoubleTuple;
  };

  template<>
  struct MEDCOUPLING_EXPORT Traits<float>
  {
    static const char ArrayTypeName[];
    static const char FieldTypeName[];
    static const char NPYStr[];
    static const char ReprStr[];
    using ArrayType = DataArrayFloat;
    using ArrayTypeCh = DataArrayFloat;
    using FieldType = MEDCouplingFieldFloat;
    using ArrayTuple = DataArrayFloatTuple;
  };
  
  template<>
  struct MEDCOUPLING_EXPORT Traits<Int32>
  {
    static const char ArrayTypeName[];
    static const char FieldTypeName[];
    static const char NPYStr[];
    static const char ReprStr[];
    static const char VTKReprStr[];
    using ArrayType = DataArrayInt32;
    using ArrayTypeCh = DataArrayInt32;
    using FieldType = MEDCouplingFieldInt32;
    using ArrayTuple = DataArrayInt32Tuple;
    using IteratorType = DataArrayInt32Iterator;
  };

  template<>
  struct MEDCOUPLING_EXPORT Traits<Int64>
  {
    static const char ArrayTypeName[];
    static const char FieldTypeName[];
    static const char NPYStr[];
    static const char ReprStr[];
    static const char VTKReprStr[];
    using ArrayType = DataArrayInt64;
    using ArrayTypeCh = DataArrayInt64;
    using FieldType = MEDCouplingFieldInt64;
    using ArrayTuple = DataArrayInt64Tuple;
    using IteratorType = DataArrayInt64Iterator;
  };

  template<>
  struct MEDCOUPLING_EXPORT Traits<char>
  {
    static const char ArrayTypeName[];
    using ArrayTypeCh = DataArrayByte;
    using ArrayType = DataArrayChar;
    using ArrayTuple = DataArrayByteTuple;
    using IteratorType = DataArrayByteIterator;
  };
}
