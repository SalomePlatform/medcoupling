// Copyright (C) 2016  CEA/DEN, EDF R&D
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

#include "MEDCouplingTraits.hxx"

using namespace MEDCoupling;

const char Traits<double>::ArrayTypeName[]="DataArrayDouble";

const char Traits<double>::FieldTypeName[]="MEDCouplingFieldDouble";

const char Traits<double>::NPYStr[]="FLOAT64";

const char Traits<double>::ReprStr[]="double";

const char Traits<float>::ArrayTypeName[]="DataArrayFloat";

const char Traits<float>::FieldTypeName[]="MEDCouplingFieldFloat";

const char Traits<float>::NPYStr[]="FLOAT32";

const char Traits<float>::ReprStr[]="float";

const char Traits<int>::ArrayTypeName[]="DataArrayInt";

const char Traits<int>::FieldTypeName[]="MEDCouplingFieldInt";

const char Traits<int>::NPYStr[]="INT32";

const char Traits<int>::ReprStr[]="int";

const char Traits<char>::ArrayTypeName[]="DataArrayChar";

const char Traits<Int64>::ArrayTypeName[]="DataArrayInt64";

const char Traits<Int64>::FieldTypeName[]="MEDCouplingFieldInt64";

const char Traits<Int64>::NPYStr[]="INT64";

const char Traits<Int64>::ReprStr[]="int64";
