// Copyright (C) 2018-2019  CEA/DEN, EDF R&D
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

#ifndef __INTERPKERNELASSERT_HXX__
#define __INTERPKERNELASSERT_HXX__

#include "InterpKernelException.hxx"

#include <sstream>

#define IKAssert(a) { bool verdict(a);           \
    if(!verdict) { std::ostringstream osszz; osszz << "Assertion \"" << #a << "\" failed into " << __FILE__ << " at line " << __LINE__ << " !"; throw INTERP_KERNEL::Exception(osszz.str()); } }

#define IKAssertMsg(a,b) { bool verdict(a);           \
    if(!verdict) { std::ostringstream osszz; osszz << "Assertion \"" << #a << "\" failed into " << __FILE__ << " at line " << __LINE__ << " with message \"" << b << "\" !"; throw INTERP_KERNEL::Exception(osszz.str()); } }

#endif
