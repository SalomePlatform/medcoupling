// Copyright (C) 2007-2016  CEA/DEN, EDF R&D
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

#ifndef __MEDFILESAFECALLER_TXX__
#define __MEDFILESAFECALLER_TXX__

#include "med.h"

#include "InterpKernelException.hxx"

#include <sstream>

// TXX extension to avoid to be installed.

// This macro calls safely MEDFile functions returning 0
#define MEDFILESAFECALLER0(a,b) { med_err retCode(a b); \
    if(retCode!=0) { std::ostringstream osszz; osszz << "Return code of MEDFile call \"" << #a << "\" is not 0 as expected ! ( Return code was "<< retCode << " at " << __FILE__ << ":" << __LINE__ << " )"; throw INTERP_KERNEL::Exception(osszz.str().c_str()); } }

#define MEDFILESAFECALLERRD0(a,b) MEDFILESAFECALLER0(a,b)

#define MEDFILESAFECALLERWR0(a,b) { med_err retCode(a b); \
    if(retCode!=0) { std::ostringstream osszz; osszz << "Return code of MEDFile call \"" << #a << "\" is not 0 as expected during writing operation ! ( Return code was "<< retCode << " at " << __FILE__ << ":" << __LINE__ << " ). Check write access on MED file ?"; throw INTERP_KERNEL::Exception(osszz.str().c_str()); } }

#endif
