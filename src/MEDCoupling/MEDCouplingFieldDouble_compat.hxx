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
// Author : Adrien Bruneton (CEA/DEN)

#ifndef __MEDCOUPLINGFIELDDOUBLE_COMPAT_HXX__
#define __MEDCOUPLINGFIELDDOUBLE_COMPAT_HXX__

#include "MEDCoupling.hxx"
#include "InterpKernelException.hxx"

namespace ParaMEDMEM
{
  /*
   * This empy shell is provided for backward compatibility with version 1 of the ICoco interface.
   * If you intend to use ICoco *with* MEDCoupling functionalities, you have to switch to the
   * version 2 of the ICoco interface (Problem_v2.h).
   *
   * The dummy implementation below only allows one to keep on compiling simple ICoco schemes (i.e.
   * the ones not using MEDCouplingFieldDouble) with the version 1 of the interface.
   *
   */
  class MEDCouplingFieldDouble {};
}

#endif
