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

#ifdef MED_ENABLE_PARMETIS
// include parmetis.h even if it is not needed here
// to avoid inclusion of c++ definitions within extern "C"
// from metis.h from parmetis.h from mpi.h(openmpi) from mpicxx.h
#include <parmetis.h>
#endif
extern "C"
{
#include "metis.h"
}

#include "MEDCouplingMemArray.hxx"
#include "MCAuto.hxx"

#include "RENUMBER_METISRenumbering.hxx"

void METISRenumbering::renumber(const int *graph, const int *index_graph, int nbCell, MEDCoupling::DataArrayInt *&iperm, MEDCoupling::DataArrayInt *&perm)
{
  MEDCoupling::MCAuto<MEDCoupling::DataArrayInt> out0(MEDCoupling::DataArrayInt::New()),out1(MEDCoupling::DataArrayInt::New());
  out0->alloc(nbCell,1); out1->alloc(nbCell,1);
  out0->fillWithZero(); out1->fillWithZero();
  int num_flag=1;
  int options=0;
  METIS_NodeND(&nbCell,(int*)index_graph,(int*)graph,&num_flag,&options,out0->getPointer(),out1->getPointer());
  iperm=out0.retn(); perm=out1.retn();
}
