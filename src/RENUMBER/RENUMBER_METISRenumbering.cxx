// Copyright (C) 2007-2023  CEA/DEN, EDF R&D
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

#ifdef MEDCOUPLING_USE_64BIT_IDS
#define ID_TYPE_SIZE 64
#else
#define ID_TYPE_SIZE 32
#endif

void METISRenumbering::renumber(const mcIdType *graph, const mcIdType *index_graph, mcIdType nbCell, MEDCoupling::DataArrayIdType *&iperm, MEDCoupling::DataArrayIdType *&perm)
{
  MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> out0(MEDCoupling::DataArrayIdType::New()),out1(MEDCoupling::DataArrayIdType::New());
  out0->alloc(nbCell,1); out1->alloc(nbCell,1);
  out0->fillWithZero(); out1->fillWithZero();
  int num_flag=1;
  int options=0;

#if ID_TYPE_SIZE == IDXTYPEWIDTH

  METIS_NodeND(&nbCell,(idx_t*)index_graph,(idx_t*)graph,&num_flag,&options,out0->getPointer(),out1->getPointer());

#else

  mcIdType indexSize = nbCell + 1, graphSize = index_graph[indexSize];
  std::vector<idx_t> indexVec( index_graph, index_graph + indexSize );
  std::vector<idx_t> graphVec( graph, graph + graphSize );
  std::vector<idx_t> out0Vec( nbCell ), out1Vec( nbCell );
  idx_t nb = static_cast<idx_t>( nbCell );
  METIS_NodeND(&nb,indexVec.data(),graphVec.data(),&num_flag,&options,out0Vec.data(),out1Vec.data());
  std::copy( out0Vec.begin(),out0Vec.end(),out0->getPointer() );
  std::copy( out1Vec.begin(),out1Vec.end(),out1->getPointer() );

#endif

  iperm=out0.retn(); perm=out1.retn();
}
