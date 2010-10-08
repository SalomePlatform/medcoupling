//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
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

#ifdef ENABLE_PARMETIS
// include parmetis.h even if it is not needed here
// to avoid inclusion of c++ definitions within extern "C"
// from metis.h from parmetis.h from mpi.h(openmpi) from mpicxx.h
#include <parmetis.h>
#endif
extern "C"
{
#include "metis.h"
}

#include "RENUMBER_METISRenumbering.hxx"

void METISRenumbering::renumber(const int* graph,const int* index_graph,int nb_cell,std::vector<int>& iperm,std::vector<int>& perm)
{
  iperm.resize(nb_cell,0);
  perm.resize(nb_cell,0);
  int num_flag=1;
  int options=0;
  METIS_NodeND(&nb_cell,(int*)index_graph,(int*)graph,&num_flag,&options,&iperm[0],&perm[0]);
}
