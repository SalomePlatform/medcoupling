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

// Creation of this C code is forced by the following.
//
// In case if Metis is a part of Parmetis V3, extern "C" {#include "metis.h"} causes
// inclusion of C++ code of MPI via parmetis.h <- mpi.h <- mpicxx.h
// that breaks compilation. To workaround this problem we create a wrapping C
// function, inclusion of whose declaration causes no problem.

#include "MEDPARTITIONER_metis.h"

#if defined(MED_ENABLE_METIS) & !defined(MED_ENABLE_PARMETIS)
  #ifndef MED_ENABLE_METIS_V5
    #include "defs.h"
  #endif // MED_ENABLE_METIS_V5
  #include "metis.h"
  #ifdef MED_ENABLE_METIS_V5
    #define idxtype idx_t
  #endif // MED_ENABLE_METIS_V5
#else
  typedef int idxtype;
#endif // defined(MED_ENABLE_METIS) & !defined(MED_ENABLE_PARMETIS)

void MEDPARTITIONER_METIS_PartGraphRecursive(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
                                             idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts,
                                             int *options, int *edgecut, idxtype *part)
{
#if defined(MED_ENABLE_METIS)
  #ifndef MED_ENABLE_METIS_V5
  METIS_PartGraphRecursive(nvtxs, xadj, adjncy, vwgt,
    adjwgt, wgtflag, numflag, nparts,
                           options, edgecut, part);
  #else
  int ncon=1;
  options[METIS_OPTION_NCUTS]=1;
  options[METIS_OPTION_NITER]=1;
  options[METIS_OPTION_UFACTOR]=1;
  METIS_PartGraphRecursive(nvtxs, &ncon, xadj, adjncy, vwgt, 0 /* vsize*/,
                           adjwgt, nparts,/* tpwgts*/ 0,/* ubvec */ 0,
                           options, edgecut, part);
  #endif
#endif
}

void MEDPARTITIONER_METIS_PartGraphKway(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
                                        idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts,
                                        int *options, int *edgecut, idxtype *part)
{
#if defined(MED_ENABLE_METIS)
  #ifndef MED_ENABLE_METIS_V5
  METIS_PartGraphKway(nvtxs, xadj, adjncy, vwgt,
    adjwgt, wgtflag, numflag, nparts,
    options, edgecut, part);
  #else
  int ncon=1;
  options[METIS_OPTION_NCUTS]=1;
  options[METIS_OPTION_NITER]=1;
  options[METIS_OPTION_UFACTOR]=1;
  METIS_PartGraphKway(nvtxs, &ncon, xadj, adjncy, vwgt, 0 /* vsize*/,
                      adjwgt, nparts, 0 , 0 /* ubvec */,
                      options, edgecut, part);
  #endif
#endif
}
