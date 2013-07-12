// Copyright (C) 2007-2013  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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
  #include "defs.h"
//  #include "struct.h"
  #include "metis.h"
#else
  typedef int idxtype;
#endif

void MEDPARTITIONER_METIS_PartGraphRecursive(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
                                             idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, 
                                             int *options, int *edgecut, idxtype *part)
{
#if defined(MED_ENABLE_METIS)
  METIS_PartGraphRecursive(nvtxs, xadj, adjncy, vwgt, 
                           adjwgt, wgtflag, numflag, nparts, 
                           options, edgecut, part);
#endif
}

void MEDPARTITIONER_METIS_PartGraphKway(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
                                        idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, 
                                        int *options, int *edgecut, idxtype *part)
{
#if defined(MED_ENABLE_METIS)
  METIS_PartGraphKway(nvtxs, xadj, adjncy, vwgt, 
                      adjwgt, wgtflag, numflag, nparts, 
                      options, edgecut, part);
#endif
}
