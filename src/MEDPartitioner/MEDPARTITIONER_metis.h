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


void MEDPARTITIONER_METIS_PartGraphRecursive(int *, int  *, int  *, int  *, int  *, int *, int *, int *, int *, int *, int  *);

void MEDPARTITIONER_METIS_PartGraphKway(int *, int  *, int  *, int  *, int  *, int *, int *, int *, int *, int *, int  *); 
