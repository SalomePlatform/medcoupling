// Copyright (C) 2025-2026  CEA, EDF
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

#include "MPITraits.hxx"

namespace MEDCoupling
{
template <>
const MPI_Datatype MPITraits<char>::MPIType = MPI_CHAR;
const MPI_Datatype MPITraits<std::int64_t>::MPIType = MPI_INT64_T;
const MPI_Datatype MPITraits<double>::MPIType = MPI_DOUBLE;
}  // namespace MEDCoupling
