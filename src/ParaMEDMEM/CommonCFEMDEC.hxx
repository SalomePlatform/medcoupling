// Copyright (C) 2026  CEA, EDF
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

#pragma once

#include "InterpolationOptions.hxx"

#include "MPIProcessorGroup.hxx"

#include "MEDCouplingUMesh.hxx"

#include "BBTreeClosest.txx"

#include "MPITraits.hxx"

#include <vector>

namespace MEDCoupling
{
/*!
 * EDF35712
 */
template <int spaceDim>
BBTreeClosestSafe<spaceDim, mcIdType>
ShareBBTreesOfAllProcs(
    MPIProcessorGroup *unionGrp, const MEDCouplingUMesh *mesh, std::vector<BBTreeClosestSafe<spaceDim, mcIdType>> &ret
);

template <class T>
std::vector<typename std::vector<T>>
all2allVector(MPIProcessorGroup *grp, const std::vector<typename std::vector<T>> &input);

template <class T>
std::vector<MCAuto<typename MPITraits<T>::ArrayType>>
all2allDA(MPIProcessorGroup *grp, const std::vector<MCAuto<typename MPITraits<T>::ArrayType>> &input);

std::vector<MCAuto<DataArrayIdType>>
ComputeNodeIdsPerProc(
    const std::vector<MCAuto<MEDCouplingUMesh>> &srcMeshes,
    const std::vector<MCAuto<DataArrayIdType>> &srcGlobalNodeIds,
    MCAuto<MEDCouplingUMesh> &wholeMesh
);

MCAuto<MEDCouplingUMesh>
ReduceMesh(
    const MEDCouplingUMesh *mesh,
    const DataArrayIdType *globalNodeIds,
    const mcIdType *bg,
    const mcIdType *end,
    MCAuto<DataArrayIdType> &globalNodeIdsOut
);

MCAuto<DataArrayDouble>
MatrixVectorMultiply(const std::vector<std::map<mcIdType, double>> &matrix, const DataArrayDouble *vectorArr);
}  // namespace MEDCoupling
