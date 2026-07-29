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

#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingFieldDiscretizationOnNodesFE.hxx"
#include "MCAuto.hxx"

#include <mpi.h>

#include <set>
#include <map>
#include <memory>
#include <vector>

namespace MEDCoupling
{
/*!
 * EDF35712
 */
class OverlapCFEMDEC : public INTERP_KERNEL::InterpolationOptions
{
   public:
    OverlapCFEMDEC(const std::set<int> &procIds, const MPI_Comm &worldComm = MPI_COMM_WORLD);
    void attachSourceField(MEDCouplingFieldDouble *srcField, DataArrayIdType *srcGlobalNodeIds);
    void attachTargetMesh(MEDCouplingUMesh *trgMesh, DataArrayIdType *trgGlobalNodeIds);
    void synchronize();
    MCAuto<MEDCouplingFieldDouble> computeTargetField();

   private:
    void basicCheck();
    void checkSameSpaceDim();
    template <int spaceDim>
    void synchronizeT(
        std::vector<MCAuto<MEDCouplingUMesh>> &srcMeshes, std::vector<MCAuto<DataArrayIdType>> &srcGlobalNodeIds
    );
    void computeMatrix(
        const std::vector<MCAuto<MEDCouplingUMesh>> &srcMeshes,
        const std::vector<MCAuto<DataArrayIdType>> &srcGlobalNodeIds
    );
    const MEDCouplingUMesh *getSourceLocalMesh() const;

   private:
    MCAuto<MEDCouplingFieldDouble> _src_field;
    MCAuto<DataArrayIdType> _src_global_node_ids;
    MCAuto<MEDCouplingUMesh> _trg_mesh;
    MCAuto<DataArrayIdType> _trg_global_node_ids;

   private:
    mcIdType _nb_nodes_src_mesh = 0;
    std::vector<MCAuto<DataArrayIdType>> _glb_nodes_in_whole_per_trg;
    std::vector<MCAuto<DataArrayIdType>> _src_rank_of_nodes_in_whole;
    std::vector<std::map<mcIdType, double>> _matrix;

   private:
    MPI_Comm _comm;
    std::unique_ptr<MPIProcessorGroup> _group;
};

}  // namespace MEDCoupling
