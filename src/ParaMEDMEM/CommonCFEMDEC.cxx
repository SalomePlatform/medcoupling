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

#include "CommonCFEMDEC.txx"

using MEDCoupling::DataArrayIdType;
using MEDCoupling::MCAuto;
using MEDCoupling::MEDCouplingUMesh;
using MEDCoupling::MPIProcessorGroup;

template BBTreeClosestSafe<1, mcIdType>
MEDCoupling::ShareBBTreesOfAllProcs<1>(
    MPIProcessorGroup *unionGrp, const MEDCouplingUMesh *mesh, std::vector<BBTreeClosestSafe<1, mcIdType>> &ret
);

template BBTreeClosestSafe<2, mcIdType>
MEDCoupling::ShareBBTreesOfAllProcs<2>(
    MPIProcessorGroup *unionGrp, const MEDCouplingUMesh *mesh, std::vector<BBTreeClosestSafe<2, mcIdType>> &ret
);

template BBTreeClosestSafe<3, mcIdType>
MEDCoupling::ShareBBTreesOfAllProcs<3>(
    MPIProcessorGroup *unionGrp, const MEDCouplingUMesh *mesh, std::vector<BBTreeClosestSafe<3, mcIdType>> &ret
);

template std::vector<typename std::vector<mcIdType>>
MEDCoupling::all2allVector<mcIdType>(MPIProcessorGroup *grp, const std::vector<typename std::vector<mcIdType>> &input);

template std::vector<typename std::vector<double>>
MEDCoupling::all2allVector<double>(MPIProcessorGroup *grp, const std::vector<typename std::vector<double>> &input);

template std::vector<MCAuto<typename MPITraits<mcIdType>::ArrayType>>
MEDCoupling::all2allDA<mcIdType>(
    MPIProcessorGroup *grp, const std::vector<MCAuto<typename MPITraits<mcIdType>::ArrayType>> &input
);

template std::vector<MCAuto<typename MPITraits<double>::ArrayType>>
MEDCoupling::all2allDA<double>(
    MPIProcessorGroup *grp, const std::vector<MCAuto<typename MPITraits<double>::ArrayType>> &input
);

MCAuto<MEDCouplingUMesh>
MEDCoupling::ReduceMesh(
    const MEDCouplingUMesh *mesh,
    const DataArrayIdType *globalNodeIds,
    const mcIdType *bg,
    const mcIdType *end,
    MCAuto<DataArrayIdType> &globalNodeIdsOut
)
{
    MCAuto<MEDCouplingUMesh> part(mesh->buildPartOfMySelf(bg, end, true));
    MCAuto<DataArrayIdType> nodeIdsFetched(part->computeFetchedNodeIds());
    MCAuto<DataArrayIdType> o2n(nodeIdsFetched->invertArrayN2O2O2N(mesh->getNumberOfNodes()));
    part->renumberNodes(o2n->begin(), nodeIdsFetched->getNumberOfTuples());
    globalNodeIdsOut = globalNodeIds->selectByTupleIdSafe(nodeIdsFetched->begin(), nodeIdsFetched->end());
    return part;
}

/*!
 * Target side : Computes wholeMesh aggregation of srcMeshes from which duplicated cells are removed. globalNodeIds are
 * used to determine common cells across source processors.
 *
 *  \param [out] wholeMesh without duplication of cells
 *  \return for each source proc node Ids in \a wholeMesh referential of its contribution
 */
std::vector<MCAuto<DataArrayIdType>>
MEDCoupling::ComputeNodeIdsPerProc(
    const std::vector<MCAuto<MEDCouplingUMesh>> &srcMeshes,
    const std::vector<MCAuto<DataArrayIdType>> &srcGlobalNodeIds,
    MCAuto<MEDCouplingUMesh> &wholeMesh
)
{
    std::size_t nbOfSrcProcs(srcGlobalNodeIds.size());
    MCAuto<DataArrayIdType> o2n(DataArrayIdType::Aggregate(FromVecAutoToVecOfConst<DataArrayIdType>(srcGlobalNodeIds)));
    MCAuto<DataArrayIdType> b(o2n->buildUniqueNotSorted());
    MCAuto<MapKeyVal<mcIdType, mcIdType>> zeMap(b->invertArrayN2O2O2NOptimized());
    o2n->transformWithIndArr(*zeMap);
    std::vector<MCAuto<DataArrayIdType>> ret(nbOfSrcProcs);
    if (!o2n->empty())
    {
        mcIdType nbOfNodesWithoutDup(o2n->getMaxAbsValueInArray() + 1);
        wholeMesh = MEDCouplingUMesh::MergeUMeshes(FromVecAutoToVecOfConst<MEDCouplingUMesh>(srcMeshes));
        wholeMesh->renumberNodes(o2n->begin(), nbOfNodesWithoutDup);
        wholeMesh->checkConsistencyLight();
        // remove ghost cells
        for (std::size_t i = 0; i < nbOfSrcProcs; ++i)
        {
            ret[i] = b->findIdForEach(srcGlobalNodeIds[i]->begin(), srcGlobalNodeIds[i]->end());
        }
        wholeMesh->zipConnectivityTraducer(0);
    }
    else
    {
        for (std::size_t i = 0; i < nbOfSrcProcs; ++i)
        {
            ret[i] = DataArrayIdType::New();
            ret[i]->alloc(0, 1);
        }
    }
    return ret;
}

namespace
{
MCAuto<DataArrayDouble>
MatrixVectorMultiplySingleCompo(const std::vector<std::map<mcIdType, double>> &matrix, const DataArrayDouble *vectorArr)
{
    std::size_t nbOfRows(matrix.size());
    MCAuto<DataArrayDouble> ret(DataArrayDouble::New());
    ret->alloc(nbOfRows, 1);
    const double *inputVec(vectorArr->begin());
    double *retPtr(ret->getPointer());
    for (std::size_t i = 0; i < nbOfRows; ++i)
    {
        const std::map<mcIdType, double> &line(matrix[i]);
        retPtr[i] = 0.0;
        for (const auto &kv : line)
        {
            retPtr[i] += kv.second * inputVec[kv.first];
        }
    }
    return ret;
}
}  // namespace

MCAuto<DataArrayDouble>
MEDCoupling::MatrixVectorMultiply(
    const std::vector<std::map<mcIdType, double>> &matrix, const DataArrayDouble *vectorArr
)
{
    std::size_t nbCompo(vectorArr->getNumberOfComponents());
    std::vector<MCAuto<DataArrayDouble>> res(nbCompo);
    for (std::size_t i = 0; i < nbCompo; ++i)
    {
        MCAuto<DataArrayDouble> curCompo(vectorArr->keepSelectedComponents({i}));
        res[i] = MatrixVectorMultiplySingleCompo(matrix, curCompo);
    }
    return MCAuto<DataArrayDouble>(DataArrayDouble::Meld(FromVecAutoToVecOfConst<DataArrayDouble>(res)));
}
