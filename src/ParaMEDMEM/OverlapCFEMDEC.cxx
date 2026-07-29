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

#include "OverlapCFEMDEC.hxx"
#include "CommInterface.hxx"
#include "CommonCFEMDEC.hxx"

using namespace MEDCoupling;

OverlapCFEMDEC::OverlapCFEMDEC(const std::set<int> &procIds, const MPI_Comm &worldComm)
{
    CommInterface comm;
    std::unique_ptr<int[]> ranks_world(new int[procIds.size()]);  // ranks of sources and targets in worldComm
    std::copy(procIds.begin(), procIds.end(), ranks_world.get());
    MPI_Group group, worldGroup;
    comm.commGroup(worldComm, &worldGroup);
    comm.groupIncl(worldGroup, (int)procIds.size(), ranks_world.get(), &group);
    comm.commCreate(worldComm, group, &_comm);
    comm.groupFree(&group);
    comm.groupFree(&worldGroup);
    if (_comm == MPI_COMM_NULL)
    {
        THROW_IK_EXCEPTION("OverlapCFEMDEC constructor : Fail to create communicator");
    }
    std::set<int> idsUnion;
    for (unsigned int i = 0; i < procIds.size(); i++) idsUnion.insert(i);
    _group.reset(new MPIProcessorGroup(comm, idsUnion, _comm));
}

void
OverlapCFEMDEC::attachSourceField(MEDCouplingFieldDouble *srcField, DataArrayIdType *srcGlobalNodeIds)
{
    _src_field.takeRef(srcField);
    _src_global_node_ids.takeRef(srcGlobalNodeIds);
}

void
OverlapCFEMDEC::attachTargetMesh(MEDCouplingUMesh *trgMesh, DataArrayIdType *trgGlobalNodeIds)
{
    _trg_mesh.takeRef(trgMesh);
    _trg_global_node_ids.takeRef(trgGlobalNodeIds);
}

template <int spaceDim>
void
OverlapCFEMDEC::synchronizeT(
    std::vector<MCAuto<MEDCouplingUMesh>> &srcMeshes, std::vector<MCAuto<DataArrayIdType>> &srcGlobalNodeIds
)
{
    std::vector<BBTreeClosest<spaceDim, mcIdType>> bbSrc, bbTrg;
    BBTreeClosest<spaceDim, mcIdType> myBBSrc(
        ShareBBTreesOfAllProcs<spaceDim>(_group.get(), getSourceLocalMesh(), bbSrc /*output*/)
    );
    ShareBBTreesOfAllProcs<spaceDim>(_group.get(), _trg_mesh, bbTrg /*output*/);
    int grpSz(_group->size());
    std::vector<std::vector<mcIdType>> tis(grpSz);
    std::vector<std::vector<double>> tds(grpSz);
    std::vector<MCAuto<DataArrayIdType>> bis(grpSz);
    std::vector<MCAuto<DataArrayDouble>> bds(grpSz);
    _glb_nodes_in_whole_per_trg.resize(grpSz);
    std::vector<std::string> ts;
    for (int iProcTrg = 0; iProcTrg < grpSz; ++iProcTrg)
    {
        std::set<const BBTreeClosest<spaceDim, mcIdType> *> blockSelectedPerTrgProc;
        std::vector<mcIdType> cellsToSendToTrg;
        const BBTreeClosest<spaceDim, mcIdType> &curBBTree(bbTrg[iProcTrg]);
        // iterate over all terminal nodes of targets
        for (const auto &leaf : curBBTree)
        {
            const BBTreeClosest<spaceDim, mcIdType> &leaf2(
                static_cast<const BBTreeClosest<spaceDim, mcIdType> &>(leaf)
            );
            // compute min of maxes over all source procs
            double zeMin(std::numeric_limits<double>::max());
            for (int iProcSrc = 0; iProcSrc < grpSz; ++iProcSrc)
            {
                bbSrc[iProcSrc].bboxMinOfMaxes(leaf2.getBBox(), zeMin);
            }
            myBBSrc.bboxSelect(leaf2.getBBox(), zeMin, blockSelectedPerTrgProc);
        }
        for (auto block : blockSelectedPerTrgProc)
        {
            const std::vector<mcIdType> &elems(block->getElements());
            cellsToSendToTrg.insert(cellsToSendToTrg.end(), elems.cbegin(), elems.cend());
        }
        MCAuto<DataArrayIdType> glbNodeIdsOfPart;
        MCAuto<MEDCouplingUMesh> part(
            MEDCoupling::ReduceMesh(
                getSourceLocalMesh(),
                _src_global_node_ids,
                cellsToSendToTrg.data(),
                cellsToSendToTrg.data() + cellsToSendToTrg.size(),
                glbNodeIdsOfPart  // output
            )
        );
        _glb_nodes_in_whole_per_trg[iProcTrg] = glbNodeIdsOfPart;
        {
            std::vector<double> td;
            std::vector<mcIdType> ti;
            ts.clear();
            part->getTinySerializationInformation(td, ti, ts);
            tis[iProcTrg] = ti;
            tds[iProcTrg] = td;
            {
                DataArrayIdType *biTmp(nullptr);
                DataArrayDouble *bdTmp(nullptr);
                part->serialize(biTmp, bdTmp);
                bis[iProcTrg] = biTmp;
                bdTmp->rearrange(1);
                bds[iProcTrg] = bdTmp;
            }
        }
    }
    //
    std::vector<std::vector<mcIdType>> tisFromSrc(all2allVector<mcIdType>(_group.get(), tis));
    std::vector<std::vector<double>> tdsFromSrc(all2allVector<double>(_group.get(), tds));
    std::vector<MCAuto<DataArrayIdType>> bisFromSrc(all2allDA<mcIdType>(_group.get(), bis));
    std::vector<MCAuto<DataArrayDouble>> bdsFromSrc(all2allDA<double>(_group.get(), bds));
    srcGlobalNodeIds = std::move(all2allDA<mcIdType>(_group.get(), _glb_nodes_in_whole_per_trg));
    //
    srcMeshes.resize(grpSz);
    for (int iProc = 0; iProc < grpSz; ++iProc)
    {
        srcMeshes[iProc] = MEDCouplingUMesh::New();
        bdsFromSrc[iProc]->rearrange(spaceDim);
        srcMeshes[iProc]->unserialization(
            tdsFromSrc[iProc], tisFromSrc[iProc], bisFromSrc[iProc], bdsFromSrc[iProc], ts
        );
    }
}

void
OverlapCFEMDEC::computeMatrix(
    const std::vector<MCAuto<MEDCouplingUMesh>> &srcMeshes, const std::vector<MCAuto<DataArrayIdType>> &srcGlobalNodeIds
)
{
    MCAuto<MEDCouplingUMesh> wholeMesh;
    _src_rank_of_nodes_in_whole = ComputeNodeIdsPerProc(srcMeshes, srcGlobalNodeIds, wholeMesh);
    this->_nb_nodes_src_mesh = wholeMesh->getNumberOfNodes();
    // compute matrix
    const double *coordsOfTrgMesh(_trg_mesh->getCoords()->begin());
    const mcIdType nbOfTrgPts(_trg_mesh->getNumberOfNodes());

    MEDCouplingFieldDiscretizationOnNodesFE::computeCrudeMatrix(
        wholeMesh, coordsOfTrgMesh, nbOfTrgPts, this->_matrix, this->getFEOptions()
    );
}

void
OverlapCFEMDEC::synchronize()
{
    checkSameSpaceDim();
    std::vector<MCAuto<MEDCouplingUMesh>> srcMeshes;
    std::vector<MCAuto<DataArrayIdType>> srcGlobalNodeIds;
    switch (_src_field->getMesh()->getSpaceDimension())
    {
        case 3:
        {
            this->synchronizeT<3>(srcMeshes, srcGlobalNodeIds);
            break;
        }
        case 2:
        {
            this->synchronizeT<2>(srcMeshes, srcGlobalNodeIds);
            break;
        }
        case 1:
        {
            this->synchronizeT<1>(srcMeshes, srcGlobalNodeIds);
            break;
        }
        default:
        {
            THROW_IK_EXCEPTION("OverlapCFEMDEC::synchronize : Manage only spaceDim 1, 2 or 3.");
        }
    }
    computeMatrix(srcMeshes, srcGlobalNodeIds);
}

MCAuto<MEDCouplingFieldDouble>
OverlapCFEMDEC::computeTargetField()
{
    int grpSz(_group->size());
    const MPI_Comm *comm(_group->getComm());
    // check same nb of compo of source field over proc
    int nbCompo(FromIdType<std::int32_t>(_src_field->getNumberOfComponents()));
    std::vector<int32_t> nbCompoOnAllProcs(grpSz);
    _group->getCommInterface().allGather(&nbCompo, 1, MPI_INT32_T, nbCompoOnAllProcs.data(), 1, MPI_INT32_T, *comm);
    for (int curNbCompo : nbCompoOnAllProcs)
    {
        if (curNbCompo != nbCompo)
        {
            THROW_IK_EXCEPTION(
                "On proc #" << _group->myRank() << " : " << "Nb of compo of source field is " << nbCompo
                            << ". Presence of nbCompo = " << curNbCompo << " in source field in MPI group !"
            );
        }
    }
    //
    std::vector<MCAuto<DataArrayDouble>> arrsToSend(grpSz);
    if (grpSz != (int)_glb_nodes_in_whole_per_trg.size())
    {
        THROW_IK_EXCEPTION("On proc #" << _group->myRank() << " : big problem detected");
    }
    for (std::size_t i = 0; i < _glb_nodes_in_whole_per_trg.size(); ++i)
    {
        MCAuto<DataArrayIdType> locIdsToSendForCurTrgProc(_src_global_node_ids->findIdForEach(
            _glb_nodes_in_whole_per_trg[i]->begin(), _glb_nodes_in_whole_per_trg[i]->end()
        ));
        arrsToSend[i] = _src_field->getArray()->selectByTupleId(
            locIdsToSendForCurTrgProc->begin(), locIdsToSendForCurTrgProc->end()
        );
    }
    MCAuto<DataArrayDouble> srcArr(DataArrayDouble::New());
    srcArr->alloc(_nb_nodes_src_mesh, nbCompo);
    std::vector<MCAuto<DataArrayDouble>> perSrcProc(all2allDA<double>(_group.get(), arrsToSend));
    for (int i = 0; i < grpSz; ++i)
    {
        perSrcProc[i]->rearrange(nbCompo);
        srcArr->setPartOfValues3(
            perSrcProc[i],
            _src_rank_of_nodes_in_whole[i]->begin(),
            _src_rank_of_nodes_in_whole[i]->end(),
            0,
            nbCompo,
            1,
            true
        );
    }
    // matrix * vector (srcArr)
    MCAuto<DataArrayDouble> res(MatrixVectorMultiply(_matrix, srcArr));
    res->copyStringInfoFrom(*_src_field->getArray());
    MCAuto<MEDCouplingFieldDouble> ret(MEDCouplingFieldDouble::New(ON_NODES_FE));
    ret->setArray(res);
    ret->setMesh(_trg_mesh);
    ret->setName(_src_field->getName());
    ret->setQuantityKind(_src_field->getQuantityKind());
    ret->checkConsistencyLight();
    return ret;
}

void
OverlapCFEMDEC::basicCheck()
{
    if (_src_field.isNull())
    {
        THROW_IK_EXCEPTION("On proc #" << _group->myRank() << " : Null source field");
    }
    try
    {
        _src_field->checkConsistencyLight();
    }
    catch (const INTERP_KERNEL::Exception &e)
    {
        THROW_IK_EXCEPTION("On proc #" << _group->myRank() << " : source field is not consistent (" << e.what() << ")");
    }
    if (_trg_mesh.isNull())
    {
        THROW_IK_EXCEPTION("On proc #" << _group->myRank() << " : Null target mesh");
    }
    try
    {
        _trg_mesh->checkConsistencyLight();
    }
    catch (const INTERP_KERNEL::Exception &e)
    {
        THROW_IK_EXCEPTION("On proc #" << _group->myRank() << " : target mesh is not consistent (" << e.what() << ")");
    }
    if (_src_global_node_ids.isNull())
    {
        THROW_IK_EXCEPTION("On proc #" << _group->myRank() << " : Null source global node ids");
    }
    if (_trg_global_node_ids.isNull())
    {
        THROW_IK_EXCEPTION("On proc #" << _group->myRank() << " : Null target global node ids");
    }
    _src_global_node_ids->checkAllocated();
    _src_global_node_ids->checkNbOfComps(1, "source global node ids array has to be single compo");
    if (_src_global_node_ids->getNumberOfTuples() != _src_field->getMesh()->getNumberOfNodes())
    {
        THROW_IK_EXCEPTION(
            "On proc #" << _group->myRank() << " : size of source global node ids array ("
                        << _src_global_node_ids->getNumberOfTuples()
                        << ") does not match the number of nodes on geometric support of source field ("
                        << _src_field->getMesh()->getNumberOfNodes() << ")"
        );
    }
    _trg_global_node_ids->checkAllocated();
    _trg_global_node_ids->checkNbOfComps(1, "target global node ids array has to be single compo");
    if (_trg_global_node_ids->getNumberOfTuples() != _trg_mesh->getNumberOfNodes())
    {
        THROW_IK_EXCEPTION(
            "On proc #" << _group->myRank() << " : size of target global node ids array ("
                        << _trg_global_node_ids->getNumberOfTuples()
                        << ") does not match the number of nodes on geometric support of target mesh ("
                        << _trg_mesh->getNumberOfNodes() << ")"
        );
    }
}

const MEDCouplingUMesh *
OverlapCFEMDEC::getSourceLocalMesh() const
{
    const MEDCouplingUMesh *ret(dynamic_cast<const MEDCouplingUMesh *>(_src_field->getMesh()));
    if (!ret)
    {
        THROW_IK_EXCEPTION("On proc #" << _group->myRank() << " : mesh behind source field is not a MEDCouplingUMesh");
    }
    return ret;
}

void
OverlapCFEMDEC::checkSameSpaceDim()
{
    basicCheck();
    if (_src_field->getMesh()->getSpaceDimension() != _trg_mesh->getSpaceDimension())
    {
        THROW_IK_EXCEPTION(
            "On proc #" << _group->myRank() << " : mismatch of space dimension : src = ("
                        << _src_field->getMesh()->getSpaceDimension() << ") trg = '(" << _trg_mesh->getSpaceDimension()
                        << ")"
        );
    }
    int unionGrpSz(_group->size());
    int spaceDim(_trg_mesh->getSpaceDimension());
    std::vector<int> spaceDimsOnAllProcs(unionGrpSz);
    _group->getCommInterface().allGather(&spaceDim, 1, MPI_INT, spaceDimsOnAllProcs.data(), 1, MPI_INT, _comm);
    for (int curSpaceDim : spaceDimsOnAllProcs)
    {
        if (curSpaceDim != spaceDim)
        {
            THROW_IK_EXCEPTION("On proc #" << _group->myRank() << " : Different spaceDim detected accross procs !");
        }
    }
}
