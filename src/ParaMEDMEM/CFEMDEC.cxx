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

#include "CFEMDEC.hxx"
#include "MPIProcessorGroup.hxx"
#include "MPITraits.hxx"
#include "MEDCouplingFieldDiscretizationOnNodesFE.hxx"

#include <functional>
#include <cstdint>
#include <cstddef>
#include <cstring>

using namespace MEDCoupling;

MPIProcessorGroup *
CFEMDECDirectAccess::getUnionGrp() const
{
    return static_cast<MPIProcessorGroup *>(_master->getUnionGrp());
}

MPIProcessorGroup *
CFEMDECDirectAccess::getSourceGrp() const
{
    return static_cast<MPIProcessorGroup *>(_master->getSourceGrp());
}

MPIProcessorGroup *
CFEMDECDirectAccess::getTargetGrp() const
{
    return static_cast<MPIProcessorGroup *>(_master->getTargetGrp());
}

MPIProcessorGroup *
CFEMDECReverseAccess::getUnionGrp() const
{
    return static_cast<MPIProcessorGroup *>(_master->getReverseUnionGrp());
}

MPIProcessorGroup *
CFEMDECReverseAccess::getSourceGrp() const
{
    return static_cast<MPIProcessorGroup *>(_master->getTargetGrp());
}

MPIProcessorGroup *
CFEMDECReverseAccess::getTargetGrp() const
{
    return static_cast<MPIProcessorGroup *>(_master->getSourceGrp());
}

void
CFEMDECOneWay::reinitializeOnNewMesh()
{
    _local_mesh_stamp_on_sync = nullptr;
}

MPIProcessorGroup *
CFEMDECOneWay::getSourceGrp() const
{
    return _to_master->getSourceGrp();
}

MPIProcessorGroup *
CFEMDECOneWay::getTargetGrp() const
{
    return _to_master->getTargetGrp();
}

MPIProcessorGroup *
CFEMDECOneWay::getUnionGrp() const
{
    return _to_master->getUnionGrp();
}

MCAuto<MEDCouplingUMesh>
CFEMDECOneWay::getLocalMesh() const
{
    return _master->getLocalMesh();
}

MCAuto<DataArrayIdType>
CFEMDECOneWay::getGlobalNodeIdsOnLocalMesh() const
{
    return _master->getGlobalNodeIdsOnLocalMesh();
}

static bool
anyChar(char *zeStart, char *zeEnd)
{
    for (char *it = zeStart; it != zeEnd; ++it)
        if (*it)
            return true;
    return false;
}

bool
CFEMDECOneWay::isMatrixToRecompute()
{
    MPI_Comm globalComm(getSourceGrp()->getWorldComm());
    if (getLocalMesh().isNull())
        THROW_IK_EXCEPTION("isMatrixToRecompute : no mesh attached !");
    char toRecomputeLocal = char((MEDCouplingUMesh *)getLocalMesh() != _local_mesh_stamp_on_sync);
    MPIProcessorGroup *grp(getUnionGrp());
    auto size(grp->size());
    std::unique_ptr<char[]> recvData(new char[size]);
    grp->getCommInterface().allGather(&toRecomputeLocal, 1, MPI_CHAR, recvData.get(), 1, MPI_CHAR, *grp->getComm());
    return anyChar(recvData.get(), recvData.get() + size);
}

template <class T>
class BCastDataArrayFunctor
{
   public:
    BCastDataArrayFunctor(std::size_t sz)
    {
        using ArrayType = typename MPITraits<T>::ArrayType;
        _arr = ArrayType::New();
        _arr->alloc(sz, 1);
    }
    T *getPointer() { return _arr->getPointer(); }
    MCAuto<typename MPITraits<T>::ArrayType> retn() { return _arr; }

   private:
    MCAuto<typename MPITraits<T>::ArrayType> _arr;
};

template <class T>
class BCastVectorFunctor
{
   public:
    BCastVectorFunctor(std::size_t sz) : _arr(sz) {}
    T *getPointer() { return _arr.data(); }
    std::vector<T> retn() { return _arr; }

   private:
    std::vector<T> _arr;
};

template <class T, class RET, class FUNCTOR>
RET
bCastTArrayFromProcInternal2(MPIProcessorGroup *grp, const T *arrToSend, mcIdType sz, int rkBroadCasting)
{
    const MPI_Comm *comm(grp->getComm());
    int sz2;
    if (grp->myRank() == rkBroadCasting)
    {
        sz2 = (int)sz;
    }
    grp->getCommInterface().broadcast(&sz2, 1, MPI_INT32_T, rkBroadCasting, *comm);
    FUNCTOR arr(sz2);
    if (grp->myRank() == rkBroadCasting)
    {
        std::copy(arrToSend, arrToSend + sz, arr.getPointer());
    }
    grp->getCommInterface().broadcast(arr.getPointer(), sz2, MPITraits<T>::MPIType, rkBroadCasting, *comm);
    return arr.retn();
}

template <class T>
MCAuto<typename MPITraits<T>::ArrayType>
bCastTArrayFromProc(MPIProcessorGroup *grp, const typename MPITraits<T>::ArrayType *arrToSend, int rkBroadCasting)
{
    const T *pt(nullptr);
    mcIdType sz(0);
    if (arrToSend)
    {
        pt = arrToSend->begin();
        sz = arrToSend->getNbOfElems();
    }
    return bCastTArrayFromProcInternal2<T, MCAuto<typename MPITraits<T>::ArrayType>, BCastDataArrayFunctor<T>>(
        grp, pt, sz, rkBroadCasting
    );
}

template <class T>
std::vector<T>
bCastTArrayFromProcInternal(MPIProcessorGroup *grp, const T *arrToSend, mcIdType sz, int rkBroadCasting)
{
    return bCastTArrayFromProcInternal2<T, std::vector<T>, BCastVectorFunctor<T>>(grp, arrToSend, sz, rkBroadCasting);
}

template <typename T>
struct VectorMCDAPolicy
{
};

template <typename T>
struct VectorMCDAPolicy<std::vector<T>>
{
    static auto size(const std::vector<T> &obj) -> decltype(obj.size()) { return obj.size(); }

    static auto begin(const std::vector<T> &obj) -> decltype(obj.begin()) { return obj.begin(); }

    static auto end(const std::vector<T> &obj) -> decltype(obj.end()) { return obj.end(); }
};

template <typename U>
struct VectorMCDAPolicy<MCAuto<U>>
{
    static auto size(const MCAuto<U> &obj) -> decltype(obj->getNbOfElems())
    {
        if (obj.isNull())
            return 0;
        return obj->getNbOfElems();
    }

    static auto begin(const MCAuto<U> &obj) -> decltype(obj->begin())
    {
        if (obj.isNull())
            return nullptr;
        return obj->begin();
    }

    static auto end(const MCAuto<U> &obj) -> decltype(obj->end())
    {
        if (obj.isNull())
            return nullptr;
        return obj->end();
    }
};

template <class T, class RET>
std::vector<RET>
all2allInternal2(MPIProcessorGroup *grp, const std::vector<RET> &input, std::function<RET(std::size_t, T *, T *)> func)
{
    using RETWrapper = VectorMCDAPolicy<RET>;
    const MPI_Comm *comm(grp->getComm());

    int size(static_cast<int>(input.size()));

    MPI_Datatype dtype = MPITraits<T>::MPIType;

    std::vector<std::int32_t> sendcounts(size);
    for (int i = 0; i < size; ++i)
    {
        sendcounts[i] = static_cast<std::int32_t>(RETWrapper::size(input[i]));
    }

    std::vector<std::int32_t> recvcounts(size, 125);
    grp->getCommInterface().allToAll(sendcounts.data(), 1, MPI_INT32_T, recvcounts.data(), 1, MPI_INT32_T, *comm);

    std::vector<std::int32_t> sdispls(size, 0);
    for (int i = 1; i < size; ++i)
    {
        sdispls[i] = sdispls[i - 1] + sendcounts[i - 1];
    }

    std::vector<T> sendbuf;
    sendbuf.reserve(sdispls.back() + sendcounts.back());
    for (const auto &v : input)
    {
        sendbuf.insert(sendbuf.end(), RETWrapper::begin(v), RETWrapper::end(v));
    }

    std::vector<int> rdispls(size, 0);
    for (int i = 1; i < size; ++i) rdispls[i] = rdispls[i - 1] + recvcounts[i - 1];

    std::vector<T> recvbuf(rdispls.back() + recvcounts.back());

    grp->getCommInterface().allToAllV(
        sendbuf.data(),
        sendcounts.data(),
        sdispls.data(),
        dtype,
        recvbuf.data(),
        recvcounts.data(),
        rdispls.data(),
        dtype,
        *comm
    );

    std::vector<RET> ret(size);
    for (int i = 0; i < size; ++i)
    {
        int nbElems(recvcounts[i]);
        ret[i] = func(nbElems, recvbuf.data() + rdispls[i], recvbuf.data() + rdispls[i] + nbElems);
    }

    return ret;
}

template <class T>
std::vector<typename std::vector<T>>
all2allVector(MPIProcessorGroup *grp, const std::vector<typename std::vector<T>> &input)
{
    return all2allInternal2<T, typename std::vector<T>>(
        grp,
        input,
        [](std::size_t nbElems, T *bg, T *end)
        {
            std::vector<T> elt(nbElems);
            std::copy(bg, end, elt.data());
            return elt;
        }
    );
}

template <class T>
std::vector<MCAuto<typename MPITraits<T>::ArrayType>>
all2allDA(MPIProcessorGroup *grp, const std::vector<MCAuto<typename MPITraits<T>::ArrayType>> &input)
{
    return all2allInternal2<T, MCAuto<typename MPITraits<T>::ArrayType>>(
        grp,
        input,
        [](std::size_t nbElems, T *bg, T *end)
        {
            using ArrayType = typename MPITraits<T>::ArrayType;
            MCAuto<ArrayType> elt(ArrayType::New());
            elt->alloc(nbElems, 1);
            std::copy(bg, end, elt->getPointer());
            return elt;
        }
    );
}

template <class T, class RET>
std::vector<RET>
gatherTArrayOnProcInternal2(
    MPIProcessorGroup *grp,
    int rkGathering,
    const T *arrToSend,
    mcIdType sz,
    std::function<RET(std::size_t, T *, T *)> func
)
{
    std::vector<RET> ret;
    const MPI_Comm *comm(grp->getComm());
    int rank(grp->myRank());
    std::vector<std::int32_t> ti_ex_2;
    if (rank == rkGathering)
    {
        ti_ex_2.resize(grp->size());
    }
    std::int64_t lenOfArr(FromIdType<std::int64_t>(sz));
    grp->getCommInterface().gather(&lenOfArr, 1, MPI_INT32_T, ti_ex_2.data(), 1, MPI_INT32_T, rkGathering, *comm);
    std::vector<T> ti_ex_3;
    std::vector<int> disps;
    if (rank == rkGathering)
    {
        std::uint64_t nbElems(0);
        std::for_each(ti_ex_2.begin(), ti_ex_2.end(), [&nbElems](std::int32_t v) { nbElems += v; });
        ti_ex_3.resize(nbElems);
        disps.resize(grp->size() + 1);
        int dispsCnt(0);
        {
            int *dispsPt(disps.data());
            std::for_each(
                ti_ex_2.begin(),
                ti_ex_2.end(),
                [&dispsCnt, &dispsPt](std::int32_t v)
                {
                    *dispsPt = dispsCnt;
                    dispsCnt += v;
                    dispsPt++;
                }
            );
            *dispsPt = (int)nbElems;
        }
    }
    grp->getCommInterface().gatherV(
        arrToSend,
        FromIdType<int>(sz),
        MPITraits<T>::MPIType,
        ti_ex_3.data(),
        ti_ex_2.data(),
        disps.data(),
        MPITraits<T>::MPIType,
        rkGathering,
        *comm
    );
    if (rank == rkGathering)
    {
        for (int i = 0; i < grp->size(); ++i)
        {
            int nbElems(disps[i + 1] - disps[i]);
            if (nbElems > 0)
            {
                ret.emplace_back(func(nbElems, ti_ex_3.data() + disps[i], ti_ex_3.data() + disps[i + 1]));
            }
        }
    }
    return ret;
}

template <class T>
std::vector<std::vector<T>>
gatherTArrayOnProcInternal(MPIProcessorGroup *grp, int rkGathering, const T *arrToSend, mcIdType sz)
{
    return gatherTArrayOnProcInternal2<T, std::vector<T>>(
        grp,
        rkGathering,
        arrToSend,
        sz,
        [](std::size_t nbElems, T *bg, T *end)
        {
            std::vector<T> elt(nbElems);
            std::copy(bg, end, elt.data());
            return elt;
        }
    );
}

template <class T>
std::vector<MCAuto<typename MPITraits<T>::ArrayType>>
gatherTArrayOnProc(MPIProcessorGroup *grp, int rkGathering, const typename MPITraits<T>::ArrayType *arrToSend)
{
    const T *pt(nullptr);
    mcIdType sz(0);
    if (arrToSend)
    {
        pt = arrToSend->begin();
        sz = arrToSend->getNbOfElems();
    }
    return gatherTArrayOnProcInternal2<T, MCAuto<typename MPITraits<T>::ArrayType>>(
        grp,
        rkGathering,
        pt,
        sz,
        [](std::size_t nbElems, T *bg, T *end)
        {
            using ArrayType = typename MPITraits<T>::ArrayType;
            MCAuto<ArrayType> elt(ArrayType::New());
            elt->alloc(nbElems, 1);
            std::copy(bg, end, elt->getPointer());
            return elt;
        }
    );
}

namespace
{
std::vector<std::int8_t>
packBits(const std::vector<bool> &bits)
{
    std::size_t n = bits.size();
    std::size_t byte_count = (n + 7) / 8;
    std::vector<std::int8_t> result(byte_count, 0);
    for (std::size_t i = 0; i < n; ++i)
    {
        if (bits[i])
        {
            std::size_t byte_index = i / 8;
            std::size_t bit_index = i % 8;

            result[byte_index] |= static_cast<std::int8_t>(1 << bit_index);
        }
    }
    return result;
}

std::vector<bool>
unpackBits(const std::int8_t *bytesPtr, std::size_t original_size)
{
    std::vector<bool> result;
    result.reserve(original_size);
    for (std::size_t i = 0; i < original_size; ++i)
    {
        std::size_t byte_index = i / 8;
        std::size_t bit_index = i % 8;

        bool bit = (bytesPtr[byte_index] >> bit_index) & 0x01;
        result.push_back(bit);
    }
    return result;
}

std::vector<std::vector<bool>>
allGatherVectBoolOnProc(MPIProcessorGroup *grp, const std::vector<bool> &structure)
{
    int nbOfProcs(grp->size());
    CommInterface ci(grp->getCommInterface());
    std::vector<std::int8_t> structure2(packBits(structure));
    std::vector<std::vector<bool>> structures(nbOfProcs);
    std::vector<mcIdType> vbPerProc(nbOfProcs);
    mcIdType szSt(structure.size());
    ci.allGather(
        &szSt, 1, MPITraits<mcIdType>::MPIType, vbPerProc.data(), 1, MPITraits<mcIdType>::MPIType, *(grp->getComm())
    );
    {
        std::vector<int> vbPerProc1(nbOfProcs), vbPerProc2(nbOfProcs + 1);
        vbPerProc2[0] = 0;
        int *vbPerProc1Ptr(vbPerProc1.data()), *vbPerProc2Ptr(vbPerProc2.data());
        std::for_each(
            vbPerProc.cbegin(),
            vbPerProc.cend(),
            [&vbPerProc1Ptr, &vbPerProc2Ptr](mcIdType v)
            {
                *vbPerProc1Ptr = int((v + 7) / 8);
                vbPerProc2Ptr[1] = vbPerProc2Ptr[0] + *vbPerProc1Ptr++;
                vbPerProc2Ptr++;
            }
        );
        std::vector<std::int8_t> data(vbPerProc2.back());
        ci.allGatherV(
            structure2.data(),
            (int)structure2.size(),
            MPI_CHAR,
            data.data(),
            vbPerProc1.data(),
            vbPerProc2.data(),
            MPI_CHAR,
            *(grp->getComm())
        );
        for (int i = 0; i < nbOfProcs; ++i)
        {
            structures[i] = unpackBits(data.data() + vbPerProc2[i], vbPerProc[i]);
        }
    }
    return structures;
}

template <int spaceDim>
void
PrintSelectionOfCells(
    int trgProcId, mcIdType nbCells, const std::set<const BBTreeClosest<spaceDim, mcIdType> *> &selectionPerTrgProc
)
{
    std::size_t sz2(0);
    for (const auto &blk : selectionPerTrgProc)
    {
        sz2 += blk->getElements().size();
    }
    std::cout << "  - to send to trg " << trgProcId << " " << sz2 << " ( in " << selectionPerTrgProc.size()
              << " blocks )" << " / " << nbCells << std::endl;
}

void
WriteVTUOfSrcMeshCandidateTrg(
    MPIProcessorGroup *unionGrp, MPIProcessorGroup *sourceGrp, const std::vector<MCAuto<MEDCouplingUMesh>> &srcMeshes
)
{
    int curProc(unionGrp->myRank()), nbProcSrc(sourceGrp->size());
    for (int i = 0; i < nbProcSrc; ++i)
    {
        std::ostringstream oss;
        oss << "file_" << curProc << "_" << i << ".vtu";
        if (srcMeshes[i]->getNumberOfCells() > 0)
        {
            srcMeshes[i]->writeVTK(oss.str());
        }
    }
}

}  // namespace

template <int spaceDim>
BBTreeClosest<spaceDim, mcIdType>
ShareBBTreesOfAllProcs(
    MPIProcessorGroup *unionGrp, const MEDCouplingUMesh *mesh, std::vector<BBTreeClosest<spaceDim, mcIdType>> &ret
)
{
    const MPI_Comm *comm(unionGrp->getComm());
    int unionGrpSz(unionGrp->size());
    MCAuto<DataArrayDouble> bbox(mesh->getBoundingBoxForBBTree());
    mcIdType nbCells(mesh->getNumberOfCells());
    const double *bboxPtr(bbox->begin());
    BBTreeClosest<spaceDim, mcIdType> myTreeBase(bboxPtr, nullptr, 0, nbCells);
    std::vector<bool> structure;
    std::vector<std::array<double, 2 * spaceDim>> bboxData;
    myTreeBase.serializeCompact(structure, bboxData);
    std::vector<std::vector<bool>> structures(allGatherVectBoolOnProc(unionGrp, structure));
    std::vector<std::vector<std::array<double, 2 * spaceDim>>> bboxes;
    {
        std::unique_ptr<double[]> result;
        std::unique_ptr<mcIdType[]> resultIndex;
        int nbProcs(unionGrp->getCommInterface().allGatherArraysTT<double>(
            *comm,
            reinterpret_cast<double *>(bboxData.data()),
            ToIdType(bboxData.size() * 2 * spaceDim),
            result,
            resultIndex
        ));
        bboxes.resize(nbProcs);
        for (int iProc = 0; iProc < nbProcs; ++iProc)
        {
            mcIdType nbOfBlocksToCpy(resultIndex[iProc + 1] - resultIndex[iProc]);
            bboxes[iProc].resize(nbOfBlocksToCpy / (2 * spaceDim));
            std::memcpy(bboxes[iProc].data(), result.get() + resultIndex[iProc], nbOfBlocksToCpy * sizeof(double));
        }
    }
    std::size_t nbOfProcs(structures.size());
    ret.resize(nbOfProcs);
    for (std::size_t iProc = 0; iProc < nbOfProcs; ++iProc)
    {
        ret[iProc] = std::move(BBTreeClosest<spaceDim, mcIdType>::DeserializeCompact(structures[iProc], bboxes[iProc]));
    }
    return myTreeBase;
}

/*!
 *  Associated method : CFEMDECOneWayTarget::dispatchMeshPartsTrgOnly
 */
template <int spaceDim>
void
CFEMDECOneWaySource::dispatchMeshPartsSrcOnly(
    MPIProcessorGroup *unionGrp,
    const MEDCouplingUMesh *mesh,
    const DataArrayIdType *glblNodeIds,
    const BBTreeClosest<spaceDim, mcIdType> &myTreeBase,
    const std::vector<BBTreeClosest<spaceDim, mcIdType>> &ret
)
{
    const MPI_Comm *comm(unionGrp->getComm());
    int unionGrpSz(unionGrp->size());
    int nbProcSrc(getSourceGrp()->size());
    int nbProcTrg(unionGrpSz - nbProcSrc);
    std::vector<std::vector<mcIdType>> tis(unionGrpSz);
    std::vector<std::vector<double>> tds(unionGrpSz);
    std::vector<MCAuto<DataArrayIdType>> bis(unionGrpSz);
    std::vector<MCAuto<DataArrayDouble>> bds(unionGrpSz);
    std::vector<MCAuto<DataArrayIdType>> globalNodeIds(unionGrpSz);
    int curProc(getSourceGrp()->myRank());
    {
        for (int iProcTrg = 0; iProcTrg < nbProcTrg; ++iProcTrg)
        {
            std::set<const BBTreeClosest<spaceDim, mcIdType> *> blockSelectedPerTrgProc;
            std::vector<mcIdType> cellsToSendToTrg;
            const BBTreeClosest<spaceDim, mcIdType> &curBBTree(ret[nbProcSrc + iProcTrg]);
            // iterate over all terminal nodes of targetProc iProcTrg
            for (const auto &leaf : curBBTree)
            {
                const BBTreeClosest<spaceDim, mcIdType> &leaf2(
                    static_cast<const BBTreeClosest<spaceDim, mcIdType> &>(leaf)
                );
                // compute min of maxes over all source procs
                double zeMin(std::numeric_limits<double>::max());
                for (int iProcSrc = 0; iProcSrc < nbProcSrc; ++iProcSrc)
                {
                    ret[iProcSrc].bboxMinOfMaxes(leaf2.getBBox(), zeMin);
                }
                myTreeBase.bboxSelect(leaf2.getBBox(), zeMin, blockSelectedPerTrgProc);
            }
            //
            // PrintSelectionOfCells<spaceDim>(iProcTrg, mesh->getNumberOfCells(), blockSelectedPerTrgProc);
            //
            for (auto block : blockSelectedPerTrgProc)
            {
                const std::vector<mcIdType> &elems(block->getElements());
                cellsToSendToTrg.insert(cellsToSendToTrg.end(), elems.cbegin(), elems.cend());
            }
            MCAuto<DataArrayIdType> glbNodeIdsOfPart;
            MCAuto<MEDCouplingUMesh> part(this->ReduceMesh(
                mesh,
                glblNodeIds,
                cellsToSendToTrg.data(),
                cellsToSendToTrg.data() + cellsToSendToTrg.size(),
                glbNodeIdsOfPart  // output
            ));
            globalNodeIds[iProcTrg + nbProcSrc] = glbNodeIdsOfPart;
            {
                std::vector<double> td;
                std::vector<mcIdType> ti;
                std::vector<std::string> ts;
                part->getTinySerializationInformation(td, ti, ts);
                tis[iProcTrg + nbProcSrc] = ti;
                tds[iProcTrg + nbProcSrc] = td;
                {
                    DataArrayIdType *biTmp(nullptr);
                    DataArrayDouble *bdTmp(nullptr);
                    part->serialize(biTmp, bdTmp);
                    bis[iProcTrg + nbProcSrc] = biTmp;
                    bdTmp->rearrange(1);
                    bds[iProcTrg + nbProcSrc] = bdTmp;
                }
            }
        }
    }
    all2allVector<mcIdType>(unionGrp, tis);
    all2allVector<double>(unionGrp, tds);
    all2allDA<mcIdType>(unionGrp, bis);
    all2allDA<double>(unionGrp, bds);
    all2allDA<mcIdType>(unionGrp, globalNodeIds);
    _glb_nodes_in_whole_per_trg = std::move(globalNodeIds);
    _glb_nodes_in_whole_per_trg.erase(
        _glb_nodes_in_whole_per_trg.begin(), _glb_nodes_in_whole_per_trg.begin() + this->getSourceGrp()->size()
    );
}

/*!
 *  Associated method : CFEMDECOneWaySource::dispatchMeshPartsSrcOnly
 */
void
CFEMDECOneWayTarget::dispatchMeshPartsTrgOnly(
    int spaceDim,
    MPIProcessorGroup *unionGrp,
    const MEDCouplingUMesh *mesh,
    std::vector<MCAuto<MEDCouplingUMesh>> &srcMeshes /*output */,
    std::vector<MCAuto<DataArrayIdType>> &srcGlobalNodeIds
)
{
    int unionGrpSz(unionGrp->size());
    int nbProcSrc(getSourceGrp()->size());
    std::vector<std::vector<mcIdType>> tis(unionGrpSz);
    std::vector<std::vector<double>> tds(unionGrpSz);
    std::vector<MCAuto<DataArrayIdType>> bis(unionGrpSz);
    std::vector<MCAuto<DataArrayDouble>> bds(unionGrpSz);
    std::vector<MCAuto<DataArrayIdType>> globalNodeIds(unionGrpSz);
    std::vector<std::string> ts;
    {
        std::vector<double> tdFake;
        std::vector<mcIdType> tiFake;
        mesh->getTinySerializationInformation(tdFake, tiFake, ts);
    }
    std::vector<std::vector<mcIdType>> tisFromSrc(all2allVector<mcIdType>(unionGrp, tis));
    std::vector<std::vector<double>> tdsFromSrc(all2allVector<double>(unionGrp, tds));
    std::vector<MCAuto<DataArrayIdType>> bisFromSrc(all2allDA<mcIdType>(unionGrp, bis));
    std::vector<MCAuto<DataArrayDouble>> bdsFromSrc(all2allDA<double>(unionGrp, bds));
    srcMeshes.resize(nbProcSrc);
    for (int i = 0; i < nbProcSrc; ++i)
    {
        srcMeshes[i] = MEDCouplingUMesh::New();
        bdsFromSrc[i]->rearrange(spaceDim);
        srcMeshes[i]->unserialization(tdsFromSrc[i], tisFromSrc[i], bisFromSrc[i], bdsFromSrc[i], ts);
    }
    srcGlobalNodeIds = all2allDA<mcIdType>(unionGrp, globalNodeIds);
    srcGlobalNodeIds.erase(srcGlobalNodeIds.begin() + nbProcSrc, srcGlobalNodeIds.end());
    // WriteVTUOfSrcMeshCandidateTrg(unionGrp, getSourceGrp(), srcMeshes); // For debug
}

template <int spaceDim>
void
DispatchMeshPartsInternal(
    CFEMDECOneWay *cfemdec,
    MPIProcessorGroup *unionGrp,
    const MEDCouplingUMesh *mesh,
    const DataArrayIdType *glblNodeIds,
    std::vector<MCAuto<MEDCouplingUMesh>> &srcMeshes,
    std::vector<MCAuto<DataArrayIdType>> &srcGlobalNodeIds
)
{
    std::vector<BBTreeClosest<spaceDim, mcIdType>> ret;
    BBTreeClosest<spaceDim, mcIdType> myTreeBase(ShareBBTreesOfAllProcs<spaceDim>(unionGrp, mesh, ret /*output*/));
    if (dynamic_cast<CFEMDECOneWaySource *>(cfemdec))
    {  // source side
        dynamic_cast<CFEMDECOneWaySource *>(cfemdec)->dispatchMeshPartsSrcOnly<spaceDim>(
            unionGrp, mesh, glblNodeIds, myTreeBase, ret
        );
    }
    else
    {  // target side
        dynamic_cast<CFEMDECOneWayTarget *>(cfemdec)->dispatchMeshPartsTrgOnly(
            spaceDim, unionGrp, mesh, srcMeshes, srcGlobalNodeIds
        );
    }
}

/*!
 * EDF34966 : dispatch relevant mesh parts
 */
void
CFEMDECOneWay::dispatchMeshParts(
    std::vector<MCAuto<MEDCouplingUMesh>> &srcMeshes, std::vector<MCAuto<DataArrayIdType>> &srcGlobalNodeIds
)
{
    MCAuto<MEDCouplingUMesh> locMesh(_master->getLocalMesh());
    MCAuto<DataArrayIdType> locGlblNodeIds(_master->getGlobalNodeIdsOnLocalMesh());
    MPIProcessorGroup *unionGrp(getUnionGrp());
    switch (locMesh->getSpaceDimension())
    {
        case 3:
        {
            DispatchMeshPartsInternal<3>(this, unionGrp, locMesh, locGlblNodeIds, srcMeshes, srcGlobalNodeIds);
            break;
        }
        case 2:
        {
            DispatchMeshPartsInternal<2>(this, unionGrp, locMesh, locGlblNodeIds, srcMeshes, srcGlobalNodeIds);
            break;
        }
        case 1:
        {
            DispatchMeshPartsInternal<1>(this, unionGrp, locMesh, locGlblNodeIds, srcMeshes, srcGlobalNodeIds);
            break;
        }
        default:
        {
            THROW_IK_EXCEPTION("CFEMDECOneWay::dispatchMeshParts : Manage only spaceDim 1, 2 or 3.")
        }
    }
}

void
CFEMDECOneWay::checkSameSpaceDim()
{
    MCAuto<MEDCouplingUMesh> locMesh(_master->getLocalMesh());
    MPIProcessorGroup *unionGrp(getUnionGrp());
    const MPI_Comm *comm(unionGrp->getComm());
    int unionGrpSz(unionGrp->size());
    int spaceDim(locMesh->getSpaceDimension());
    std::vector<int> spaceDimsOnAllProcs(unionGrpSz);
    unionGrp->getCommInterface().allGather(&spaceDim, 1, MPI_INT, spaceDimsOnAllProcs.data(), 1, MPI_INT, *comm);
    for (int curSpaceDim : spaceDimsOnAllProcs)
    {
        if (curSpaceDim != spaceDim)
        {
            THROW_IK_EXCEPTION("Different spaceDim detected accross procs !");
        }
    }
}

void
CFEMDECOneWay::synchronize()
{
    checkSameSpaceDim();
    reinitializeOnNewMesh();
    std::vector<MCAuto<MEDCouplingUMesh>> srcMeshes;
    std::vector<MCAuto<DataArrayIdType>> srcGlobalNodeIds;
    dispatchMeshParts(srcMeshes, srcGlobalNodeIds);
    computeMatrix(srcMeshes, srcGlobalNodeIds);
    // keep track of local mesh pointer to determine if matrix computation is needed
    _local_mesh_stamp_on_sync = _master->getLocalMesh();
}

void
CFEMDECOneWay::smartSynchronize()
{
    if (isMatrixToRecompute())
        synchronize();
}

static std::size_t
FromVectorStringToCharSize(const std::vector<std::string> &vs)
{
    std::size_t sizeOfStr(0);
    std::for_each(
        vs.cbegin(), vs.cend(), [&sizeOfStr](const std::string &elt) { sizeOfStr += (int)elt.size() + 1; }
    );  // +1 for C nul char at the end
    return sizeOfStr;
}

static std::unique_ptr<char[]>
FromVectorStringToChar(const std::vector<std::string> &vs)
{
    std::size_t retSz(FromVectorStringToCharSize(vs));
    std::unique_ptr<char[]> ret(new char[retSz]);
    std::size_t pos(0);
    std::for_each(
        vs.cbegin(),
        vs.cend(),
        [&pos, &ret](const std::string &elt)
        {
            std::size_t sz(elt.size());
            std::copy(elt.cbegin(), elt.cbegin() + sz, ret.get() + pos);
            ret[pos + sz] = '\0';
            pos += sz + 1;
        }
    );
    return ret;
}

static std::vector<std::string>
FromCharToVectorString(const char *st, std::size_t len, int szOfOutVec)
{
    std::vector<std::string> ret(szOfOutVec);
    std::size_t pos(0);
    for (auto iComp = 0; iComp < szOfOutVec; ++iComp)
    {
        std::size_t pos2(pos);
        while (pos2 < len)
        {
            if (st[pos2] == '\0')
                break;
            pos2++;
        }
        ret[iComp] = std::string(st + pos, st + pos2);
        pos = pos2 + 1;
    }
    return ret;
}

/*!
 * Synchronized with CFEMDECOneWayTarget::receiveFromSource
 */
void
CFEMDECOneWaySource::sendToTarget(MEDCouplingFieldDouble *srcFieldOnLocal)
{
    smartSynchronize();
    checkMesh(srcFieldOnLocal);
    // management of number of components
    int nbCompo((int)srcFieldOnLocal->getNumberOfComponents());
    std::vector<std::string> comps(srcFieldOnLocal->getArray()->getInfoOnComponents());
    int sizeOfStr(int(FromVectorStringToCharSize(comps)));
    int nbSrcProcs(getSourceGrp()->size());
    {
        int tmp[2] = {nbCompo, sizeOfStr};
        std::unique_ptr<int[]> nbCompos(new int[2 * nbSrcProcs]);
        _to_master->getSourceGrp()->getCommInterface().allGather(
            tmp, 2, MPI_INT32_T, nbCompos.get(), 2, MPI_INT32_T, *getSourceGrp()->getComm()
        );
        for (auto iSrcProc = 0; iSrcProc < nbSrcProcs; ++iSrcProc)
        {
            if (nbCompos[2 * iSrcProc + 0] != nbCompo)
            {
                THROW_IK_EXCEPTION(
                    "Nb of compos is not the same for all source procs ! ( " << nbCompos[2 * iSrcProc + 0]
                                                                             << " != " << nbCompo << " )"
                );
            }
            if (nbCompos[2 * iSrcProc + 1] != sizeOfStr)
            {
                THROW_IK_EXCEPTION(
                    "Size of strings in components is not the same for all source procs ! ( "
                    << nbCompos[2 * iSrcProc + 1] << " != " << sizeOfStr << " )"
                );
            }
        }
    }
    // send nbCompo to all target procs. First send nbCompo to proc #0 of targets
    // send compostrings to all target procs. First send compostrings to proc #0 of targets
    int rkDest(getUnionGrp()->translateRank(getTargetGrp(), 0));
    if (getSourceGrp()->myRank() == 0)
    {
        {
            int tmp[2] = {nbCompo, sizeOfStr};
            getUnionGrp()->getCommInterface().send(tmp, 2, MPI_INT32_T, rkDest, 1234, *getUnionGrp()->getComm());
        }
        {
            std::unique_ptr<char[]> stringsToSend(FromVectorStringToChar(comps));
            getUnionGrp()->getCommInterface().send(
                stringsToSend.get(), sizeOfStr, MPI_INT8_T, rkDest, 1235, *getUnionGrp()->getComm()
            );
        }
    }
    MPIProcessorGroup *unionGrp(getUnionGrp());
    int srcGrpSz(getSourceGrp()->size());
    std::vector<MCAuto<DataArrayDouble>> arrsToSend(unionGrp->size());
    MCAuto<DataArrayIdType> globalNodeIds(getGlobalNodeIdsOnLocalMesh());
    const DataArrayDouble *array(srcFieldOnLocal->getArray());

    for (std::size_t i = 0; i < _glb_nodes_in_whole_per_trg.size(); ++i)
    {
        MCAuto<DataArrayIdType> locIdsToSendForCurTrgProc(
            globalNodeIds->findIdForEach(_glb_nodes_in_whole_per_trg[i]->begin(), _glb_nodes_in_whole_per_trg[i]->end())
        );
        arrsToSend[nbSrcProcs + i] =
            array->selectByTupleId(locIdsToSendForCurTrgProc->begin(), locIdsToSendForCurTrgProc->end());
    }
    all2allDA<double>(unionGrp, arrsToSend);
}

/*
 * check that input field is consistent with assigned local mesh.
 */
void
CFEMDECOneWaySource::checkMesh(MEDCouplingFieldDouble *field)
{
    if (!field)
        THROW_IK_EXCEPTION("Input field is nullptr !");
    if (!field->getMesh())
        THROW_IK_EXCEPTION("Input field has a nullptr geometrical support !");
    if (field->getMesh() != (MEDCouplingUMesh *)_master->getLocalMesh())
        THROW_IK_EXCEPTION("Input field lies on a geometrical support different than local mesh specified !");
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

static MCAuto<DataArrayDouble>
MatrixVectorMultiply(const std::vector<std::map<mcIdType, double>> &matrix, const DataArrayDouble *vectorArr)
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

/*!
 * Target side : Computes wholeMesh aggregation of srcMeshes from which duplicated cells are removed. globalNodeIds are
 * used to determine common cells across source processors.
 *
 *  \param [out] wholeMesh without duplication of cells
 *  \return for each source proc node Ids in \a wholeMesh referential of its contribution
 */
std::vector<MCAuto<DataArrayIdType>>
ComputeNodeIdsPerProc(
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
    mcIdType nbOfNodesWithoutDup(o2n->getMaxAbsValueInArray() + 1);
    wholeMesh = MEDCouplingUMesh::MergeUMeshes(FromVecAutoToVecOfConst<MEDCouplingUMesh>(srcMeshes));
    wholeMesh->renumberNodes(o2n->begin(), nbOfNodesWithoutDup);
    wholeMesh->checkConsistencyLight();
    // remove ghost cells
    std::vector<MCAuto<DataArrayIdType>> ret(nbOfSrcProcs);
    for (std::size_t i = 0; i < nbOfSrcProcs; ++i)
    {
        ret[i] = b->findIdForEach(srcGlobalNodeIds[i]->begin(), srcGlobalNodeIds[i]->end());
    }
    wholeMesh->zipConnectivityTraducer(0);
    return ret;
}
}  // namespace

/*!
 * Synchronized with CFEMDECOneWaySource::sendToTarget
 */
MCAuto<MEDCouplingFieldDouble>
CFEMDECOneWayTarget::receiveFromSource()
{
    constexpr char NB_COMPO_POS = 0, SIZE_OF_STR_POS = 1;
    smartSynchronize();
    int tmp[2] = {0 /*nbCompo*/, 0 /*sizeOfStr*/};
    std::unique_ptr<char[]> stringsToRecv;
    if (_to_master->getTargetGrp()->myRank() == 0)
    {
        int rkOrig(getUnionGrp()->translateRank(_to_master->getSourceGrp(), 0));
        MPI_Status status;
        getUnionGrp()->getCommInterface().recv(tmp, 2, MPI_INT32_T, rkOrig, 1234, *getUnionGrp()->getComm(), &status);
        stringsToRecv.reset(new char[tmp[SIZE_OF_STR_POS]]);
        getUnionGrp()->getCommInterface().recv(
            stringsToRecv.get(), tmp[SIZE_OF_STR_POS], MPI_INT8_T, rkOrig, 1235, *getUnionGrp()->getComm(), &status
        );
    }
    // broadcast nb compo over target group procs
    _to_master->getTargetGrp()->getCommInterface().broadcast(
        tmp, 2, MPI_INT32_T, 0, *_to_master->getTargetGrp()->getComm()
    );
    if (_to_master->getTargetGrp()->myRank() != 0)
    {
        stringsToRecv.reset(new char[tmp[SIZE_OF_STR_POS]]);
    }
    // broadcast strings of components over target group procs
    _to_master->getTargetGrp()->getCommInterface().broadcast(
        stringsToRecv.get(), tmp[SIZE_OF_STR_POS], MPI_INT8_T, 0, *_to_master->getTargetGrp()->getComm()
    );
    std::vector<std::string> compsInfo(
        FromCharToVectorString(stringsToRecv.get(), tmp[SIZE_OF_STR_POS], tmp[NB_COMPO_POS])
    );
    //
    int rkGathering(getUnionGrp()->translateRank(_to_master->getTargetGrp(), 0));
    std::vector<MCAuto<DataArrayDouble>> arrs;
    // broadcast of arrs on all procs of targets and assign srcValues
    MCAuto<DataArrayDouble> srcArr = DataArrayDouble::New();
    srcArr->alloc(_nb_nodes_src_mesh, tmp[0]);
    std::vector<MCAuto<DataArrayDouble>> dummy(getUnionGrp()->size());
    std::vector<MCAuto<DataArrayDouble>> perSrcProc(all2allDA<double>(getUnionGrp(), dummy));
    int szSrcGrp(_to_master->getSourceGrp()->size());
    for (int i = 0; i < szSrcGrp; ++i)
    {
        perSrcProc[i]->rearrange(tmp[NB_COMPO_POS]);
        srcArr->setPartOfValues3(
            perSrcProc[i],
            _src_rank_of_nodes_in_whole[i]->begin(),
            _src_rank_of_nodes_in_whole[i]->end(),
            0,
            tmp[NB_COMPO_POS],
            1,
            true
        );
    }
    // matrix * vector (srcArr)
    MCAuto<DataArrayDouble> res(MatrixVectorMultiply(_matrix, srcArr));
    res->setInfoOnComponents(compsInfo);
    //
    MCAuto<MEDCouplingFieldDouble> ret(MEDCouplingFieldDouble::New(ON_NODES_FE));
    ret->setArray(res);
    ret->setMesh(_master->getLocalMesh());
    ret->setName("Field");
    ret->checkConsistencyLight();
    return ret;
}

void
CFEMDECOneWayTarget::reinitializeOnNewMesh()
{
    CFEMDECOneWay::reinitializeOnNewMesh();
    _nb_nodes_src_mesh = 0;
}

void
CFEMDECOneWayTarget::computeMatrix(
    const std::vector<MCAuto<MEDCouplingUMesh>> &srcMeshes, const std::vector<MCAuto<DataArrayIdType>> &srcGlobalNodeIds
)
{
    MCAuto<MEDCouplingUMesh> wholeMesh;
    _src_rank_of_nodes_in_whole = ComputeNodeIdsPerProc(srcMeshes, srcGlobalNodeIds, wholeMesh);
    this->_nb_nodes_src_mesh = wholeMesh->getNumberOfNodes();
    // send srcNodeIds to source groups
    int nbProcUnion(getUnionGrp()->size());
    // compute matrix
    const double *coordsOfTrgMesh(_master->getLocalMesh()->getCoords()->begin());
    const mcIdType nbOfTrgPts(_master->getLocalMesh()->getNumberOfNodes());

    MEDCouplingFieldDiscretizationOnNodesFE::computeCrudeMatrix(
        wholeMesh, coordsOfTrgMesh, nbOfTrgPts, this->_matrix, this->getFEOptions()
    );
}

MCAuto<MEDCouplingUMesh>
CFEMDECOneWaySource::ReduceMesh(
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

CFEMDEC::CFEMDEC(ProcessorGroup &source_group, ProcessorGroup &target_group)
    : DisjointDECAbstract(source_group, target_group)
{
}

CFEMDEC::CFEMDEC(const std::set<int> &src_ids, const std::set<int> &trg_ids, const MPI_Comm &world_comm)
    : DisjointDECAbstract(src_ids, trg_ids, world_comm)
{
}

ProcessorGroup *
CFEMDEC::getReverseUnionGrp() const
{
    if (!_reverse_union_group.get())
    {
        // warning to order. Use fuseNotOrdered and not fuse for collectives
        _reverse_union_group.reset(getTargetGrp()->fuseNotOrdered(*getSourceGrp()));
    }
    return _reverse_union_group.get();
}

void
CFEMDEC::attachLocalMesh(MEDCouplingUMesh *mesh, DataArrayIdType *globalNodeIds)
{
    _local_mesh = MCAuto<MEDCouplingUMesh>::TakeRef(mesh);
    _global_node_ids = MCAuto<DataArrayIdType>::TakeRef(globalNodeIds);
    getEngine()->reinitializeOnNewMesh();
}

void
CFEMDEC::sendToTarget(MEDCouplingFieldDouble *srcFieldOnLocal)
{
    getEngine()->sendToTarget(srcFieldOnLocal);
}

void
CFEMDEC::sendToSource(MEDCouplingFieldDouble *trgFieldOnLocal)
{
    getReverseEngine()->sendToTarget(trgFieldOnLocal);
}

MCAuto<MEDCouplingFieldDouble>
CFEMDEC::receiveFromSource()
{
    return getEngine()->receiveFromSource();
}

MCAuto<MEDCouplingFieldDouble>
CFEMDEC::receiveFromTarget()
{
    return getReverseEngine()->receiveFromSource();
}

template <class T>
CFEMDECOneWay *
CFEMDEC::getEngineInternal()
{
    if (!_engine)
    {
        _engine = getClass(std::make_unique<T>(this));
        _engine->copyOptions(*this);
    }
    return _engine.get();
}

CFEMDECOneWay *
CFEMDEC::getEngine()
{
    return getEngineInternal<CFEMDECDirectAccess>();
}

template <class T>
CFEMDECOneWay *
CFEMDEC::getReverseEngineInternal()
{
    if (!_reverse_engine)
    {
        _reverse_engine = getReversedClass(std::make_unique<T>(this));
        _reverse_engine->copyOptions(*this);
    }
    return _reverse_engine.get();
}

CFEMDECOneWay *
CFEMDEC::getReverseEngine()
{
    return getReverseEngineInternal<CFEMDECReverseAccess>();
}

template <class SRC, class TRG>
std::unique_ptr<CFEMDECOneWay>
CFEMDEC::getClassInternal(std::unique_ptr<CFEMDECAccessToMaster> &&accessToMaster)
{
    if (getSourceGrp()->containsMyRank())
    {
        return std::make_unique<SRC>(std::move(accessToMaster), this);
    }
    else
    {
        return std::make_unique<TRG>(std::move(accessToMaster), this);
    }
}

std::unique_ptr<CFEMDECOneWay>
CFEMDEC::getClass(std::unique_ptr<CFEMDECAccessToMaster> &&accessToMaster)
{
    return this->getClassInternal<CFEMDECOneWaySource, CFEMDECOneWayTarget>(std::move(accessToMaster));
}

std::unique_ptr<CFEMDECOneWay>
CFEMDEC::getReversedClass(std::unique_ptr<CFEMDECAccessToMaster> &&accessToMaster)
{
    return this->getClassInternal<CFEMDECOneWayTarget, CFEMDECOneWaySource>(std::move(accessToMaster));
}
