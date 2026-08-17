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

#include "CommonCFEMDEC.hxx"
#include "MPITraits.hxx"
#include "MEDCouplingMemArray.hxx"

#include <cstdint>
#include <cstddef>
#include <cstring>
#include <memory>

using namespace MEDCoupling;

namespace
{

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
}  // namespace

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

    if (size != grp->size())
    {
        THROW_IK_EXCEPTION(
            "All2All : size of in input (" << size << ") is not equal to size of MPI group (" << grp->size() << ")"
        );
    }

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
MEDCoupling::all2allVector(MPIProcessorGroup *grp, const std::vector<typename std::vector<T>> &input)
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
MEDCoupling::all2allDA(MPIProcessorGroup *grp, const std::vector<MCAuto<typename MPITraits<T>::ArrayType>> &input)
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
BBTreeClosestSafe<spaceDim, mcIdType>
MEDCoupling::ShareBBTreesOfAllProcs(
    MPIProcessorGroup *unionGrp, const MEDCouplingUMesh *mesh, std::vector<BBTreeClosestSafe<spaceDim, mcIdType>> &ret
)
{
    const MPI_Comm *comm(unionGrp->getComm());
    MCAuto<DataArrayDouble> bbox(mesh->getBoundingBoxForBBTree());
    mcIdType nbCells(mesh->getNumberOfCells());
    const double *bboxPtr(bbox->begin());
    BBTreeClosestSafe<spaceDim, mcIdType> myTreeBase(bboxPtr, nullptr, 0, nbCells);
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
        ret[iProc] =
            std::move(BBTreeClosestSafe<spaceDim, mcIdType>::DeserializeCompact(structures[iProc], bboxes[iProc]));
    }
    return myTreeBase;
}
