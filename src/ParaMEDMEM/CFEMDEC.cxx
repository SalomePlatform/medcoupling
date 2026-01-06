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

using namespace MEDCoupling;

MPIProcessorGroup *
CFEMDECAccessToMaster::getUnionGrp() const
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

MCAuto<MEDCouplingUMesh>
CFEMDECOneWay::getLocalMesh() const
{
    return _master->getLocalMesh();
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
    MPIProcessorGroup *grp(_to_master->getUnionGrp());
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

void
CFEMDECOneWay::synchronize()
{
    int rkGathering(_to_master->getUnionGrp()->translateRank(_to_master->getTargetGrp(), 0));
    // MCAuto< DataArray>
    std::vector<double> td;
    std::vector<mcIdType> ti;
    std::vector<std::string> ts;
    MCAuto<MEDCouplingUMesh> locMesh(_master->getLocalMesh());
    MCAuto<DataArrayIdType> bi;
    MCAuto<DataArrayDouble> bd;
    int spaceDim(locMesh->getSpaceDimension());
    if (_to_master->getSourceGrp()->containsMyRank())
    {
        locMesh->getTinySerializationInformation(td, ti, ts);
        {
            DataArrayIdType *biTmp(nullptr);
            DataArrayDouble *bdTmp(nullptr);
            locMesh->serialize(biTmp, bdTmp);
            bi = biTmp;
            bd = bdTmp;
        }
    }
    else
    {
        std::vector<double> tdFake;
        std::vector<mcIdType> tiFake;
        locMesh->getTinySerializationInformation(tdFake, tiFake, ts);
    }
    // retrieve all srcMeshes on proc0 of target ( see EDF31187 )
    std::vector<std::vector<mcIdType>> tiGathered(
        gatherTArrayOnProcInternal<mcIdType>(_to_master->getUnionGrp(), rkGathering, ti.data(), ti.size())
    );
    std::vector<std::vector<double>> tdGathered(
        gatherTArrayOnProcInternal<double>(_to_master->getUnionGrp(), rkGathering, td.data(), td.size())
    );
    std::vector<MCAuto<DataArrayIdType>> biGathered(
        gatherTArrayOnProc<mcIdType>(_to_master->getUnionGrp(), rkGathering, bi)
    );
    std::vector<MCAuto<DataArrayDouble>> bdGathered(
        gatherTArrayOnProc<double>(_to_master->getUnionGrp(), rkGathering, bd)
    );
    std::vector<MCAuto<MEDCouplingUMesh>> srcMeshes;
    if (_to_master->getUnionGrp()->myRank() == rkGathering)
    {
        std::size_t nbPartsExp(_to_master->getSourceGrp()->size());
        if (tiGathered.size() != nbPartsExp || tdGathered.size() != nbPartsExp || biGathered.size() != nbPartsExp ||
            bdGathered.size() != nbPartsExp)
            THROW_IK_EXCEPTION("All gathered parts must have a size equal to src group ( " << nbPartsExp << " ) !");
        srcMeshes.resize(nbPartsExp);
        for (std::size_t i = 0; i < nbPartsExp; i++)
        {
            // assume that local spacedim is equal to remote spacedim
            bdGathered[i]->rearrange(spaceDim);
            MCAuto<MEDCouplingUMesh> mesh(MEDCouplingUMesh::New());
            mesh->unserialization(tdGathered[i], tiGathered[i], biGathered[i], bdGathered[i], ts);
            srcMeshes[i] = mesh;
        }
    }
    // global node ids
    MCAuto<DataArrayIdType> gni;
    if (_to_master->getSourceGrp()->containsMyRank())
    {
        gni = _master->getGlobalNodeIdsOnLocalMesh();
    }
    std::vector<MCAuto<DataArrayIdType>> srcGlobalNodeIds(
        gatherTArrayOnProc<mcIdType>(_to_master->getUnionGrp(), rkGathering, gni)
    );
    reinitializeOnNewMesh();
    computeMatrix(srcMeshes, srcGlobalNodeIds);
    _local_mesh_stamp_on_sync = locMesh;
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
    int nbSrcProcs(_to_master->getSourceGrp()->size());
    {
        int tmp[2] = {nbCompo, sizeOfStr};
        std::unique_ptr<int[]> nbCompos(new int[2 * nbSrcProcs]);
        _to_master->getSourceGrp()->getCommInterface().allGather(
            tmp, 2, MPI_INT32_T, nbCompos.get(), 2, MPI_INT32_T, *_to_master->getSourceGrp()->getComm()
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
    int rkDest(_to_master->getUnionGrp()->translateRank(_to_master->getTargetGrp(), 0));
    if (_to_master->getSourceGrp()->myRank() == 0)
    {
        {
            int tmp[2] = {nbCompo, sizeOfStr};
            _to_master->getUnionGrp()->getCommInterface().send(
                tmp, 2, MPI_INT32_T, rkDest, 1234, *_to_master->getUnionGrp()->getComm()
            );
        }
        {
            std::unique_ptr<char[]> stringsToSend(FromVectorStringToChar(comps));
            _to_master->getUnionGrp()->getCommInterface().send(
                stringsToSend.get(), sizeOfStr, MPI_INT8_T, rkDest, 1235, *_to_master->getUnionGrp()->getComm()
            );
        }
    }
    gatherTArrayOnProc<double>(_to_master->getUnionGrp(), rkDest, srcFieldOnLocal->getArray());
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

static MCAuto<DataArrayDouble>
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
        int rkOrig(_to_master->getUnionGrp()->translateRank(_to_master->getSourceGrp(), 0));
        MPI_Status status;
        _to_master->getUnionGrp()->getCommInterface().recv(
            tmp, 2, MPI_INT32_T, rkOrig, 1234, *_to_master->getUnionGrp()->getComm(), &status
        );
        stringsToRecv.reset(new char[tmp[SIZE_OF_STR_POS]]);
        _to_master->getUnionGrp()->getCommInterface().recv(
            stringsToRecv.get(),
            tmp[SIZE_OF_STR_POS],
            MPI_INT8_T,
            rkOrig,
            1235,
            *_to_master->getUnionGrp()->getComm(),
            &status
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
    int rkGathering(_to_master->getUnionGrp()->translateRank(_to_master->getTargetGrp(), 0));
    std::vector<MCAuto<DataArrayDouble>> arrs;
    {
        MCAuto<DataArrayDouble> arr;
        arrs = gatherTArrayOnProc<double>(_to_master->getUnionGrp(), rkGathering, arr);
    }
    // broadcast of arrs on all procs of targets and assign srcValues
    MCAuto<DataArrayDouble> srcArr = DataArrayDouble::New();
    srcArr->alloc(_nb_nodes_src_mesh, tmp[0]);
    //
    int szSrcGrp(_to_master->getSourceGrp()->size());
    if (_to_master->getTargetGrp()->myRank() != 0)
        arrs.resize(szSrcGrp);
    for (int i = 0; i < szSrcGrp; ++i)
    {
        MCAuto<DataArrayDouble> arrFieldSrc(bCastTArrayFromProc<double>(_to_master->getTargetGrp(), arrs[i], 0));
        arrFieldSrc->rearrange(tmp[NB_COMPO_POS]);
        srcArr->setPartOfValuesBase3(
            arrFieldSrc,
            _src_rank_of_nodes_in_whole[i]->begin(),
            _src_rank_of_nodes_in_whole[i]->end(),
            0,
            tmp[NB_COMPO_POS],
            1
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

static std::vector<MCAuto<DataArrayIdType>>
ComputeNodeIdsPerProc(
    const std::vector<MCAuto<MEDCouplingUMesh>> &srcMeshes,
    const std::vector<MCAuto<DataArrayIdType>> &srcGlobalNodeIds,
    MCAuto<MEDCouplingUMesh> &wholeMesh
)
{
    MCAuto<DataArrayIdType> o2n(DataArrayIdType::Aggregate(FromVecAutoToVecOfConst<DataArrayIdType>(srcGlobalNodeIds)));
    MCAuto<DataArrayIdType> b(o2n->buildUniqueNotSorted());
    MCAuto<MapKeyVal<mcIdType, mcIdType>> zeMap((b->invertArrayN2O2O2NOptimized()));
    o2n->transformWithIndArr(*zeMap);
    mcIdType nbOfNodesWithoutDup(o2n->getMaxAbsValueInArray() + 1);
    wholeMesh = MEDCouplingUMesh::MergeUMeshes(FromVecAutoToVecOfConst<MEDCouplingUMesh>(srcMeshes));
    wholeMesh->renumberNodes(o2n->begin(), nbOfNodesWithoutDup);
    wholeMesh->checkConsistencyLight();
    // remove ghost cells
    std::size_t nbOfSrcProcs(srcGlobalNodeIds.size());
    std::vector<MCAuto<DataArrayIdType>> ret(nbOfSrcProcs);
    for (std::size_t i = 0; i < nbOfSrcProcs; ++i)
    {
        ret[i] = b->findIdForEach(srcGlobalNodeIds[i]->begin(), srcGlobalNodeIds[i]->end());
    }
    wholeMesh->zipConnectivityTraducer(0);
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
    bool isProc0Target(_to_master->getTargetGrp()->myRank() == 0);
    MCAuto<MEDCouplingUMesh> wholeMesh;
    if (isProc0Target)
        _src_rank_of_nodes_in_whole = ComputeNodeIdsPerProc(srcMeshes, srcGlobalNodeIds, wholeMesh);
    // broadcast wholeMesh over all target procs
    std::vector<double> td;
    std::vector<mcIdType> ti;
    std::vector<std::string> ts;
    MCAuto<DataArrayIdType> bi;
    MCAuto<DataArrayDouble> bd;
    if (isProc0Target)
    {
        wholeMesh->getTinySerializationInformation(td, ti, ts);
        {
            DataArrayIdType *biTmp(nullptr);
            DataArrayDouble *bdTmp(nullptr);
            wholeMesh->serialize(biTmp, bdTmp);
            bi = biTmp;
            bd = bdTmp;
        }
    }
    else
    {
        std::vector<double> tdFake;
        std::vector<mcIdType> tiFake;
        _master->getLocalMesh()->getTinySerializationInformation(tdFake, tiFake, ts);
    }
    MPIProcessorGroup *trgGrp(_to_master->getTargetGrp());
    ti = bCastTArrayFromProcInternal<mcIdType>(trgGrp, ti.data(), ti.size(), 0);
    td = bCastTArrayFromProcInternal<double>(trgGrp, td.data(), td.size(), 0);
    bi = bCastTArrayFromProc<mcIdType>(trgGrp, bi, 0);
    bd = bCastTArrayFromProc<double>(trgGrp, bd, 0);
    bd->rearrange(_master->getLocalMesh()->getSpaceDimension());
    if (!isProc0Target)
    {
        wholeMesh = MEDCouplingUMesh::New();
        wholeMesh->unserialization(td, ti, bi, bd, ts);
    }
    this->_nb_nodes_src_mesh = wholeMesh->getNumberOfNodes();
    // broadcast _src_rank_of_nodes_in_whole over all target procs
    int nbSrcProcs(_to_master->getSourceGrp()->size());
    if (!isProc0Target)
    {
        _src_rank_of_nodes_in_whole.resize(nbSrcProcs);
    }
    for (int i = 0; i < nbSrcProcs; ++i)
    {
        MCAuto<DataArrayIdType> srcRk(_src_rank_of_nodes_in_whole[i]);
        MCAuto<DataArrayIdType> srcRk2(bCastTArrayFromProc<mcIdType>(trgGrp, srcRk, 0));
        if (!isProc0Target)
            _src_rank_of_nodes_in_whole[i] = srcRk2;
    }
    // compute matrix
    const double *coordsOfTrgMesh(_master->getLocalMesh()->getCoords()->begin());
    const mcIdType nbOfTrgPts(_master->getLocalMesh()->getNumberOfNodes());

    MEDCouplingFieldDiscretizationOnNodesFE::computeCrudeMatrix(
        wholeMesh, coordsOfTrgMesh, nbOfTrgPts, this->_matrix, this->getFEOptions()
    );
}

CFEMDEC::CFEMDEC(ProcessorGroup &source_group, ProcessorGroup &target_group)
    : DisjointDECAbstract(source_group, target_group)
{
}

CFEMDEC::CFEMDEC(const std::set<int> &src_ids, const std::set<int> &trg_ids, const MPI_Comm &world_comm)
    : DisjointDECAbstract(src_ids, trg_ids, world_comm)
{
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
