
// Copyright (C) 2007-2026  CEA, EDF
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

#include "OverlapManagerDEC.hxx"
#include "MPIProcessorGroup.hxx"
#include "CommInterface.hxx"

#include <numeric>
#include <algorithm>
#include <iostream>

using namespace MEDCoupling;

constexpr bool debug = false;

struct PackItem
{
    mcIdType index;
    int origin;
    int local_pos;
};

static MPI_Datatype
mpiPackItemType()
{
    MPI_Datatype MPI_PACKITEM;
    int blocklen[3] = {1, 1, 1};
    MPI_Aint disp[3];
    MPI_Aint base;
    PackItem dummy;

    MPI_Get_address(&dummy, &base);
    MPI_Get_address(&dummy.index, &disp[0]);
    MPI_Get_address(&dummy.origin, &disp[1]);
    MPI_Get_address(&dummy.local_pos, &disp[2]);

    disp[0] -= base;
    disp[1] -= base;
    disp[2] -= base;

    MPI_Datatype types[3] = {MPI_ID_TYPE, MPI_INT, MPI_INT};
    MPI_Type_create_struct(3, blocklen, disp, types, &MPI_PACKITEM);
    MPI_Type_commit(&MPI_PACKITEM);
    return MPI_PACKITEM;
};

struct OwnerItem
{
    int local_pos;
    int owner;
    int multiplicity;
};

static MPI_Datatype
mpiOwnerItemType()
{
    MPI_Datatype MPI_OWNERITEM;
    int blocklen[3] = {1, 1, 1};
    MPI_Aint disp[3];
    MPI_Aint base;
    OwnerItem dummy;

    MPI_Get_address(&dummy, &base);
    MPI_Get_address(&dummy.local_pos, &disp[0]);
    MPI_Get_address(&dummy.owner, &disp[1]);
    MPI_Get_address(&dummy.multiplicity, &disp[2]);

    disp[0] -= base;
    disp[1] -= base;
    disp[2] -= base;

    MPI_Datatype types[3] = {MPI_INT, MPI_INT, MPI_INT};
    MPI_Type_create_struct(3, blocklen, disp, types, &MPI_OWNERITEM);
    MPI_Type_commit(&MPI_OWNERITEM);
    return MPI_OWNERITEM;
};

// Debug fonction
namespace DBG
{

// surcharge de l'opérateur <<
std::ostream &
operator<<(std::ostream &os, const PackItem &p)
{
    os << "{index=" << p.index << ", origin=" << p.origin << ", local_pos=" << p.local_pos << " }";
    return os;
}

// surcharge de l'opérateur <<
std::ostream &
operator<<(std::ostream &os, const OwnerItem &p)
{
    os << "{local_pos=" << p.local_pos << ", owner=" << p.owner << ", multiplicity=" << p.multiplicity << " }";
    return os;
}

template <typename T>
void
print(const std::vector<T> &vec, int rank)
{
    if (rank < 0)
    {
        std::cout << "vec(" << vec.size() << ")=[";
    }
    else
    {
        std::cout << "<rank=" << rank << "> vec(" << vec.size() << ")=[";
    }
    for (auto &v : vec)
    {
        std::cout << v << ", ";
    }
    std::cout << "]" << std::endl;
}

template <typename T>
void
print(const std::vector<std::vector<T>> &vec, int rank)
{
    std::cout << "<rank=" << rank << "> vec(" << vec.size() << ")=[";
    for (auto &v : vec)
    {
        std::cout << ", ";
        print(v, -1);
    }
    std::cout << "]" << std::endl;
}

template <typename K, typename V>
void
print(const std::unordered_map<K, V> &map, int rank)
{
    std::cout << "<rank=" << rank << "> map(" << map.size() << ")=[";

    for (const auto &kv : map)
    {
        std::cout << "{";
        std::cout << kv.first << ": " << kv.second;
        std::cout << "}, ";
    }

    std::cout << "]" << std::endl;
}

template <typename K, typename T>
void
print(const K N, const T *&vec, int rank)
{
    std::cout << "<rank=" << rank << "> vec(" << N << ")=[";
    for (int i = 0; i < N; i++)
    {
        std::cout << vec[i] << ", ";
    }
    std::cout << "]" << std::endl;
}
}  // namespace DBG

OverlapManagerDEC::OverlapManagerDEC()
    : _global_ids(nullptr),
      _own_originalField(false),
      _originalField(nullptr),
      _own_restrictedField(true),
      _restrictedField(nullptr),
      _local_group(nullptr),
      _mpi_group(nullptr),
      _hasLocalOverlap(false),
      _hasGlobalOverlap(false)
{
}

OverlapManagerDEC::~OverlapManagerDEC() { release(); }

void
OverlapManagerDEC::releaseRestrictedField()
{
    if (_own_restrictedField)
    {
        delete _restrictedField;
    }
    _restrictedField = nullptr;
    _own_restrictedField = true;
}

void
OverlapManagerDEC::releaseOriginalField()
{
    if (_own_originalField)
    {
        delete _originalField;
    }
    _originalField = nullptr;
    _own_originalField = false;
}

void
OverlapManagerDEC::releaseMPIGroup()
{
    if (_mpi_group != nullptr)
    {
        delete _mpi_group;
    }
    _mpi_group = nullptr;
};

void
OverlapManagerDEC::release()
{
    this->releaseMPIGroup();
    _local_group = nullptr;
    this->releaseRestrictedField();
    this->releaseOriginalField();
    _originalField = nullptr;
    _global_ids = nullptr;
    _owner.clear();
    _joints_recv.clear();
    _joints_send.clear();
    _kept.clear();
    _removed.clear();
    _o2r.clear();
}

void
OverlapManagerDEC::attachLocalField(const ParaFIELD *field, const MPI_Comm &comm)
{
    // cleaning between two call
    this->releaseMPIGroup();
    this->releaseRestrictedField();
    this->releaseOriginalField();

    _local_group = field->getProcGroup();
    // needd MPI group to avoid dead-lock
    _mpi_group = new MPIProcessorGroup(_local_group->getCommInterface(), _local_group->getProcIDs(), comm);

    _originalField = field;

    // compute ownership
    this->computeOwnership();

    if (_hasGlobalOverlap)
    {
        _restrictedField = this->restrict(_originalField);
        _own_restrictedField = true;
    }
};

/*********************************************************************
 * To knwon if we have to recompute the ownership
 *********************************************************************/
bool
OverlapManagerDEC::recomputeOwnership() const
{
    bool sameNum = false;

    if (_global_ids != nullptr && _originalField->returnGlobalNumbering() != nullptr)
    {
        const auto Nr = _originalField->returnGlobalNumbering()->getNbOfElems();
        const mcIdType *gln_ptr = _originalField->returnGlobalNumbering()->getConstPointer();

        const auto No = _global_ids->getNbOfElems();
        const mcIdType *glo_ptr = _global_ids->getConstPointer();

        sameNum = Nr == No;
        if (Nr == No)
        {
            sameNum = std::equal(gln_ptr, gln_ptr + Nr, glo_ptr);
        }
    }

    CommInterface ci = _local_group->getCommInterface();

    const MPI_Comm *comm = _mpi_group->getComm();

    if ((int)_joints_recv.size() != _mpi_group->size())
    {
        sameNum = false;
    }

    bool sameNumGl;

    ci.allReduce(&sameNum, &sameNumGl, 1, MPI_C_BOOL, MPI_LAND, *comm);

    if constexpr (debug)
    {
        const int rank = _mpi_group->myRank();
        std::cout << "<rank=" << rank << ">  recomputeOwnership : " << !sameNumGl << std::endl;
    }

    return !sameNumGl;
};

/*********************************************************************
 * Return the ParaFIELD to use. The original one if no overlapping is detected
 * else a the restricted field on the non-overlapped mesh.
 *********************************************************************/

const ParaFIELD *
OverlapManagerDEC::getRestrictedField(bool &originalOrComputed) const
{
    if (_hasGlobalOverlap)
    {
        originalOrComputed = true;
        return _restrictedField;
    }
    originalOrComputed = false;
    return _originalField;
};

void
OverlapManagerDEC::doNotOwnRestrictedFieldAnymore()
{
    _own_restrictedField = false;
}

void
OverlapManagerDEC::ownOriginalField()
{
    _own_originalField = true;
}

/*********************************************************************
 * Déterminer propriétaire de chaque cellule (simple modulo)
 *********************************************************************/
void
OverlapManagerDEC::computeOwnership()
{
    CommInterface commInterface = _local_group->getCommInterface();

    const MPI_Comm *comm = _mpi_group->getComm();

    const int rank = _mpi_group->myRank();
    const int size = _mpi_group->size();

    if constexpr (debug)
    {
        std::cout << "MPI: " << rank << " / " << size << std::endl;
    }

    // No overlap
    if (_originalField->returnGlobalNumbering() == nullptr || size <= 1)
    {
        if constexpr (debug)
        {
            std::cout << "<rank=" << rank << ">  No overlap " << std::endl;
        }

        _hasLocalOverlap = false;
        _hasGlobalOverlap = false;
        return;
    }

    if (!this->recomputeOwnership())
    {
        return;
    }

    _global_ids = _originalField->returnGlobalNumbering();

    const auto N = _global_ids->getNbOfElems();
    const mcIdType *gids_ptr = _global_ids->getConstPointer();

    if constexpr (debug)
    {
        std::cout << "<rank=" << rank << ">  global_idx " << std::endl;
        DBG::print(N, gids_ptr, rank);
    }

    // ===================================================================
    // 1) Pack local items for hash owners
    // ===================================================================
    std::vector<int> sendCounts(size, 0);
    for (int i = 0; i < N; i++)
    {
        sendCounts[gids_ptr[i] % size]++;
    }

    std::vector<int> sendDisp(size, 0);
    for (int i = 1; i < size; i++)
    {
        sendDisp[i] = sendDisp[i - 1] + sendCounts[i - 1];
    }

    const auto totalSend = sendDisp[size - 1] + sendCounts[size - 1];
    std::vector<PackItem> sendBuf(totalSend);
    std::vector<int> counter(size, 0);

    for (int i = 0; i < N; i++)
    {
        int h = int(gids_ptr[i] % size);
        int pos = sendDisp[h] + counter[h]++;
        sendBuf[pos] = {gids_ptr[i], rank, i};
    }

    if constexpr (debug)
    {
        std::cout << "<rank=" << rank << ">  sendCounts " << std::endl;
        DBG::print(sendCounts, rank);
        std::cout << "<rank=" << rank << ">  counter " << std::endl;
        DBG::print(counter, rank);
        std::cout << "<rank=" << rank << ">  sendBuf " << std::endl;
        DBG::print(sendBuf, rank);
    }

    // ===================================================================
    // 2) Exchange counts
    // ===================================================================
    std::vector<int> recvCounts(size);
    commInterface.allToAll(sendCounts.data(), 1, MPI_INT, recvCounts.data(), 1, MPI_INT, *comm);

    if constexpr (debug)
    {
        std::cout << "<rank=" << rank << ">  recvCounts " << std::endl;
        DBG::print(recvCounts, rank);
    }

    std::vector<int> recvDisp(size, 0);
    int totalRecv = 0;
    for (int i = 0; i < size; i++)
    {
        recvDisp[i] = totalRecv;
        totalRecv += recvCounts[i];
    }

    std::vector<PackItem> recvBuf(totalRecv);

    // ===================================================================
    // 3) Single Alltoallv : PackItem (3 elements: mcIdType,int,int)
    // ===================================================================
    // => Here we do it safely: create a contiguous MPI type.
    MPI_Datatype MPI_PACKITEM = mpiPackItemType();

    commInterface.allToAllV(
        sendBuf.data(),
        sendCounts.data(),
        sendDisp.data(),
        MPI_PACKITEM,
        recvBuf.data(),
        recvCounts.data(),
        recvDisp.data(),
        MPI_PACKITEM,
        *comm
    );

    MPI_Type_free(&MPI_PACKITEM);

    if constexpr (debug)
    {
        std::cout << "<rank=" << rank << ">  recvBuf " << std::endl;
        DBG::print(recvBuf, rank);
    }

    // ===================================================================
    // 4) Compute minimal owner per index
    // ===================================================================
    std::unordered_map<mcIdType, int> minOwner, multiplicity;
    minOwner.reserve(recvBuf.size());
    multiplicity.reserve(recvBuf.size());

    for (auto &it : recvBuf)
    {
        auto p = minOwner.find(it.index);
        if (p == minOwner.end())
        {
            minOwner[it.index] = it.origin;
        }
        else
        {
            p->second = std::min(p->second, it.origin);
        }
        auto m = multiplicity.find(it.index);
        if (m == multiplicity.end())
        {
            multiplicity[it.index] = 1;
        }
        else
        {
            multiplicity[it.index]++;
        }
    }

    if constexpr (debug)
    {
        std::cout << "<rank=" << rank << ">  minOwner " << std::endl;
        DBG::print(minOwner, rank);
        std::cout << "<rank=" << rank << ">  multiplicity " << std::endl;
        DBG::print(multiplicity, rank);
    }

    // ===================================================================
    // 5) Build OwnerItem to return to original ranks (same Alltoallv)
    // ===================================================================
    std::vector<OwnerItem> retBuf(totalRecv);

    for (int i = 0; i < totalRecv; i++)
    {
        const auto index = recvBuf[i].index;
        retBuf[i] = {recvBuf[i].local_pos, minOwner[index], multiplicity[index]};
    }

    if constexpr (debug)
    {
        std::cout << "<rank=" << rank << ">  retBuf " << std::endl;
        DBG::print(retBuf, rank);
    }

    // ===================================================================
    // 6) Single Alltoallv: return OwnerItem to origin ranks
    // ===================================================================
    // => Here we do it safely: create a contiguous MPI type.
    MPI_Datatype MPI_OWNERITEM = mpiOwnerItemType();

    std::vector<OwnerItem> recvOwners(totalSend);

    commInterface.allToAllV(
        retBuf.data(),
        recvCounts.data(),
        recvDisp.data(),
        MPI_OWNERITEM,
        recvOwners.data(),
        sendCounts.data(),
        sendDisp.data(),
        MPI_OWNERITEM,
        *comm
    );

    MPI_Type_free(&MPI_OWNERITEM);

    if constexpr (debug)
    {
        std::cout << "<rank=" << rank << ">  recvOwners " << std::endl;
        DBG::print(recvOwners, rank);
    }

    // ===================================================================
    // 8) Reconstruct ordered result
    // ===================================================================
    _owner.resize(N);
    _shared.clear();
    _shared.reserve(N / 4);

    for (auto &o : recvOwners)
    {
        _owner[o.local_pos] = o.owner;
        if (o.owner != rank)
        {
            _hasLocalOverlap = true;
        }

        if (o.owner == rank && o.multiplicity > 1)
        {
            _shared.push_back(o.local_pos);
        }
    }

    commInterface.allReduce(&_hasLocalOverlap, &_hasGlobalOverlap, 1, MPI_C_BOOL, MPI_LOR, *comm);

    if constexpr (debug)
    {
        std::cout << "<rank=" << rank << ">  _owner " << std::endl;
        DBG::print(_owner, rank);
        std::cout << "<rank=" << rank << ">  _shared " << std::endl;
        DBG::print(_shared, rank);
    }

    _kept.clear();
    _kept.reserve(N);
    _removed.clear();
    _removed.reserve(N / 2);
    _o2r.clear();
    _o2r.resize(N, -1);

    mcIdType id_res = 0;
    for (mcIdType i = 0; i < N; i++)
    {
        if (_owner[i] == rank)
        {
            _kept.push_back(i);
            _o2r[i] = id_res;
            id_res++;
        }
        else
            _removed.push_back(i);
    }

    if constexpr (debug)
    {
        std::cout << "<rank=" << rank << ">  _removed " << std::endl;
        DBG::print(_removed, rank);
        std::cout << "<rank=" << rank << ">  _kept " << std::endl;
        DBG::print(_kept, rank);
    }

    // ===================================================================
    // 9) Create connection beetween cells
    // ===================================================================

    std::unordered_map<mcIdType, mcIdType> sharedGlobNum;
    for (const auto &lid : _shared)
    {
        const mcIdType g_id = gids_ptr[lid];
        sharedGlobNum[g_id] = lid;
    }

    if constexpr (debug)
    {
        std::cout << "<rank=" << rank << ">  sharedGlobNum " << std::endl;
        DBG::print(sharedGlobNum, rank);
    }

    /*===========================================*/
    /* 9.1) build _joints_recv: removed -> owner        */
    /*===========================================*/
    _joints_recv.clear();
    _joints_recv.resize(size);

    for (const auto &lid : _removed)
    {
        const mcIdType gid = gids_ptr[lid];
        const int owner = _owner[lid];
        _joints_recv[owner].push_back(lid);
    }

    /*===========================================*/
    /* 9.2) exchange request counts                */
    /*===========================================*/
    std::vector<int> sendCounts2(size, 0), recvCounts2(size, 0);

    for (int r = 0; r < size; r++)
    {
        sendCounts2[r] = (int)_joints_recv[r].size();
    }

    int totalSend2 = std::accumulate(sendCounts2.begin(), sendCounts2.end(), 0);

    if constexpr (debug)
    {
        std::cout << "<rank=" << rank << ">  sendCounts " << std::endl;
        DBG::print(sendCounts2, rank);
    }

    commInterface.allToAll(sendCounts2.data(), 1, MPI_INT, recvCounts2.data(), 1, MPI_INT, *comm);

    int totalRecv2 = std::accumulate(recvCounts2.begin(), recvCounts2.end(), 0);

    std::vector<int> sdisp2(size, 0), rdisp2(size, 0);
    for (int i = 1; i < size; i++)
    {
        sdisp2[i] = sdisp2[i - 1] + sendCounts2[i - 1];
        rdisp2[i] = rdisp2[i - 1] + recvCounts2[i - 1];
    }

    if constexpr (debug)
    {
        std::cout << "<rank=" << rank << ">  sdisp2 " << std::endl;
        DBG::print(sdisp2, rank);
        std::cout << "<rank=" << rank << ">  rdisp2 " << std::endl;
        DBG::print(rdisp2, rank);
    }

    /*===========================================*/
    /* 9.3) exchange requested global IDs          */
    /*===========================================*/
    std::vector<mcIdType> sendBuf2;
    sendBuf2.reserve(totalSend2);

    for (int r = 0; r < size; r++)
    {
        for (const auto &lid : _joints_recv[r])
        {
            const mcIdType gid = gids_ptr[lid];
            sendBuf2.push_back(gid);
        }
    }

    if constexpr (debug)
    {
        std::cout << "<rank=" << rank << ">  sendBuf2 " << std::endl;
        DBG::print(sendBuf2, rank);
    }

    std::vector<mcIdType> recvBuf2(totalRecv2);

    commInterface.allToAllV(
        sendBuf2.data(),
        sendCounts2.data(),
        sdisp2.data(),
        MPI_ID_TYPE,
        recvBuf2.data(),
        recvCounts2.data(),
        rdisp2.data(),
        MPI_ID_TYPE,
        *comm
    );

    if constexpr (debug)
    {
        std::cout << "<rank=" << rank << ">  recvBuf2" << std::endl;
        DBG::print(recvBuf2, rank);
    }

    _joints_send.clear();
    _joints_send.resize(size);

    mcIdType offset = 0;
    for (int r = 0; r < size; r++)
    {
        const mcIdType nb_val = recvCounts2[r];
        for (mcIdType i = 0; i < nb_val; i++)
        {
            const mcIdType gid = recvBuf2[offset];
            _joints_send[r].push_back(sharedGlobNum[gid]);
            offset++;
        }
    }

    if constexpr (debug)
    {
        std::cout << "<rank=" << rank << ">  _joints_send" << std::endl;
        DBG::print(_joints_send, rank);
        std::cout << "<rank=" << rank << ">  _joints_recv" << std::endl;
        DBG::print(_joints_recv, rank);
    }
}

/*********************************************************************
 * Restreindre un maillage
 *********************************************************************/
ParaMESH *OverlapManagerDEC::restrict(const MEDCouplingMesh *mesh) const
{
    ParaMESH *paramesh = new ParaMESH(
        static_cast<MEDCouplingPointSet *>(const_cast<MEDCouplingMesh *>(mesh)), *_local_group, mesh->getName()
    );
    return paramesh;
}

/*********************************************************************
 * Restreindre un ParaFIELD
 *********************************************************************/
ParaFIELD *OverlapManagerDEC::restrict(const ParaFIELD *pfield) const
{
    MEDCouplingFieldDouble *field = _originalField->getField();

    MCAuto<MEDCouplingFieldDouble> field_res(field->buildSubPart(_kept.data(), _kept.data() + _kept.size()));
    field_res->zipCoords();

    if constexpr (debug)
    {
        int rank = _mpi_group->myRank();

        std::cout << "<rank=" << rank << "> Number of cells: before = " << field->getMesh()->getNumberOfCells()
                  << ", after = " << field_res->getMesh()->getNumberOfCells() << std::endl;
        std::cout << "<rank=" << rank << "> Number of nodes: before = " << field->getMesh()->getNumberOfNodes()
                  << ", after = " << field_res->getMesh()->getNumberOfNodes() << std::endl;
        std::cout << "<rank=" << rank << "> Number of values: before = " << field->getNumberOfValues()
                  << ", after = " << field_res->getNumberOfValues() << std::endl;
    }

    ParaMESH *pmesh_res = this->restrict(field_res->getMesh());

    ParaFIELD *pfield_res = new ParaFIELD(field_res, pmesh_res, *_local_group);
    pfield_res->setOwnSupport(true);

    return pfield_res;
}

/*********************************************************************
 * Pré- traitement : Update array of restrictedArray
 * since originalField can be update without call to attachLocalField
 *********************************************************************/
void
OverlapManagerDEC::updateLocalArray()
{
    if (!_hasGlobalOverlap)
    {
        return;
    }

    const auto arr_res = _restrictedField->getField()->getArray();
    const auto arr_orig = _originalField->getField()->getArray();

    const mcIdType nbCmp = _originalField->getField()->getNumberOfComponents();

    if constexpr (debug)
    {
        int rank = _mpi_group->myRank();
        DBG::print(_kept, rank);
        std::cout << "<rank=" << rank << "> ArrayRes: " << arr_res->getNumberOfTuples() << " x "
                  << arr_res->getNumberOfComponents() << std::endl;
        std::cout << "<rank=" << rank << "> ArrayOri: " << arr_orig->getNumberOfTuples() << " x "
                  << arr_orig->getNumberOfComponents() << std::endl;
    }

    // Check
    const char msg[] = "OverlapManagerDEC::updateLocalArray";
    arr_orig->checkNbOfTuples(_owner.size(), msg);
    arr_res->checkNbOfTuplesAndComp(_kept.size(), nbCmp, msg);
    double *arr_res_ptr(arr_res->getPointer());
    const double *arr_orig_ptr(arr_orig->begin());

    const mcIdType nbTuple = arr_res->getNumberOfTuples();

    for (mcIdType i = 0; i < nbTuple; i++)
    {
        const mcIdType tup_ori = _kept[i];
        for (mcIdType j = 0; j < nbCmp; j++)
        {
            arr_res_ptr[i * nbCmp + j] = arr_orig_ptr[nbCmp * tup_ori + j];
        }
    }

    // Mauvais sens. Je ne sais pas si ça existe
    // arr_res->setPartOfValues3(arr_orig, _kept.data(), _kept.data() + _kept.size(), 0, nbCmp, 1);
};

/*********************************************************************
 * Post traitement : Update ghost values
 *********************************************************************/
void
OverlapManagerDEC::synchronizeGhosts()
{
    if constexpr (debug)
    {
        std::cout << "Overlap: " << _hasGlobalOverlap << std::endl;
    }

    if (!_hasGlobalOverlap)
    {
        return;
    }

    CommInterface ci = _local_group->getCommInterface();
    const MPI_Comm *comm = _mpi_group->getComm();

    int rank = _mpi_group->myRank();
    int size = _mpi_group->size();

    if constexpr (debug)
    {
        std::cout << "MPI: " << rank << " / " << size << std::endl;
    }

    ParaFIELD *pf = _restrictedField;
    auto *res_arr = pf->getField()->getArray();
    const int ncpt = (int)res_arr->getNumberOfComponents();

    if constexpr (debug)
    {
        const auto N = res_arr->getNbOfElems();
        const double *data = res_arr->getConstPointer();

        if constexpr (debug)
        {
            std::cout << "<rank=" << rank << ">  _restrictedField " << std::endl;
            DBG::print(N, data, rank);
        }
    }

    /*===========================================*/
    /* 1) prepare counts                         */
    /*===========================================*/
    std::vector<int> sendCounts(size, 0), recvCounts(size, 0);

    for (int r = 0; r < size; r++)
    {
        sendCounts[r] = (int)_joints_send[r].size() * ncpt;
        recvCounts[r] = (int)_joints_recv[r].size() * ncpt;
    }

    if constexpr (debug)
    {
        std::cout << "<rank=" << rank << ">  sendCounts " << std::endl;
        DBG::print(sendCounts, rank);
        std::cout << "<rank=" << rank << ">  recvCounts " << std::endl;
        DBG::print(recvCounts, rank);
    }

    const int totalRecv = std::accumulate(recvCounts.begin(), recvCounts.end(), 0);

    std::vector<int> sdisp(size, 0), rdisp(size, 0);
    for (int i = 1; i < size; i++)
    {
        sdisp[i] = sdisp[i - 1] + sendCounts[i - 1];
        rdisp[i] = rdisp[i - 1] + recvCounts[i - 1];
    }

    if constexpr (debug)
    {
        std::cout << "<rank=" << rank << ">  sdisp " << std::endl;
        DBG::print(sdisp, rank);
        std::cout << "<rank=" << rank << ">  rdisp " << std::endl;
        DBG::print(rdisp, rank);
    }

    /*===========================================*/
    /* 2) Prepare replies: owners send values     */
    /*===========================================*/
    std::vector<double> buffSend;
    buffSend.reserve(totalRecv * ncpt);

    for (const auto &js : _joints_send)
    {
        for (const auto &lid : js)
        {
            const mcIdType lid_res = _o2r[lid];
            for (int c = 0; c < ncpt; c++)
            {
                buffSend.push_back(res_arr->getIJ(lid_res, c));
            }
        }
    }

    if constexpr (debug)
    {
        std::cout << "<rank=" << rank << ">  buffSend " << std::endl;
        DBG::print(buffSend, rank);
    }

    std::vector<double> buffRecv(totalRecv * ncpt);

    ci.allToAllV(
        buffSend.data(),
        sendCounts.data(),
        sdisp.data(),
        MPI_DOUBLE,
        buffRecv.data(),
        recvCounts.data(),
        rdisp.data(),
        MPI_DOUBLE,
        *comm
    );

    if constexpr (debug)
    {
        std::cout << "<rank=" << rank << ">  buffRecv " << std::endl;
        DBG::print(buffRecv, rank);
    }

    /*===========================================*/
    /* 3) Reconstruct original field            */
    /*===========================================*/

    auto *orig_arr = _originalField->getField()->getArray();

    auto orig_arr_ptr(orig_arr->getPointer());
    const auto *res_arr_ptr(res_arr->begin());
    /* copy owned values */
    for (const auto &lid : _kept)
    {
        const mcIdType lid_res = _o2r[lid];
        for (int c = 0; c < ncpt; c++)
        {
            orig_arr_ptr[lid * ncpt + c] = res_arr_ptr[lid_res * ncpt + c];
        }
    }

    /* copy ghost values */
    mcIdType offset = 0;
    for (const auto &jr : _joints_recv)
    {
        for (const auto &lid : jr)
        {
            for (int c = 0; c < ncpt; c++)
            {
                orig_arr_ptr[lid * ncpt + c] = buffRecv[offset++];
            }
        }
    }

    if constexpr (debug)
    {
        const auto No = orig_arr->getNbOfElems();
        const double *datao = orig_arr->getConstPointer();

        if constexpr (debug)
        {
            std::cout << "<rank=" << rank << ">  _originalField " << std::endl;
            DBG::print(No, datao, rank);
        }
    }
};
