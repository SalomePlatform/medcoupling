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

#pragma once

#include "DisjointDECAbstract.hxx"
#include "InterpolationOptions.hxx"

#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MCAuto.hxx"

#include <map>
#include <memory>
#include <vector>

namespace MEDCoupling
{

class CFEMDEC;
class MPIProcessorGroup;

class CFEMDECAccessToMaster
{
   public:
    CFEMDECAccessToMaster(CFEMDEC *master) : _master(master) {}
    MPIProcessorGroup *getUnionGrp() const;
    virtual MPIProcessorGroup *getSourceGrp() const = 0;
    virtual MPIProcessorGroup *getTargetGrp() const = 0;

   protected:
    CFEMDEC *_master = nullptr;
};

class CFEMDECDirectAccess : public CFEMDECAccessToMaster
{
   public:
    CFEMDECDirectAccess(CFEMDEC *master) : CFEMDECAccessToMaster(master) {}
    MPIProcessorGroup *getSourceGrp() const override;
    MPIProcessorGroup *getTargetGrp() const override;
};

class CFEMDECReverseAccess : public CFEMDECAccessToMaster
{
   public:
    CFEMDECReverseAccess(CFEMDEC *master) : CFEMDECAccessToMaster(master) {}
    MPIProcessorGroup *getSourceGrp() const override;
    MPIProcessorGroup *getTargetGrp() const override;
};

class CFEMDECOneWay : public INTERP_KERNEL::InterpolationOptions
{
   public:
    CFEMDECOneWay(std::unique_ptr<CFEMDECAccessToMaster> &&accessToMaster, CFEMDEC *master)
        : _to_master(std::move(accessToMaster)), _master(master)
    {
    }
    virtual void reinitializeOnNewMesh();
    //
    MPIProcessorGroup *getSourceGrp() const;
    MPIProcessorGroup *getTargetGrp() const;
    //
    void smartSynchronize();
    //
    virtual void sendToTarget(MEDCouplingFieldDouble *srcFieldOnLocal) = 0;
    virtual MCAuto<MEDCouplingFieldDouble> receiveFromSource() = 0;
    virtual void computeMatrix(
        const std::vector<MCAuto<MEDCouplingUMesh> > &srcMeshes,
        const std::vector<MCAuto<DataArrayIdType> > &srcGlobalNodeIds
    ) = 0;

   private:
    MCAuto<MEDCouplingUMesh> getLocalMesh() const;
    bool isMatrixToRecompute();
    void synchronize();

   protected:
    std::unique_ptr<CFEMDECAccessToMaster> _to_master;
    CFEMDEC *_master = nullptr;
    MEDCouplingUMesh *_local_mesh_stamp_on_sync = nullptr;
};

class CFEMDECOneWaySource : public CFEMDECOneWay
{
   public:
    CFEMDECOneWaySource(std::unique_ptr<CFEMDECAccessToMaster> &&accessToMaster, CFEMDEC *master)
        : CFEMDECOneWay(std::move(accessToMaster), master)
    {
    }
    void sendToTarget(MEDCouplingFieldDouble *srcFieldOnLocal) override;
    MCAuto<MEDCouplingFieldDouble> receiveFromSource() override { return MCAuto<MEDCouplingFieldDouble>(); }
    void computeMatrix(
        const std::vector<MCAuto<MEDCouplingUMesh> > &srcMeshes,
        const std::vector<MCAuto<DataArrayIdType> > &srcGlobalNodeIds
    ) override
    {
    }

   private:
    void checkMesh(MEDCouplingFieldDouble *field);
};

class CFEMDECOneWayTarget : public CFEMDECOneWay
{
   public:
    CFEMDECOneWayTarget(std::unique_ptr<CFEMDECAccessToMaster> &&accessToMaster, CFEMDEC *master)
        : CFEMDECOneWay(std::move(accessToMaster), master)
    {
    }
    void reinitializeOnNewMesh();
    void sendToTarget(MEDCouplingFieldDouble *srcFieldOnLocal) override {}
    MCAuto<MEDCouplingFieldDouble> receiveFromSource() override;
    void computeMatrix(
        const std::vector<MCAuto<MEDCouplingUMesh> > &srcMeshes,
        const std::vector<MCAuto<DataArrayIdType> > &srcGlobalNodeIds
    ) override;

   private:
    mcIdType _nb_nodes_src_mesh = 0;
    std::vector<MCAuto<DataArrayIdType> > _src_rank_of_nodes_in_whole;
    std::vector<std::map<mcIdType, double> > _matrix;
};

class CFEMDEC : public DisjointDECAbstract, public INTERP_KERNEL::InterpolationOptions
{
   public:
    CFEMDEC() {}
    CFEMDEC(ProcessorGroup &source_group, ProcessorGroup &target_group);
    CFEMDEC(const std::set<int> &src_ids, const std::set<int> &trg_ids, const MPI_Comm &world_comm = MPI_COMM_WORLD);
    //
    void attachLocalMesh(MEDCouplingUMesh *mesh, DataArrayIdType *globalNodeIds);
    MCAuto<MEDCouplingUMesh> getLocalMesh() const { return _local_mesh; }
    MCAuto<DataArrayIdType> getGlobalNodeIdsOnLocalMesh() const { return _global_node_ids; }
    void sendToTarget(MEDCouplingFieldDouble *srcFieldOnLocal);
    void sendToSource(MEDCouplingFieldDouble *trgFieldOnLocal);
    MCAuto<MEDCouplingFieldDouble> receiveFromSource();
    MCAuto<MEDCouplingFieldDouble> receiveFromTarget();
    //
    void synchronize() override {}
    void sendRecvData(bool way = true) override {}

   private:
    CFEMDECOneWay *getEngine();
    CFEMDECOneWay *getReverseEngine();
    std::unique_ptr<CFEMDECOneWay> getClass(std::unique_ptr<CFEMDECAccessToMaster> &&accessToMaster);
    std::unique_ptr<CFEMDECOneWay> getReversedClass(std::unique_ptr<CFEMDECAccessToMaster> &&accessToMaster);

   private:
    template <class T>
    CFEMDECOneWay *getEngineInternal();
    template <class T>
    CFEMDECOneWay *getReverseEngineInternal();
    template <class SRC, class TRG>
    std::unique_ptr<CFEMDECOneWay> getClassInternal(std::unique_ptr<CFEMDECAccessToMaster> &&accessToMaster);

   private:
    MCAuto<MEDCouplingUMesh> _local_mesh;
    MCAuto<DataArrayIdType> _global_node_ids;
    std::unique_ptr<CFEMDECOneWay> _engine;
    std::unique_ptr<CFEMDECOneWay> _reverse_engine;
};

}  // namespace MEDCoupling
