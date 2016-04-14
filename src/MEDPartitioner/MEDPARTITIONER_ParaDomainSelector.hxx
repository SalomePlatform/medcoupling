// Copyright (C) 2007-2016  CEA/DEN, EDF R&D
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

#ifndef __MEDPARTITIONER_PARADOMAINSELECTOR_HXX__
#define __MEDPARTITIONER_PARADOMAINSELECTOR_HXX__

#include "MEDPARTITIONER.hxx"

#include <memory>
#include <vector>

namespace MEDCoupling
{
  class MEDCouplingUMesh;
}

namespace MEDPARTITIONER
{
  class Graph;
  class JointExchangeData;

  /*!
   * \brief Communication helper in parallel mode
   */
  class MEDPARTITIONER_EXPORT ParaDomainSelector
  {
  public:
    ParaDomainSelector(bool mesure_memory=false);
    ~ParaDomainSelector();

    //processor rank
    int rank() const { return _rank; }
    //number of processors
    int nbProcs() const { return _world_size; }
    //true if is running on different hosts
    bool isOnDifferentHosts() const;
    //true if the domain with domainIndex is to be loaded on this proc
    bool isMyDomain(int domainIndex) const;
    //processor id where the domain with domainIndex resides
    int getProcessorID(int domainIndex) const;
    //Set nb of required domains. (Used to sort joints via jointId())
    void setNbDomains(int nb) { _nb_result_domains = nb; }
    //identifier for a joint
    int jointId( int local_domain, int distant_domain ) const;
  
    int getNbTotalCells() { return _cell_shift_by_domain.back(); }
    int getNbTotalNodes() { return _node_shift_by_domain.back(); };
    int getNbTotalFaces() { return _face_shift_by_domain.back(); };

    //Collect nb of entities on procs
    void gatherNbOf(const std::vector<MEDCoupling::MEDCouplingUMesh*>& domain_meshes);
  
    //distribution of the graph vertices among the processors
    int* getProcVtxdist() const;

    //nb of nodes on processors with lower rank
    int getProcNodeShift() const;
    //nb of cells in domains with lower index
    int getDomainCellShift(int domainIndex) const;
    //nb of nodes in domains with lower index
    int getDomainNodeShift(int domainIndex) const;

    //Gather graphs from all processors into one
    std::auto_ptr<Graph> gatherGraph(const Graph* graph) const;

    //Set nb of cell/cell pairs in a joint between domains
    void setNbCellPairs( int nb_cell_pairs, int dist_domain, int loc_domain );
    //Gather size of each proc/proc joint
    void gatherNbCellPairs();
    //nb of cell/cell pairs in a joint between domains on different procs
    int getNbCellPairs( int dist_domain, int loc_domain ) const;

    //get the first global id of sub-entity for the joint
    int getFisrtGlobalIdOfSubentity( int loc_domain, int dist_domain ) const;
    //Send-receive local ids of joint faces
    int* exchangeSubentityIds( int loc_domain, int dist_domain,
                               const std::vector<int>& loc_ids_here ) const;
    //time passed from construction in seconds
    double getPassedTime() const;

    //Evaluate current memory usage and return the maximal one in KB
    int evaluateMemory() const;

    void sendMesh(const MEDCoupling::MEDCouplingUMesh& mesh, int target) const;
    void recvMesh(MEDCoupling::MEDCouplingUMesh*& mesh, int source) const;
  private:
    int _rank; //my rank
    int _world_size; //nb of processors
    int _nb_result_domains; //required nb of domains

    std::vector< int > _nb_cell_pairs_by_joint;
    std::vector< int > _nb_vert_of_procs; //graph vertices
    std::vector< int > _cell_shift_by_domain;
    std::vector< int > _node_shift_by_domain;
    std::vector< int > _face_shift_by_domain;

    double _init_time;
    bool _mesure_memory;
    mutable int _init_memory;
    mutable int _max_memory;
  };
}
#endif
