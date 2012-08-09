// Copyright (C) 2007-2012  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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

#ifndef __MEDPARTITIONER_MESHCOLLECTION_HXX__
#define __MEDPARTITIONER_MESHCOLLECTION_HXX__

#include "MEDPARTITIONER.hxx"
#include "MEDPARTITIONER_Graph.hxx"

#include "MEDCouplingUMesh.hxx"

#include <map>
#include <vector>
#include <string>

#include "BBTree.txx"

namespace ParaMEDMEM
{
  class MEDCouplingUMesh;
  class DataArrayInt;
}

namespace MEDPARTITIONER
{
  class Topology;
  class MeshCollectionDriver;
  class ParaDomainSelector;
  class SkyLineArray;
  class ConnectZone;
  class JointFinder;
  
  typedef enum{MedAscii, MedXml, Undefined} DriverType;
  typedef std::multimap<std::pair<int,int>, std::pair<int,int> > NodeMapping ;
  typedef std::vector<std::pair<int,int> >  NodeList;
  
  class MEDPARTITIONER_EXPORT MeshCollection
  {
  public:
    MeshCollection();
    //Constructing from an existing mesh and a new topology
    MeshCollection(MeshCollection&, Topology*, bool family_splitting=false, bool create_empty_groups=false);
    //Constructing the mesh collection from a file
    MeshCollection(const std::string& filename);
    //Constructing the mesh collection from a file
    MeshCollection(const std::string& filename, ParaDomainSelector& domainSelector);
    //Constructing the mesh collection from a file
    MeshCollection(const std::string& filename, const std::string& meshname);
    ~MeshCollection();
    bool isParallelMode() const { return _domain_selector; }
    
    //writing to a distributed file
    void write(const std::string& filename);
    
    //getting the driver
    MeshCollectionDriver *retrieveDriver();
    MeshCollectionDriver *getDriver() const;
    void setDriverType(MEDPARTITIONER::DriverType type) { _driver_type=type; }

    //creation of the cell graph
    void buildCellGraph(MEDPARTITIONER::SkyLineArray* & array,int *& edgeweights );
    //creation and partition of the associated graph
    Topology* createPartition(int nbdomain, Graph::splitter_type type = Graph::METIS,
                              const std::string& ="", int* edgeweights=0, int* verticesweights=0);

    //creation of a user specified partition
    Topology* createPartition(const int* partition);

    //getting mesh dimension
    int getMeshDimension() const;
    int getNbOfLocalMeshes() const;
    int getNbOfGlobalMeshes() const { return _mesh.size(); }
    int getNbOfLocalCells() const;
    int getNbOfLocalFaces() const;
    
    //getting a reference to mesh vector
    std::vector<ParaMEDMEM::MEDCouplingUMesh*>& getMesh();
    std::vector<ParaMEDMEM::MEDCouplingUMesh*>& getFaceMesh();
    std::vector<std::vector<ParaMEDMEM::MEDCouplingUMesh*> >& getGroupMeshes();

    ParaMEDMEM::MEDCouplingUMesh* getMesh(int idomain) const;
    ParaMEDMEM::MEDCouplingUMesh* getFaceMesh(int idomain);
    std::vector<ParaMEDMEM::MEDCouplingUMesh*>& getGroupMeshes(int idomain);

    std::vector<ParaMEDMEM::DataArrayInt*>& getCellFamilyIds() { return _cell_family_ids; }
    std::vector<ParaMEDMEM::DataArrayInt*>& getFaceFamilyIds() { return _face_family_ids; }
    
    std::map<std::string, ParaMEDMEM::DataArrayInt*>& getMapDataArrayInt() { return _map_dataarray_int; }
    std::map<std::string, ParaMEDMEM::DataArrayDouble*>& getMapDataArrayDouble() { return _map_dataarray_double; }

    std::map<std::string,int>& getFamilyInfo() { return _family_info; }
    std::map<std::string, std::vector<std::string> >& getGroupInfo() { return _group_info; }

    ParaMEDMEM::DataArrayDouble* getField(std::string descriptionField, int iold);
    std::vector<std::string>&  getFieldDescriptions() { return _field_descriptions; }
    void prepareFieldDescriptions();
    void filterFaceOnCell();
     
    //getting a reference to connect zones vector
    std::vector<MEDPARTITIONER::ConnectZone*>& getCZ();

    //getting a pointer to topology
    Topology* getTopology() const ;
    ParaDomainSelector* getParaDomainSelector() const { return _domain_selector; }
    //setting a new topology
    void setTopology(Topology* topology);

    //getting/setting the name of the global mesh (as opposed 
    //to the name of a subdomain \a nn, which is name_nn) 
    std::string getName() const { return _name; }
    void setName(const std::string& name) { _name=name; }
    void setDomainNames(const std::string& name);

    //getting/setting the description of the global mesh
    std::string getDescription() const { return _description; }
    void setDescription(const std::string& name) { _description=name; }

    //creates the node mapping between an old collection and the present one
    void createNodeMapping(MeshCollection& initialCollection, 
                           std::multimap<std::pair<int,int>,std::pair<int,int> >& nodeMapping);
    
    void castCellMeshes(MeshCollection& initialCollection, 
                        std::vector<std::vector<std::vector<int> > >& new2oldIds);
    
    //creates faces on the new collection
    void castFaceMeshes(MeshCollection& initialCollection,
                        const std::multimap<std::pair<int,int>, std::pair<int,int> >& nodeMapping,
                        std::vector<std::vector<std::vector<int> > >& new2oldIds);

  private:
    void castIntField(std::vector<ParaMEDMEM::MEDCouplingUMesh*>& meshesCastFrom,
                       std::vector<ParaMEDMEM::MEDCouplingUMesh*>& meshesCastTo,
                       std::vector<ParaMEDMEM::DataArrayInt*>& arrayFrom,
                       std::string nameArrayTo);

    void castAllFields(MeshCollection& initialCollection,
                       std::string nameArrayTo);

    void findCommonDistantNodes(std::vector<std::vector<std::multimap<int,int> > >& commonDistantNodes);

    
    void remapIntField(int inew, int iold, 
                        const ParaMEDMEM::MEDCouplingUMesh& sourceMesh,
                        const ParaMEDMEM::MEDCouplingUMesh& targetMesh,
                        const int* fromArray,
                        std::string nameArrayTo,
                       const BBTree<3,int>* tree);

    void remapDoubleField(int inew, int iold,
                           ParaMEDMEM::DataArrayDouble* fromArray,
                           std::string nameArrayTo,
                           std::string descriptionField);
  private:

    //link to mesh_collection topology
    Topology* _topology;
    
    //control over topology
    bool _owns_topology;
    
    //Driver for read/write operations
    MeshCollectionDriver* _driver;
    
    //Parallelizer - mark of parallel execution mode
    ParaDomainSelector* _domain_selector;
    
    //links to meshes
    std::vector<ParaMEDMEM::MEDCouplingUMesh*> _mesh;
    std::vector<ParaMEDMEM::MEDCouplingUMesh*> _face_mesh;
    
    //index of a non empty mesh within _mesh (in parallel mode all of meshes can be empty)
    int _i_non_empty_mesh;
    
    //links to connectzones
    std::vector<MEDPARTITIONER::ConnectZone*> _connect_zones;

    //family ids storages
    std::vector<ParaMEDMEM::DataArrayInt*> _cell_family_ids;
    std::vector<ParaMEDMEM::DataArrayInt*> _face_family_ids;
    
    //DataArrayInt* storages
    std::map<std::string, ParaMEDMEM::DataArrayInt*> _map_dataarray_int;
    //DataArrayDouble* storages
    std::map<std::string, ParaMEDMEM::DataArrayDouble*> _map_dataarray_double;
    
    //fields to be partitioned
    std::vector<std::string> _field_descriptions;
    
    //group family conversion
    std::map<std::string, int> _family_info;
    std::map<std::string, std::vector<std::string> > _group_info;
  
    //list of groups that are not to be splitted
    std::vector<std::string> _indivisible_regions;

    //name of global mesh
    std::string _name;

    //description of global mesh
    std::string _description;

    //specifies the driver associated to the collection
    DriverType _driver_type;

    //flag specifying that the splitter should create boundary constituent entity
    //so that they are written in joints
    bool _subdomain_boundary_creates;

    //flag specifying that families must be preserved by the splitting
    bool _family_splitting;

    //flag specifying that groups must be created on all domains, even if they are empty
    bool _create_empty_groups;

    JointFinder* _joint_finder;
  };
}
#endif
