//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
#ifndef __MEDPARTITIONER_MESHCOLLECTION_HXX__
#define __MEDPARTITIONER_MESHCOLLECTION_HXX__

#include "MEDPARTITIONER.hxx"

#include "MEDPARTITIONER_Graph.hxx"
#include "MEDCouplingUMesh.hxx"

//#include "MEDPARTITIONER_FaceModel.hxx"
//#include "boost/shared_ptr.hpp"
#include <vector>
#include <map>
#include <string>
namespace ParaMEDMEM
{
  class MEDCouplingUMesh;
  class DataArrayInt;
}


namespace MEDPARTITIONER
{

  class Topology;
  class MESHCollectionDriver;
  class ParaDomainSelector;
  class MEDSKYLINEARRAY;
  class CONNECTZONE;
  class JointFinder;
  typedef enum{MedAscii, MedXML, Undefined} DriverType;

  typedef std::multimap<std::pair<int,int>, std::pair<int,int> > NodeMapping ;
  typedef std::vector<std::pair<int,int> >  NodeList;
  
  class MEDPARTITIONER_EXPORT MESHCollection
  {

  public:

    //Default constructor
    MESHCollection();

    //Constructing from an existing mesh and a new topology
    MESHCollection(MESHCollection&, Topology*, bool family_splitting=false, bool create_empty_groups=false);

    //Constructing the mesh collection from a file
    MESHCollection(const std::string& filename);

    //Constructing the mesh collection from a file
    MESHCollection(const std::string& filename, ParaDomainSelector& domainSelector);

    //Constructing the mesh collection from a file
    MESHCollection(const std::string& filename, const std::string& meshname);

    ~MESHCollection();
        
    bool isParallelMode() const { return _domain_selector; }
        
    //writing to a distributed file
    void write(const std::string& filename);

    //getting the driver
    MESHCollectionDriver* retrieveDriver();
    MESHCollectionDriver* getDriver() const;
    void setDriverType(MEDPARTITIONER::DriverType type) {_driver_type=type;}

    //creation of the cell graph
    void buildCellGraph(MEDPARTITIONER::MEDSKYLINEARRAY* & array,int *& edgeweights );

    //creation and partition of the associated graph
    Topology* createPartition(int nbdomain, Graph::splitter_type type = Graph::METIS,
                              const std::string& ="", int* edgeweights=0, int* verticesweights=0);

    //creation of a user specified partition
    Topology* createPartition(const int* partition);
    
     //getting mesh dimension
     int getMeshDimension() const ;

     //getting a reference to mesh vector
     std::vector<ParaMEDMEM::MEDCouplingUMesh*>& getMesh();
     std::vector<ParaMEDMEM::MEDCouplingUMesh*>& getFaceMesh();
     std::vector<std::vector<ParaMEDMEM::MEDCouplingUMesh*> >& getGroupMeshes();

     ParaMEDMEM::MEDCouplingUMesh* getMesh(int idomain) const;
     ParaMEDMEM::MEDCouplingUMesh* getFaceMesh(int idomain);
     std::vector<ParaMEDMEM::MEDCouplingUMesh*>& getGroupMeshes(int idomain);

     std::vector<ParaMEDMEM::DataArrayInt*>& getCellFamilyIds() {return _cellFamilyIds;}
     std::vector<ParaMEDMEM::DataArrayInt*>& getFaceFamilyIds() {return _faceFamilyIds;}
     
     std::map<std::string, ParaMEDMEM::DataArrayInt*>& getMapDataArrayInt() {return _mapDataArrayInt;}
     std::map<std::string, ParaMEDMEM::DataArrayDouble*>& getMapDataArrayDouble() {return _mapDataArrayDouble;}

     std::map<std::string,int>& getFamilyInfo() {return _familyInfo;}
     std::map<std::string, std::vector<std::string> >& getGroupInfo() {return _groupInfo;}

     ParaMEDMEM::DataArrayDouble* getField(std::string descriptionField, int iold);
     std::vector<std::string>&  getFieldDescriptions() {return _fieldDescriptions;}
     void prepareFieldDescriptions();
     void filterFaceOnCell();
     
    //getting a reference to connect zones vector
    std::vector<MEDPARTITIONER::CONNECTZONE*>& getCZ();

    //getting a pointer to topology
    Topology* getTopology() const ;
    ParaDomainSelector* getParaDomainSelector() const{return _domain_selector;}
    //settig a new topology
    void setTopology(Topology* topology);

    //getting/setting the name of the global mesh (as opposed 
    //to the name of a subdomain \a nn, which is name_nn) 
    std::string getName() const {return _name;}
    void setName(const std::string& name){_name=name;}
    void setDomainNames(const std::string& name);

    //getting/setting the description of the global mesh
    std::string getDescription() const {return _description;}
    void setDescription(const std::string& name) { _description=name;}

    //creates the node mapping between an old collection and the present one
    void createNodeMapping(MESHCollection& initialCollection, 
                           std::multimap<std::pair<int,int>,std::pair<int,int> >& nodeMapping);

    void castCellMeshes(MESHCollection& initialCollection, 
                        std::vector<std::vector<std::vector<int> > >& new2oldIds);
    
    //creates faces on the new collection
    void castFaceMeshes(MESHCollection& initialCollection,
                    const std::multimap<std::pair<int,int>, std::pair<int,int> >& nodeMapping,
                    std::vector<std::vector<std::vector<int> > >& new2oldIds);

  private:

    void castIntField2(std::vector<ParaMEDMEM::MEDCouplingUMesh*>& meshesCastFrom,
                       std::vector<ParaMEDMEM::MEDCouplingUMesh*>& meshesCastTo,
                       std::vector<ParaMEDMEM::DataArrayInt*>& arrayFrom,
                       std::string nameArrayTo);

    void castAllFields(MESHCollection& initialCollection,
                       std::string nameArrayTo);

    void findCommonDistantNodes(std::vector<std::vector<std::multimap<int,int> > >& commonDistantNodes);

    void remapIntField(const ParaMEDMEM::MEDCouplingUMesh& sourceMesh,
                       const ParaMEDMEM::MEDCouplingUMesh& targetMesh,
                       const int* fromArray,
                       int* toArray);
    
    void remapIntField2(int inew, int iold, 
                       const ParaMEDMEM::MEDCouplingUMesh& sourceMesh,
                       const ParaMEDMEM::MEDCouplingUMesh& targetMesh,
                       const int* fromArray,
                       std::string nameArrayTo);

    void remapDoubleField3(int inew, int iold,
                           ParaMEDMEM::DataArrayDouble* fromArray,
                           std::string nameArrayTo,
                           std::string descriptionField);

    //!link to mesh_collection topology
    Topology* _topology;

    //!control over topology
    bool _owns_topology;

    //!link to graph
    //Graph* _cell_graph;

    //! Driver for read/write operations
    MESHCollectionDriver* _driver;

    //! Parallelizer - mark of parallel execution mode
    ParaDomainSelector* _domain_selector;

    //!links to meshes
    std::vector<ParaMEDMEM::MEDCouplingUMesh*> _mesh;
    std::vector<ParaMEDMEM::MEDCouplingUMesh*> _faceMesh;
   
    //!index of a non empty mesh within _mesh (in parallel mode all of meshes can be empty)
    int _i_non_empty_mesh;

    //!links to connectzones
    std::vector<MEDPARTITIONER::CONNECTZONE*> _connect_zones;

    //!family ids storages
    std::vector<ParaMEDMEM::DataArrayInt*> _cellFamilyIds;
    std::vector<ParaMEDMEM::DataArrayInt*> _faceFamilyIds;
    
    //!DataArrayInt* storages
    std::map<std::string, ParaMEDMEM::DataArrayInt*> _mapDataArrayInt;
    //!DataArrayDouble* storages
    std::map<std::string, ParaMEDMEM::DataArrayDouble*> _mapDataArrayDouble;
    
    std::vector<std::string> _fieldDescriptions; //fields to be partitioned
    
    //!group family conversion
    std::map<std::string, int> _familyInfo;
    std::map<std::string, std::vector<std::string> > _groupInfo;
  
    //!list of groups that are not to be splitted
    std::vector<std::string> _indivisible_regions;

    //!name of global mesh
    std::string _name;

    //!description of global mesh
    std::string _description;

    //! specifies the driver associated to the collection
    DriverType _driver_type;

    /*! flag specifying that the splitter should create boundary constituent entity
      so that they are written in joints*/
    bool _subdomain_boundary_creates;

        /*! flag specifying that families must be preserved by the
                splitting*/
    bool _family_splitting;

        /*! flag specifying that groups must be created on all domains, 
                even if they are empty*/
    bool _create_empty_groups;

    JointFinder* _joint_finder;
  };

}//of namespace

#endif /*MESHCOLLECTION_*/
