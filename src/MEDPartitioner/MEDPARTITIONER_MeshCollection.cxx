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

#include "MEDPARTITIONER_MeshCollection.hxx"
#include "MEDPARTITIONER_MeshCollectionDriver.hxx"
#include "MEDPARTITIONER_MeshCollectionMedXmlDriver.hxx"
#include "MEDPARTITIONER_MeshCollectionMedAsciiDriver.hxx"
#include "MEDPARTITIONER_ParaDomainSelector.hxx"
#include "MEDPARTITIONER_Topology.hxx"
#include "MEDPARTITIONER_ParallelTopology.hxx"

#ifdef HAVE_MPI2
#include "MEDPARTITIONER_JointFinder.hxx"
#endif

#include "MEDPARTITIONER_Graph.hxx"
#include "MEDPARTITIONER_UserGraph.hxx"
#include "MEDPARTITIONER_Utils.hxx" 

#include "MEDLoaderBase.hxx"
#include "MEDLoader.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingNormalizedUnstructuredMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "PointLocator3DIntersectorP0P0.hxx"

#include "MEDCouplingAutoRefCountObjectPtr.hxx"
#include "BBTree.txx"

#ifdef HAVE_MPI2
#include <mpi.h>
#endif

#if defined(MED_ENABLE_PARMETIS) || defined(MED_ENABLE_METIS)
#include "MEDPARTITIONER_MetisGraph.hxx"
#endif

#ifdef MED_ENABLE_SCOTCH
#include "MEDPARTITIONER_ScotchGraph.hxx"
#endif

#include <set>
#include <vector>
#include <string>
#include <limits>
#include <iostream>
#include <fstream>

MEDPARTITIONER::MeshCollection::MeshCollection()
  : _topology(0),
    _owns_topology(false),
    _driver(0),
    _domain_selector( 0 ),
    _i_non_empty_mesh(-1),
    _driver_type(MEDPARTITIONER::MedXml),
    _subdomain_boundary_creates(false),
    _family_splitting(false),
    _create_empty_groups(false),
    _joint_finder(0)
{
}

/*!constructor creating a new mesh collection (mesh series + topology)
 *from an old collection and a new topology
 * 
 * On output, the constructor has built the meshes corresponding to the new mesh collection.
 * The new topology has been updated so that face and node mappings are included.
 * The families have been cast to their projections in the new topology.
 * 
 * \param initial_collection collection from which the data (coordinates, connectivity) are taken
 * \param topology topology containing the cell mappings
 */

MEDPARTITIONER::MeshCollection::MeshCollection(MeshCollection& initialCollection, 
                                               Topology* topology, 
                                               bool family_splitting, 
                                               bool create_empty_groups)
  : _topology(topology),
    _owns_topology(false),
    _driver(0),
    _domain_selector( initialCollection._domain_selector ),
    _i_non_empty_mesh(-1),
    _name(initialCollection._name),
    _driver_type(MEDPARTITIONER::MedXml),
    _subdomain_boundary_creates(false),
    _family_splitting(family_splitting),
    _create_empty_groups(create_empty_groups),
    _joint_finder(0)
{
  std::vector<std::vector<std::vector<int> > > new2oldIds(initialCollection.getTopology()->nbDomain());
  castCellMeshes(initialCollection, new2oldIds);

  //defining the name for the collection and the underlying meshes
  setName(initialCollection.getName());

  /////////////////
  //treating faces
  /////////////////

#ifdef HAVE_MPI2
  if (MyGlobals::_Verbose>0 && MyGlobals::_World_Size>1)
    MPI_Barrier(MPI_COMM_WORLD); //synchronize verbose messages
#endif
  if (MyGlobals::_Is0verbose)
    std::cout<<"treating faces"<<std::endl;
  NodeMapping nodeMapping;
  //nodeMapping contains the mapping between old nodes and new nodes
  // (iolddomain,ioldnode)->(inewdomain,inewnode)
  createNodeMapping(initialCollection, nodeMapping);
  std::vector<std::vector<std::vector<int> > > new2oldFaceIds;
  castFaceMeshes(initialCollection, nodeMapping, new2oldFaceIds);

  ////////////////////
  //treating families
  ////////////////////
#ifdef HAVE_MPI2
  if (MyGlobals::_Verbose>0 && MyGlobals::_World_Size>1)
    MPI_Barrier(MPI_COMM_WORLD); //synchronize verbose messages
#endif
  if (MyGlobals::_Is0verbose)
    {
      if (isParallelMode())
        std::cout << "ParallelMode on " << topology->nbDomain() << " Domains" << std::endl;
      else
        std::cout << "NOT ParallelMode on " << topology->nbDomain() << " Domains" << std::endl;
    }
  if (MyGlobals::_Is0verbose>10)
    std::cout<<"treating cell and face families"<<std::endl;
  
  castIntField(initialCollection.getMesh(),
                this->getMesh(),
                initialCollection.getCellFamilyIds(),
                "cellFamily");
  castIntField(initialCollection.getFaceMesh(), 
                this->getFaceMesh(),
                initialCollection.getFaceFamilyIds(),
                "faceFamily");

  //treating groups
#ifdef HAVE_MPI2
  if (MyGlobals::_Verbose>0 && MyGlobals::_World_Size>1)
    MPI_Barrier(MPI_COMM_WORLD); //synchronize verbose messages
#endif
  if (MyGlobals::_Is0verbose)
    std::cout << "treating groups" << std::endl;
  _family_info=initialCollection.getFamilyInfo();
  _group_info=initialCollection.getGroupInfo();
  
#ifdef HAVE_MPI2
  if (MyGlobals::_Verbose>0 && MyGlobals::_World_Size>1)
    MPI_Barrier(MPI_COMM_WORLD); //synchronize verbose messages
#endif
  if (MyGlobals::_Is0verbose)
    std::cout << "treating fields" << std::endl;
  castAllFields(initialCollection,"cellFieldDouble");
  if (_i_non_empty_mesh<0)
    {
      for (int i=0; i<_mesh.size(); i++)
        {
          if (_mesh[i])
            {
              _i_non_empty_mesh=i; //first existing one local
              break;
            }
        }
    }

}

/*!
  Creates the meshes using the topology underlying he mesh collection and the mesh data 
  coming from the ancient collection
  \param initialCollection collection from which the data is extracted to create the new meshes
*/

void MEDPARTITIONER::MeshCollection::castCellMeshes(MeshCollection& initialCollection,
                                                    std::vector<std::vector<std::vector<int> > >& new2oldIds)
{
  if (MyGlobals::_Verbose>10)
    std::cout << "proc " << MyGlobals::_Rank << " : castCellMeshes" << std::endl;
  if (_topology==0)
    throw INTERP_KERNEL::Exception("Topology has not been defined on call to castCellMeshes");
  
  int nbNewDomain=_topology->nbDomain();
  int nbOldDomain=initialCollection.getTopology()->nbDomain();
  
  _mesh.resize(nbNewDomain);
  int rank=MyGlobals::_Rank;
  //splitting the initial domains into smaller bits
  std::vector<std::vector<ParaMEDMEM::MEDCouplingUMesh*> > splitMeshes;
  splitMeshes.resize(nbNewDomain);
  for (int inew=0; inew<nbNewDomain; inew++)
    {
      splitMeshes[inew].resize(nbOldDomain, (ParaMEDMEM::MEDCouplingUMesh*)0);
    }

  for (int iold=0; iold<nbOldDomain; iold++)
    {
      if (!isParallelMode() || initialCollection._domain_selector->isMyDomain(iold))
        {
          int size=(initialCollection._mesh)[iold]->getNumberOfCells();
          std::vector<int> globalids(size);
          initialCollection.getTopology()->getCellList(iold, &globalids[0]);
          std::vector<int> ilocalnew(size); //local
          std::vector<int> ipnew(size);     //idomain old
          _topology->convertGlobalCellList(&globalids[0],size,&ilocalnew[0],&ipnew[0]);
      
          new2oldIds[iold].resize(nbNewDomain);
          for (int i=0; i<(int)ilocalnew.size(); i++)
            {
              new2oldIds[iold][ipnew[i]].push_back(i);
            }
          for (int inew=0; inew<nbNewDomain; inew++)
            {
              splitMeshes[inew][iold]=(ParaMEDMEM::MEDCouplingUMesh*)
                (initialCollection.getMesh())[iold]->buildPartOfMySelf(&new2oldIds[iold][inew][0],
                                                                       &new2oldIds[iold][inew][0]+new2oldIds[iold][inew].size(),
                                                                       true);
              if (MyGlobals::_Verbose>400)
                std::cout<< "proc " << rank << " : a splitMesh iold inew NbCells " << iold << " " << inew << " "
                         << splitMeshes[inew][iold]->getNumberOfCells() << std::endl;
            }
        }
    }
#ifdef HAVE_MPI2
  if (isParallelMode())
    {
      //if (MyGlobals::_Verbose>300) std::cout<<"proc "<<rank<<" : castCellMeshes send/receive"<<std::endl;
      for (int iold=0; iold<nbOldDomain; iold++)
        for(int inew=0; inew<nbNewDomain; inew++)
          {
            if (initialCollection._domain_selector->isMyDomain(iold) && _domain_selector->isMyDomain(inew))
              continue;

            if(initialCollection._domain_selector->isMyDomain(iold))
              _domain_selector->sendMesh(*(splitMeshes[inew][iold]),_domain_selector->getProcessorID(inew));

            if (_domain_selector->isMyDomain(inew))
              _domain_selector->recvMesh(splitMeshes[inew][iold],_domain_selector->getProcessorID(iold));

          }
    }
#endif

  //fusing the split meshes
  if (MyGlobals::_Verbose>200)
    std::cout << "proc " << rank << " : castCellMeshes fusing" << std::endl;
  for (int inew=0; inew<nbNewDomain ;inew++)
    {
      std::vector<const ParaMEDMEM::MEDCouplingUMesh*> meshes;
    
      for (int i=0; i<(int)splitMeshes[inew].size(); i++)
        if (splitMeshes[inew][i]!=0) 
          if (splitMeshes[inew][i]->getNumberOfCells()>0)
            meshes.push_back(splitMeshes[inew][i]);

      if (!isParallelMode()||_domain_selector->isMyDomain(inew))
        {
          if (meshes.size()==0) 
            {
              _mesh[inew]=CreateEmptyMEDCouplingUMesh();
              std::cout << "WARNING : castCellMeshes fusing : no meshes try another number of processors" << std::endl;
            }
          else
            {
              _mesh[inew]=ParaMEDMEM::MEDCouplingUMesh::MergeUMeshes(meshes);
              bool areNodesMerged;
              int nbNodesMerged;
              if (meshes.size()>1)
                {
                  ParaMEDMEM::DataArrayInt* array=_mesh[inew]->mergeNodes(1e-12,areNodesMerged,nbNodesMerged);
                  array->decrRef(); // array is not used in this case
                }
              _mesh[inew]->zipCoords();
                
            }
        }
      for (int i=0;i<(int)splitMeshes[inew].size();i++)
        if (splitMeshes[inew][i]!=0)
          splitMeshes[inew][i]->decrRef();
    }  
  if (MyGlobals::_Verbose>300)
    std::cout << "proc " << rank << " : castCellMeshes end fusing" << std::endl;
}

/*!
  \param initialCollection source mesh collection 
  \param nodeMapping structure containing the correspondency between nodes in the initial collection and the node(s) in the new collection
*/
void MEDPARTITIONER::MeshCollection::createNodeMapping( MeshCollection& initialCollection, NodeMapping& nodeMapping)
{
  using std::vector;
  using std::make_pair;
  //  NodeMapping reverseNodeMapping;
  for (int iold=0; iold<initialCollection.getTopology()->nbDomain();iold++)
    {

      double* bbox;
      BBTree<3>* tree; 
      if (!isParallelMode() || (_domain_selector->isMyDomain(iold)))
        {
          //      std::map<pair<double,pair<double, double> >, int > nodeClassifier;
          int nvertices=initialCollection.getMesh(iold)->getNumberOfNodes();
          bbox=new double[nvertices*2*3];
          ParaMEDMEM::DataArrayDouble* coords = initialCollection.getMesh(iold)->getCoords();
          double* coordsPtr=coords->getPointer();

          for (int i=0; i<nvertices*3;i++)
            {
              bbox[i*2]=coordsPtr[i]-1e-8;
              bbox[i*2+1]=coordsPtr[i]+1e-8;
            }
          tree=new BBTree<3>(bbox,0,0,nvertices,1e-9);
        }

      for (int inew=0; inew<_topology->nbDomain(); inew++)
        {
#ifdef HAVE_MPI2
          //sending meshes for parallel computation
          if (isParallelMode() && _domain_selector->isMyDomain(inew) && !_domain_selector->isMyDomain(iold))  
            _domain_selector->sendMesh(*(getMesh(inew)), _domain_selector->getProcessorID(iold));
          else if (isParallelMode() && !_domain_selector->isMyDomain(inew)&& _domain_selector->isMyDomain(iold))
            {
              ParaMEDMEM::MEDCouplingUMesh* mesh;
              _domain_selector->recvMesh(mesh, _domain_selector->getProcessorID(inew));
              ParaMEDMEM::DataArrayDouble* coords = mesh->getCoords();
              for (int inode=0; inode<mesh->getNumberOfNodes();inode++)
                {
                  double* coordsPtr=coords->getPointer()+inode*3;
                  vector<int> elems;
                  tree->getElementsAroundPoint(coordsPtr,elems);
                  if (elems.size()==0) continue;
                  nodeMapping.insert(make_pair(make_pair(iold,elems[0]),make_pair(inew,inode)));
                }
              mesh->decrRef();
            }
          else if (!isParallelMode() || (_domain_selector->isMyDomain(inew) && _domain_selector->isMyDomain(iold)))
#else
          if (!isParallelMode() || (_domain_selector->isMyDomain(inew) && _domain_selector->isMyDomain(iold)))
#endif
            {
              ParaMEDMEM::DataArrayDouble* coords = getMesh(inew)->getCoords();
              for (int inode=0; inode<_mesh[inew]->getNumberOfNodes();inode++)
                {
                  double* coordsPtr=coords->getPointer()+inode*3;
                  vector<int> elems;
                  tree->getElementsAroundPoint(coordsPtr,elems);
                  if (elems.size()==0) continue;
                  nodeMapping.insert(make_pair(make_pair(iold,elems[0]),make_pair(inew,inode)));
                }
            }
        }
      if (!isParallelMode() || (_domain_selector->isMyDomain(iold)))
        {
          delete tree;
          delete[] bbox;
        }
    } 

}

void getNodeIds(ParaMEDMEM::MEDCouplingUMesh& meshOne, ParaMEDMEM::MEDCouplingUMesh& meshTwo, std::vector<int>& nodeIds)
{
  using std::vector;
  if (!&meshOne || !&meshTwo) return;  //empty or not existing
  double* bbox;
  BBTree<3>* tree;
  int nv1=meshOne.getNumberOfNodes();
  bbox=new double[nv1*6];
  ParaMEDMEM::DataArrayDouble* coords=meshOne.getCoords();
  double* coordsPtr=coords->getPointer();
  for (int i=0; i<nv1*3; i++)
    {
      bbox[i*2]=coordsPtr[i]-1e-8;
      bbox[i*2+1]=coordsPtr[i]+1e-8;
    }
  tree=new BBTree<3>(bbox,0,0,nv1,1e-9);
  
  int nv2=meshTwo.getNumberOfNodes();
  nodeIds.resize(nv2,-1);
  coords=meshTwo.getCoords();
  for (int inode=0; inode<nv2; inode++)
    {
      double* coordsPtr2=coords->getPointer()+inode*3;
      vector<int> elems;
      tree->getElementsAroundPoint(coordsPtr2,elems);
      if (elems.size()==0) continue;
      nodeIds[inode]=elems[0];
    }
  delete tree;
  delete[] bbox;
}

/*!
  creates the face meshes on the new domains from the faces on the old domain and the node mapping
  faces at the interface are duplicated
*/
void MEDPARTITIONER::MeshCollection::castFaceMeshes(MeshCollection& initialCollection,
                                                    const std::multimap<std::pair<int,int>, std::pair<int,int> >& nodeMapping,
                                                    std::vector<std::vector<std::vector<int> > >& new2oldIds)
{

  //splitMeshes structure will contain the partition of 
  //the old faces on the new ones
  //splitMeshes[4][2] contains the faces from old domain 2
  //that have to be added to domain 4
  
  using std::vector;
  using std::map;
  using std::multimap;
  using std::pair;
  using std::make_pair;
  
  if (MyGlobals::_Verbose>10)
    std::cout << "proc " << MyGlobals::_Rank << " : castFaceMeshes" << std::endl;
  if (_topology==0)
    throw INTERP_KERNEL::Exception("Topology has not been defined on call to castFaceMeshes");
  
  int nbNewDomain=_topology->nbDomain();
  int nbOldDomain=initialCollection.getTopology()->nbDomain();
  
  vector<ParaMEDMEM::MEDCouplingUMesh*>& meshesCastFrom=initialCollection.getFaceMesh();
  vector<ParaMEDMEM::MEDCouplingUMesh*>& meshesCastTo=this->getFaceMesh();
  
  vector< vector<ParaMEDMEM::MEDCouplingUMesh*> > splitMeshes;
  
  splitMeshes.resize(nbNewDomain);
  for (int inew=0; inew<nbNewDomain; inew++)
    {
      splitMeshes[inew].resize(nbOldDomain, (ParaMEDMEM::MEDCouplingUMesh*)0);
    }
  new2oldIds.resize(nbOldDomain);
  for (int iold=0; iold<nbOldDomain; iold++) new2oldIds[iold].resize(nbNewDomain);
  
  //init null pointer for empty meshes
  for (int inew=0; inew<nbNewDomain; inew++)
    {
      for (int iold=0; iold<nbOldDomain; iold++)
        {
          splitMeshes[inew][iold]=0; //null for empty meshes
          new2oldIds[iold][inew].clear();
        }
    }

  //loop over the old domains to analyse the faces and decide 
  //on which new domain they belong
  for (int iold=0; iold<nbOldDomain; iold++)
    {
      if (isParallelMode() && !initialCollection._domain_selector->isMyDomain(iold)) continue;
      if (MyGlobals::_Verbose>400)
        std::cout<<"proc "<<MyGlobals::_Rank<<" : castFaceMeshes iodDomain "<<iold<<std::endl;
      //initial face mesh known : in my domain
      if (meshesCastFrom[iold] != 0)
        {
          for (int ielem=0; ielem<meshesCastFrom[iold]->getNumberOfCells(); ielem++)
            {
              vector<int> nodes;
              meshesCastFrom[iold]->getNodeIdsOfCell(ielem,nodes);
          
              map <int,int> faces;

              //analysis of element ielem
              //counters are set for the element
              //for each source node, the mapping is interrogated and the domain counters 
              //are incremented for each target node
              //the face is considered as going to target domains if the counter of the domain 
              //is equal to the number of nodes
              for (int inode=0; inode<(int)nodes.size(); inode++)
                {
                  typedef multimap<pair<int,int>,pair<int,int> >::const_iterator MI;
                  int mynode=nodes[inode];

                  pair <MI,MI> myRange = nodeMapping.equal_range(make_pair(iold,mynode));
                  for (MI iter=myRange.first; iter!=myRange.second; iter++)
                    {
                      int inew=iter->second.first;
                      if (faces.find(inew)==faces.end())
                        faces[inew]=1;
                      else
                        faces[inew]++;
                    }
                }
          
              for (map<int,int>::iterator iter=faces.begin(); iter!=faces.end(); iter++)
                {
                  if (iter->second==(int)nodes.size())
                    //cvw eligible but may be have to be face of a cell of this->getMesh()[inew]?
                    //it is not sure here...
                    //done before writeMedfile on option?... see filterFaceOnCell()
                    new2oldIds[iold][iter->first].push_back(ielem);
                }
            }
      
          //creating the splitMeshes from the face ids
          for (int inew=0; inew<nbNewDomain; inew++)
            {
              if (meshesCastFrom[iold]->getNumberOfCells() > 0)
                {
                  splitMeshes[inew][iold]=
                    (ParaMEDMEM::MEDCouplingUMesh*) 
                    ( meshesCastFrom[iold]->buildPartOfMySelf(&new2oldIds[iold][inew][0],
                                                              &new2oldIds[iold][inew][0]+new2oldIds[iold][inew].size(),
                                                              true) 
                      );
                  splitMeshes[inew][iold]->zipCoords();
                }
              else
                {
                  splitMeshes[inew][iold]=CreateEmptyMEDCouplingUMesh();
                }
            }
        }
      else
        {
          std::cout<<"proc "<<MyGlobals::_Rank<<" : castFaceMeshes empty mesh from iodDomain "<<iold<<std::endl;
        }
    }
  
#ifdef HAVE_MPI2
  //send/receive stuff
  if (isParallelMode())
    {
      ParaMEDMEM::MEDCouplingUMesh *empty=CreateEmptyMEDCouplingUMesh();
      for (int iold=0; iold<nbOldDomain; iold++)
        for (int inew=0; inew<nbNewDomain; inew++)
          {
            if (initialCollection._domain_selector->isMyDomain(iold) && !_domain_selector->isMyDomain(inew))
              {
                if (splitMeshes[inew][iold] != 0)
                  {
                    _domain_selector->sendMesh(*(splitMeshes[inew][iold]), _domain_selector->getProcessorID(inew));
                    //std::cout << "proc " << MyGlobals::_Rank << " : castFaceMeshes send " <<inew<<" "<<iold<<" "<<splitMeshes[inew][iold]->getNumberOfCells()<<std::endl;
                  }
                else
                  {
                    _domain_selector->sendMesh(*(empty), _domain_selector->getProcessorID(inew));
                    //std::cout << "proc " << MyGlobals::_Rank << " : castFaceMeshes sen0 " <<inew<<" "<<iold<<std::endl;
                  }
              }
            if (!initialCollection._domain_selector->isMyDomain(iold) && _domain_selector->isMyDomain(inew))
              _domain_selector->recvMesh(splitMeshes[inew][iold], _domain_selector->getProcessorID(iold));
              int nb=0;
              if (splitMeshes[inew][iold])
                nb=splitMeshes[inew][iold]->getNumberOfCells();
              //std::cout << "proc " << MyGlobals::_Rank << " : castFaceMeshes recv "<<inew<<" "<<iold<<" "<<nb<<std::endl;//" "<<splitMeshes[inew][iold]->getNumberOfCells()<<std::endl;
          }
      empty->decrRef();
    }
#endif

  //fusing the split meshes
  if (MyGlobals::_Verbose>200)
    std::cout << "proc " << MyGlobals::_Rank << " : castFaceMeshes fusing" << std::endl;
  meshesCastTo.resize(nbNewDomain);
  for (int inew=0; inew<nbNewDomain; inew++)
    {
      vector<const ParaMEDMEM::MEDCouplingUMesh*> myMeshes;
      for (int iold=0; iold<nbOldDomain; iold++)
        {
          ParaMEDMEM::MEDCouplingUMesh *umesh=splitMeshes[inew][iold];
          if (umesh!=0)
            if (umesh->getNumberOfCells()>0)
                myMeshes.push_back(umesh);
        }
      
      if (myMeshes.size()>0)
        {
          meshesCastTo[inew]=ParaMEDMEM::MEDCouplingUMesh::MergeUMeshes(myMeshes);
        }
      else
        {
          ParaMEDMEM::MEDCouplingUMesh *empty=CreateEmptyMEDCouplingUMesh();
          meshesCastTo[inew]=empty;
        }
      for (int iold=0; iold<nbOldDomain; iold++)
        if (splitMeshes[inew][iold]!=0)
          splitMeshes[inew][iold]->decrRef();
    }
  if (MyGlobals::_Verbose>300)
    std::cout << "proc " << MyGlobals::_Rank << " : castFaceMeshes end fusing" << std::endl;
}



void MEDPARTITIONER::MeshCollection::castIntField(std::vector<ParaMEDMEM::MEDCouplingUMesh*>& meshesCastFrom,
                                   std::vector<ParaMEDMEM::MEDCouplingUMesh*>& meshesCastTo,
                                   std::vector<ParaMEDMEM::DataArrayInt*>& arrayFrom,
                                   std::string nameArrayTo)
{
  using std::vector;
  
  int ioldMax=meshesCastFrom.size();
  int inewMax=meshesCastTo.size();


  //preparing bounding box trees for accelerating source-target node identifications
  if (MyGlobals::_Verbose>99)
    std::cout<<"making accelerating structures"<<std::endl;
  std::vector<BBTree<3,int>* > acceleratingStructures(ioldMax);
  std::vector<ParaMEDMEM::DataArrayDouble*>bbox(ioldMax);
  for (int iold =0; iold< ioldMax; iold++)
    if (isParallelMode() && _domain_selector->isMyDomain(iold))
      {
        ParaMEDMEM::DataArrayDouble* sourceCoords=meshesCastFrom[iold]->getBarycenterAndOwner();
        bbox[iold]=sourceCoords->computeBBoxPerTuple(1.e-6);
        acceleratingStructures[iold]=new BBTree<3,int> (bbox[iold]->getConstPointer(),0,0,bbox[iold]->getNumberOfTuples());
      }

  // send-recv operations
#ifdef HAVE_MPI2
  for (int inew=0; inew<inewMax; inew++)
    {
      for (int iold=0; iold<ioldMax; iold++)
        {
          //sending arrays for distant domains
          if (isParallelMode() && _domain_selector->isMyDomain(iold) && !_domain_selector->isMyDomain(inew))
            {
              //send mesh
              _domain_selector->sendMesh(*meshesCastFrom[iold],_domain_selector->getProcessorID(inew));
              //send vector
              int size=arrayFrom[iold]->getNumberOfTuples(); //cvw may be -1!
              vector<int>sendIds;
              if (MyGlobals::_Verbose>400) std::cout<<"proc "<<_domain_selector->rank()<<" : castIntField SendIntVec size "<<size<<std::endl;
              if (size>0) //no empty
                {
                  sendIds.resize(size);
                  std::copy(arrayFrom[iold]->getPointer(),arrayFrom[iold]->getPointer()+size,&sendIds[0]);
                }
              else //empty
                {
                  size=0;
                  sendIds.resize(size);
                }
              SendIntVec(sendIds, _domain_selector->getProcessorID(inew));
            }
          //receiving arrays from distant domains
          if (isParallelMode() && !_domain_selector->isMyDomain(iold) && _domain_selector->isMyDomain(inew))
            {
              //receive mesh
              vector<int> recvIds;
              ParaMEDMEM::MEDCouplingUMesh* recvMesh;
              _domain_selector->recvMesh(recvMesh,_domain_selector->getProcessorID(iold));
              //receive vector
              if (MyGlobals::_Verbose>400) std::cout<<"proc "<<_domain_selector->rank()<<" : castIntField recIntVec "<<std::endl;
              RecvIntVec(recvIds, _domain_selector->getProcessorID(iold));
              remapIntField(inew,iold,*recvMesh,*meshesCastTo[inew],&recvIds[0],nameArrayTo,0);
              recvMesh->decrRef(); //cww is it correct?
            }
        }
    }
#endif
  
  //local contributions and aggregation
  for (int inew=0; inew<inewMax; inew++)    
    {
      for (int iold=0; iold<ioldMax; iold++)
        if (!isParallelMode() || ( _domain_selector->isMyDomain(iold) && _domain_selector->isMyDomain(inew)))
          {
            remapIntField(inew,iold,*meshesCastFrom[iold],*meshesCastTo[inew],arrayFrom[iold]->getConstPointer(),nameArrayTo,acceleratingStructures[iold]);
          }
    }
  for (int iold=0; iold<ioldMax;iold++)
    if (isParallelMode() && _domain_selector->isMyDomain(iold)) 
      {
        bbox[iold]->decrRef();
        delete acceleratingStructures[iold];
      }
}

void MEDPARTITIONER::MeshCollection::remapIntField(int inew, int iold,
                                                    const ParaMEDMEM::MEDCouplingUMesh& sourceMesh,
                                                    const ParaMEDMEM::MEDCouplingUMesh& targetMesh,
                                                    const int* fromArray,
                                                    std::string nameArrayTo,
                                                    const BBTree<3,int>* myTree)
{

  if (sourceMesh.getNumberOfCells()<=0) return; //empty mesh could exist
  ParaMEDMEM::DataArrayDouble* targetCoords=targetMesh.getBarycenterAndOwner();
  const double*  tc=targetCoords->getConstPointer();
  int targetSize=targetMesh.getNumberOfCells();
  int sourceSize=sourceMesh.getNumberOfCells();
if (MyGlobals::_Verbose>200)
std::cout<<"remap vers target de taille "<<targetSize<<std::endl;
  std::vector<int> ccI;
  std::string str,cle;
  str=nameArrayTo+"_toArray";
  cle=Cle1ToStr(str,inew);
  int* toArray;
  
  const BBTree<3>* tree;
  bool cleantree=false;
  ParaMEDMEM::DataArrayDouble* sourceBBox=0;
  if (myTree==0)
    {
      sourceBBox=sourceMesh.getBarycenterAndOwner()->computeBBoxPerTuple(1e-8);
      tree=new BBTree<3> (sourceBBox->getConstPointer(),0,0, sourceBBox->getNumberOfTuples(),1e-10);
      cleantree=true;
    }
  else tree=myTree;
  //first time iold : create and initiate 
  if (_map_dataarray_int.find(cle)==_map_dataarray_int.end())
    {
      if (MyGlobals::_Is0verbose>100)
        std::cout << "create " << cle << " size " << targetSize << std::endl;
      ParaMEDMEM::DataArrayInt* p=ParaMEDMEM::DataArrayInt::New();
      p->alloc(targetSize,1);
      p->fillWithZero();
      toArray=p->getPointer();
      _map_dataarray_int[cle]=p;
    }
  else //other times iold: refind and complete
    {
      toArray=_map_dataarray_int.find(cle)->second->getPointer();
    }
  
  for (int itargetnode=0; itargetnode<targetSize; itargetnode++)    
    {
      std::vector<int> intersectingElems;
      tree->getElementsAroundPoint(tc+itargetnode*3,intersectingElems); // to be changed in 2D
      if (intersectingElems.size()!=0)
        {
          int isourcenode=intersectingElems[0];
          toArray[itargetnode]=fromArray[isourcenode];
          ccI.push_back(itargetnode);
          ccI.push_back(isourcenode);
        }
    }
  if (MyGlobals::_Verbose>200)
    std::cout<<"nb points trouves"<<ccI.size()/2<<std::endl;
  //memories intersection for future same job on fields (if no existing cle=no intersection)
  str=Cle2ToStr(nameArrayTo+"_ccI",inew,iold);
  if (MyGlobals::_Verbose>700)
    std::cout << "proc " << MyGlobals::_Rank << " : map memorize '" << str << "'\n";

  _map_dataarray_int[str]=CreateDataArrayIntFromVector(ccI, 2);
  
  targetCoords->decrRef();
  if (cleantree) delete tree;
  if (sourceBBox !=0) sourceBBox->decrRef();
  
 }

void MEDPARTITIONER::MeshCollection::castAllFields(MeshCollection& initialCollection, std::string nameArrayTo)
{
  if (nameArrayTo!="cellFieldDouble") 
    throw INTERP_KERNEL::Exception("Error castAllField only on cellFieldDouble");

  std::string nameTo="typeData=6"; //resume the type of field casted 
  // send-recv operations
  int ioldMax=initialCollection.getMesh().size();
  int inewMax=this->getMesh().size();
  int iFieldMax=initialCollection.getFieldDescriptions().size();
  if (MyGlobals::_Verbose>10)
    std::cout << "castAllFields with:\n" << ReprVectorOfString(initialCollection.getFieldDescriptions()) << std::endl;
  //see collection.prepareFieldDescriptions()
  for (int ifield=0; ifield<iFieldMax; ifield++)
    {
      std::string descriptionField=initialCollection.getFieldDescriptions()[ifield];
      if (descriptionField.find(nameTo)==std::string::npos)
        continue; //only nameTo accepted in Fields name description
#ifdef HAVE_MPI2
      for (int inew=0; inew<inewMax; inew++)
        {
          for (int iold=0; iold<ioldMax; iold++)
            {
              //sending arrays for distant domains
              if (isParallelMode() && _domain_selector->isMyDomain(iold) && !_domain_selector->isMyDomain(inew))
                {
                  int target=_domain_selector->getProcessorID(inew);
                  ParaMEDMEM::DataArrayDouble* field=initialCollection.getField(descriptionField,iold);
                  if (MyGlobals::_Verbose>10) 
                    std::cout << "proc " << _domain_selector->rank() << " : castAllFields sendDouble" << std::endl;
                  SendDataArrayDouble(field, target);
                }
              //receiving arrays from distant domains
              if (isParallelMode() && !_domain_selector->isMyDomain(iold) && _domain_selector->isMyDomain(inew))
                {
                  int source=_domain_selector->getProcessorID(iold);
                  //receive vector
                  if (MyGlobals::_Verbose>10) 
                    std::cout << "proc " << _domain_selector->rank() << " : castAllFields recvDouble" << std::endl;
                  ParaMEDMEM::DataArrayDouble* field=RecvDataArrayDouble(source);
                  remapDoubleField(inew,iold,field,nameArrayTo,descriptionField);
                }
            }
        }
#endif
      //local contributions and aggregation
      for (int inew=0; inew<inewMax; inew++)
        {
          for (int iold=0; iold<ioldMax; iold++)
            if (!isParallelMode() || ( _domain_selector->isMyDomain(iold) && _domain_selector->isMyDomain(inew)))
              {
                ParaMEDMEM::DataArrayDouble* field=initialCollection.getField(descriptionField,iold);
                remapDoubleField(inew,iold,field,nameArrayTo,descriptionField);
              }
        }
    }
}

void MEDPARTITIONER::MeshCollection::remapDoubleField(int inew, int iold, 
                                                       ParaMEDMEM::DataArrayDouble* fromArray,
                                                       std::string nameArrayTo,
                                                       std::string descriptionField)
//here we use 'cellFamily_ccI inew iold' set in remapIntField
{
  if (nameArrayTo!="cellFieldDouble") 
    throw INTERP_KERNEL::Exception("Error remapDoubleField only on cellFieldDouble");
  std::string key=Cle2ToStr("cellFamily_ccI",inew,iold);
  
  std::map<std::string,ParaMEDMEM::DataArrayInt*>::iterator it1;
  it1=_map_dataarray_int.find(key);
  if (it1==_map_dataarray_int.end())
    {
      std::cerr << "proc " << MyGlobals::_Rank << " : remapDoubleField key '" << key << "' not found" << std::endl;
      std::cerr << " trying remap of field double on cells : " << descriptionField << std::endl;
      return;
    }
  //create ccI in remapIntField
  ParaMEDMEM::DataArrayInt *ccI=it1->second;
  if (MyGlobals::_Verbose>300)
    std::cout << "proc " << MyGlobals::_Rank << " : remapDoubleField " << key << " size " << ccI->getNbOfElems() << std::endl;
  
  int nbcell=this->getMesh()[inew]->getNumberOfCells();
  int nbcomp=fromArray->getNumberOfComponents();
  int nbPtGauss=StrToInt(ExtractFromDescription(descriptionField, "nbPtGauss="));
  
  std::string tag="inewFieldDouble="+IntToStr(inew);
  key=descriptionField+SerializeFromString(tag);
  int fromArrayNbOfElem=fromArray->getNbOfElems();
  int fromArrayNbOfComp=fromArray->getNumberOfComponents();
  int fromArrayNbOfCell=fromArrayNbOfElem/fromArrayNbOfComp/nbPtGauss;
  
  if (MyGlobals::_Verbose>1000)
    {
      std::cout<<"proc " << MyGlobals::_Rank << " nbcell " << nbcell << " nbcomp " << nbcomp << " nbPtGauss " << nbPtGauss <<
        " fromArray nbOfElems " << fromArrayNbOfElem <<
        " nbTuples " << fromArray->getNumberOfTuples() <<
        " nbcells " << fromArrayNbOfCell <<
        " nbComponents " << fromArray->getNumberOfComponents() << std::endl;
    }

  ParaMEDMEM::DataArrayDouble* field=0;
  std::map<std::string,ParaMEDMEM::DataArrayDouble*>::iterator it2;
  it2=_map_dataarray_double.find(key);
  if (it2==_map_dataarray_double.end())
    {
      if (MyGlobals::_Verbose>300)
        std::cout << "proc "<< MyGlobals::_Rank << " : remapDoubleField key '" << key << "' not found and create it" << std::endl;
      field=ParaMEDMEM::DataArrayDouble::New();
      _map_dataarray_double[key]=field;
      field->alloc(nbcell*nbPtGauss,nbcomp);
      field->fillWithZero();
    }
  else
    {
      field=it2->second;
      if (field->getNumberOfTuples()!=nbcell*nbPtGauss || field->getNumberOfComponents()!=nbcomp)
        {
          std::cerr << "proc " << MyGlobals::_Rank << " : remapDoubleField3 pb of size " <<
            " trying remap of field double on cells : \n" << descriptionField << std::endl;
          return;
        }
    }
  
  if (nbPtGauss==1)
    {
      field->setPartOfValuesAdv(fromArray,ccI);
    }
  else
    {
      //replaced by setPartOfValuesAdv if nbPtGauss==1
      int iMax=ccI->getNbOfElems();
      int* pccI=ccI->getPointer();
      double* pField=field->getPointer();
      double* pFrom=fromArray->getPointer();
      int itarget, isource, delta=nbPtGauss*nbcomp;
      for (int i=0; i<iMax; i=i+2)  //cell
        {
          itarget=pccI[i];
          isource=pccI[i+1];
          if ((itarget<0) || (itarget>=nbcell) || (isource<0) || (isource>=fromArrayNbOfCell))
            throw INTERP_KERNEL::Exception("Error field override");
          int ita=itarget*delta;
          int iso=isource*delta;
          for (int k=0; k<delta; k++) pField[ita+k]=pFrom[iso+k]; //components and gausspoints
        }
    }
}

/*! constructing the MESH collection from a distributed file
 *
 * \param filename name of the master file containing the list of all the MED files
 */
MEDPARTITIONER::MeshCollection::MeshCollection(const std::string& filename)
  : _topology(0),
    _owns_topology(true),
    _driver(0),
    _domain_selector( 0 ),
    _i_non_empty_mesh(-1),
    _driver_type(MEDPARTITIONER::Undefined),
    _subdomain_boundary_creates(false),
    _family_splitting(false),
    _create_empty_groups(false),
    _joint_finder(0)
{
  try
    {
      _driver=new MeshCollectionMedXmlDriver(this);
      _driver->read (filename.c_str());
      _driver_type = MedXml;
    }
  catch(...) 
    {  // Handle all exceptions
      if ( _driver ) delete _driver; _driver=0;
      try
        {
          _driver=new MeshCollectionMedAsciiDriver(this);
          _driver->read (filename.c_str());
          _driver_type=MedAscii;
        }
      catch(...)
        {
          delete _driver;
          _driver=0;
          throw INTERP_KERNEL::Exception("file does not comply with any recognized format");
        }
    }
  for ( int idomain = 0; idomain < (int)_mesh.size(); ++idomain )
    if ( _mesh[idomain] && _mesh[idomain]->getNumberOfNodes() > 0 )
      _i_non_empty_mesh = idomain;
}

/*! Constructing the MESH collection from selected domains of a distributed file
 * 
 * \param filename  - name of the master file containing the list of all the MED files
 * \param domainSelector - selector of domains to load
 */
MEDPARTITIONER::MeshCollection::MeshCollection(const std::string& filename, ParaDomainSelector& domainSelector)
  : _topology(0),
    _owns_topology(true),
    _driver(0),
    _domain_selector( &domainSelector ),
    _i_non_empty_mesh(-1),
    _driver_type(MEDPARTITIONER::Undefined),
    _subdomain_boundary_creates(false),
    _family_splitting(false),
    _create_empty_groups(false),
    _joint_finder(0)
{
  std::string myfile=filename;
  std::size_t found=myfile.find(".xml");
  if (found!=std::string::npos) //file .xml
    {
      try
        {
          _driver=new MeshCollectionMedXmlDriver(this);
          _driver->read ( filename.c_str(), _domain_selector );
          _driver_type = MedXml;
        }
      catch(...)
        {  // Handle all exceptions
          delete _driver;
          throw INTERP_KERNEL::Exception("file .xml does not comply with any recognized format");
        }
    }
  else 
    {
      found=myfile.find(".med");
      if (found!=std::string::npos) //file .med single
        {
          //make a temporary file .xml and retry MedXmlDriver
          std::string xml="\
<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n \
<root>\n \
  <version maj=\"2\" min=\"3\" ver=\"1\"/>\n \
  <description what=\"\" when=\"\"/>\n \
  <content>\n \
    <mesh name=\"$meshName\"/>\n \
  </content>\n \
  <splitting>\n \
    <subdomain number=\"1\"/>\n \
    <global_numbering present=\"no\"/>\n \
  </splitting>\n \
  <files>\n \
    <subfile id=\"1\">\n \
      <name>$fileName</name>\n \
      <machine>localhost</machine>\n \
    </subfile>\n \
  </files>\n \
  <mapping>\n \
    <mesh name=\"$meshName\">\n \
      <chunk subdomain=\"1\">\n \
        <name>$meshName</name>\n \
      </chunk>\n \
    </mesh>\n \
  </mapping>\n \
</root>\n";
          std::vector<std::string> meshNames=MEDLoader::GetMeshNames(myfile.c_str());
          xml.replace(xml.find("$fileName"),9,myfile);
          xml.replace(xml.find("$meshName"),9,meshNames[0]);
          xml.replace(xml.find("$meshName"),9,meshNames[0]);
          xml.replace(xml.find("$meshName"),9,meshNames[0]);
          std::string nameFileXml(myfile);
          nameFileXml.replace(nameFileXml.find(".med"),4,".xml");
          std::string nameFileXmlDN,nameFileXmlBN;
          MEDLoaderBase::getDirAndBaseName(nameFileXml,nameFileXmlDN,nameFileXmlBN);
          nameFileXml=MEDLoaderBase::joinPath(nameFileXmlDN,"medpartitioner_"+nameFileXmlBN);
          if (_domain_selector->rank()==0) //only on to write it
            {
              std::ofstream f(nameFileXml.c_str());
              f<<xml;
              f.close();
            }
#ifdef HAVE_MPI2
           if (MyGlobals::_World_Size>1)
             MPI_Barrier(MPI_COMM_WORLD); //wait for creation of nameFileXml
#endif
          try
            {
              _driver=new MeshCollectionMedXmlDriver(this);
              _driver->read ( nameFileXml.c_str(), _domain_selector );
              _driver_type = MedXml;
            }
          catch(...)
            {  // Handle all exceptions
              delete _driver; _driver=0;
              throw INTERP_KERNEL::Exception("file medpartitioner_xxx.xml does not comply with any recognized format");
            }
        }
      else //no extension
        {
          try
            {
              _driver=new MeshCollectionMedAsciiDriver(this);
              _driver->read ( filename.c_str(), _domain_selector );
              _driver_type=MedAscii;
            }
          catch(...)
            {
              delete _driver;
              _driver=0;
              throw INTERP_KERNEL::Exception("file name with no extension does not comply with any recognized format");
            }
        }
    }
  // find non-empty domain mesh
  for ( int idomain = 0; idomain < (int)_mesh.size(); ++idomain )
    if ( _mesh[idomain] && _mesh[idomain]->getNumberOfNodes() > 0 )
      _i_non_empty_mesh = idomain;
  
  try
    {
      //check for all proc/file compatibility of _field_descriptions
#ifdef HAVE_MPI2
      _field_descriptions=AllgathervVectorOfString(MyGlobals::_Field_Descriptions);
#else
      _field_descriptions=MyGlobals::_Field_Descriptions;
#endif
    }
  catch(INTERP_KERNEL::Exception& e)
    {
      std::cerr << "proc " << MyGlobals::_Rank << " : INTERP_KERNEL_Exception : " << e.what() << std::endl;
      throw INTERP_KERNEL::Exception("Something wrong verifying coherency of files med ands fields");
    }
#ifdef HAVE_MPI2
  try
    {
      //check for all proc/file compatibility of _family_info
      std::vector<std::string> v2=AllgathervVectorOfString(VectorizeFromMapOfStringInt(_family_info));
      _family_info=DevectorizeToMapOfStringInt(v2);
    }
  catch(INTERP_KERNEL::Exception& e)
    {
      std::cerr << "proc " << MyGlobals::_Rank << " : INTERP_KERNEL_Exception : " << e.what() << std::endl;
      throw INTERP_KERNEL::Exception("Something wrong merging all familyInfo");
    }

  try
    {
      //check for all proc/file compatibility of _group_info
      std::vector<std::string> v2=AllgathervVectorOfString(VectorizeFromMapOfStringVectorOfString(_group_info));
      _group_info=DeleteDuplicatesInMapOfStringVectorOfString(DevectorizeToMapOfStringVectorOfString(v2));
    }
  catch(INTERP_KERNEL::Exception& e)
    {
      std::cerr << "proc " << MyGlobals::_Rank << " : INTERP_KERNEL_Exception : " << e.what() << std::endl;
      throw INTERP_KERNEL::Exception("Something wrong merging all groupInfo");
    }
#endif
}

/*! constructing the MESH collection from a sequential MED-file
 * 
 * \param filename MED file
 * \param meshname name of the mesh that is to be read
 */
MEDPARTITIONER::MeshCollection::MeshCollection(const std::string& filename, const std::string& meshname)
  : _topology(0),
    _owns_topology(true),
    _driver(0),
    _domain_selector( 0 ),
    _i_non_empty_mesh(-1),
    _name(meshname),
    _driver_type(MEDPARTITIONER::MedXml),
    _subdomain_boundary_creates(false),
    _family_splitting(false),
    _create_empty_groups(false),
    _joint_finder(0)
{
  try // avoid memory leak in case of inexistent filename
    {
      retrieveDriver()->readSeq (filename.c_str(),meshname.c_str());
    }
  catch (...)
    {
      delete _driver;
      _driver=0;
      throw INTERP_KERNEL::Exception("problem reading .med files");
    }
  if ( _mesh[0] && _mesh[0]->getNumberOfNodes() > 0 )
    _i_non_empty_mesh = 0;
}

MEDPARTITIONER::MeshCollection::~MeshCollection()
{
  for (int i=0; i<(int)_mesh.size();i++)
    if (_mesh[i]!=0) _mesh[i]->decrRef(); 
  
  for (int i=0; i<(int)_cell_family_ids.size();i++)
    if (_cell_family_ids[i]!=0)
      _cell_family_ids[i]->decrRef();
  
  for (int i=0; i<(int)_face_mesh.size();i++)
    if (_face_mesh[i]!=0)
      _face_mesh[i]->decrRef();
  
  for (int i=0; i<(int)_face_family_ids.size();i++)
    if (_face_family_ids[i]!=0)
      _face_family_ids[i]->decrRef();
  
  for (std::map<std::string, ParaMEDMEM::DataArrayInt*>::iterator it=_map_dataarray_int.begin() ; it!=_map_dataarray_int.end(); it++ )
    if ((*it).second!=0)
      (*it).second->decrRef();
  
  for (std::map<std::string, ParaMEDMEM::DataArrayDouble*>::iterator it=_map_dataarray_double.begin() ; it!=_map_dataarray_double.end(); it++ )
    if ((*it).second!=0)
      (*it).second->decrRef();
  
  delete _driver;
  if (_topology!=0 && _owns_topology)
    delete _topology;
#ifdef HAVE_MPI2
  delete _joint_finder;
#endif
}

/*! constructing the MESH collection from a file
 * 
 * The method creates as many MED-files as there are domains in the 
 * collection. It also creates a master file that lists all the MED files.
 * The MED files created in ths manner contain joints that describe the 
 * connectivity between subdomains.
 * 
 * \param filename name of the master file that will contain the list of the MED files
 * 
 */
void MEDPARTITIONER::MeshCollection::write(const std::string& filename)
{
  //building the connect zones necessary for writing joints
  //   if (_topology->nbDomain()>1)
  //     buildConnectZones();
  //suppresses link with driver so that it can be changed for writing
  delete _driver;
  _driver=0;
  retrieveDriver()->write (filename.c_str(), _domain_selector);
}

/*! creates or gets the link to the collection driver
 */
MEDPARTITIONER::MeshCollectionDriver* MEDPARTITIONER::MeshCollection::retrieveDriver()
{
  if (_driver==0)
    {
      switch(_driver_type)
        {
        case MedXml:
          _driver=new MeshCollectionMedXmlDriver(this);
          break;
        case MedAscii:
          _driver=new MeshCollectionMedAsciiDriver(this);
          break;
        default:
          throw INTERP_KERNEL::Exception("Unrecognized driver");
        }
    }
  return _driver;
}


/*! gets an existing driver
 * 
 */
MEDPARTITIONER::MeshCollectionDriver* MEDPARTITIONER::MeshCollection::getDriver() const
{
  return _driver;
}

/*!
 * retrieves the mesh dimension
*/
int MEDPARTITIONER::MeshCollection::getMeshDimension() const
{
  return _i_non_empty_mesh < 0 ? -1 : _mesh[_i_non_empty_mesh]->getMeshDimension();
}

int MEDPARTITIONER::MeshCollection::getNbOfLocalMeshes() const
{
  int nb=0;
  for (int i=0; i<_mesh.size(); i++)
    {
      if (_mesh[i]) nb++;
    }
  return nb;
}

int MEDPARTITIONER::MeshCollection::getNbOfLocalCells() const
{
  int nb=0;
  for (int i=0; i<_mesh.size(); i++)
    {
      if (_mesh[i]) nb=nb+_mesh[i]->getNumberOfCells();
    }
  return nb;
}

int MEDPARTITIONER::MeshCollection::getNbOfLocalFaces() const
{
  int nb=0;
  for (int i=0; i<_face_mesh.size(); i++)
    {
      if (_face_mesh[i]) nb=nb+_face_mesh[i]->getNumberOfCells();
    }
  return nb;
}

std::vector<ParaMEDMEM::MEDCouplingUMesh*>& MEDPARTITIONER::MeshCollection::getMesh()
{
  return _mesh;
}

std::vector<ParaMEDMEM::MEDCouplingUMesh*>& MEDPARTITIONER::MeshCollection::getFaceMesh()
{
  return _face_mesh;
}

ParaMEDMEM::MEDCouplingUMesh* MEDPARTITIONER::MeshCollection::getMesh(int idomain) const
{
  return _mesh[idomain];
}

ParaMEDMEM::MEDCouplingUMesh* MEDPARTITIONER::MeshCollection::getFaceMesh(int idomain)
{
  return _face_mesh[idomain];
}

std::vector<MEDPARTITIONER::ConnectZone*>& MEDPARTITIONER::MeshCollection::getCZ()
{
  return _connect_zones;
}

MEDPARTITIONER::Topology* MEDPARTITIONER::MeshCollection::getTopology() const
{
  return _topology;
}

void MEDPARTITIONER::MeshCollection::setTopology(Topology* topo)
{
  if (_topology!=0)
    {
      throw INTERP_KERNEL::Exception("topology is already set");
    }
  else
    _topology = topo;
}

/*! Method creating the cell graph
 * 
 * \param array returns the pointer to the structure that contains the graph 
 * \param edgeweight returns the pointer to the table that contains the edgeweights
 *        (only used if indivisible regions are required)
 */
void MEDPARTITIONER::MeshCollection::buildCellGraph(MEDPARTITIONER::SkyLineArray* & array, int *& edgeweights )
{
  using std::multimap;
  using std::vector;
  using std::make_pair;
  using std::pair;
  
  std::multimap< int, int > node2cell;
  std::map< pair<int,int>, int > cell2cellcounter;
  std::multimap<int,int> cell2cell;

  std::vector<std::vector<std::multimap<int,int> > > commonDistantNodes;
  int nbdomain=_topology->nbDomain();
#ifdef HAVE_MPI2
  if (isParallelMode())
    {
      _joint_finder=new JointFinder(*this);
      _joint_finder->findCommonDistantNodes();
      commonDistantNodes=_joint_finder->getDistantNodeCell();
    }

  if (MyGlobals::_Verbose>500)
    _joint_finder->print();
#endif

  if (MyGlobals::_Verbose>50)
    std::cout<<"getting nodal connectivity"<<std::endl;
  //looking for reverse nodal connectivity i global numbering
  for (int idomain=0; idomain<nbdomain; idomain++)
    {
      if (isParallelMode() && !_domain_selector->isMyDomain(idomain))
        continue;
    
      ParaMEDMEM::DataArrayInt* index=ParaMEDMEM::DataArrayInt::New();
      ParaMEDMEM::DataArrayInt* revConn=ParaMEDMEM::DataArrayInt::New();
      int nbNodes=_mesh[idomain]->getNumberOfNodes();
      _mesh[idomain]->getReverseNodalConnectivity(revConn,index);
      //problem saturation over 1 000 000 nodes for 1 proc
      if (MyGlobals::_Verbose>100)
        std::cout << "proc " << MyGlobals::_Rank << " : getReverseNodalConnectivity done on " << nbNodes << " nodes" << std::endl;
      int* index_ptr=index->getPointer();
      int* revConnPtr=revConn->getPointer();
      for (int i=0; i<nbNodes; i++)
        {
          for (int icell=index_ptr[i]; icell<index_ptr[i+1]; icell++)
            {
              int globalNode=_topology->convertNodeToGlobal(idomain,i);
              int globalCell=_topology->convertCellToGlobal(idomain,revConnPtr[icell]);
              node2cell.insert(make_pair(globalNode, globalCell));
            }
        }
      revConn->decrRef();
      index->decrRef();
#ifdef HAVE_MPI2
      for (int iother=0; iother<nbdomain; iother++)
        {
          std::multimap<int,int>::iterator it;
          int isource=idomain;
          int itarget=iother;
          for (it=_joint_finder->getDistantNodeCell()[isource][itarget].begin(); 
               it!=_joint_finder->getDistantNodeCell()[isource][itarget].end(); it++)
            {
              int globalNode=_topology->convertNodeToGlobal(idomain,(*it).first);
              int globalCell=(*it).second;
              node2cell.insert(make_pair(globalNode, globalCell));
            }
        }
#endif
    }  //endfor idomain

  //creating graph arcs (cell to cell relations)
  //arcs are stored in terms of (index,value) notation
  // 0 3 5 6 6
  // 1 2 3 2 3 3 
  // means 6 arcs (0,1), (0,2), (0,3), (1,2), (1,3), (2,3)
  // in present version arcs are not doubled but reflexive (1,1) arcs are present for each cell
 
  //warning here one node have less than or equal effective number of cell with it
  //but cell could have more than effective nodes
  //because other equals nodes in other domain (with other global inode)
  if (MyGlobals::_Verbose>50)
    std::cout<< "proc " << MyGlobals::_Rank << " : creating graph arcs on nbNodes " << _topology->nbNodes() << std::endl;
 
  for (int inode=0;inode<_topology->nbNodes();inode++)
    {
      typedef multimap<int,int>::const_iterator MI;
      std::pair <MI,MI> nodeRange=node2cell.equal_range(inode);
      for (MI cell1=nodeRange.first;cell1!=nodeRange.second;cell1++)
        for (MI cell2=nodeRange.first;cell2!=cell1;cell2++)
          {
            int icell1=cell1->second;
            int icell2=cell2->second;
            if (icell1>icell2) {int tmp=icell1; icell1=icell2; icell2=tmp;}
            std::map<pair<int,int>,int>::iterator it=cell2cellcounter.find(make_pair(icell1,icell2));
            if (it==cell2cellcounter.end()) cell2cellcounter.insert(make_pair(make_pair(icell1,icell2),1));
            else (it->second)++;
          }
    }      
  // for (int icell1=0; icell1<_topology->nbCells(); icell1++)  //on all nodes
  //   {
  //     typedef multimap<int,int>::const_iterator MI;
  //     std::pair <MI,MI> nodeRange=cell2node.equal_range(icell1);
  //     for (MI node1=nodeRange.first; node1!=nodeRange.second; node1++)  //on nodes with icell
  //       {
  //         std::pair<MI,MI> cellRange=node2cell.equal_range(node1->second);
  //         for (MI cell2=cellRange.first; cell2!=cellRange.second; cell2++)  //on one of these cell
  //           {
  //             int icell2=cell2->second;
  //             std::map<pair<int,int>,int>::iterator it=cell2cellcounter.find(make_pair(icell1,icell2));
  //             if (it==cell2cellcounter.end()) cell2cellcounter.insert(make_pair(make_pair(icell1,icell2),1));
  //             else (it->second)++;
  //           }
  //       }
  //   }


  //converting the counter to a multimap structure
  for (std::map<pair<int,int>,int>::const_iterator it=cell2cellcounter.begin();
       it!=cell2cellcounter.end();
       it++)
    if (it->second>=3) 
      {
        cell2cell.insert(std::make_pair(it->first.first,it->first.second)); //should be adapted for 2D!
        cell2cell.insert(std::make_pair(it->first.second, it->first.first));
      }

  
  if (MyGlobals::_Verbose>50)
    std::cout << "proc " << MyGlobals::_Rank << " : create skylinearray" << std::endl;
  //filling up index and value to create skylinearray structure
  std::vector <int> index,value;
  index.push_back(0);
  int idep=0;
  
  for (int idomain=0; idomain<nbdomain; idomain++)
    {
      if (isParallelMode() && !_domain_selector->isMyDomain(idomain)) continue;
      int nbCells=_mesh[idomain]->getNumberOfCells();
      for (int icell=0; icell<nbCells; icell++)
        {
          int size=0;
          int globalCell=_topology->convertCellToGlobal(idomain,icell);
          multimap<int,int>::iterator it;
          pair<multimap<int,int>::iterator,multimap<int,int>::iterator> ret;
          ret=cell2cell.equal_range(globalCell);
          for (it=ret.first; it!=ret.second; ++it)
            {
              int ival=(*it).second; //no adding one existing yet
              for (int i=idep ; i<idep+size ; i++)
                {
                  if (value[i]==ival)
                    {
                      ival= -1; break;
                    }
                }
              if (ival!= -1)
                {
                  value.push_back(ival);
                  size++;
                }
            }
          idep=index[index.size()-1]+size;
          index.push_back(idep);
        }
    }
  
  array=new MEDPARTITIONER::SkyLineArray(index,value);

  if (MyGlobals::_Verbose>100)
    {
      std::cout << "\nproc " << _domain_selector->rank() << " : end MeshCollection::buildCellGraph " <<
        index.size()-1 << " " << value.size() << std::endl;
      int max=index.size()>15?15:index.size();
      if (index.size()>1)
        {
          for (int i=0; i<max; ++i)
            std::cout<<index[i]<<" ";
          std::cout << "... " << index[index.size()-1] << std::endl;
          for (int i=0; i<max; ++i)
            std::cout<< value[i] << " ";
          int ll=index[index.size()-1]-1;
          std::cout << "... (" << ll << ") " << value[ll-1] << " " << value[ll] << std::endl;
        }
    }
  
}


/*! Creates the partition corresponding to the cell graph and the partition number
 * 
 * \param nbdomain number of subdomains for the newly created graph
 * 
 * returns a topology based on the new graph
 */
MEDPARTITIONER::Topology* MEDPARTITIONER::MeshCollection::createPartition(int nbdomain,
                                                                          Graph::splitter_type split, 
                                                                          const std::string& options_string,
                                                                          int *user_edge_weights,
                                                                          int *user_vertices_weights)
{
  if (MyGlobals::_Verbose>10)
    std::cout << "proc " << MyGlobals::_Rank << " : MeshCollection::createPartition : Building cell graph" << std::endl;
  
  if (nbdomain <1)
    throw INTERP_KERNEL::Exception("Number of subdomains must be > 0");
  MEDPARTITIONER::SkyLineArray* array=0;
  int* edgeweights=0;
  buildCellGraph(array,edgeweights);
  
  Graph* cellGraph;
  switch (split)
    {
    case Graph::METIS:
#if defined(MED_ENABLE_PARMETIS) || defined(MED_ENABLE_METIS)
      if (MyGlobals::_Verbose>10)
        std::cout << "METISGraph" << std::endl;
      cellGraph=new METISGraph(array,edgeweights);
#else
      throw INTERP_KERNEL::Exception("MeshCollection::createPartition : PARMETIS/METIS is not available. Check your products, please.");
#endif
      break;
    case Graph::SCOTCH:
#ifdef MED_ENABLE_SCOTCH
      if (MyGlobals::_Verbose>10)
        std::cout << "SCOTCHGraph" << std::endl;
      cellGraph=new SCOTCHGraph(array,edgeweights);
#else
      throw INTERP_KERNEL::Exception("MeshCollection::createPartition : SCOTCH is not available. Check your products, please.");
#endif
      break;
    }

  //!user-defined weights
  if (user_edge_weights!=0) 
    cellGraph->setEdgesWeights(user_edge_weights);
  if (user_vertices_weights!=0)
    cellGraph->setVerticesWeights(user_vertices_weights);

  if (MyGlobals::_Is0verbose>10)
    std::cout << "partitioning graph on " << nbdomain << " domains" << std::endl;
  cellGraph->partGraph(nbdomain, options_string, _domain_selector);

  if (MyGlobals::_Is0verbose>10)
    std::cout << "building new topology" << std::endl;
  //cellGraph is a shared pointer 
  Topology *topology=0;
  topology=new ParallelTopology (cellGraph, getTopology(), nbdomain, getMeshDimension());
  //cleaning
  delete [] edgeweights;
  delete cellGraph;
  if (MyGlobals::_Verbose>11)
    std::cout << "proc " << MyGlobals::_Rank << " : end MeshCollection::createPartition" << std::endl;
  return topology;
}

/*! Creates a topology for a partition specified by the user
 * 
 * \param table user-specified partition (for each cell contains the domain number from 0 to n-1)
 * 
 * returns a topology based on the new partition
 */
MEDPARTITIONER::Topology* MEDPARTITIONER::MeshCollection::createPartition(const int* partition)
{
  MEDPARTITIONER::SkyLineArray* array=0;
  int* edgeweights=0;

  buildCellGraph(array,edgeweights);
  Graph* cellGraph;
  std::set<int> domains;
  for (int i=0; i<_topology->nbCells(); i++)
    {
      domains.insert(partition[i]);
    }
  cellGraph=new UserGraph(array, partition, _topology->nbCells());
  
  //cellGraph is a shared pointer 
  Topology *topology=0;
  int nbdomain=domains.size();
  topology=new ParallelTopology (cellGraph, getTopology(), nbdomain, getMeshDimension());
  //  if (array!=0) delete array;
  delete cellGraph;
  return topology;
}

void MEDPARTITIONER::MeshCollection::setDomainNames(const std::string& name)
{
  for (int i=0; i<_topology->nbDomain(); i++)
    {
      std::ostringstream oss;
      oss<<name<<"_"<<i;
      if (!isParallelMode() || _domain_selector->isMyDomain(i))
        _mesh[i]->setName(oss.str().c_str());
    }
}

ParaMEDMEM::DataArrayDouble *MEDPARTITIONER::MeshCollection::getField(std::string descriptionField, int iold)
//getField look for and read it if not done, and assume decrRef() in ~MeshCollection;
//something like MEDCouplingFieldDouble *f2=MEDLoader::ReadFieldCell(name.c_str(),f1->getMesh()->getName(),0,f1->getName(),0,1);
{
  int rank=MyGlobals::_Rank;
  std::string tag="ioldFieldDouble="+IntToStr(iold);
  std::string descriptionIold=descriptionField+SerializeFromString(tag);
  if (_map_dataarray_double.find(descriptionIold)!=_map_dataarray_double.end())
    {
      if (MyGlobals::_Verbose>300)
        std::cout << "proc " << rank << " : YET READ getField : " << descriptionIold << std::endl;
      ParaMEDMEM::DataArrayDouble* res=_map_dataarray_double[descriptionIold];
      return res;
    }
  if (MyGlobals::_Verbose>200)
    std::cout << "proc " << rank << " : TO BE READ getField : " << descriptionIold << std::endl;
  std::string description, fileName, meshName, fieldName;
  int typeField, DT, IT, entity;
  fileName=MyGlobals::_File_Names[iold];
  if (MyGlobals::_Verbose>10) 
    std::cout << "proc " << MyGlobals::_Rank << " : in " << fileName << " " << iold << " " << descriptionIold << std::endl;
  FieldShortDescriptionToData(descriptionIold, fieldName, typeField, entity, DT, IT);
  meshName=MyGlobals::_Mesh_Names[iold];
  
  ParaMEDMEM::MEDCouplingFieldDouble* f2=MEDLoader::ReadField((ParaMEDMEM::TypeOfField) typeField,
                                                              fileName.c_str(), meshName.c_str(), 0, fieldName.c_str(), DT, IT);
  
  ParaMEDMEM::DataArrayDouble* res=f2->getArray();
  //to know names of components
  std::vector<std::string> browse=BrowseFieldDouble(f2);
  std::string localFieldInformation=descriptionIold+SerializeFromVectorOfString(browse);
  if (MyGlobals::_Verbose>10)
    std::cout << "proc " << MyGlobals::_Rank << " : localFieldInformation : " << localFieldInformation << std::endl;
  MyGlobals::_General_Informations.push_back(localFieldInformation);
  res->incrRef();
  f2->decrRef();
  _map_dataarray_double[descriptionIold]=res; 
  return res;
}

void MEDPARTITIONER::MeshCollection::prepareFieldDescriptions()
//to have unique valid fields names/pointers/descriptions for partitionning
//filter _field_descriptions to be in all procs compliant and equal
{
  int nbfiles=MyGlobals::_File_Names.size(); //nb domains
  std::vector<std::string> r2;
  //from allgatherv then vector(procs) of serialised vector(fields) of vector(description) data
  for (int i=0; i<(int)_field_descriptions.size(); i++)
    {
      std::vector<std::string> r1=DeserializeToVectorOfString(_field_descriptions[i]);
      for (int ii=0; ii<(int)r1.size(); ii++)
        r2.push_back(r1[ii]);
    }
  //here vector(procs*fields) of serialised vector(description) data
  _field_descriptions=r2;
  int nbfields=_field_descriptions.size(); //on all domains
  if ((nbfields%nbfiles)!=0)
    {
      if (MyGlobals::_Rank==0)
        {
          std::cerr<< "\nERROR : incoherent number of fields references in all files .med\n" << std::endl
              << "fileMedNames :" << std::endl
              << ReprVectorOfString(MyGlobals::_File_Names)
              << "field_descriptions :" << std::endl
              << ReprVectorOfString(MyGlobals::_Field_Descriptions);
        }
      throw INTERP_KERNEL::Exception("incoherent number of fields references in all files .med\n");
    }
  _field_descriptions.resize(nbfields/nbfiles);
  for (int i=0; i<(int)_field_descriptions.size(); i++)
    {
      std::string str=_field_descriptions[i];
      str=EraseTagSerialized(str,"idomain=");
      str=EraseTagSerialized(str,"fileName=");
      _field_descriptions[i]=str;
    }
}

//returns true if inodes of a face are in inodes of a cell
bool isFaceOncell(std::vector< int >& inodesFace, std::vector< int >&  inodesCell)
{
  int ires=0;
  int nbok=inodesFace.size();
  for (int i=0; i<nbok; i++)
    {
      int ii=inodesFace[i];
      if (ii<0)
        std::cout << "isFaceOncell problem inodeface<0" << std::endl;
      for (int j=0; j<(int)inodesCell.size(); j++)
        {
          if (ii==inodesCell[j])
            {
              ires=ires+1;
              break; //inode of face found
            }
        }
      if (ires<i+1)
        break; //inode of face not found do not continue...
    }
  return (ires==nbok);
}

void MEDPARTITIONER::MeshCollection::filterFaceOnCell()
{
  for (int inew=0; inew<_topology->nbDomain(); inew++)
    {
      if (isParallelMode() && _domain_selector->isMyDomain(inew))
        {
          if (MyGlobals::_Verbose>200) 
            std::cout << "proc " << MyGlobals::_Rank << " : filterFaceOnCell on inewDomain " << inew << " nbOfFaces " << _face_mesh[inew]->getNumberOfCells() << std::endl;
          ParaMEDMEM::MEDCouplingUMesh* mcel=_mesh[inew];
          ParaMEDMEM::MEDCouplingUMesh* mfac=_face_mesh[inew];
      
          //to have cellnode=f(facenode)... inodeCell=nodeIds[inodeFace]
          std::vector<int> nodeIds;
          getNodeIds(*mcel, *mfac, nodeIds);
          if (nodeIds.size()==0)
            continue;  //one empty mesh nothing to do
      
          ParaMEDMEM::DataArrayInt *revNodalCel=ParaMEDMEM::DataArrayInt::New();
          ParaMEDMEM::DataArrayInt *revNodalIndxCel=ParaMEDMEM::DataArrayInt::New();
          mcel->getReverseNodalConnectivity(revNodalCel,revNodalIndxCel);
          int *revC=revNodalCel->getPointer();
          int *revIndxC=revNodalIndxCel->getPointer();
      
          std::vector< int > faceOnCell;
          std::vector< int > faceNotOnCell;
          int nbface=mfac->getNumberOfCells();
          for (int iface=0; iface<nbface; iface++)
            {
              bool ok;
              std::vector< int > inodesFace;
              mfac->getNodeIdsOfCell(iface, inodesFace);
              int nbnodFace=inodesFace.size();
              //set inodesFace in mcel
              for (int i=0; i<nbnodFace; i++) inodesFace[i]=nodeIds[inodesFace[i]];
              int inod=inodesFace[0];
              if (inod<0)
                std::cout << "filterFaceOnCell problem 1" << std::endl;
              int nbcell=revIndxC[inod+1]-revIndxC[inod];
              for (int j=0; j<nbcell; j++) //look for each cell with inod
                {
                  int icel=revC[revIndxC[inod]+j];
                  std::vector< int > inodesCell;
                  mcel->getNodeIdsOfCell(icel, inodesCell);
                  ok=isFaceOncell(inodesFace, inodesCell);
                  if (ok) break;
                }
              if (ok)
                {
                  faceOnCell.push_back(iface);
                }
              else
                {
                  faceNotOnCell.push_back(iface);
                  if (MyGlobals::_Is0verbose>300)
                    std::cout << "face NOT on cell " << iface << " " << faceOnCell.size()-1 << std::endl;
                }
            }
      
          revNodalCel->decrRef();
          revNodalIndxCel->decrRef();
      
          std::string keyy;
          keyy=Cle1ToStr("filterFaceOnCell",inew);
          _map_dataarray_int[keyy]=CreateDataArrayIntFromVector(faceOnCell);
          keyy=Cle1ToStr("filterNotFaceOnCell",inew);
          _map_dataarray_int[keyy]=CreateDataArrayIntFromVector(faceNotOnCell);
        }
    }
}
