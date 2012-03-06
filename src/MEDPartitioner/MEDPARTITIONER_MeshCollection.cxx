// Copyright (C) 2007-2011  CEA/DEN, EDF R&D
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
#include "MEDPARTITIONER_ParallelTopology.hxx"
#include "MEDPARTITIONER_ParaDomainSelector.hxx"
#include "MEDPARTITIONER_JointFinder.hxx"
#include "MEDPARTITIONER_Graph.hxx"
#include "MEDPARTITIONER_UserGraph.hxx"
#include "MEDPARTITIONER_Utils.hxx" 

#include "MEDLoader.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingNormalizedUnstructuredMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "PointLocator3DIntersectorP0P0.hxx"
#include "BBTree.txx"

#ifdef HAVE_MPI2
#include <mpi.h>
#endif

#ifdef ENABLE_METIS
#include "MEDPARTITIONER_MetisGraph.hxx"
#endif
#ifdef ENABLE_SCOTCH
#include "MEDPARTITIONER_ScotchGraph.hxx"
#endif

#include <vector>
#include <string>
#include <limits>
#include <set>
#include <iostream>
#include <fstream>

using namespace std;
using namespace ParaMEDMEM;
using namespace MEDPARTITIONER;

MeshCollection::MeshCollection()
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
  cout<<"coucou"<<endl;
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

MeshCollection::MeshCollection(MeshCollection& initialCollection, 
                               Topology* topology, 
                               bool family_splitting, 
                               bool create_empty_groups)
  : _name(initialCollection._name),
    _topology(topology),
    _owns_topology(false),
    _driver(0),
    _domain_selector( initialCollection._domain_selector ),
    _i_non_empty_mesh(-1),
    _driver_type(MEDPARTITIONER::MedXml),
    _subdomain_boundary_creates(false),
    _family_splitting(family_splitting),
    _create_empty_groups(create_empty_groups),
    _joint_finder(0)
{
  std::vector<std::vector<std::vector<int> > > new2oldIds(initialCollection.getTopology()->nbDomain());
  if (MyGlobals::_Verbose>10) std::cout<<"proc "<<MyGlobals::_Rank<<" : castCellMeshes"<<std::endl;
  castCellMeshes(initialCollection, new2oldIds);

  //defining the name for the collection and the underlying meshes
  setName(initialCollection.getName());

  /////////////////
  //treating faces
  /////////////////

  if (MyGlobals::_Is0verbose) std::cout<<"treating faces"<<std::endl;
  NodeMapping nodeMapping;
  //nodeMapping contains the mapping between old nodes and new nodes
  // (iolddomain,ioldnode)->(inewdomain,inewnode)
  createNodeMapping(initialCollection, nodeMapping);
  //cvw std::cout<<"castMeshes"<<std::endl;
  std::vector<std::vector<std::vector<int> > > new2oldFaceIds;
  castFaceMeshes(initialCollection, nodeMapping, new2oldFaceIds);

  ////////////////////
  //treating families
  ////////////////////

  if (MyGlobals::_Is0verbose) 
    if (isParallelMode()) std::cout<<"ParallelMode on "<<topology->nbDomain()<<" Domains"<<std::endl;
    else std::cout<<"NOT ParallelMode on "<<topology->nbDomain()<<" Domains"<<std::endl;
    
  if (MyGlobals::_Is0verbose>10) std::cout<<"treating cell and face families"<<std::endl;
  
  castIntField2(initialCollection.getMesh(),
                this->getMesh(),
                initialCollection.getCellFamilyIds(),
                "cellFamily");
  castIntField2(initialCollection.getFaceMesh(), 
                this->getFaceMesh(),
                initialCollection.getFaceFamilyIds(),
                "faceFamily");

  //////////////////
  //treating groups
  //////////////////
  if (MyGlobals::_Is0verbose) std::cout<<"treating groups"<<std::endl;
  _familyInfo=initialCollection.getFamilyInfo();
  _groupInfo=initialCollection.getGroupInfo();
  
  //////////////////
  //treating fields
  //////////////////
  if (MyGlobals::_Is0verbose) std::cout<<"treating fields"<<std::endl;
  //cvwat08
  castAllFields(initialCollection,"cellFieldDouble");
  //castAllFields(initialCollection,"faceFieldsIds");
}

/*!
  Creates the meshes using the topology underlying he mesh collection and the mesh data 
  coming from the ancient collection
  \param initialCollection collection from which the data is extracted to create the new meshes
*/

void MeshCollection::castCellMeshes(
                                    MeshCollection& initialCollection,
                                    std::vector<std::vector<std::vector<int> > >& new2oldIds)
{
  if (_topology==0) throw INTERP_KERNEL::Exception(LOCALIZED("Topology has not been defined on call to castCellMeshes"));
  
  int nbNewDomain=_topology->nbDomain();
  int nbOldDomain=initialCollection.getTopology()->nbDomain();
  
  _mesh.resize(nbNewDomain);
  int rank=MyGlobals::_Rank;
  //if (MyGlobals::_Verbose>10) std::cout<<"proc "<<rank<<" : castCellMeshes splitting"<<std::endl;
  //splitting the initial domains into smaller bits
  std::vector<std::vector<ParaMEDMEM::MEDCouplingUMesh*> > splitMeshes;
  splitMeshes.resize(nbNewDomain);
  for (int inew=0; inew<nbNewDomain; inew++)
    {
      splitMeshes[inew].resize(nbOldDomain, (ParaMEDMEM::MEDCouplingUMesh*)0);
      /*std::fill( &(splitMeshes[inew][0]), 
        &(splitMeshes[inew][0])+splitMeshes[inew].size(), 
        (ParaMEDMEM::MEDCouplingUMesh*)0 );*/
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
          //cvw work locally
          _topology->convertGlobalCellList(&globalids[0],size,&ilocalnew[0],&ipnew[0]);
      
          new2oldIds[iold].resize(nbNewDomain);
          for (int i=0; i<ilocalnew.size(); i++)
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
                std::cout<<"proc "<<rank<<" : a splitMesh iold inew NbCells "<<iold<<" "<<inew<<" "
                         <<splitMeshes[inew][iold]->getNumberOfCells()<<std::endl;
            }
        }
    }

  if (isParallelMode())
    {
      //if (MyGlobals::_Verbose>300) std::cout<<"proc "<<rank<<" : castCellMeshes send/receive"<<std::endl;
      for (int iold=0; iold<nbOldDomain; iold++)
        for(int inew=0; inew<nbNewDomain; inew++)
          {
            if (initialCollection._domain_selector->isMyDomain(iold) && _domain_selector->isMyDomain(inew)) continue;

            if(initialCollection._domain_selector->isMyDomain(iold))
              _domain_selector->sendMesh(*(splitMeshes[inew][iold]),_domain_selector->getProcessorID(inew));

            if (_domain_selector->isMyDomain(inew))
              _domain_selector->recvMesh(splitMeshes[inew][iold],_domain_selector->getProcessorID(iold));

          }
    }
  
  //fusing the split meshes
  if (MyGlobals::_Verbose>200) std::cout<<"proc "<<rank<<" : castCellMeshes fusing"<<std::endl;
  for (int inew=0; inew<nbNewDomain ;inew++)
    {
      vector<const ParaMEDMEM::MEDCouplingUMesh*> meshes;
    
      for (int i=0; i< splitMeshes[inew].size();i++)
        if (splitMeshes[inew][i]!=0) 
          if (splitMeshes[inew][i]->getNumberOfCells()>0)
            meshes.push_back(splitMeshes[inew][i]);

      if (!isParallelMode()||_domain_selector->isMyDomain(inew))
        {
          if (meshes.size()==0) 
            {
              _mesh[inew]=CreateEmptyMEDCouplingUMesh();
              //throw INTERP_KERNEL::Exception(LOCALIZED("castCellMeshes fusing : no meshes"));
              cout<<"WARNING : castCellMeshes fusing : no meshes try another number of processors"<<endl;
            }
          else
            {
              _mesh[inew]=ParaMEDMEM::MEDCouplingUMesh::MergeUMeshes(meshes);
              bool areNodesMerged;
              int nbNodesMerged;
              ParaMEDMEM::DataArrayInt* array=_mesh[inew]->mergeNodes(1e-12,areNodesMerged,nbNodesMerged);
              array->decrRef(); // array is not used in this case
              _mesh[inew]->zipCoords();
            }
        }
      for (int i=0; i< splitMeshes[inew].size(); i++)
        if (splitMeshes[inew][i]!=0) splitMeshes[inew][i]->decrRef();
    }  
  if (MyGlobals::_Verbose>300) std::cout<<"proc "<<rank<<" : castCellMeshes end fusing"<<std::endl;
}

/*!
  \param initialCollection source mesh collection 
  \param nodeMapping structure containing the correspondency between nodes in the initial collection and the node(s) in the new collection
*/
void MeshCollection::createNodeMapping( MeshCollection& initialCollection, NodeMapping& nodeMapping)
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
          bbox=new double[nvertices*6];
          ParaMEDMEM::DataArrayDouble* coords = initialCollection.getMesh(iold)->getCoords();
          double* coordsPtr=coords->getPointer();

          for (int i=0; i<nvertices*3;i++)
            {
              bbox[i*2]=coordsPtr[i]-1e-6;
              bbox[i*2+1]=coordsPtr[i]+1e-6;
            }
          tree=new BBTree<3>(bbox,0,0,nvertices,1e-9);
        }

      for (int inew=0; inew<_topology->nbDomain(); inew++) //cvwat12
        {
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

//getNodeIds(meshCell, meshFace, nodeIds)
//inodeCell=nodeIds[inodeFace]
//(put the biggest mesh in One)
//if no corresponding node then inodeCell = -1
void getNodeIds(ParaMEDMEM::MEDCouplingUMesh& meshOne, ParaMEDMEM::MEDCouplingUMesh& meshTwo, vector<int>& nodeIds)
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
      bbox[i*2]=coordsPtr[i]-1e-6;
      bbox[i*2+1]=coordsPtr[i]+1e-6;
    }
  tree=new BBTree<3>(bbox,0,0,nv1,1e-9);
  
  int nv2=meshTwo.getNumberOfNodes();
  nodeIds.resize(nv2,-1);
  coords=meshTwo.getCoords();
  for (int inode=0; inode<nv2; inode++)
    {
      double* coordsPtr=coords->getPointer()+inode*3;
      vector<int> elems;
      tree->getElementsAroundPoint(coordsPtr,elems);
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
void MeshCollection::castFaceMeshes(MeshCollection& initialCollection,
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
  
  vector<ParaMEDMEM::MEDCouplingUMesh*>& meshesCastFrom=initialCollection.getFaceMesh();
  vector<ParaMEDMEM::MEDCouplingUMesh*>& meshesCastTo=this->getFaceMesh();
  
  vector< vector<ParaMEDMEM::MEDCouplingUMesh*> > splitMeshes;
  int newSize=_topology->nbDomain();
  int fromSize=meshesCastFrom.size();
  
  splitMeshes.resize(newSize);
  for (int inew=0;inew<newSize;inew++) splitMeshes[inew].resize(fromSize);
  new2oldIds.resize(fromSize);
  for (int iold=0;iold<fromSize;iold++) new2oldIds[iold].resize(newSize);
  
  //init null pointer for empty meshes
  for (int inew=0;inew<newSize;inew++)
    {
      for (int iold=0;iold<fromSize;iold++)
        {
          splitMeshes[inew][iold]=0; //null for empty meshes
          new2oldIds[iold][inew].clear();
        }
    }

  //loop over the old domains to analyse the faces and decide 
  //on which new domain they belong
  
  for (int iold=0; iold<fromSize;iold++)
    {
      if (isParallelMode() && !_domain_selector->isMyDomain(iold)) continue;
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
              for (int inode=0; inode<nodes.size(); inode++)
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
                  if (iter->second==nodes.size())
                    //cvw eligible but may be have to be face of a cell of this->getMesh()[inew]?
                    //it is not sure here...
                    //done before writeMedfile on option?... see filterFaceOnCell()
                    new2oldIds[iold][iter->first].push_back(ielem);
                }
            }
      
          //creating the splitMeshes from the face ids
          for (int inew=0;inew<_topology->nbDomain();inew++)
            {
              if (meshesCastFrom[iold]->getNumberOfCells() > 0) //cvw
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
                  //std::cout<<"one empty mesh from "<<iold<<std::endl; //cvw
                  splitMeshes[inew][iold]=CreateEmptyMEDCouplingUMesh();
                }
            }
        }
      else
        {
          std::cout<<"proc "<<MyGlobals::_Rank<<" : castFaceMeshes empty mesh from iodDomain "<<iold<<std::endl;
        }
    }
  
  // send/receive stuff
  if (isParallelMode())
    {
      ParaMEDMEM::MEDCouplingUMesh *empty=CreateEmptyMEDCouplingUMesh();
      for (int iold=0; iold<fromSize; iold++)
        for (int inew=0; inew<newSize; inew++)
          {
            /*std::cout<<"iold inew "<<iold<<" "<<inew<<" "<<
              _domain_selector->isMyDomain(iold)<<" "<<
              _domain_selector->isMyDomain(inew)<<std::endl;*/
            if (_domain_selector->isMyDomain(iold) && !_domain_selector->isMyDomain(inew))
              if (splitMeshes[inew][iold] != 0) {
                //cvw std::cout<<"send NOT empty mesh "<<splitMeshes[inew][iold]->getName()<<" "<<inew<<"<-"<<iold<<std::endl;
                _domain_selector->sendMesh(*(splitMeshes[inew][iold]), _domain_selector->getProcessorID(inew));
              }
              else {
                //std::cout<<"send empty mesh "<<inew<<"<-"<<iold<<std::endl;
                _domain_selector->sendMesh(*(empty), _domain_selector->getProcessorID(inew));
              }
            if (!_domain_selector->isMyDomain(iold) && _domain_selector->isMyDomain(inew))
              _domain_selector->recvMesh(splitMeshes[inew][iold], _domain_selector->getProcessorID(iold));
          }
      empty->decrRef();
    }
  
  //recollecting the bits of splitMeshes to fuse them into one
  if (MyGlobals::_Verbose>300) std::cout<<"proc "<<MyGlobals::_Rank<<" : fuse splitMeshes"<<std::endl;
  meshesCastTo.resize(newSize);
  for (int inew=0; inew < newSize; inew++)
    {
      vector<const ParaMEDMEM::MEDCouplingUMesh*> myMeshes;
      for (int iold=0; iold<fromSize; iold++)
        {
          ParaMEDMEM::MEDCouplingUMesh *umesh=splitMeshes[inew][iold];
          if (umesh!=0)
            if (umesh->getNumberOfCells()>0) 
              {
                myMeshes.push_back(umesh);
              }
          //else {
          //  std::cout<<"one empty mesh "<<inew<<" "<<iold<<std::endl;
          //}
        }
      
      if (myMeshes.size()>0)
        {
          meshesCastTo[inew]=ParaMEDMEM::MEDCouplingUMesh::MergeUMeshes(myMeshes);
        }
      else
        {
          //std::cout<<"one empty meshes to merge "<<inew<<std::endl;
          //ParaMEDMEM::MEDCouplingUMesh *empty=ParaMEDMEM::MEDCouplingUMesh::New(); //empty one
          //empty->setName("emptyMesh");
          //empty->setMeshDimension(3);
          //empty->allocateCells(0);
          ParaMEDMEM::MEDCouplingUMesh *empty=CreateEmptyMEDCouplingUMesh();
          meshesCastTo[inew]=empty;
        }
      //      meshesCastTo[inew]->zipCoords();
      for (int iold=0; iold<fromSize; iold++)
        if (splitMeshes[inew][iold]!=0) splitMeshes[inew][iold]->decrRef();
    }
  //if (MyGlobals::_Verbose>1) std::cout<<"proc "<<MyGlobals::_Rank<<" : end fuse"<<std::endl;
}

void MeshCollection::remapIntField(const ParaMEDMEM::MEDCouplingUMesh& sourceMesh,
                                   const ParaMEDMEM::MEDCouplingUMesh& targetMesh,
                                   const int* fromArray,
                                   int* toArray)
{
  using std::vector;
  if (sourceMesh.getNumberOfCells()<=0) return; //empty mesh could exist
  //cvw std::cout<<"remapIntField "<<sourceMesh.getNumberOfCells()<<" "<<targetMesh.getNumberOfCells()<<std::endl;
  ParaMEDMEM::DataArrayDouble* sourceCoords=sourceMesh.getBarycenterAndOwner();
  ParaMEDMEM::DataArrayDouble* targetCoords=targetMesh.getBarycenterAndOwner();
   
  ParaMEDMEM::MEDCouplingUMesh* tmpMesh=ParaMEDMEM::MEDCouplingUMesh::New();
  tmpMesh->setCoords(sourceCoords);
  vector<int> c;
  vector<int> cI;
  tmpMesh->getNodeIdsNearPoints(targetCoords->getConstPointer(),targetMesh.getNumberOfCells(),1e-10,c,cI);
  if (cI.size()!= targetMesh.getNumberOfCells()+1) 
    throw INTERP_KERNEL::Exception(LOCALIZED("Error in source/target projection"));
  for (int itargetnode=0; itargetnode<targetMesh.getNumberOfCells();itargetnode++)    
    {
      if (cI[itargetnode]==cI[itargetnode+1]) continue;
      int isourcenode=c[cI[itargetnode]];
      toArray[itargetnode]=fromArray[isourcenode];
    } 
  sourceCoords->decrRef();
  targetCoords->decrRef();
  tmpMesh->decrRef();
}

void MeshCollection::castIntField2(std::vector<ParaMEDMEM::MEDCouplingUMesh*>& meshesCastFrom,
                                   std::vector<ParaMEDMEM::MEDCouplingUMesh*>& meshesCastTo,
                                   std::vector<ParaMEDMEM::DataArrayInt*>& arrayFrom,
                                   std::string nameArrayTo)
{
  using std::vector;
  
  int ioldMax=meshesCastFrom.size();
  int inewMax=meshesCastTo.size();
  // send-recv operations
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
              //std::cout<<"proc "<<_domain_selector->rank()<<" : castIntField SendIntVec "<<size<<std::endl;
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
              remapIntField2(inew,iold,*recvMesh,*meshesCastTo[inew],&recvIds[0],nameArrayTo);
              recvMesh->decrRef(); //cww is it correct?
            }
        }
    }
  
  //local contributions and aggregation
  for (int inew=0; inew<inewMax; inew++)    
    {
      for (int iold=0; iold<ioldMax; iold++)
        if (!isParallelMode() || ( _domain_selector->isMyDomain(iold) && _domain_selector->isMyDomain(inew)))
          {
            remapIntField2(inew,iold,*meshesCastFrom[iold],*meshesCastTo[inew],arrayFrom[iold]->getConstPointer(),nameArrayTo);
          }
    }
}

void MeshCollection::remapIntField2(int inew,int iold,
                                    const ParaMEDMEM::MEDCouplingUMesh& sourceMesh,
                                    const ParaMEDMEM::MEDCouplingUMesh& targetMesh,
                                    const int* fromArray,
                                    string nameArrayTo)
//here we store ccI for next use in future call of castAllFields and remapDoubleField2
{
  //cout<<"remapIntField2 "<<Cle2ToStr(nameArrayTo,inew,iold)<<endl;
  if (sourceMesh.getNumberOfCells()<=0) return; //empty mesh could exist
  //cvw std::cout<<"remapIntField "<<sourceMesh.getNumberOfCells()<<" "<<targetMesh.getNumberOfCells()<<std::endl;
  ParaMEDMEM::DataArrayDouble* sourceCoords=sourceMesh.getBarycenterAndOwner();
  ParaMEDMEM::DataArrayDouble* targetCoords=targetMesh.getBarycenterAndOwner();
   
  ParaMEDMEM::MEDCouplingUMesh* tmpMesh=ParaMEDMEM::MEDCouplingUMesh::New();
  tmpMesh->setCoords(sourceCoords);
  vector<int> c;
  vector<int> cI;
  vector<int> ccI; //memorize intersection target<-source(inew,iold)
  string str,cle;
  str=nameArrayTo+"_toArray";
  cle=Cle1ToStr(str,inew);
  int* toArray;
  int targetSize=targetMesh.getNumberOfCells();
  //first time iold : create and initiate 
  if (_mapDataArrayInt.find(cle)==_mapDataArrayInt.end())
    {
      if (MyGlobals::_Is0verbose>100) cout<<"create "<<cle<<" size "<<targetSize<<endl;
      ParaMEDMEM::DataArrayInt* p=DataArrayInt::New();
      p->alloc(targetSize,1);
      p->fillWithZero();
      toArray=p->getPointer();
      _mapDataArrayInt[cle]=p;
    }
  else //other times iold: refind and complete
    {
      toArray=_mapDataArrayInt.find(cle)->second->getPointer();
    }
  tmpMesh->getNodeIdsNearPoints(targetCoords->getConstPointer(),targetSize,1e-10,c,cI);
  if (cI.size()!=targetSize+1) 
    throw INTERP_KERNEL::Exception(LOCALIZED("Error in source/target projection"));
  for (int itargetnode=0; itargetnode<targetSize; itargetnode++)    
    {
      if (cI[itargetnode]==cI[itargetnode+1]) continue;
      int isourcenode=c[cI[itargetnode]];
      toArray[itargetnode]=fromArray[isourcenode];
      //memorize intersection 
      ccI.push_back(itargetnode); //next we'll do toArray[ccI[i]]=fromArray[ccI[i+1]]
      ccI.push_back(isourcenode);
    }
  //ccI.push_back(sourceMesh.getNumberOfCells()); //additionnal information at end??
  
  //memories intersection for future same job on fields (if no existing cle=no intersection)
  str=Cle2ToStr(nameArrayTo+"_ccI",inew,iold);
  if (MyGlobals::_Verbose>700) cout<<"proc "<<MyGlobals::_Rank<<" : map memorize '"<<str<<"'\n";
  _mapDataArrayInt[str]=CreateDataArrayIntFromVector(ccI, 2);
  sourceCoords->decrRef();
  targetCoords->decrRef();
  tmpMesh->decrRef();
}

void MeshCollection::castAllFields(MeshCollection& initialCollection, string nameArrayTo) //cvwat08
{
  using std::vector;
  if (nameArrayTo!="cellFieldDouble") 
    throw INTERP_KERNEL::Exception(LOCALIZED("Error castAllField only on cellFieldDouble"));

  string nameTo="typeData=6"; //resume the type of field casted 
  // send-recv operations
  int ioldMax=initialCollection.getMesh().size();
  int inewMax=this->getMesh().size();
  int iFieldMax=initialCollection.getFieldDescriptions().size();
  if (MyGlobals::_Verbose>10) cout<<"castAllFields with:\n"<<ReprVectorOfString(initialCollection.getFieldDescriptions())<<endl;
  //std::vector<std::string> initialCollection.getFieldDescriptions()
  //getFieldDescriptions() is a string description of field coherent and tested and set BEFORE.
  //see collection.prepareFieldDescriptions()
  for (int ifield=0; ifield<iFieldMax; ifield++)
    {
      string descriptionField=initialCollection.getFieldDescriptions()[ifield];
      if (descriptionField.find(nameTo)==string::npos) continue; //only nameTo accepted in Fields name description
      for (int inew=0; inew<inewMax; inew++)
        {
          for (int iold=0; iold<ioldMax; iold++)
            {
              //sending arrays for distant domains
              if (isParallelMode() && _domain_selector->isMyDomain(iold) && !_domain_selector->isMyDomain(inew))
                {
                  int target=_domain_selector->getProcessorID(inew);
                  ParaMEDMEM::DataArrayDouble* field=initialCollection.getField(descriptionField,iold); //cvwat14
                  //getField look for and read it if not done, and assume decrRef() in ~MeshCollection;
                  if (MyGlobals::_Verbose>10) 
                    std::cout<<"proc "<<_domain_selector->rank()<<" : castAllFields sendDouble"<<std::endl;
                  SendDataArrayDouble(field, target);
                }
              //receiving arrays from distant domains
              if (isParallelMode() && !_domain_selector->isMyDomain(iold) && _domain_selector->isMyDomain(inew))
                {
                  int source=_domain_selector->getProcessorID(iold);
                  //receive vector
                  if (MyGlobals::_Verbose>10) 
                    std::cout<<"proc "<<_domain_selector->rank()<<" : castAllFields recvDouble"<<std::endl;
                  ParaMEDMEM::DataArrayDouble* field=RecvDataArrayDouble(source);
                  remapDoubleField3(inew,iold,field,nameArrayTo,descriptionField);
                }
            }
        }
  
      //local contributions and aggregation
      for (int inew=0; inew<inewMax; inew++)
        {
          for (int iold=0; iold<ioldMax; iold++)
            if (!isParallelMode() || ( _domain_selector->isMyDomain(iold) && _domain_selector->isMyDomain(inew)))
              {
                ParaMEDMEM::DataArrayDouble* field=initialCollection.getField(descriptionField,iold); //cvwat14
                remapDoubleField3(inew,iold,field,nameArrayTo,descriptionField);
              }
        }
    }
}

void MeshCollection::remapDoubleField3(int inew, int iold, 
                                       ParaMEDMEM::DataArrayDouble* fromArray,
                                       string nameArrayTo,
                                       string descriptionField)
//here we use 'cellFamily_ccI inew iold' set in remapIntField2
{
  using std::vector;
  if (nameArrayTo!="cellFieldDouble") 
    throw INTERP_KERNEL::Exception(LOCALIZED("Error remapDoubleField3 only on cellFieldDouble"));
  string cle=Cle2ToStr("cellFamily_ccI",inew,iold);
  
  map<string,ParaMEDMEM::DataArrayInt*>::iterator it1;
  it1=_mapDataArrayInt.find(cle);
  if (it1==_mapDataArrayInt.end())
    {
      cerr<<"proc "<<MyGlobals::_Rank<<" : remapDoubleField3 cle '"<<cle<<"' not found"<<endl;
      cerr<<" trying remap of field double on cells : "<<descriptionField<<endl;
      return;
    }
  //create ccI in remapIntField2
  ParaMEDMEM::DataArrayInt* ccI=it1->second;
  if (MyGlobals::_Verbose>300) cout<<"proc "<<MyGlobals::_Rank<<" : remapDoubleField3 "<<cle<<" size "<<ccI->getNbOfElems()<<endl;
  //cout<<descriptionField<<endl;
  
  int nbcell=this->getMesh()[inew]->getNumberOfCells(); //number of cell of mesh
  int nbcomp=fromArray->getNumberOfComponents();
  int nbPtGauss=StrToInt(ExtractFromDescription(descriptionField, "nbPtGauss="));
  
  //int nbk=fromArray->getNumberOfTuples();
  
  //cle=reprGenericDescription(descriptionField)+" "+IntToStr(inew);
  string tag="inewFieldDouble="+IntToStr(inew);
  cle=descriptionField+SerializeFromString(tag);
  //cout<<"descriptionField in remapDoubleField3 : "<<descriptionField<<endl;
  int fromArrayNbOfElem=fromArray->getNbOfElems();
  int fromArrayNbOfComp=fromArray->getNumberOfComponents();
  int fromArrayNbOfCell=fromArrayNbOfElem/fromArrayNbOfComp/nbPtGauss;
  
  if (MyGlobals::_Verbose>1000)
    {
      cout<<"proc "<<MyGlobals::_Rank<<" nbcell "<<nbcell<<" nbcomp "<<nbcomp<<" nbPtGauss "<<nbPtGauss<<
        " fromArray nbOfElems "<<fromArrayNbOfElem<<
        " nbTuples "<<fromArray->getNumberOfTuples()<<
        " nbcells "<<fromArrayNbOfCell<<
        " nbComponents "<<fromArray->getNumberOfComponents()<<endl;
    }

  ParaMEDMEM::DataArrayDouble* field=0;
  map<string,ParaMEDMEM::DataArrayDouble*>::iterator it2;
  it2=_mapDataArrayDouble.find(cle);
  if (it2==_mapDataArrayDouble.end())
    {
      if (MyGlobals::_Verbose>300) cout<<"proc "<<MyGlobals::_Rank<<" : remapDoubleField3 cle '"<<cle<<"' not found and create it"<<endl;
      field=DataArrayDouble::New();
      _mapDataArrayDouble[cle]=field;
      field->alloc(nbcell*nbPtGauss,nbcomp);
      field->fillWithZero();
    }
  else
    {
      field=it2->second;
      if (field->getNumberOfTuples()!=nbcell*nbPtGauss || field->getNumberOfComponents()!=nbcomp)
        {
          cerr<<"proc "<<MyGlobals::_Rank<<" : remapDoubleField3 pb of size "<<
            " trying remap of field double on cells : \n"<<descriptionField<<endl;
          return;
        }
    }
  //cout<<"proc "<<MyGlobals::_Rank<<" : remapDoubleField3 "<<cle<<" size "<<ccI->getNbOfElems()<<endl;
  
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
            throw INTERP_KERNEL::Exception(LOCALIZED("Error field override"));
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
MeshCollection::MeshCollection(const std::string& filename)
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
  // char filenamechar[256];
  // strcpy(filenamechar,filename.c_str());
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
          if ( _driver ) delete _driver; _driver=0;
          throw INTERP_KERNEL::Exception(LOCALIZED("file does not comply with any recognized format"));
        }
    }
  for ( int idomain = 0; idomain < _mesh.size(); ++idomain )
    if ( _mesh[idomain] && _mesh[idomain]->getNumberOfNodes() > 0 )
      _i_non_empty_mesh = idomain;
}

/*! Constructing the MESH collection from selected domains of a distributed file
 * 
 * \param filename  - name of the master file containing the list of all the MED files
 * \param domainSelector - selector of domains to load
 */
MeshCollection::MeshCollection(const std::string& filename, ParaDomainSelector& domainSelector) //cvwat01
  : _topology(0),
    _owns_topology(true),
    _driver(0),
    //cvw _domain_selector( domainSelector.nbProcs() > 1 ? & domainSelector : 0 ),
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
          _driver=new MeshCollectionMedXmlDriver(this); //cvwat02
          _driver->read ( (char*)filename.c_str(), _domain_selector );
          _driver_type = MedXml;
        }
      catch(...)
        {  // Handle all exceptions
          if ( _driver ) delete _driver; _driver=0;
          throw INTERP_KERNEL::Exception(LOCALIZED("file .xml does not comply with any recognized format"));
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
          //std::cout<<xml<<std::endl;
          std::string nameFileXml=myfile;
          nameFileXml.replace(nameFileXml.find(".med"),4,".xml");
          nameFileXml="medpartitioner_"+nameFileXml;
          if (_domain_selector->rank()==0) //only on to write it
            {
              std::ofstream f(nameFileXml.c_str());
              f<<xml;
              f.close();
            }
#ifdef HAVE_MPI2
          MPI_Barrier(MPI_COMM_WORLD); //wait for creation of nameFileXml
#endif
          try
            {
              _driver=new MeshCollectionMedXmlDriver(this);
              _driver->read ( (char*)nameFileXml.c_str(), _domain_selector );
              _driver_type = MedXml;
            }
          catch(...)
            {  // Handle all exceptions
              if ( _driver ) delete _driver; _driver=0;
              throw INTERP_KERNEL::Exception(LOCALIZED("file medpartitioner_xxx.xml does not comply with any recognized format"));
            }
        }
      else //no extension
        {
          try
            {
              _driver=new MeshCollectionMedAsciiDriver(this);
              _driver->read ( (char*)filename.c_str(), _domain_selector );
              _driver_type=MedAscii;
            }
          catch(...)
            {
              if ( _driver ) delete _driver; _driver=0;
              throw INTERP_KERNEL::Exception(LOCALIZED("file name with no extension does not comply with any recognized format"));
            }
        }
    }
  
  /*done in MeshCollectionMedXmlDriver read
    if ( isParallelMode() )
    // to know nb of cells on each proc to compute global cell ids from locally global
    _domain_selector->gatherNbOf( getMesh() );*/

  // find non-empty domain mesh
  for ( int idomain = 0; idomain < _mesh.size(); ++idomain )
    if ( _mesh[idomain] && _mesh[idomain]->getNumberOfNodes() > 0 )
      _i_non_empty_mesh = idomain;
  
  try
    {
      //check for all proc/file compatibility of _fieldDescriptions
      //*MyGlobals::_File_Names=AllgathervVectorOfString(*MyGlobals::_File_Names);
      _fieldDescriptions=AllgathervVectorOfString(MyGlobals::_Field_Descriptions); //cvwat07
    }
  catch(INTERP_KERNEL::Exception& e)
    {
      cerr<<"proc "<<MyGlobals::_Rank<<" : INTERP_KERNEL_Exception : "<<e.what()<<endl;
      throw INTERP_KERNEL::Exception(LOCALIZED("Something wrong verifying coherency of files med ands fields"));
    }
  
  try
    {
      //check for all proc/file compatibility of _familyInfo //cvwat05
      vector<string> v2=AllgathervVectorOfString(VectorizeFromMapOfStringInt(_familyInfo));
      _familyInfo=DevectorizeToMapOfStringInt(v2);
    }
  catch(INTERP_KERNEL::Exception& e)
    {
      cerr<<"proc "<<MyGlobals::_Rank<<" : INTERP_KERNEL_Exception : "<<e.what()<<endl;
      throw INTERP_KERNEL::Exception(LOCALIZED("Something wrong merging all familyInfo"));
    }

  try
    {
      //check for all proc/file compatibility of _groupInfo
      vector<string> v2=AllgathervVectorOfString(
                                                 VectorizeFromMapOfStringVectorOfString(_groupInfo));
      _groupInfo=DeleteDuplicatesInMapOfStringVectorOfString(
                                                             DevectorizeToMapOfStringVectorOfString(v2));
    }
  catch(INTERP_KERNEL::Exception& e)
    {
      cerr<<"proc "<<MyGlobals::_Rank<<" : INTERP_KERNEL_Exception : "<<e.what()<<endl;
      throw INTERP_KERNEL::Exception(LOCALIZED("Something wrong merging all groupInfo"));
    }

  
  //std::vector< std::string > _meshes=MEDLoader::GetMeshNames(filename);
  //std::vector< std::string > _fields=MEDLoader::GetAllFieldNamesOnMesh(filename,meshname[0]);
  //cout<<"number of fields "<<_fields.size()<<endl;
  
}

/*! constructing the MESH collection from a sequential MED-file
 * 
 * \param filename MED file
 * \param meshname name of the mesh that is to be read
 */
MeshCollection::MeshCollection(const std::string& filename, const std::string& meshname)
  : _name(meshname),
    _topology(0),
    _owns_topology(true),
    _driver(0),
    _domain_selector( 0 ),
    _i_non_empty_mesh(-1),
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
      if ( _driver ) delete _driver; _driver=0;
      throw INTERP_KERNEL::Exception(LOCALIZED("problem reading .med files"));
    }
  if ( _mesh[0] && _mesh[0]->getNumberOfNodes() > 0 )
    _i_non_empty_mesh = 0;
}

MeshCollection::~MeshCollection()
{
  for (int i=0; i<_mesh.size();i++)
    if (_mesh[i]!=0) _mesh[i]->decrRef(); 
  
  for (int i=0; i<_cellFamilyIds.size();i++)
    if (_cellFamilyIds[i]!=0) _cellFamilyIds[i]->decrRef();
  
  for (int i=0; i<_faceMesh.size();i++)
    if (_faceMesh[i]!=0) _faceMesh[i]->decrRef();
  
  for (int i=0; i<_faceFamilyIds.size();i++)
    if (_faceFamilyIds[i]!=0) _faceFamilyIds[i]->decrRef();
  
  for (map<string, ParaMEDMEM::DataArrayInt*>::iterator it=_mapDataArrayInt.begin() ; it!=_mapDataArrayInt.end(); it++ )
    if ((*it).second!=0) (*it).second->decrRef();
  
  for (map<string, ParaMEDMEM::DataArrayDouble*>::iterator it=_mapDataArrayDouble.begin() ; it!=_mapDataArrayDouble.end(); it++ )
    if ((*it).second!=0) (*it).second->decrRef();
  
  if (_driver !=0) {delete _driver; _driver=0;}
  if (_topology!=0 && _owns_topology) {delete _topology; _topology=0;}
  
  if (_joint_finder!=0) {delete _joint_finder; _joint_finder=0;}
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
void MeshCollection::write(const std::string& filename)
{
  //building the connect zones necessary for writing joints
  //   if (_topology->nbDomain()>1)
  //     buildConnectZones();
  //suppresses link with driver so that it can be changed for writing
  if (_driver!=0) delete _driver;
  _driver=0;

  //char filenamechar[256];
  //  strcpy(filenamechar,filename.c_str());
  retrieveDriver()->write (filename.c_str(), _domain_selector);
}

/*! creates or gets the link to the collection driver
 */
MeshCollectionDriver* MeshCollection::retrieveDriver()
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
          throw INTERP_KERNEL::Exception(LOCALIZED("Unrecognized driver"));
        }
    }
  return _driver;
}


/*! gets an existing driver
 * 
 */
MeshCollectionDriver* MeshCollection::getDriver() const {
  return _driver;
}

// /*! retrieves the mesh dimension*/
int MeshCollection::getMeshDimension() const {
  return _i_non_empty_mesh < 0 ? -1 : _mesh[_i_non_empty_mesh]->getMeshDimension();
}

std::vector<ParaMEDMEM::MEDCouplingUMesh*>& MeshCollection::getMesh() {
  return _mesh;
}

std::vector<ParaMEDMEM::MEDCouplingUMesh*>& MeshCollection::getFaceMesh() {
  return _faceMesh;
}

ParaMEDMEM::MEDCouplingUMesh* MeshCollection::getMesh(int idomain) const {
  return _mesh[idomain];
}

ParaMEDMEM::MEDCouplingUMesh* MeshCollection::getFaceMesh(int idomain) {
  return _faceMesh[idomain];
}

std::vector<MEDPARTITIONER::ConnectZone*>& MeshCollection::getCZ() {
  return _connect_zones;
}

Topology* MeshCollection::getTopology() const {
  return _topology;
}

void MeshCollection::setTopology(Topology* topo) {
  if (_topology!=0)
    {
      throw INTERP_KERNEL::Exception(LOCALIZED("topology is already set"));
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
void MeshCollection::buildCellGraph(MEDPARTITIONER::SkyLineArray* & array, int *& edgeweights ) //cvwat09
{
  using std::multimap;
  using std::vector;
  using std::make_pair;
  using std::pair;
  
  multimap< int, int > node2cell;
  multimap< int, int > cell2cell;
  multimap< int, int > cell2node;

  vector<vector<multimap<int,int> > > commonDistantNodes;
  int nbdomain=_topology->nbDomain();
  if (isParallelMode())
    {
      _joint_finder=new JointFinder(*this);
      _joint_finder->findCommonDistantNodes();
      commonDistantNodes=_joint_finder->getDistantNodeCell();
    }

  if (MyGlobals::_Verbose>500) _joint_finder->print();
  
  //looking for reverse nodal connectivity i global numbering
  for (int idomain=0; idomain<nbdomain; idomain++)
    {
      if (isParallelMode() && !_domain_selector->isMyDomain(idomain)) continue;
    
      /*obsolete
        int offsetCell=0, offsetNode=0;
        if (isParallelMode())
        {
        offsetCell=_domain_selector->getDomainCellShift(idomain);
        offsetNode=_domain_selector->getDomainNodeShift(idomain);
        }*/
    
      ParaMEDMEM::DataArrayInt* index=ParaMEDMEM::DataArrayInt::New();
      ParaMEDMEM::DataArrayInt* revConn=ParaMEDMEM::DataArrayInt::New();
      int nbNodes=_mesh[idomain]->getNumberOfNodes();
      //cout<<"proc "<<MyGlobals::_Rank<<" idomain "<<idomain<<" nbNodes "<<nbNodes<<" offsetCell "<<offsetCell<<" offsetNode "<<offsetNode<<endl;
      _mesh[idomain]->getReverseNodalConnectivity(revConn,index);
      //problem saturation over 1 000 000 nodes for 1 proc
      if (MyGlobals::_Verbose>100) cout<<"proc "<<MyGlobals::_Rank<<" getReverseNodalConnectivity done on "<<nbNodes<<" nodes"<<endl;
      int* index_ptr=index->getPointer();
      int* revConnPtr=revConn->getPointer();
      //if (MyGlobals::_Verbose>100) cout<<"proc "<<MyGlobals::_Rank<<" create node2cell on local nodes with global numerotation idomain|inode|icell\n";
      for (int i=0; i<nbNodes; i++)
        {
          for (int icell=index_ptr[i]; icell<index_ptr[i+1]; icell++)
            {
              /*cvw local
                node2cell.insert(make_pair(i, revConnPtr[icell]));
                cout<<" "<<idomain<<"|"<<i<<"|"<<revConnPtr[icell];
                cell2node.insert(make_pair(revConnPtr[icell], i));
              */
              int globalNode=_topology->convertNodeToGlobal(idomain,i);
              int globalCell=_topology->convertCellToGlobal(idomain,revConnPtr[icell]);
              node2cell.insert(make_pair(globalNode, globalCell));
              //cvw cout<<" "<<idomain<<"|"<<i<<"#"<< globalNode<<"|"<<revConnPtr[icell]<<"#"<<globalCell;
              cell2node.insert(make_pair(globalCell, globalNode));
            }
        }
      revConn->decrRef();
      index->decrRef();
      //vector<vector<multimap<int,int> > > dNC=getDistantNodeCell()
      for (int iother=0; iother<nbdomain; iother++)
        {
          std::multimap<int,int>::iterator it;
          int isource=idomain;
          int itarget=iother;
          for (it=_joint_finder->_distant_node_cell[isource][itarget].begin(); 
               it!=_joint_finder->_distant_node_cell[isource][itarget].end(); it++)
            {
              int globalNode=_topology->convertNodeToGlobal(idomain,(*it).first);
              int globalCell=(*it).second;
              node2cell.insert(make_pair(globalNode, globalCell));
              //cout<<"processor "<<MyGlobals::_Rank<<" : isource "<<isource<<" itarget "<<itarget<<
              //      " "<<(*it).first<<"~"<<globalNode<<"~"<<globalCell<<endl;
              cell2node.insert(make_pair(globalCell, globalNode));
            }
        }
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
  if (MyGlobals::_Verbose>100) cout<<"proc "<<MyGlobals::_Rank<<" creating graph arcs on nbNodes "<<_topology->nbNodes()<<endl;
  for (int inode=0; inode<_topology->nbNodes(); inode++)  //on all nodes
    {
      typedef multimap<int,int>::const_iterator MI;
      pair <MI,MI> myRange=node2cell.equal_range(inode);
      for (MI cell1=myRange.first; cell1!=myRange.second; cell1++)  //on cells with inode
        {
          pair <MI,MI> myNodes1=cell2node.equal_range(cell1->second);  //nodes of one cell
          for (MI cell2=myRange.first; cell2!=myRange.second; cell2++)  //on one of these cell
            {
              if (cell1->second!=cell2->second)  //in cells of my domain
                {
                  pair <MI,MI> myNodes2=cell2node.equal_range(cell2->second); //on nodes of this cell
                  int nbcn=0; //number of common nodes between cells: at least 3 for cells
                  for (MI it1=myNodes1.first; it1!=myNodes1.second; ++it1)
                    {
                      for (MI it2=myNodes2.first; it2!=myNodes2.second; ++it2)
                        {
                          if ((*it1).second==(*it2).second)
                            {
                              nbcn=nbcn+1 ; break;
                            }
                        }
                    }
                  if (nbcn>=3)  //cvw TODO if 2d cells set at 2
                    cell2cell.insert(make_pair(cell1->second,cell2->second));
                  //note here there is some global numerotation of cell which come from other domain (not mydomain)
                  //cout<<" "<<MyGlobals::_Rank<<"!"<<cell1->second<<"!"<<cell2->second; //cvw
                }
            }
        }
    }
  
  if (MyGlobals::_Verbose>100) cout<<"proc "<<MyGlobals::_Rank<<" create skylinearray"<<endl;
  //filling up index and value to create skylinearray structure
  vector <int> index,value;
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
          //cout<<" "<<MyGlobals::_Rank<<"$"<<icell<<"$"<<globalCell;
          for (it=ret.first; it!=ret.second; ++it)
            {
              //cout<<" "<<MyGlobals::_Rank<<"$"<<icell<<"$"<<globalCell<<"$"<<(*it).second<<endl;
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
                  //cout<<"|"<<ival;
                  size++;
                }
            }
          //cout<<" ";
          idep=index[index.size()-1]+size;
          index.push_back(idep);
        }
    }
  
  array=new MEDPARTITIONER::SkyLineArray(index,value);

  if (MyGlobals::_Verbose>100)
    {
      std::cout<<"\nproc "<<_domain_selector->rank()<<" : end MeshCollection::buildCellGraph "<<
        index.size()-1<<" "<<value.size()<<std::endl;
      if (index.size()>1)
        {
          for (int i=0; i<10; ++i) cout<<index[i]<<" ";
          cout<<"... "<<index[index.size()-1]<<endl;
          for (int i=0; i<15; ++i) cout<<value[i]<<" ";
          int ll=index[index.size()-1]-1;
          cout<<"... ("<<ll<<") "<<value[ll-1]<<" "<<value[ll]<<endl;
        }
    }
  
}


/*! Creates the partition corresponding to the cell graph and the partition number
 * 
 * \param nbdomain number of subdomains for the newly created graph
 * 
 * returns a topology based on the new graph
 */
Topology* MeshCollection::createPartition(int nbdomain, //cvwat06
                                          Graph::splitter_type split, 
                                          const std::string& options_string,
                                          int* user_edge_weights,
                                          int* user_vertices_weights)
{
  using std::cout;
  using std::endl;
  
  if (MyGlobals::_Verbose>10) cout<<"proc "<<MyGlobals::_Rank<<" : MeshCollection::createPartition : Building cell graph"<<endl;
  
  if (nbdomain <1) throw INTERP_KERNEL::Exception(LOCALIZED("Number of subdomains must be > 0"));
  MEDPARTITIONER::SkyLineArray* array=0;
  int* edgeweights=0;

  //cout<<"Building cell graph... ";
  //   if ( _domain_selector )
  //     buildCellGraphParallel(array,edgeweights);
  //   else
  buildCellGraph(array,edgeweights); //cvwat09
  //MPI_Barrier(MPI_COMM_WORLD);
  //cout<<"proc "<<MyGlobals::_Rank<<" :end barrier CellGraph done"<<endl;
  Graph* cellGraph;
  switch (split)
    {
    case Graph::METIS:
#ifdef ENABLE_METIS
      if (MyGlobals::_Verbose>10) cout<<"METISGraph"<<endl;
      cellGraph=(Graph*)(new METISGraph(array,edgeweights));
#else
      throw INTERP_KERNEL::Exception(LOCALIZED("METIS Graph is not available. Check your products, please."));
#endif
      break;
    case Graph::SCOTCH:
#ifdef ENABLE_SCOTCH
      if (MyGlobals::_Verbose>10) cout<<"SCOTCHGraph"<<endl;
      cellGraph=(Graph*)(new SCOTCHGraph(array,edgeweights));
#else
      throw INTERP_KERNEL::Exception(LOCALIZED("SCOTCH Graph is not available. Check your products, please."));
#endif
      break;
    }

  //!user-defined weights
  if (user_edge_weights!=0) 
    cellGraph->setEdgesWeights(user_edge_weights);
  if (user_vertices_weights!=0)
    cellGraph->setVerticesWeights(user_vertices_weights);

  if (MyGlobals::_Is0verbose>10) cout<<"partitioning graph on "<<nbdomain<<" domains"<<endl;
  cellGraph->partGraph(nbdomain, options_string, _domain_selector);  //cvwat10

  if (MyGlobals::_Is0verbose>10) cout<<"building new topology"<<endl;
  //cellGraph is a shared pointer 
  Topology* topology=new ParallelTopology (cellGraph, getTopology(), nbdomain, getMeshDimension());

  //cleaning
  if (edgeweights!=0) delete[] edgeweights;
  // if (array!=0) delete array;
  delete cellGraph;
  if (MyGlobals::_Verbose>11) cout<<"proc "<<MyGlobals::_Rank<<" : end MeshCollection::createPartition"<<endl;
  return topology;
}

/*! Creates a topology for a partition specified by the user
 * 
 * \param table user-specified partition (for each cell contains the domain number from 0 to n-1)
 * 
 * returns a topology based on the new partition
 */
Topology* MeshCollection::createPartition(const int* partition)
{
  using std::set;
  
  MEDPARTITIONER::SkyLineArray* array=0;
  int* edgeweights=0;

  buildCellGraph(array,edgeweights);
  Graph* cellGraph;
  set<int> domains;
  for (int i=0; i<_topology->nbCells(); i++)
    {
      domains.insert(partition[i]);
    }
  int nbdomain=domains.size();
  
  cellGraph=(Graph*)(new UserGraph(array, partition, _topology->nbCells()));
  
  //cellGraph is a shared pointer 
  Topology* topology = new ParallelTopology (cellGraph, getTopology(), nbdomain, getMeshDimension());
  
  //  if (array!=0) delete array;
  delete cellGraph;
  return topology;
}

void MeshCollection::setDomainNames(const std::string& name)
{
  for (int i=0; i<_topology->nbDomain(); i++)
    {
      std::ostringstream oss;
      oss<<name<<"_"<<i;
      if (!isParallelMode() || _domain_selector->isMyDomain(i))
        _mesh[i]->setName(oss.str().c_str());
    }
}

ParaMEDMEM::DataArrayDouble* MeshCollection::getField(std::string descriptionField, int iold)
//getField look for and read it if not done, and assume decrRef() in ~MeshCollection;
//something like MEDCouplingFieldDouble *f2=MEDLoader::ReadFieldCell(name.c_str(),f1->getMesh()->getName(),0,f1->getName(),0,1);
{
  int rank=MyGlobals::_Rank;
  string tag="ioldFieldDouble="+IntToStr(iold);
  string descriptionIold=descriptionField+SerializeFromString(tag);
  if (_mapDataArrayDouble.find(descriptionIold)!=_mapDataArrayDouble.end())
    {
      if (MyGlobals::_Verbose>300) cout<<"proc "<<rank<<" : YET READ getField : "<<descriptionIold<<endl;
      DataArrayDouble* res=_mapDataArrayDouble[descriptionIold];
      //cout<<res->reprZip()<<endl;
      return res;
    }
  if (MyGlobals::_Verbose>200) cout<<"proc "<<rank<<" : TO BE READ getField : "<<descriptionIold<<endl;
  string description, fileName, meshName, fieldName;
  int idomain, typeField, DT, IT, entity;
  idomain=iold;
  fileName=MyGlobals::_File_Names[iold];
  if (MyGlobals::_Verbose>10) 
    cout<<"proc "<<MyGlobals::_Rank<<" : in "<<fileName<<" "<<iold<<" "<<descriptionIold<<endl;
  //cout<<"\n\n"<<"ON_CELLS "<<ON_CELLS<<" ON_NODES "<<ON_NODES<<" ON_GAUSS_PT "<<ON_GAUSS_PT<<" ON_GAUSS_NE "<<ON_GAUSS_NE<<endl;;
  //typeField ON_CELLS 0 ON_NODES 1 ON_GAUSS_PT 2 ON_GAUSS_NE 3;
  //FieldDescriptionToData(descriptionField, &idomain, &fileName, &meshName, &fieldName, &typeField, &DT, &IT);
  FieldShortDescriptionToData(descriptionIold, fieldName, typeField, entity, DT, IT);
  meshName=MyGlobals::_Mesh_Names[iold];
  
  //MEDCouplingFieldDouble* f2=MEDLoader::ReadFieldCell(
  //         fileName.c_str(), meshName.c_str(), meshDimRelToMax, fieldName.c_str(), DT, IT);
  MEDCouplingFieldDouble* f2=MEDLoader::ReadField((ParaMEDMEM::TypeOfField) typeField,
                                                  fileName.c_str(), meshName.c_str(), 0, fieldName.c_str(), DT, IT);
  
  DataArrayDouble* res=f2->getArray();
  //to know names of components
  vector <string> browse=BrowseFieldDouble(f2);
  //done yet 
  //double time=f2->getTime(IT,DT);
  //browse.push_back("time="+DoubleToStr(time));
  string localFieldInformation=descriptionIold+SerializeFromVectorOfString(browse);
  if (MyGlobals::_Verbose>10) cout<<"proc "<<MyGlobals::_Rank<<" : localFieldInformation : "<<localFieldInformation<<endl;
  MyGlobals::_General_Informations.push_back(localFieldInformation);
  res->incrRef();  //free field, keep res
  f2->decrRef();
  _mapDataArrayDouble[descriptionIold]=res; 
  
  //duplicate it! because f2->decRef!!
  //DataArrayDouble* res=f2->getArray()->deepCpy();
  //f2->decrRef();
  //cout<<res->reprZip()<<endl;
  //have to put it in map for next needs.. decRef later...~MeshCollection
  return res;
}

void MeshCollection::prepareFieldDescriptions()
//to have unique valid fields names/pointers/descriptions for partitionning
//filter _fieldDescriptions to be in all procs compliant and equal
{
  int nbfiles=MyGlobals::_File_Names.size(); //nb domains
  vector<string> r2;
  //from allgatherv then vector(procs) of serialised vector(fields) of vector(description) data
  for (int i=0; i<_fieldDescriptions.size(); i++)
    {
      vector<string> r1=DeserializeToVectorOfString(_fieldDescriptions[i]);
      for (int i=0; i<r1.size(); i++) r2.push_back(r1[i]);
    }
  //here vector(procs*fields) of serialised vector(description) data
  _fieldDescriptions=r2;
  int nbfields=_fieldDescriptions.size(); //on all domains
  if ((nbfields%nbfiles)!=0)
    {
      if (MyGlobals::_Rank==0)
        {
          cerr<<"\nERROR : incoherent number of fields references in all files .med\n"<<endl
              <<"fileMedNames :"<<endl
              <<ReprVectorOfString(MyGlobals::_File_Names)
              <<"fieldDescriptions :"<<endl
              <<ReprVectorOfString(MyGlobals::_Field_Descriptions); //cvwat07
        }
      throw INTERP_KERNEL::Exception(LOCALIZED("incoherent number of fields references in all files .med\n"));
    }
  _fieldDescriptions.resize(nbfields/nbfiles);
  for (int i=0; i<_fieldDescriptions.size(); i++)
    {
      string str=_fieldDescriptions[i];
      str=EraseTagSerialized(str,"idomain=");
      str=EraseTagSerialized(str,"fileName=");
      _fieldDescriptions[i]=str;
    }
}

//returns true if inodes of a face are in inodes of a cell
bool isFaceOncell(vector< int >& inodesFace,vector< int >&  inodesCell)
{
  int ires=0;
  int nbok=inodesFace.size();
  for (int i=0; i<nbok; i++)
    {
      int ii=inodesFace[i];
      if (ii<0) cout<<"isFaceOncell problem inodeface<0"<<endl;
      for (int j=0; j<inodesCell.size(); j++)
        {
          if (ii==inodesCell[j])
            {
              ires=ires+1; break; //inode of face found
            }
        }
      if (ires<i+1) break; //inode of face not found do not continue...
    }
  return (ires==nbok);
}

void MeshCollection::filterFaceOnCell()
{
  //meshesCells=_mesh;
  //meshesFaces=_faceMesh;
  for (int inew=0; inew<_topology->nbDomain(); inew++)
    {
      if (isParallelMode() && _domain_selector->isMyDomain(inew))
        {
          if (MyGlobals::_Verbose>200) 
            std::cout<<"proc "<<MyGlobals::_Rank<<" : filterFaceOnCell on inewDomain "<<inew<<
              " nbOfFaces "<<_faceMesh[inew]->getNumberOfCells()<<endl;
          ParaMEDMEM::MEDCouplingUMesh* mcel=_mesh[inew];
          ParaMEDMEM::MEDCouplingUMesh* mfac=_faceMesh[inew];
      
          //to have cellnode=f(facenode)... inodeCell=nodeIds[inodeFace]
          vector<int> nodeIds;
          //cout<<"proc "<<MyGlobals::_Rank<<" : nodeIds beg "<<inew<<" "<<mcel<<" "<<mfac<<endl;
          getNodeIds(*mcel, *mfac, nodeIds);
          if (nodeIds.size()==0) continue;  //one empty mesh nothing to do
      
          DataArrayInt *revNodalCel=DataArrayInt::New();
          DataArrayInt *revNodalIndxCel=DataArrayInt::New();
          mcel->getReverseNodalConnectivity(revNodalCel,revNodalIndxCel);
          int *revC=revNodalCel->getPointer();
          int *revIndxC=revNodalIndxCel->getPointer();
      
          vector< int > faceOnCell;
          vector< int > faceNotOnCell;
          int nbface=mfac->getNumberOfCells();
          for (int iface=0; iface<nbface; iface++)
            {
              bool ok;
              vector< int > inodesFace;
              mfac->getNodeIdsOfCell(iface, inodesFace);
              int nbnodFace=inodesFace.size();
              //set inodesFace in mcel
              for (int i=0; i<nbnodFace; i++) inodesFace[i]=nodeIds[inodesFace[i]];
              int inod=inodesFace[0];
              if (inod<0) cout<<"filterFaceOnCell problem 1"<<endl;
              int nbcell=revIndxC[inod+1]-revIndxC[inod];
              for (int j=0; j<nbcell; j++) //look for each cell with inod
                {
                  int icel=revC[revIndxC[inod]+j];
                  vector< int > inodesCell;
                  mcel->getNodeIdsOfCell(icel, inodesCell);
                  ok=isFaceOncell(inodesFace, inodesCell);
                  if (ok) break;
                }
              if (ok)
                {
                  faceOnCell.push_back(iface);
                  //if (MyGlobals::_Is0verbose) cout<<"face on cell "<<iface<<" "<<faceOnCell.size()-1<<endl;
                }
              else
                {
                  faceNotOnCell.push_back(iface);
                  if (MyGlobals::_Is0verbose) cout<<"face NOT on cell "<<iface<<" "<<faceOnCell.size()-1<<endl;
                }
            }
      
          revNodalCel->decrRef();
          revNodalIndxCel->decrRef();
      
          string cle;
          cle=Cle1ToStr("filterFaceOnCell",inew);
          _mapDataArrayInt[cle]=CreateDataArrayIntFromVector(faceOnCell);
          cle=Cle1ToStr("filterNotFaceOnCell",inew);
          _mapDataArrayInt[cle]=CreateDataArrayIntFromVector(faceNotOnCell);
      
          /*ParaMEDMEM::DataArrayInt* index=ParaMEDMEM::DataArrayInt::New();
            ParaMEDMEM::DataArrayInt* revConn=ParaMEDMEM::DataArrayInt::New();
            _mesh[idomain]->getReverseNodalConnectivity(revConn,index);
            int* index_ptr=index->getPointer();*/
      
          /*if (MyGlobals::_Is0verbose)
            {
            cout<<"proc "<<MyGlobals::_Rank<<" : nodeIds end "<<inew<<" "<<nodeIds.size()<<endl;
            for (int i=0; i<nodeIds.size(); i++) cout<<" "<<nodeIds[i];
            cout<<endl;
            }*/

        }
    }
}

/*
void MeshCollection::buildBoundaryOnCellMeshes()
//no used... yet
{
  //cout<<"buildBoundaryOnCellMeshes"<<endl;
  //meshesCells=_mesh;
  //meshesFaces=_faceMesh;
  for (int inew=0; inew<_topology->nbDomain(); inew++)
  {
    if (isParallelMode() && _domain_selector->isMyDomain(inew))
    {
      if (MyGlobals::_Verbose>1) std::cout<<"proc "<<MyGlobals::_Rank<<" : filterFaceOnCell on "<<inew<<" "<<_faceMesh[inew]->getNumberOfCells()<<endl;
      ParaMEDMEM::MEDCouplingUMesh* mcel=_mesh[inew];
      //ParaMEDMEM::MEDCouplingUMesh& mfac=_faceMesh[inew];
      
      DataArrayInt *desc=DataArrayInt::New();
      DataArrayInt *descIndx=DataArrayInt::New();
      DataArrayInt *revDesc=DataArrayInt::New();
      DataArrayInt *revDescIndx=DataArrayInt::New();
      //
      MEDCouplingUMesh *meshDM1=mcel->buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx);
      revDesc->decrRef();
      desc->decrRef();
      descIndx->decrRef();
      int nbOfCells=meshDM1->getNumberOfCells();
      const int *revDescIndxC=revDescIndx->getConstPointer();
      std::vector<int> boundaryCells;
      for(int i=0; i<nbOfCells; i++)
        if(revDescIndxC[i+1]-revDescIndxC[i]==1)
          boundaryCells.push_back(i);
      revDescIndx->decrRef();
      bool keepCoords=true;
      MEDCouplingUMesh *ret=(MEDCouplingUMesh *)meshDM1->buildPartOfMySelf(&boundaryCells[0],&boundaryCells[0]+boundaryCells.size(),keepCoords);
      meshDM1->decrRef();
      //don't know what to do with result yet..
      //_faceMesh[inew]->decrRef();
      //_faceMesh[inew]=ret;
    }
  }
}
*/

