// Copyright (C) 2007-2024  CEA, EDF
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

#include "MEDPARTITIONER_JointFinder.hxx"

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <utility>
#include <map>
#include <vector>

#include "MCType.hxx"
#include "MCIdType.hxx"
#include "MEDPARTITIONER_MeshCollection.hxx"
#include "MEDPARTITIONER_Topology.hxx"
#include "MEDPARTITIONER_ParaDomainSelector.hxx"
#include "MEDPARTITIONER_Utils.hxx"

#include "MEDCouplingUMesh.hxx"

/*!
 * Method contributing to the distant cell graph
 */
MEDPARTITIONER::JointFinder::JointFinder(const MeshCollection& mc):_mesh_collection(mc),_domain_selector(mc.getParaDomainSelector()),_topology(mc.getTopology())
{
}

MEDPARTITIONER::JointFinder::~JointFinder()
= default;

void MEDPARTITIONER::JointFinder::findCommonDistantNodes()
{
  int const nbdomain=_topology->nbDomain();
  _distant_node_cell.resize(nbdomain);
  _node_node.resize(nbdomain);
  for (int i=0; i<nbdomain; i++) 
    {
      _distant_node_cell[i].resize(nbdomain);
      _node_node[i].resize(nbdomain);
    }
  int const nbproc=_domain_selector->nbProcs();
  std::vector<BBTreeOfDim* > bbtree(nbdomain,(BBTreeOfDim*) nullptr);
  std::vector<double* > bbxi(nbdomain,(double*) nullptr);
  std::vector<MEDCoupling::DataArrayIdType*> rev(nbdomain,(MEDCoupling::DataArrayIdType*) nullptr);
  std::vector<MEDCoupling::DataArrayIdType*> revIndx(nbdomain,(MEDCoupling::DataArrayIdType*) nullptr);
  //int meshDim=-1;
  int spaceDim=-1;

  //init rev and revIndx and bbtree for my domain (of me:proc n)
  for (int mydomain=0; mydomain<nbdomain; mydomain++)
    {
      if(!_domain_selector->isMyDomain(mydomain))
        continue;
      const MEDCoupling::MEDCouplingUMesh* myMesh=_mesh_collection.getMesh(mydomain);
      //meshDim = myMesh->getMeshDimension();
      spaceDim= myMesh->getSpaceDimension();
      rev[mydomain] = MEDCoupling::DataArrayIdType::New();
      revIndx[mydomain] = MEDCoupling::DataArrayIdType::New();
      myMesh->getReverseNodalConnectivity(rev[mydomain],revIndx[mydomain]);
      auto* bbx=new double[2*spaceDim*myMesh->getNumberOfNodes()];
      for (int i=0; i<myMesh->getNumberOfNodes()*spaceDim; i++)
        {
          const double* coords=myMesh->getCoords()->getConstPointer();
          bbx[2*i]=(coords[i])-1e-12;
          bbx[2*i+1]=bbx[2*i]+2e-12;
        }
      bbtree[mydomain]=new BBTreeOfDim( spaceDim, bbx,nullptr,0,myMesh->getNumberOfNodes(),-1e-12);
      //keep bbx because need it in getIntersectingElems
      //no delete [] bbx yet
      bbxi[mydomain]=bbx;
    }

  //send my domains to other proc an receive other domains from other proc
  for (int isource=0; isource<nbdomain; isource++)
    {
      for (int itarget=0; itarget<nbdomain; itarget++)
        {
          const MEDCoupling::MEDCouplingUMesh* sourceMesh=_mesh_collection.getMesh(isource);
          if (_domain_selector->isMyDomain(isource)&&_domain_selector->isMyDomain(itarget))
            continue;
          if (_domain_selector->isMyDomain(isource))
            {
              //preparing data for treatment on target proc
              int const targetProc = _domain_selector->getProcessorID(itarget);

              std::vector<double> vec(spaceDim*sourceMesh->getNumberOfNodes());
              std::copy(sourceMesh->getCoords()->getConstPointer(),sourceMesh->getCoords()->getConstPointer()+sourceMesh->getNumberOfNodes()*spaceDim,&vec[0]);
              SendDoubleVec(vec,targetProc);

              //retrieving target data for storage in commonDistantNodes array
              std::vector<mcIdType> localCorrespondency;
              RecvIntVec(localCorrespondency, targetProc);
              for (std::size_t i=0; i<localCorrespondency.size()/2; i++)
                {
                  _distant_node_cell[isource][itarget].insert(std::make_pair(localCorrespondency[2*i],localCorrespondency[2*i+1]));
                }
    
            }

          if (_domain_selector->isMyDomain(itarget))
            {
              //receiving data from source proc
              int const sourceProc = isource%nbproc;
              std::vector<double> recvVec;
              RecvDoubleVec(recvVec,sourceProc);
              std::map<mcIdType,mcIdType> commonNodes; // (local nodes, distant nodes) list
              for (mcIdType inode=0; inode<ToIdType(recvVec.size()/spaceDim); inode++)
                {
                  auto* bbox=new double[2*spaceDim];
                  for (int i=0; i<spaceDim; i++)
                    {
                      bbox[2*i]=recvVec[inode*spaceDim+i]-1e-12;
                      bbox[2*i+1]=bbox[2*i]+2e-12;
                    }
                  std::vector<mcIdType> inodes;
                  bbtree[itarget]->getIntersectingElems(bbox,inodes);
                  delete [] bbox;
      
                  if (inodes.size()>0) 
                    {
                      commonNodes.insert(std::make_pair(inodes[0],inode));
                    }
          
                }
              std::vector<mcIdType> nodeCellCorrespondency;
              for (auto & commonNode : commonNodes)
                {
                  _node_node[itarget][isource].push_back(std::make_pair(commonNode.first, commonNode.second));//storing node pairs in a vector
                  const mcIdType* revIndxPtr=revIndx[itarget]->getConstPointer();
                  const mcIdType* revPtr=rev[itarget]->getConstPointer();
                  for (mcIdType icell=revIndxPtr[commonNode.first]; icell<revIndxPtr[commonNode.first+1]; icell++)
                    {
                      nodeCellCorrespondency.push_back(commonNode.second); //
                      mcIdType const globalCell=_topology->convertCellToGlobal(itarget,revPtr[icell]);
                      nodeCellCorrespondency.push_back(globalCell);
                    }
                }
              SendIntVec(nodeCellCorrespondency, sourceProc); //itarget proc send to other (otherLocalNode-itargetGlobalCell)
            }
        }
    }
    
  //free  rev(nbdomain) revIndx(nbdomain) bbtree(nbdomain) bbxi(nbdomain)
  for (int i=0; i<nbdomain; i++)
    {
      if (rev[i]!=nullptr)
        rev[i]->decrRef();
      if (revIndx[i]!=nullptr)
        revIndx[i]->decrRef();
      if (bbtree[i]!=nullptr)
        delete bbtree[i];
      if (bbxi[i]!=nullptr)
        delete [] bbxi[i];
    }

  if (MyGlobals::_Verbose>100) 
    std::cout << "proc " << _domain_selector->rank() << " : end JointFinder::findCommonDistantNodes" << std::endl;
}

std::vector<std::vector<std::multimap<mcIdType,mcIdType> > >& MEDPARTITIONER::JointFinder::getDistantNodeCell()
{
  return _distant_node_cell;
}

std::vector<std::vector<std::vector<std::pair<mcIdType,mcIdType> > > >& MEDPARTITIONER::JointFinder::getNodeNode()
{
  return _node_node;
}

void MEDPARTITIONER::JointFinder::print()
//it is for debug on small arrays under mpi 2,3 cpus
{
  int const nbdomain=_topology->nbDomain();
  //MPI_Barrier(MPI_COMM_WORLD);
  if (MyGlobals::_Is0verbose>0) 
    std::cout << "\nJointFinder print node-node (nn)iproc|itarget|isource|i|inodefirst-inodesecond\n\n" <<
      "JointFinder print distantNode=cell (nc)iproc|itarget|isource|inode=icell\n\n";
  for (int isource=0; isource<nbdomain; isource++)
    {
      for (int itarget=0; itarget<nbdomain; itarget++)
        {
          for (std::size_t i=0; i<_node_node[itarget][isource].size(); i++)
            std::cout << " nn" << _domain_selector->rank() << itarget << "|" << isource << "|" << i << "|" <<
              _node_node[itarget][isource][i].first << "-" <<
              _node_node[itarget][isource][i].second;
        }
    }
  std::cout<<std::endl;
  //MPI_Barrier(MPI_COMM_WORLD);
  for (int isource=0; isource<nbdomain; isource++)
    {
      for (int itarget=0; itarget<nbdomain; itarget++)
        {
          std::multimap<mcIdType,mcIdType>::iterator it;
          for (it=_distant_node_cell[isource][itarget].begin() ; it!=_distant_node_cell[isource][itarget].end(); it++)
            {
              std::cout << " nc" << _domain_selector->rank() << "|" << itarget << "|" << isource << "|" << (*it).first << "=" << (*it).second;
            }
        }
    }
  std::cout << std::endl;
  //MPI_Barrier(MPI_COMM_WORLD);
}
