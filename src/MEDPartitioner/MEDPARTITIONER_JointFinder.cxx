#include "MEDPARTITIONER_JointFinder.hxx"
#include "MEDPARTITIONER_MESHCollection.hxx"
#include "MEDPARTITIONER_Topology.hxx"
#include "MEDPARTITIONER_ParaDomainSelector.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDPARTITIONER_utils.hxx"
#include "BBTree.txx"


/*! Method contributing to the distant cell graph
 */
using namespace MEDPARTITIONER;

JointFinder::JointFinder(const MESHCollection& mc):_mesh_collection(mc), _topology(mc.getTopology()),_domain_selector(mc.getParaDomainSelector()) 
{
}

JointFinder::~JointFinder()
{
  //if (MyGlobals::_is0verbose>100) cout<<"TODO ~JointFinder"<<endl;
}

void JointFinder::findCommonDistantNodes()
{
  int nbdomain=_topology->nbDomain();
  _distant_node_cell.resize(nbdomain);
  _node_node.resize(nbdomain);
  for (int i=0; i<nbdomain; i++) 
  {
    _distant_node_cell[i].resize(nbdomain);
    _node_node[i].resize(nbdomain);
  }
  int nbproc=_domain_selector->nbProcs();
  std::vector<BBTree<3>* > bbtree(nbdomain,(BBTree<3>*) 0); 
  std::vector<ParaMEDMEM::DataArrayInt*> rev(nbdomain,(DataArrayInt*) 0);
  std::vector<ParaMEDMEM::DataArrayInt*> revIndx(nbdomain,(DataArrayInt*) 0);
  int meshDim;
  int spaceDim;
  
  //init rev and revIndx and bbtree for my domain (of me:proc n)
  for (int mydomain=0; mydomain<nbdomain; mydomain++)
  {
    if(!_domain_selector->isMyDomain(mydomain)) continue;
    const ParaMEDMEM::MEDCouplingUMesh* myMesh=_mesh_collection.getMesh(mydomain);
    meshDim=myMesh->getMeshDimension();
    spaceDim= myMesh->getSpaceDimension();
    rev[mydomain] = ParaMEDMEM::DataArrayInt::New();
    revIndx[mydomain] = ParaMEDMEM::DataArrayInt::New();
    myMesh->getReverseNodalConnectivity(rev[mydomain],revIndx[mydomain]);
    double* bbx=new double[2*spaceDim*myMesh->getNumberOfNodes()];
    for (int i=0; i<myMesh->getNumberOfNodes()*spaceDim; i++)
    {
      const double* coords=myMesh->getCoords()->getConstPointer();
      bbx[2*i]=(coords[i])-1e-12;
      bbx[2*i+1]=bbx[2*i]+2e-12;
    }
    bbtree[mydomain]=new BBTree<3> (bbx,0,0,myMesh->getNumberOfNodes(),-1e-12);
    delete[] bbx;
  }

  //send my domains to other proc an receive other domains from other proc
  for (int isource=0; isource<nbdomain; isource++)
  {
    for (int itarget=0; itarget<nbdomain; itarget++)
    {
      const ParaMEDMEM::MEDCouplingUMesh* sourceMesh=_mesh_collection.getMesh(isource);
      if (_domain_selector->isMyDomain(isource)&&_domain_selector->isMyDomain(itarget)) continue;
      if (_domain_selector->isMyDomain(isource))
      {
        //preparing data for treatment on target proc
        int targetProc = _domain_selector->getProcessorID(itarget);
    
        std::vector<double> vec(spaceDim*sourceMesh->getNumberOfNodes());
        //cvw cout<<"\nproc "<<_domain_selector->rank()<<" : numberOfNodes "<<sourceMesh->getNumberOfNodes()<<endl;
        std::copy(sourceMesh->getCoords()->getConstPointer(),sourceMesh->getCoords()->getConstPointer()+sourceMesh->getNumberOfNodes()*spaceDim,&vec[0]);
        sendDoubleVec(vec,targetProc);
    
        //retrieving target data for storage in commonDistantNodes array
        std::vector<int> localCorrespondency;
        recvIntVec(localCorrespondency, targetProc);
        //cvw cout<<"\nproc "<<_domain_selector->rank()<<" : nodeCellCorrespondency ";
        for (int i=0; i<localCorrespondency.size()/2; i++)
        {
          _distant_node_cell[isource][itarget].insert(std::make_pair(localCorrespondency[2*i],localCorrespondency[2*i+1]));
          //cvw cout<<" "<<localCorrespondency[2*i]<<"/"<<localCorrespondency[2*i+1];
        }
    
      }
      if (_domain_selector->isMyDomain(itarget))
      {
        //receiving data from source proc
        int sourceProc = isource%nbproc;
        std::vector<double> recvVec;
        recvDoubleVec(recvVec,sourceProc);
        std::map<int,int> commonNodes; // (local nodes, distant nodes) list
        //cvw cout<<"\nproc "<<_domain_selector->rank()<<" : commonNodes ";
        for (int inode=0; inode<(recvVec.size()/meshDim); inode++)
        {
          double* bbox=new double[2*spaceDim];
          for (int i=0; i<spaceDim; i++)
          {
            bbox[2*i]=recvVec[inode*spaceDim+i]-1e-12;
            bbox[2*i+1]=bbox[2*i]+2e-12;
          }
          std::vector<int> inodes;
          bbtree[itarget]->getIntersectingElems(bbox,inodes);
          delete[] bbox;
      
          if (inodes.size()>0) 
          {
            commonNodes.insert(std::make_pair(inodes[0],inode));
            //cvw cout<<" "<<inodes[0]<<"/"<<inode;
          }
          
        }
        std::vector<int> nodeCellCorrespondency;
        for (std::map<int,int>::iterator iter=commonNodes.begin(); iter!=commonNodes.end(); iter++)
        {
          _node_node[itarget][isource].push_back(std::make_pair(iter->first, iter->second));//storing node pairs in a vector
          const int* revIndxPtr=revIndx[itarget]->getConstPointer();
          const int* revPtr=rev[itarget]->getConstPointer();
          for (int icell=revIndxPtr[iter->first]; icell<revIndxPtr[iter->first+1]; icell++)
          {
            nodeCellCorrespondency.push_back(iter->second); //
            int globalCell=_topology->convertCellToGlobal(itarget,revPtr[icell]);
            nodeCellCorrespondency.push_back(globalCell);
            //nodeCellCorrespondency.push_back(revPtr[icell]); //need to set at global numerotation 
            //cout<<"processor "<<MyGlobals::_rank<<" : isource "<<isource<<" itarget "<<itarget<<
            //      " node "<<iter->second<<" cellLoc "<<revPtr[icell]<<" cellGlob "<<globalCell<<endl;
          }
        }
        //std::cout<<"proc "<<_domain_selector->rank()<<" : JointFinder sendIntVec "<<_domain_selector->rank()<<std::endl;  //cvwdebug
        sendIntVec(nodeCellCorrespondency, sourceProc); //itarget proc send to other (otherLocalNode-itargetGlobalCell)
        }
      }
    }
    //free  rev(nbdomain) revIndx(nbdomain) bbtree(nbdomain)
    for (int i=0; i<nbdomain; i++)
    {
      if (rev[i]!=0) rev[i]->decrRef();
      if (revIndx[i]!=0) revIndx[i]->decrRef();
      if (bbtree[i]!=0) delete bbtree[i];
    }

    if (MyGlobals::_verbose>100) 
      std::cout<<"proc "<<_domain_selector->rank()<<" : end JointFinder::findCommonDistantNodes"<<std::endl;
}

std::vector<std::vector<std::multimap<int,int> > > & JointFinder::getDistantNodeCell()
{
  return _distant_node_cell;
}

std::vector<std::vector<std::vector<std::pair<int,int> > > >& JointFinder::getNodeNode()
{
  return _node_node;
}

void JointFinder::print()
//it is for debug on small arrays under mpi 2,3 cpus
{
   int nbdomain=_topology->nbDomain();
   //MPI_Barrier(MPI_COMM_WORLD);
   if (MyGlobals::_is0verbose>0) 
     cout<<"\nJointFinder print node-node (nn)iproc|itarget|isource|i|inodefirst-inodesecond\n\n"<<
           "JointFinder print distantNode=cell (nc)iproc|itarget|isource|inode=icell\n\n";
   for (int isource=0; isource<nbdomain; isource++)
   {
     for (int itarget=0; itarget<nbdomain; itarget++)
     {
       for (int i=0; i<_node_node[itarget][isource].size(); i++)
         cout<<" nn"<<_domain_selector->rank()<<itarget<<"|"<<isource<<"|"<<i<<"|"<<
               _node_node[itarget][isource][i].first<<"-"<<
               _node_node[itarget][isource][i].second;
     }
   }
   cout<<endl;
   //MPI_Barrier(MPI_COMM_WORLD);
   
   //cout<<"proc "<<_domain_selector->rank()<<" : JointFinder _distant_node_cell itarget/isource/inode=icell"<<endl;
   for (int isource=0; isource<nbdomain; isource++)
   {
     for (int itarget=0; itarget<nbdomain; itarget++)
     {
       std::multimap<int,int>::iterator it;
       for (it=_distant_node_cell[isource][itarget].begin() ; it!=_distant_node_cell[isource][itarget].end(); it++)
       {
         cout<<" nc"<<_domain_selector->rank()<<"|"<<itarget<<"|"<<isource<<"|"<<(*it).first<<"="<<(*it).second;
       }
     }
   }
   cout<<endl;
   //MPI_Barrier(MPI_COMM_WORLD);
}
