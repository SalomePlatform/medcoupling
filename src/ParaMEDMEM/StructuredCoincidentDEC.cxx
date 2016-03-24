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

#include <mpi.h>
#include "CommInterface.hxx"
#include "Topology.hxx"
#include "BlockTopology.hxx"
#include "ComponentTopology.hxx"
#include "ParaFIELD.hxx"
#include "MPIProcessorGroup.hxx"
#include "StructuredCoincidentDEC.hxx"
#include "InterpKernelUtilities.hxx"

#include <iostream>

using namespace std;

namespace MEDCoupling
{

  /*!
    \anchor StructuredCoincidentDEC-det
    \class StructuredCoincidentDEC

    This class aims at \ref interpolation "remapping fields" that have identical
    structured supports (=the same underlying mesh) but different parallel topologies
    (=different sub-domains in the structured mesh). It can be used to couple
    together multi-physics codes that operate on the same domain
    with different partitioning. This can be useful for example if one of
    the computation is much faster than the other. It can also be used 
    to couple together codes that share an interface that was generated
    in the same manner (with identical global ids). 
    Also, this \ref para-dec "DEC" can be used for fields that have component topologies,
    i.e., components that are scattered over several processors.

    The remapping between the two supports is based on identity of global
    ids, instead of geometrical considerations (as it is the case for
    \ref InterpKernelDEC-det "InterpKernelDEC").
    Therefore, beware that this \ref para-dec "DEC" can not be used
    for coincident meshes if they do *not* have the exact same numbering.

    With this %DEC no projection, and no interpolation of the field data is done, contrary
    to what happens in \ref InterpKernelDEC-det "InterpKernelDEC". It is just
    a matter of allocating the values from one side to the other, using directly the cell
    identifiers.

    As all the other DECs, its usage requires two phases :
    - a setup phase during which the topologies are exchanged so that
    the target side knows from which processors it should expect 
    the data.
    - a send/recv phase during which the field data is actually transferred.

    This example illustrates the sending of a field with 
    the \c StructuredCoincidentDEC : 
    \code
    ...
    StructuredCoincidentDEC dec(groupA, groupB);
    dec.attachLocalField(field);
    dec.synchronize();
    if (groupA.containsMyRank())
    dec.recvData();
    else if (groupB.containsMyRank())
    dec.sendData();
    ...
    \endcode

    Creating a ParaFIELD to be attached to the %DEC is done in exactly the same way as for
    the other DECs, if only the partitioning of the support mesh differs.
    In the case where the
    fields have also different *component* topologies, creating the ParaFIELD
    requires some more effort. See the \ref para-over "parallelism" section for more details.
  */


  StructuredCoincidentDEC::StructuredCoincidentDEC():_topo_source(0),_topo_target(0),
                                                     _send_counts(0),_recv_counts(0),
                                                     _send_displs(0),_recv_displs(0),
                                                     _recv_buffer(0),_send_buffer(0)
  {  
  }


  StructuredCoincidentDEC::~StructuredCoincidentDEC()
  {
    delete [] _send_buffer;
    delete [] _recv_buffer;
    delete []_send_displs;
    delete [] _recv_displs;
    delete [] _send_counts;
    delete [] _recv_counts;
    if (! _source_group->containsMyRank())
      delete _topo_source;
    if(!_target_group->containsMyRank())
      delete _topo_target;
  }

  StructuredCoincidentDEC::StructuredCoincidentDEC(ProcessorGroup& local_group, ProcessorGroup& distant_group):
      DisjointDEC(local_group,distant_group),
      _topo_source(0),_topo_target(0),
      _send_counts(0),_recv_counts(0),
      _send_displs(0),_recv_displs(0),
      _recv_buffer(0),_send_buffer(0)
  {
  }

  /*! Synchronization process for exchanging topologies
   */
  void StructuredCoincidentDEC::synchronizeTopology()
  {
    if (_source_group->containsMyRank())
      _topo_source = dynamic_cast<BlockTopology*>(_local_field->getTopology());
    if (_target_group->containsMyRank())
      _topo_target = dynamic_cast<BlockTopology*>(_local_field->getTopology());
  
    // Transmitting source topology to target code 
    broadcastTopology(_topo_source,1000);
    // Transmitting target topology to source code
    broadcastTopology(_topo_target,2000);
    if (_topo_source->getNbElements() != _topo_target->getNbElements())
      throw INTERP_KERNEL::Exception("Incompatible dimensions for target and source topologies");

  }

  /*! Creates the arrays necessary for the data transfer
   * and fills the send array with the values of the 
   * source field
   *  */
  void StructuredCoincidentDEC::prepareSourceDE()
  {
    ////////////////////////////////////
    //Step 1 : _buffer array creation 
  
    if (!_topo_source->getProcGroup()->containsMyRank())
      return;
    MPIProcessorGroup* group=new MPIProcessorGroup(_topo_source->getProcGroup()->getCommInterface());
  
    int myranksource = _topo_source->getProcGroup()->myRank();
  
    vector <int>* target_arrays=new vector<int>[_topo_target->getProcGroup()->size()];
  
    //cout<<" topotarget size"<<  _topo_target->getProcGroup()->size()<<endl;
  
    int nb_local = _topo_source-> getNbLocalElements();
    for (int ielem=0; ielem< nb_local ; ielem++)
      {
        //  cout <<"source local :"<<myranksource<<","<<ielem<<endl; 
        int global = _topo_source->localToGlobal(make_pair(myranksource, ielem));
        //  cout << "global "<<global<<endl;
        pair<int,int> target_local =_topo_target->globalToLocal(global);
        //  cout << "target local : "<<target_local.first<<","<<target_local.second<<endl; 
        target_arrays[target_local.first].push_back(target_local.second); 
      }  
  
    int union_size=group->size();
  
    _send_counts=new int[union_size];
    _send_displs=new int[union_size];
    _recv_counts=new int[union_size];
    _recv_displs=new int[union_size];
     
    for (int i=0; i< union_size; i++)
      {
        _send_counts[i]=0;
        _recv_counts[i]=0;
        _recv_displs[i]=0;
      }
    _send_displs[0]=0;
  
    for (int iproc=0; iproc < _topo_target->getProcGroup()->size(); iproc++)
      {
        //converts the rank in target to the rank in union communicator
        int unionrank=group->translateRank(_topo_target->getProcGroup(),iproc);
        _send_counts[unionrank]=target_arrays[iproc].size();
      }
  
    for (int iproc=1; iproc<group->size();iproc++)
      _send_displs[iproc]=_send_displs[iproc-1]+_send_counts[iproc-1];
  
    _send_buffer = new double [nb_local ];

    /////////////////////////////////////////////////////////////
    //Step 2 : filling the _buffers with the source field values 

    int* counter=new int [_topo_target->getProcGroup()->size()];
    counter[0]=0;  
    for (int i=1; i<_topo_target->getProcGroup()->size(); i++)
      counter[i]=counter[i-1]+target_arrays[i-1].size();
    
      
    const double* value = _local_field->getField()->getArray()->getPointer();
    //cout << "Nb local " << nb_local<<endl;
    for (int ielem=0; ielem<nb_local ; ielem++)
      {
        int global = _topo_source->localToGlobal(make_pair(myranksource, ielem));
        pair<int,int> target_local =_topo_target->globalToLocal(global);
        //cout <<"global : "<< global<<" local :"<<target_local.first<<" "<<target_local.second;
        //cout <<"counter[]"<<counter[target_local.first]<<endl;
        _send_buffer[counter[target_local.first]++]=value[ielem];
    
      }
    delete[] target_arrays;
    delete[] counter;
    delete group;
  }

  /*!
   *  Creates the buffers for receiving the fields on the target side
   */
  void StructuredCoincidentDEC::prepareTargetDE()
  {
    if (!_topo_target->getProcGroup()->containsMyRank())
      return;
    MPIProcessorGroup* group=new MPIProcessorGroup(_topo_source->getProcGroup()->getCommInterface());
  
    int myranktarget = _topo_target->getProcGroup()->myRank();
  
    vector < vector <int> > source_arrays(_topo_source->getProcGroup()->size());
    int nb_local = _topo_target-> getNbLocalElements();
    for (int ielem=0; ielem< nb_local ; ielem++)
      {
        //  cout <<"TS target local :"<<myranktarget<<","<<ielem<<endl; 
        int global = _topo_target->localToGlobal(make_pair(myranktarget, ielem));
        //cout << "TS global "<<global<<endl;
        pair<int,int> source_local =_topo_source->globalToLocal(global);
        //  cout << "TS source local : "<<source_local.first<<","<<source_local.second<<endl; 
        source_arrays[source_local.first].push_back(source_local.second); 
      }  
    int union_size=group->size();
    _recv_counts=new int[union_size];
    _recv_displs=new int[union_size];
    _send_counts=new int[union_size];
    _send_displs=new int[union_size];
    
    for (int i=0; i< union_size; i++)
      {
        _send_counts[i]=0;
        _recv_counts[i]=0;
        _recv_displs[i]=0;
      }
    for (int iproc=0; iproc < _topo_source->getProcGroup()->size(); iproc++)
      {
        //converts the rank in target to the rank in union communicator
        int unionrank=group->translateRank(_topo_source->getProcGroup(),iproc);
        _recv_counts[unionrank]=source_arrays[iproc].size();
      }
    for (int i=1; i<union_size; i++)
      _recv_displs[i]=_recv_displs[i-1]+_recv_counts[i-1];
    _recv_buffer=new double[nb_local];
    
    delete group;
  }

 
  /*!
   * Synchronizing a topology so that all the 
   * group possesses it.
   * 
   * \param topo Topology that is transmitted. It is read on processes where it already exists, and it is created and filled on others.
   * \param tag Communication tag associated with this operation.
   */
  void StructuredCoincidentDEC::broadcastTopology(BlockTopology*& topo, int tag)
  {
    MPI_Status status;
  
    int* serializer=0;
    int size;
  
    MPIProcessorGroup* group=new MPIProcessorGroup(*_comm_interface);
  
    // The master proc creates a send buffer containing
    // a serialized topology
    int rank_master;
  
    if (topo!=0 && topo->getProcGroup()->myRank()==0)
      {
        MESSAGE ("Master rank");
        topo->serialize(serializer, size);
        rank_master = group->translateRank(topo->getProcGroup(),0);
        MESSAGE("Master rank world number is "<<rank_master);
        MESSAGE("World Size is "<<group->size());
        for (int i=0; i< group->size(); i++)
          {
            if (i!= rank_master)
              _comm_interface->send(&rank_master,1,MPI_INT, i,tag+i,*(group->getComm()));
          }
      }
    else
      {
        MESSAGE(" rank "<<group->myRank()<< " waiting ...");
        _comm_interface->recv(&rank_master, 1,MPI_INT, MPI_ANY_SOURCE, tag+group->myRank(), *(group->getComm()),&status);
        MESSAGE(" rank "<<group->myRank()<< "received master rank"<<rank_master);
      }
    // The topology is broadcasted to all processsors in the group
    _comm_interface->broadcast(&size, 1,MPI_INT,rank_master,*(group->getComm()));
  
    int* buffer=new int[size];
    if (topo!=0 && topo->getProcGroup()->myRank()==0)
      copy(serializer, serializer+size, buffer); 
    _comm_interface->broadcast(buffer,size,MPI_INT,rank_master,*(group->getComm()));
  
    // Processors which did not possess the source topology 
    // unserialize it
  
    BlockTopology* topotemp=new BlockTopology();
    topotemp->unserialize(buffer, *_comm_interface);
  
    if (topo==0) 
      topo=topotemp;
    else 
      delete topotemp;
  
    // Memory cleaning
    delete[] buffer;
    if (serializer!=0)
      delete[] serializer;
    MESSAGE (" rank "<<group->myRank()<< " unserialize is over");
    delete group;
  }



  void StructuredCoincidentDEC::recvData()
  {
    //MPI_COMM_WORLD is used instead of group because there is no
    //mechanism for creating the union group yet
    MESSAGE("recvData");
    for (int i=0; i< 4; i++)
      cout << _recv_counts[i]<<" ";
    cout <<endl;
    for (int i=0; i< 4; i++)
      cout << _recv_displs[i]<<" ";
    cout <<endl;
  
    cout<<"start AllToAll"<<endl;
    MPI_Comm comm = *(dynamic_cast<MPIProcessorGroup*>(_union_group)->getComm());
    _comm_interface->allToAllV(_send_buffer, _send_counts, _send_displs, MPI_DOUBLE, 
                               _recv_buffer, _recv_counts, _recv_displs, MPI_DOUBLE,comm);
    cout<<"end AllToAll"<<endl;

    int nb_local = _topo_target->getNbLocalElements();
    //double* value=new double[nb_local];
    double* value=const_cast<double*>(_local_field->getField()->getArray()->getPointer());
  
    int myranktarget=_topo_target->getProcGroup()->myRank();
    vector<int> counters(_topo_source->getProcGroup()->size());
    counters[0]=0;
    for (int i=0; i<_topo_source->getProcGroup()->size()-1; i++)
      {
        MPIProcessorGroup* group=new MPIProcessorGroup(*_comm_interface);
        int worldrank=group->translateRank(_topo_source->getProcGroup(),i);
        counters[i+1]=counters[i]+_recv_counts[worldrank];
        delete group;
      }
  
    for (int ielem=0; ielem<nb_local ; ielem++)
      {
        int global = _topo_target->localToGlobal(make_pair(myranktarget, ielem));
        pair<int,int> source_local =_topo_source->globalToLocal(global);
        value[ielem]=_recv_buffer[counters[source_local.first]++];
      }
  
  
    //_local_field->getField()->setValue(value);
  }

  void StructuredCoincidentDEC::sendData()
  {
    MESSAGE ("sendData");
    for (int i=0; i< 4; i++)
      cout << _send_counts[i]<<" ";
    cout <<endl;
    for (int i=0; i< 4; i++)
      cout << _send_displs[i]<<" ";
    cout <<endl;
    cout <<"start AllToAll"<<endl;
    MPI_Comm comm = *(dynamic_cast<MPIProcessorGroup*>(_union_group)->getComm());
    _comm_interface->allToAllV(_send_buffer, _send_counts, _send_displs, MPI_DOUBLE, 
                               _recv_buffer, _recv_counts, _recv_displs, MPI_DOUBLE,comm);
    cout<<"end AllToAll"<<endl;
  }

  /*! Prepares a DEC for data exchange

    This method broadcasts the topologies from source to target 
    so that the target side can analyse from which processors it
    is expected to receive data. 
  */
  
  void StructuredCoincidentDEC::synchronize()
  {
    if (_source_group->containsMyRank())
      {
        synchronizeTopology();
        prepareSourceDE();
      }
    else if (_target_group->containsMyRank())
      {
        synchronizeTopology();
        prepareTargetDE();
      }
  }
}

