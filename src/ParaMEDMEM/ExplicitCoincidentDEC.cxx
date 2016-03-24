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
#include "ExplicitCoincidentDEC.hxx"
#include "ExplicitMapping.hxx"
#include "InterpKernelUtilities.hxx"

using namespace std;

namespace MEDCoupling
{

  /*!
    \anchor ExplicitCoincidentDEC-det
    \class ExplicitCoincidentDEC


    This class aims at \ref interpolation "remapping fields" that have identical
    supports (=the same underlying mesh) but different parallel topologies
    (=different sub-domains in the mesh). It can be used to couple
    together multi-physics codes that operate on the same domain
    with different partitioning.

    It is very similar to what the \ref StructuredCoincidentDEC-det "StructuredCoincidentDEC"
    does, except that it works with an arbitrary user-defined topology.

    The remapping between the two supports is based on identity of global
    ids, instead of geometrical considerations (as it is the case for
    \ref InterpKernelDEC-det "InterpKernelDEC").
    Therefore, beware that this \ref para-dec "DEC" can not be used
    for coincident meshes if they do *not* have the exact same numbering.

    With this \ref para-dec "DEC" no projection, and no interpolation of the field data is done, contrary
    to what happens in \ref InterpKernelDEC-det "InterpKernelDEC". It is just
    a matter of allocating the values from one side to the other, using directly the cell
    identifiers.

    As all the other DECs, its usage requires two phases :
    - a setup phase during which the topologies are exchanged so that
    the target side knows from which processors it should expect
    the data.
    - a send/recv phase during which the field data is actually transferred.

    This example illustrates the sending of a field with
    the \c ExplicitCoincidentDEC :
    \code
    ...
    ExplicitCoincidentDEC dec(groupA, groupB);
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


  /*! Constructor
   */
  ExplicitCoincidentDEC::ExplicitCoincidentDEC():
      _toposource(0),_topotarget(0),
      _targetgroup(0), _sourcegroup(0),
      _sendcounts(0), _recvcounts(0),
      _senddispls(0), _recvdispls(0),
      _recvbuffer(0), _sendbuffer(0),
      _distant_elems(), _explicit_mapping()
  {  
  }

  ExplicitCoincidentDEC::~ExplicitCoincidentDEC()
  {
  }

  /*! Synchronization process for exchanging topologies
   */
  void ExplicitCoincidentDEC::synchronize()
  {
    if (_source_group->containsMyRank())
      {
        _toposource = dynamic_cast<ExplicitTopology*>(_local_field->getTopology());
        _sourcegroup= _toposource->getProcGroup()->createProcGroup();
        _targetgroup=_toposource->getProcGroup()->createComplementProcGroup();
      }
    if (_target_group->containsMyRank())
      {
        _topotarget = dynamic_cast<ExplicitTopology*>(_local_field->getTopology());
        _sourcegroup= _topotarget->getProcGroup()->createComplementProcGroup();
        _targetgroup=_topotarget->getProcGroup()->createProcGroup();
      }
  
    // Exchanging
  
    // Transmitting source topology to target code 
    broadcastTopology(_toposource,_topotarget,1000);
    transferMappingToSource();
  }

  /*! Creates the arrays necessary for the data transfer
   * and fills the send array with the values of the 
   * source field
   *  */
  void ExplicitCoincidentDEC::prepareSourceDE()
  {
    ////////////////////////////////////
    //Step 1 : buffer array creation 
  
    if (!_toposource->getProcGroup()->containsMyRank())
      return;
    MPIProcessorGroup* group=new MPIProcessorGroup(_sourcegroup->getCommInterface());
  
    // Warning : the size of the target side is implicitly deduced
    //from the size of MPI_COMM_WORLD
    int target_size = _toposource->getProcGroup()->getCommInterface().worldSize()- _toposource->getProcGroup()->size()  ;
  
    vector<int>* target_arrays=new vector<int>[target_size];
  
    int nb_local = _toposource-> getNbLocalElements();

    int union_size=group->size();
  
    _sendcounts=new int[union_size];
    _senddispls=new int[union_size];
    _recvcounts=new int[union_size];
    _recvdispls=new int[union_size];
  
    for (int i=0; i< union_size; i++)
      {
        _sendcounts[i]=0;
        _recvcounts[i]=0;
        _recvdispls[i]=0;
      }
    _senddispls[0]=0;
 
    int* counts=_explicit_mapping.getCounts();
    for (int i=0; i<group->size(); i++)
      _sendcounts[i]=counts[i];
  
    for (int iproc=1; iproc<group->size();iproc++)
      _senddispls[iproc]=_senddispls[iproc-1]+_sendcounts[iproc-1];
  
    _sendbuffer = new double [nb_local * _toposource->getNbComponents()];
  
    /////////////////////////////////////////////////////////////
    //Step 2 : filling the buffers with the source field values 
  
    int* counter=new int [target_size];
    counter[0]=0;  
    for (int i=1; i<target_size; i++)
      counter[i]=counter[i-1]+target_arrays[i-1].size();
  
  
    const double* value = _local_field->getField()->getArray()->getPointer();
  
    int* bufferindex= _explicit_mapping.getBufferIndex();
  
    for (int ielem=0; ielem<nb_local; ielem++)
      {
        int ncomp = _toposource->getNbComponents();
        for (int icomp=0; icomp<ncomp; icomp++)
          {
            _sendbuffer[ielem*ncomp+icomp]=value[bufferindex[ielem]*ncomp+icomp];
          }  
      }
    delete[] target_arrays;
    delete[] counter;
  }

  /*!
   *  Creates the buffers for receiving the fields on the target side
   */
  void ExplicitCoincidentDEC::prepareTargetDE()
  {
    if (!_topotarget->getProcGroup()->containsMyRank())
      return;
    MPIProcessorGroup* group=new MPIProcessorGroup(_topotarget->getProcGroup()->getCommInterface());

    vector < vector <int> > source_arrays(_sourcegroup->size());
    int nb_local = _topotarget-> getNbLocalElements();
    for (int ielem=0; ielem< nb_local ; ielem++)
      {
        //pair<int,int> source_local =_distant_elems[ielem];
        pair <int,int> source_local=_explicit_mapping.getDistantNumbering(ielem);
        source_arrays[source_local.first].push_back(source_local.second); 
      }  
    int union_size=group->size();
    _recvcounts=new int[union_size];
    _recvdispls=new int[union_size];
    _sendcounts=new int[union_size];
    _senddispls=new int[union_size];
    
    for (int i=0; i< union_size; i++)
      {
        _sendcounts[i]=0;
        _recvcounts[i]=0;
        _recvdispls[i]=0;
      }
    for (int iproc=0; iproc < _sourcegroup->size(); iproc++)
      {
        //converts the rank in target to the rank in union communicator
        int unionrank=group->translateRank(_sourcegroup,iproc);
        _recvcounts[unionrank]=source_arrays[iproc].size()*_topotarget->getNbComponents();
      }
    for (int i=1; i<union_size; i++)
      _recvdispls[i]=_recvdispls[i-1]+_recvcounts[i-1];
    _recvbuffer=new double[nb_local*_topotarget->getNbComponents()];
    
  }

 
  /*!
   * Synchronizing a topology so that all the groups get it.
   * 
   * \param toposend Topology that is transmitted. It is read on processes where it already exists, and it is created and filled on others.
   * \param toporecv Topology which is received.
   * \param tag Communication tag associated with this operation.
   */
  void ExplicitCoincidentDEC::broadcastTopology(const ExplicitTopology* toposend, ExplicitTopology* toporecv, int tag)
  {
    MPI_Status status;
  
    int* serializer=0;
    int size;
  
    MPIProcessorGroup* group=new MPIProcessorGroup(*_comm_interface);
  
    // The send processors serialize the send topology
    // and send the buffers to the recv procs
    if (toposend !=0 && toposend->getProcGroup()->containsMyRank())
      {
        toposend->serialize(serializer, size);
        for (int iproc=0; iproc< group->size(); iproc++)
          {
            int itarget=iproc;
            if (!toposend->getProcGroup()->contains(itarget))
              {
                _comm_interface->send(&size,1,MPI_INT, itarget,tag+itarget,*(group->getComm()));
                _comm_interface->send(serializer, size, MPI_INT, itarget, tag+itarget,*(group->getComm()));          
              }
          }
      }
    else
      {
        vector <int> size2(group->size());
        int myworldrank=group->myRank();
        for (int iproc=0; iproc<group->size();iproc++)
          {
            int isource = iproc;
            if (!toporecv->getProcGroup()->contains(isource))
              {
                int nbelem;
                _comm_interface->recv(&nbelem, 1, MPI_INT, isource, tag+myworldrank, *(group->getComm()), &status);
                int* buffer = new int[nbelem];
                _comm_interface->recv(buffer, nbelem, MPI_INT, isource,tag+myworldrank, *(group->getComm()), &status);        
      
                ExplicitTopology* topotemp=new ExplicitTopology();
                topotemp->unserialize(buffer, *_comm_interface);
                delete[] buffer;
        
                for (int ielem=0; ielem<toporecv->getNbLocalElements(); ielem++)
                  {
                    int global = toporecv->localToGlobal(ielem);
                    int sendlocal=topotemp->globalToLocal(global);
                    if (sendlocal!=-1)
                      {
                        size2[iproc]++;
                        _explicit_mapping.pushBackElem(make_pair(iproc,sendlocal));
                      }
                  }
                delete topotemp;
              }
          }  
      }  
    MESSAGE (" rank "<<group->myRank()<< " broadcastTopology is over");
  }

  void ExplicitCoincidentDEC::transferMappingToSource()
  {

    MPIProcessorGroup* group=new MPIProcessorGroup(*_comm_interface);
  
    // sending source->target mapping which is stored by target
    //in _distant_elems from target to source
    if (_topotarget!=0 && _topotarget->getProcGroup()->containsMyRank())
      {
        int world_size = _topotarget->getProcGroup()->getCommInterface().worldSize()  ;
        int* nb_transfer_union=new int[world_size];
        int* dummy_recv=new int[world_size];
        for (int i=0; i<world_size; i++)
          nb_transfer_union[i]=0;
        //converts the rank in target to the rank in union communicator
    
        for (int i=0; i<  _explicit_mapping.nbDistantDomains(); i++)
          {
            int unionrank=group->translateRank(_sourcegroup,_explicit_mapping.getDistantDomain(i));
            nb_transfer_union[unionrank]=_explicit_mapping.getNbDistantElems(i);
          }
        _comm_interface->allToAll(nb_transfer_union, 1, MPI_INT, dummy_recv, 1, MPI_INT, MPI_COMM_WORLD);
      
        int* sendbuffer= _explicit_mapping.serialize(_topotarget->getProcGroup()->myRank());
      
        int* sendcounts= new int [world_size];
        int* senddispls = new int [world_size];
        for (int i=0; i< world_size; i++)
          {
            sendcounts[i]=2*nb_transfer_union[i];
            if (i==0)
              senddispls[i]=0;
            else
              senddispls[i]=senddispls[i-1]+sendcounts[i-1];
          }
        int* recvcounts=new int[world_size];
        int* recvdispls=new int[world_size];
        int *dummyrecv=0;
        for (int i=0; i <world_size; i++)
          {
            recvcounts[i]=0;
            recvdispls[i]=0;
          }
        _comm_interface->allToAllV(sendbuffer, sendcounts, senddispls, MPI_INT, dummyrecv, recvcounts, senddispls, MPI_INT, MPI_COMM_WORLD);
      
      }
    //receiving in the source subdomains the mapping sent by targets
    else
      {
        int world_size = _toposource->getProcGroup()->getCommInterface().worldSize()  ;
        int* nb_transfer_union=new int[world_size];
        int* dummy_send=new int[world_size];
        for (int i=0; i<world_size; i++)
          dummy_send[i]=0;
        _comm_interface->allToAll(dummy_send, 1, MPI_INT, nb_transfer_union, 1, MPI_INT, MPI_COMM_WORLD);
      
        int total_size=0;
        for (int i=0; i< world_size; i++)
          total_size+=nb_transfer_union[i];
        int nbtarget = _targetgroup->size();
        int* targetranks = new int[ nbtarget];
        for (int i=0; i<nbtarget; i++)
          targetranks[i]=group->translateRank(_targetgroup,i);
        int* mappingbuffer= new int [total_size*2];
        int* sendcounts= new int [world_size];
        int* senddispls = new int [world_size];
        int* recvcounts=new int[world_size];
        int* recvdispls=new int[world_size];
        for (int i=0; i< world_size; i++)
          {
            recvcounts[i]=2*nb_transfer_union[i];
            if (i==0)
              recvdispls[i]=0;
            else
              recvdispls[i]=recvdispls[i-1]+recvcounts[i-1];
          }

        int *dummysend=0;
        for (int i=0; i <world_size; i++)
          {
            sendcounts[i]=0;
            senddispls[i]=0;
          }
        _comm_interface->allToAllV(dummysend, sendcounts, senddispls, MPI_INT, mappingbuffer, recvcounts, recvdispls, MPI_INT, MPI_COMM_WORLD);
        _explicit_mapping.unserialize(world_size,nb_transfer_union,nbtarget, targetranks, mappingbuffer);
      }
  }

  void ExplicitCoincidentDEC::recvData()
  {
    //MPI_COMM_WORLD is used instead of group because there is no
    //mechanism for creating the union group yet
    MESSAGE("recvData");

    cout<<"start AllToAll"<<endl;
    _comm_interface->allToAllV(_sendbuffer, _sendcounts, _senddispls, MPI_DOUBLE, 
                               _recvbuffer, _recvcounts, _recvdispls, MPI_DOUBLE,MPI_COMM_WORLD);
    cout<<"end AllToAll"<<endl;
    int nb_local = _topotarget->getNbLocalElements();
    double* value=new double[nb_local*_topotarget->getNbComponents()];

    vector<int> counters(_sourcegroup->size());
    counters[0]=0;
    for (int i=0; i<_sourcegroup->size()-1; i++)
      {
        MPIProcessorGroup* group=new MPIProcessorGroup(*_comm_interface);
        int worldrank=group->translateRank(_sourcegroup,i);
        counters[i+1]=counters[i]+_recvcounts[worldrank];
      }
  
    for (int ielem=0; ielem<nb_local ; ielem++)
      {
        pair<int,int> distant_numbering=_explicit_mapping.getDistantNumbering(ielem);
        int iproc=distant_numbering.first; 
        int ncomp =  _topotarget->getNbComponents();
        for (int icomp=0; icomp< ncomp; icomp++)
          value[ielem*ncomp+icomp]=_recvbuffer[counters[iproc]*ncomp+icomp];
        counters[iproc]++;
      }  
    _local_field->getField()->getArray()->useArray(value,true,CPP_DEALLOC,nb_local,_topotarget->getNbComponents());
  }

  void ExplicitCoincidentDEC::sendData()
  {
    MESSAGE ("sendData");
    for (int i=0; i< 4; i++)
      cout << _sendcounts[i]<<" ";
    cout <<endl;
    for (int i=0; i< 4; i++)
      cout << _senddispls[i]<<" ";
    cout <<endl;
    //MPI_COMM_WORLD is used instead of group because there is no
    //mechanism for creating the union group yet
    cout <<"start AllToAll"<<endl;
    _comm_interface->allToAllV(_sendbuffer, _sendcounts, _senddispls, MPI_DOUBLE, 
                               _recvbuffer, _recvcounts, _recvdispls, MPI_DOUBLE,MPI_COMM_WORLD);
  }
}

