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

#include "CommInterface.hxx" 
#include "ProcessorGroup.hxx"
#include "MPIProcessorGroup.hxx"
#include "MPIAccessDEC.hxx"
#include "MxN_Mapping.hxx"

using namespace std;

namespace MEDCoupling
{

  MxN_Mapping::MxN_Mapping(const ProcessorGroup& source_group, const ProcessorGroup& target_group,const DECOptions& dec_options)
    : DECOptions(dec_options),
      _union_group(source_group.fuse(target_group)),
      _nb_comps(0), _sending_ids(), _recv_ids()
  {
    _access_DEC = new MPIAccessDEC(source_group,target_group,getAsynchronous());
    _access_DEC->setTimeInterpolator(getTimeInterpolationMethod());
    _send_proc_offsets.resize(_union_group->size()+1,0);
    _recv_proc_offsets.resize(_union_group->size()+1,0);
  
  }

  MxN_Mapping::~MxN_Mapping()
  {
    delete _union_group;
    delete _access_DEC;
  }


  /*!
    Method registering a new element for correspondence with a distant element
    \param distant_proc proc rank of the distant processor (in terms of the union group)
    \param distant_element id of the element on the distant processor
  */
  void MxN_Mapping::addElementFromSource(int distant_proc, int distant_element)
  {
    _sending_ids.push_back(make_pair(distant_proc,distant_element));
    for (int i=distant_proc; i<_union_group->size(); i++)
      _send_proc_offsets[i+1]++;
  }

  void MxN_Mapping::initialize()
  {
    _sending_ids.clear();
    std::fill(_send_proc_offsets.begin(),_send_proc_offsets.end(),0);
  }

  void MxN_Mapping::prepareSendRecv()
  {
    CommInterface comm_interface=_union_group->getCommInterface();
    // sending count pattern
    int* nbsend=new int[_union_group->size()];
    int* nbrecv=new int[_union_group->size()];
    for (int i=0; i<_union_group->size(); i++)
      {
        nbsend[i]=_send_proc_offsets[i+1]-_send_proc_offsets[i];
      }
  
    MPIProcessorGroup* group = static_cast<MPIProcessorGroup*>(_union_group);
    const MPI_Comm* comm=group->getComm();
    comm_interface.allToAll(nbsend, 1, MPI_INT,
                            nbrecv, 1, MPI_INT,
                            *comm);
         
    for (int i=0; i<_union_group->size(); i++)
      {
        for (int j=i+1;j<_union_group->size()+1; j++)
          _recv_proc_offsets[j]+=nbrecv[i];
    
      } 

    delete[] nbsend;
    delete[] nbrecv;

    _recv_ids.resize(_recv_proc_offsets[_union_group->size()]);
    int* isendbuf=0;
    int* irecvbuf=0;
    if (_sending_ids.size()>0)
      isendbuf = new int[_sending_ids.size()];
    if (_recv_ids.size()>0)  
      irecvbuf = new int[_recv_ids.size()];
    int* sendcounts = new int[_union_group->size()];
    int* senddispls=new int[_union_group->size()];
    int* recvcounts=new int[_union_group->size()];
    int* recvdispls=new int[_union_group->size()];
    for (int i=0; i< _union_group->size(); i++)
      {
        sendcounts[i]=_send_proc_offsets[i+1]-_send_proc_offsets[i];
        senddispls[i]=_send_proc_offsets[i];
        recvcounts[i]=_recv_proc_offsets[i+1]-_recv_proc_offsets[i];
        recvdispls[i]=_recv_proc_offsets[i];
      }
    vector<int> offsets = _send_proc_offsets;
    for (int i=0; i<(int)_sending_ids.size();i++)
      {
        int iproc = _sending_ids[i].first;
        isendbuf[offsets[iproc]]=_sending_ids[i].second;
        offsets[iproc]++;
      }
    comm_interface.allToAllV(isendbuf, sendcounts, senddispls, MPI_INT,
                             irecvbuf, recvcounts, recvdispls, MPI_INT,
                             *comm);
                           
    for (int i=0; i< _recv_proc_offsets[_union_group->size()]; i++)
      _recv_ids[i]=irecvbuf[i];                           
 
    if (_sending_ids.size()>0)
      delete[] isendbuf;
    if (_recv_ids.size()>0)  
      delete[] irecvbuf;
    delete[] sendcounts;
    delete[]recvcounts;
    delete[]senddispls;
    delete[] recvdispls;
  }

  /*! Exchanging field data between two groups of processes
   * 
   * \param field MEDCoupling field containing the values to be sent
   * 
   * The ids that were defined by addElementFromSource method
   * are sent.
   */ 
  void MxN_Mapping::sendRecv(double* sendfield, MEDCouplingFieldDouble& field) const 
  {
    CommInterface comm_interface=_union_group->getCommInterface();
    const MPIProcessorGroup* group = static_cast<const MPIProcessorGroup*>(_union_group);
 
    int nbcomp=field.getArray()->getNumberOfComponents();
    double* sendbuf=0;
    double* recvbuf=0;
    if (_sending_ids.size() >0)
      sendbuf = new double[_sending_ids.size()*nbcomp];
    if (_recv_ids.size()>0)
      recvbuf = new double[_recv_ids.size()*nbcomp];
    
    int* sendcounts = new int[_union_group->size()];
    int* senddispls=new int[_union_group->size()];
    int* recvcounts=new int[_union_group->size()];
    int* recvdispls=new int[_union_group->size()];
  
    for (int i=0; i< _union_group->size(); i++)
      {
        sendcounts[i]=nbcomp*(_send_proc_offsets[i+1]-_send_proc_offsets[i]);
        senddispls[i]=nbcomp*(_send_proc_offsets[i]);
        recvcounts[i]=nbcomp*(_recv_proc_offsets[i+1]-_recv_proc_offsets[i]);
        recvdispls[i]=nbcomp*(_recv_proc_offsets[i]);
      }
    //building the buffer of the elements to be sent
    vector<int> offsets = _send_proc_offsets;

    for (int i=0; i<(int)_sending_ids.size();i++)
      { 
        int iproc = _sending_ids[i].first;
        for (int icomp=0; icomp<nbcomp; icomp++)
          sendbuf[offsets[iproc]*nbcomp+icomp]=sendfield[i*nbcomp+icomp];
        offsets[iproc]++;
      }
  
    //communication phase
    switch (getAllToAllMethod())
      {
      case Native:
        {
          const MPI_Comm* comm = group->getComm();
          comm_interface.allToAllV(sendbuf, sendcounts, senddispls, MPI_DOUBLE,
                                   recvbuf, recvcounts, recvdispls, MPI_DOUBLE,
                                   *comm);
        }
        break;
      case PointToPoint:
        _access_DEC->allToAllv(sendbuf, sendcounts, senddispls, MPI_DOUBLE,
                              recvbuf, recvcounts, recvdispls, MPI_DOUBLE);
        break;
      }
  
    //setting the received values in the field
    DataArrayDouble *fieldArr=field.getArray();
    double* recvptr=recvbuf;                         
    for (int i=0; i< _recv_proc_offsets[_union_group->size()]; i++)
      {
        for (int icomp=0; icomp<nbcomp; icomp++)
          {
            double temp = fieldArr->getIJ(_recv_ids[i],icomp);
            fieldArr->setIJ(_recv_ids[i],icomp,temp+*recvptr);
            recvptr++;
          }
      }   
    if (sendbuf!=0 && getAllToAllMethod()== Native)
      delete[] sendbuf;
    if (recvbuf !=0)
      delete[] recvbuf;
    delete[] sendcounts;
    delete[] recvcounts;
    delete[] senddispls; 
    delete[] recvdispls;
  
  }

  /*! Exchanging field data between two groups of processes
   * 
   * \param field MEDCoupling field containing the values to be sent
   * 
   * The ids that were defined by addElementFromSource method
   * are sent.
   */ 
  void MxN_Mapping::reverseSendRecv(double* recvfield, MEDCouplingFieldDouble& field) const 
  {
    CommInterface comm_interface=_union_group->getCommInterface();
    const MPIProcessorGroup* group = static_cast<const MPIProcessorGroup*>(_union_group);

    int nbcomp=field.getArray()->getNumberOfComponents();
    double* sendbuf=0;
    double* recvbuf=0;
    if (_recv_ids.size() >0)
      sendbuf = new double[_recv_ids.size()*nbcomp];
    if (_sending_ids.size()>0)
      recvbuf = new double[_sending_ids.size()*nbcomp];

    int* sendcounts = new int[_union_group->size()];
    int* senddispls=new int[_union_group->size()];
    int* recvcounts=new int[_union_group->size()];
    int* recvdispls=new int[_union_group->size()];

    for (int i=0; i< _union_group->size(); i++)
      {
        sendcounts[i]=nbcomp*(_recv_proc_offsets[i+1]-_recv_proc_offsets[i]);
        senddispls[i]=nbcomp*(_recv_proc_offsets[i]);
        recvcounts[i]=nbcomp*(_send_proc_offsets[i+1]-_send_proc_offsets[i]);
        recvdispls[i]=nbcomp*(_send_proc_offsets[i]);
      }
    //building the buffer of the elements to be sent
    vector<int> offsets = _recv_proc_offsets;
    DataArrayDouble *fieldArr=field.getArray();
    for (int iproc=0; iproc<_union_group->size();iproc++)
      for (int i=_recv_proc_offsets[iproc]; i<_recv_proc_offsets[iproc+1]; i++)
        {
          for (int icomp=0; icomp<nbcomp; icomp++)
            sendbuf[i*nbcomp+icomp]=fieldArr->getIJ(_recv_ids[i],icomp);
        }

    //communication phase
    switch (getAllToAllMethod())
      {
      case Native:
        {
          const MPI_Comm* comm = group->getComm();
          comm_interface.allToAllV(sendbuf, sendcounts, senddispls, MPI_DOUBLE,
                                   recvbuf, recvcounts, recvdispls, MPI_DOUBLE,
                                   *comm);
        }
        break;
      case PointToPoint:
        _access_DEC->allToAllv(sendbuf, sendcounts, senddispls, MPI_DOUBLE,
                               recvbuf, recvcounts, recvdispls, MPI_DOUBLE);
        break;
      }

    //setting the received values in the field
    double* recvptr=recvbuf;                         
    for (int i=0; i< _send_proc_offsets[_union_group->size()]; i++)
      {
        for (int icomp=0; icomp<nbcomp; icomp++)
          {
            recvfield[i*nbcomp+icomp]=*recvptr;
            recvptr++;
          }
      }
    if (sendbuf!=0 && getAllToAllMethod() == Native)
      delete[] sendbuf;
    if (recvbuf!=0)
      delete[] recvbuf;
    delete[] sendcounts;
    delete[] recvcounts;
    delete[] senddispls; 
    delete[] recvdispls;
  }

  ostream & operator<< (ostream & f ,const AllToAllMethod & alltoallmethod )
  {
    switch (alltoallmethod)
      {
      case Native :
        f << " Native ";
        break;
      case PointToPoint :
        f << " PointToPoint ";
        break;
      default :
        f << " UnknownAllToAllMethod ";
        break;
      }
    return f;
  }
}
