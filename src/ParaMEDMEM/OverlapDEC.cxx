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

#include "OverlapDEC.hxx"
#include "CommInterface.hxx"
#include "ParaFIELD.hxx"
#include "MPIProcessorGroup.hxx"
#include "OverlapElementLocator.hxx"
#include "OverlapInterpolationMatrix.hxx"
/*!
    \defgroup overlapdec OverlapDEC
    The OverlapDEC enables the \ref InterpKerRemapGlobal "conservative remapping" of fields between two parallel codes. This remapping is based on the computation of intersection volumes on a \b same \b processor \b group. On this processor group is defined a field template defined on A and an another on B. The computation is possible for 3D meshes, 2D meshes, 3D-surface meshes, 1D meshes and 2D-curve meshes. Dimensions must be similar for code A and code B.
    The main difference with \ref interpkerneldec is that this DEC manage 2 field-templates on each process in processor group A and B called in code source and target.
    Furthermore all process in process group cooperates in global interpolation matrix computation. In this sense InterpKernelDEC is a specialization of OverlapDEC.

    \section ParaMEDMEMOverlapDECAlgorithmDescription Algorithm Description

    Let's consider the following use case that is ran in ParaMEDMEMTest_OverlapDEC.cxx to describes the different steps in computation. Processor group contains 3 processors.
    \anchor ParaMEDMEMOverlapDECImgTest1
    \image html OverlapDEC1.png "Example showing the use case to explain the different steps."

    \subsection ParaMEDMEMOverlapDECAlgoStep1 Step 1 : Bounding box exchange and global interaction between procs computation.

    In order to reduce as fas as possible number of exchange among processors every procs computes a bounding box for A and B. Then a AllToAll is performed so that
    every procs can compute the \b global interactions between procs.
    This computation leads every procs to compute the same global TODO list expressed as a list of pair. A pair (x,y) means that proc \b x fieldtemplate A can interacts (because bounding boxes interacts) with fieltemplate B of proc \b y.
    In the \ref ParaMEDMEMOverlapDECImgTest1 "example above" this computation leads to the following a \b global TODO list :

    \b (0,0),(0,1),(1,0),(1,2),(2,0),(2,1),(2,2)

    Here pair (0,2) does not appear because bounding box of fieldtemplateA of proc#2 does not interact with proc#0 of fieldtemplate B.

    Stage performed by ParaMEDMEM::OverlapElementLocator::computeBoundingBoxes.

    \subsection ParaMEDMEMOverlapDECAlgoStep2 Step 2 : Computation of local TODO list

    Starting from global interaction previously computed in \ref ParaMEDMEMOverlapDECAlgoStep1 "Step 1", each proc computes the TODO list per proc.
    The following rules is chosen : a pair (x,y) can be treated by either proc #x or proc #y, in order to reduce the amount of data transfert accross
    procs. The algorithm chosen for loadbalancing is the following : Each proc as an empty \b local TODO list at the beginning. The for each pair (k,m) in
    \b global TODO list if proc#k has less temporary local list than proc#m pair (k,m) is added to temparary local TODO list of proc#k.
    If proc#m has less temporary local TODO list than proc#k pair (k,m) is added to temparary local TODO list of proc#m.
    If proc#k and proc#m have the same amount of temporary local TODO list pair (k,m) is added to temparary local TODO list of proc#k.

    In the \ref ParaMEDMEMOverlapDECImgTest1 "example above" this computation leads to the following a local TODO list :

    - proc#0 : (0,0)
    - proc#1 : (0,1),(1,0)
    - proc#2 : (1,2),(2,0),(2,1),(2,2)

    Stage performed by ParaMEDMEM::OverlapElementLocator::computeBoundingBoxes too.
    
    The algorithm described here is not perfect for this use case. It will be enhanced soon, I hope.

    At this stage each proc knows precisely its \b local TODO list (at interpolation sense). The \b local TODO list of other procs than local
    is kept for future computation.

    \subsection ParaMEDMEMOverlapDECAlgoStep3 Step 3 : Matrix echange between procs

    At this step knows its \b local TODO list, the aim now is to exchange field-templates between procs. Each proc computes knowing TODO list per
    proc computed in \ref ParaMEDMEMOverlapDECAlgoStep2 "Step 2" the exchange TODO list :

    In the \ref ParaMEDMEMOverlapDECImgTest1 "example above" the exchange TODO list gives the following results :

    Sending TODO list per proc :

    - proc #0 : Send fieldtemplate A to Proc#1, Send fieldtemplate B to Proc#1, Send fieldtemplate B to Proc#2
    - Proc #1 : Send fieldtemplate A to Proc#2, Send fieldtemplate B to Proc#2
    - Proc #2 : No send.

    Receiving TODO list per proc :

    - proc #0 : No receiving
    - proc #1 : receiving fieldtemplate A from Proc#0,  receiving fieldtemplate B from Proc#0
    - proc #2 : receiving fieldtemplate B from Proc#0, receiving fieldtemplate A from Proc#1,  receiving fieldtemplate B from Proc#1

    To avoid as far as possible big amount of exchange between procs only relevant part of mesh. For a proc#k sending fieldtemplate A to fieldtemplate B
    of proc #m. In this case proc#k computes part of mesh A in boundingbox B of proc#m. It implies that the corresponding cellIds or nodeIds of
    corresponding part are sent to proc #m too.

    Let's consider the couple (k,m) in TODO list. This couple is treated by either k or m as seen \ref ParaMEDMEMOverlapDECAlgoStep2 "here in Step2".

    As it will be dealt in Step 6, at the end for final matrix-vector computation the result matrix of the couple (k,m) anywhere it is computed (proc #k or proc #m)
    it will be stored in \b proc#m.

    - If proc #k is in charge (performs matrix computation) of this couple (k,m) target ids (cells or nodes) of mesh in proc #m are renumbered, because proc #m has stripped 
    its target mesh to avoid big amount of data to transfer. In this case as it is finally proc #m in charge finally of the matrix, proc #k must keep preciously the
    source ids needed to be sent to proc#m. No problem will appear for matrix assembling in proc m, for source ids because no restriction done.
    Concerning source ids to be sent for matrix-vector computation, proc k will known precisely which source ids field values to send to proc #m.
    This is incarnated by OverlapMapping::keepTracksOfTargetIds in proc m.

    - If proc #m is in charge (performs matrix computation) of this couple (k,m) source ids (cells or nodes) of mesh in proc #k are renumbered, because proc #k has stripped
    its source mesh to avoid big amount of data to transfer. In this case as it is finally proc #m in charge finally of the matrix, proc #m receive the source ids
    from remote proc #k so the matrix is directly OK, no need of renumbering will be needed in \ref ParaMEDMEMOverlapDECAlgoStep5 "Step 5". But proc k must
    keep tracks of sent ids to proc m for matrix-vector computation.
    This is incarnated by OverlapMapping::keepTracksOfSourceIds in proc k.

    This step is performed in ParaMEDMEM::OverlapElementLocator::exchangeMeshes method.

    \subsection ParaMEDMEMOverlapDECAlgoStep4 Step 4 : The interpolation matrix computation

    After mesh exchange in \ref ParaMEDMEMOverlapDECAlgoStep3 "Step3" each proc has all information to perform its \b local TODO list computed in
    \ref ParaMEDMEMOverlapDECAlgoStep2 "Step2". This step is so potentially CPU costly. That's why the \b local TODO list per proc is expected to
    be as well balanced as possible.

    The interpolation is performed as \ref ParaMEDMEM::MEDCouplingRemapper "Remapper" does.

    This operation is performed by OverlapInterpolationMatrix::addContribution method.

    \subsection ParaMEDMEMOverlapDECAlgoStep5 Step 5 : Global matrix construction.
    
    After having performed the TODO list at the end of \ref ParaMEDMEMOverlapDECAlgoStep4 "Step4" it is needed to assemble the final matrix.
    
    The final aim is to have a distributed matrix \f$ M_k \f$ on each proc#k. In order to reduce data exchange during matrix product process.
    \f$ M_k \f$ is built using sizeof(Proc group) \c std::vector< \c std::map<int,double> \c >.

    For a proc#k, it is necessary to fetch info of all matrix built in \ref ParaMEDMEMOverlapDECAlgoStep4 "Step4" where the first element in pair
    is equal to k.

    After this step, the matrix repartition is the following after call ParaMEDMEM::OverlapMapping::prepare :

    - proc#0 : (0,0),(1,0),(2,0)
    - proc#1 : (0,1),(2,1)
    - proc#2 : (1,2),(2,2)

    Tuple (2,1) computed on proc 2 is stored at the end of "prepare" in proc 1. So it is an example of item 0 in \ref ParaMEDMEMOverlapDECAlgoStep2 "Step2".
    Tuple (0,1) computed on proc 1 and stored in proc 1 too. So it is an example of item 1 in \ref ParaMEDMEMOverlapDECAlgoStep2 "Step2".

    In ParaMEDMEM::OverlapMapping::_proc_ids_to_send_vector_st will contain :

    - Proc#0 : 0,1
    - Proc#1 : 0,2
    - Proc#2 : 0,1,2

    In ParaMEDMEM::OverlapMapping::_proc_ids_to_recv_vector_st will contain :

    - Proc#0 : 0,1,2
    - Proc#1 : 0,2
    - Proc#2 : 1,2

    The method in charge to perform this is : ParaMEDMEM::OverlapMapping::prepare.
*/
namespace ParaMEDMEM
{
  OverlapDEC::OverlapDEC(const std::set<int>& procIds, const MPI_Comm& world_comm):_own_group(true),_interpolation_matrix(0),
                                                                                   _source_field(0),_own_source_field(false),
                                                                                   _target_field(0),_own_target_field(false)
  {
    ParaMEDMEM::CommInterface comm;
    int *ranks_world=new int[procIds.size()]; // ranks of sources and targets in world_comm
    std::copy(procIds.begin(),procIds.end(),ranks_world);
    MPI_Group group,world_group;
    comm.commGroup(world_comm,&world_group);
    comm.groupIncl(world_group,procIds.size(),ranks_world,&group);
    delete [] ranks_world;
    MPI_Comm theComm;
    comm.commCreate(world_comm,group,&theComm);
    comm.groupFree(&group);
    if(theComm==MPI_COMM_NULL)
      {
        _group=0;
        return ;
      }
    std::set<int> idsUnion;
    for(std::size_t i=0;i<procIds.size();i++)
      idsUnion.insert(i);
    _group=new MPIProcessorGroup(comm,idsUnion,theComm);
  }

  OverlapDEC::~OverlapDEC()
  {
    if(_own_group)
      delete _group;
    if(_own_source_field)
      delete _source_field;
    if(_own_target_field)
      delete _target_field;
    delete _interpolation_matrix;
  }

  void OverlapDEC::sendRecvData(bool way)
  {
    if(way)
      sendData();
    else
      recvData();
  }

  void OverlapDEC::sendData()
  {
    _interpolation_matrix->multiply();
  }

  void OverlapDEC::recvData()
  {
    throw INTERP_KERNEL::Exception("Not implemented yet !!!!");
    //_interpolation_matrix->transposeMultiply();
  }
  
  void OverlapDEC::synchronize()
  {
    if(!isInGroup())
      return ;
    delete _interpolation_matrix;
    _interpolation_matrix=new OverlapInterpolationMatrix(_source_field,_target_field,*_group,*this,*this);
    OverlapElementLocator locator(_source_field,_target_field,*_group);
    locator.copyOptions(*this);
    locator.exchangeMeshes(*_interpolation_matrix);
    std::vector< std::pair<int,int> > jobs=locator.getToDoList();
    std::string srcMeth=locator.getSourceMethod();
    std::string trgMeth=locator.getTargetMethod();
    for(std::vector< std::pair<int,int> >::const_iterator it=jobs.begin();it!=jobs.end();it++)
      {
        const MEDCouplingPointSet *src=locator.getSourceMesh((*it).first);
        const DataArrayInt *srcIds=locator.getSourceIds((*it).first);
        const MEDCouplingPointSet *trg=locator.getTargetMesh((*it).second);
        const DataArrayInt *trgIds=locator.getTargetIds((*it).second);
        _interpolation_matrix->addContribution(src,srcIds,srcMeth,(*it).first,trg,trgIds,trgMeth,(*it).second);
      }
    _interpolation_matrix->prepare(locator.getProcsInInteraction());
    _interpolation_matrix->computeDeno();
  }

  void OverlapDEC::attachSourceLocalField(ParaFIELD *field, bool ownPt)
  {
    if(!isInGroup())
      return ;
    if(_own_source_field)
      delete _source_field;
    _source_field=field;
    _own_source_field=ownPt;
  }

  void OverlapDEC::attachTargetLocalField(ParaFIELD *field, bool ownPt)
  {
    if(!isInGroup())
      return ;
    if(_own_target_field)
      delete _target_field;
    _target_field=field;
    _own_target_field=ownPt;
  }

  bool OverlapDEC::isInGroup() const
  {
    if(!_group)
      return false;
    return _group->containsMyRank();
  }
}
