// Copyright (C) 2007-2014  CEA/DEN, EDF R&D
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
// Author : Anthony Geay (CEA/DEN)

#include "OverlapMapping.hxx"
#include "MPIProcessorGroup.hxx"

#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingAutoRefCountObjectPtr.hxx"

#include "InterpKernelAutoPtr.hxx"

#include <numeric>
#include <algorithm>

using namespace ParaMEDMEM;

OverlapMapping::OverlapMapping(const ProcessorGroup& group):_group(group)
{
}

/*!
 * This method keeps tracks of source ids to know in step 6 of main algorithm, which tuple ids to send away.
 * This method incarnates item#1 of step2 algorithm.
 */
void OverlapMapping::keepTracksOfSourceIds(int procId, DataArrayInt *ids)
{
  ids->incrRef();
  _src_ids_st2.push_back(ids);
  _src_proc_st2.push_back(procId);
}

/*!
 * This method keeps tracks of target ids to know in step 6 of main algorithm.
 * This method incarnates item#0 of step2 algorithm.
 */
void OverlapMapping::keepTracksOfTargetIds(int procId, DataArrayInt *ids)
{
  ids->incrRef();
  _trg_ids_st2.push_back(ids);
  _trg_proc_st2.push_back(procId);
}

/*!
 * This method stores from a matrix in format Target(rows)/Source(cols) for a source procId 'srcProcId' and for a target procId 'trgProcId'.
 * All ids (source and target) are in format of local ids. 
 */
void OverlapMapping::addContributionST(const std::vector< std::map<int,double> >& matrixST, const DataArrayInt *srcIds, int srcProcId, const DataArrayInt *trgIds, int trgProcId)
{
  _matrixes_st.push_back(matrixST);
  _source_proc_id_st.push_back(srcProcId);
  _target_proc_id_st.push_back(trgProcId);
  if(srcIds)
    {//item#1 of step2 algorithm in proc m. Only to know in advanced nb of recv ids [ (0,1) computed on proc1 and Matrix-Vector on proc1 ]
      _nb_of_src_ids_proc_st2.push_back(srcIds->getNumberOfTuples());
      _src_ids_proc_st2.push_back(srcProcId);
    }
  else
    {//item#0 of step2 algorithm in proc k
      std::set<int> s;
      for(std::vector< std::map<int,double> >::const_iterator it1=matrixST.begin();it1!=matrixST.end();it1++)
        for(std::map<int,double>::const_iterator it2=(*it1).begin();it2!=(*it1).end();it2++)
          s.insert((*it2).first);
      _src_ids_zip_st2.resize(_src_ids_zip_st2.size()+1);
      _src_ids_zip_st2.back().insert(_src_ids_zip_st2.back().end(),s.begin(),s.end());
      _src_ids_zip_proc_st2.push_back(trgProcId);
    }
}

/*!
 * 'procsInInteraction' gives the global view of interaction between procs.
 * In 'procsInInteraction' for a proc with id i, is in interaction with procs listed in procsInInteraction[i].
 *
 * This method is in charge to send matrixes in AlltoAll mode.
 * After the call of this method 'this' contains the matrixST for all source elements of the current proc
 */
void OverlapMapping::prepare(const std::vector< std::vector<int> >& procsInInteraction, int nbOfTrgElems)
{
  CommInterface commInterface=_group.getCommInterface();
  const MPIProcessorGroup *group=static_cast<const MPIProcessorGroup*>(&_group);
  const MPI_Comm *comm=group->getComm();
  int grpSize=_group.size();
  INTERP_KERNEL::AutoPtr<int> nbsend=new int[grpSize];
  INTERP_KERNEL::AutoPtr<int> nbsend2=new int[grpSize];
  INTERP_KERNEL::AutoPtr<int> nbsend3=new int[grpSize];
  std::fill<int *>(nbsend,nbsend+grpSize,0);
  int myProcId=_group.myRank();
  _proc_ids_to_recv_vector_st.clear();
  int curProc=0;
  for(std::vector< std::vector<int> >::const_iterator it1=procsInInteraction.begin();it1!=procsInInteraction.end();it1++,curProc++)
    if(std::find((*it1).begin(),(*it1).end(),myProcId)!=(*it1).end())
      _proc_ids_to_recv_vector_st.push_back(curProc);
  _proc_ids_to_send_vector_st=procsInInteraction[myProcId];
  for(std::size_t i=0;i<_matrixes_st.size();i++)
    if(_source_proc_id_st[i]==myProcId)
      nbsend[_target_proc_id_st[i]]=_matrixes_st[i].size();
  INTERP_KERNEL::AutoPtr<int> nbrecv=new int[grpSize];
  commInterface.allToAll(nbsend,1,MPI_INT,nbrecv,1,MPI_INT,*comm);
  //exchanging matrix
  //first exchanging offsets+ids_source
  INTERP_KERNEL::AutoPtr<int> nbrecv1=new int[grpSize];
  INTERP_KERNEL::AutoPtr<int> nbrecv2=new int[grpSize];
  //
  int *tmp=0;
  serializeMatrixStep0ST(nbrecv,
                         tmp,nbsend2,nbsend3,
                         nbrecv1,nbrecv2);
  INTERP_KERNEL::AutoPtr<int> bigArr=tmp;
  INTERP_KERNEL::AutoPtr<int> bigArrRecv=new int[nbrecv2[grpSize-1]+nbrecv1[grpSize-1]];
  commInterface.allToAllV(bigArr,nbsend2,nbsend3,MPI_INT,
                          bigArrRecv,nbrecv1,nbrecv2,MPI_INT,
                          *comm);// sending ids of sparse matrix (n+1 elems)
  //second phase echange target ids
  std::fill<int *>(nbsend2,nbsend2+grpSize,0);
  INTERP_KERNEL::AutoPtr<int> nbrecv3=new int[grpSize];
  INTERP_KERNEL::AutoPtr<int> nbrecv4=new int[grpSize];
  double *tmp2=0;
  int lgthOfArr=serializeMatrixStep1ST(nbrecv,bigArrRecv,nbrecv1,nbrecv2,
                                       tmp,tmp2,
                                       nbsend2,nbsend3,nbrecv3,nbrecv4);
  INTERP_KERNEL::AutoPtr<int> bigArr2=tmp;
  INTERP_KERNEL::AutoPtr<double> bigArrD2=tmp2;
  INTERP_KERNEL::AutoPtr<int> bigArrRecv2=new int[lgthOfArr];
  INTERP_KERNEL::AutoPtr<double> bigArrDRecv2=new double[lgthOfArr];
  commInterface.allToAllV(bigArr2,nbsend2,nbsend3,MPI_INT,
                          bigArrRecv2,nbrecv3,nbrecv4,MPI_INT,
                          *comm);
  commInterface.allToAllV(bigArrD2,nbsend2,nbsend3,MPI_DOUBLE,
                          bigArrDRecv2,nbrecv3,nbrecv4,MPI_DOUBLE,
                          *comm);
  //finishing
  unserializationST(nbOfTrgElems,nbrecv,bigArrRecv,nbrecv1,nbrecv2,
                    bigArrRecv2,bigArrDRecv2,nbrecv3,nbrecv4);
  //updating _src_ids_zip_st2 and _src_ids_zip_st2 with received matrix.
  updateZipSourceIdsForFuture();
  //finish to fill _the_matrix_st with already in place matrix in _matrixes_st
  finishToFillFinalMatrixST();
  //printTheMatrix();
}

/*!
 * Compute denominators.
 */
void OverlapMapping::computeDenoGlobConstraint()
{
  _the_deno_st.clear();
  std::size_t sz1=_the_matrix_st.size();
  _the_deno_st.resize(sz1);
  for(std::size_t i=0;i<sz1;i++)
    {
      std::size_t sz2=_the_matrix_st[i].size();
      _the_deno_st[i].resize(sz2);
      for(std::size_t j=0;j<sz2;j++)
        {
          double sum=0;
          std::map<int,double>& mToFill=_the_deno_st[i][j];
          const std::map<int,double>& m=_the_matrix_st[i][j];
          for(std::map<int,double>::const_iterator it=m.begin();it!=m.end();it++)
            sum+=(*it).second;
          for(std::map<int,double>::const_iterator it=m.begin();it!=m.end();it++)
            mToFill[(*it).first]=sum;
        }
    }
}

/*!
 * Compute denominators.
 */
void OverlapMapping::computeDenoConservativeVolumic(int nbOfTuplesTrg)
{
  CommInterface commInterface=_group.getCommInterface();
  int myProcId=_group.myRank();
  //
  _the_deno_st.clear();
  std::size_t sz1=_the_matrix_st.size();
  _the_deno_st.resize(sz1);
  std::vector<double> deno(nbOfTuplesTrg);
  for(std::size_t i=0;i<sz1;i++)
    {
      const std::vector< std::map<int,double> >& mat=_the_matrix_st[i];
      int curSrcId=_the_matrix_st_source_proc_id[i];
      std::vector<int>::iterator isItem1=std::find(_trg_proc_st2.begin(),_trg_proc_st2.end(),curSrcId);
      int rowId=0;
      if(isItem1==_trg_proc_st2.end() || curSrcId==myProcId)//item1 of step2 main algo. Simple, because rowId of mat are directly target ids.
        {
          for(std::vector< std::map<int,double> >::const_iterator it1=mat.begin();it1!=mat.end();it1++,rowId++)
            for(std::map<int,double>::const_iterator it2=(*it1).begin();it2!=(*it1).end();it2++)
              deno[rowId]+=(*it2).second;
        }
      else
        {//item0 of step2 main algo. More complicated.
          std::vector<int>::iterator fnd=isItem1;//std::find(_trg_proc_st2.begin(),_trg_proc_st2.end(),curSrcId);
          int locId=std::distance(_trg_proc_st2.begin(),fnd);
          const DataArrayInt *trgIds=_trg_ids_st2[locId];
          const int *trgIds2=trgIds->getConstPointer();
          for(std::vector< std::map<int,double> >::const_iterator it1=mat.begin();it1!=mat.end();it1++,rowId++)
            for(std::map<int,double>::const_iterator it2=(*it1).begin();it2!=(*it1).end();it2++)
              deno[trgIds2[rowId]]+=(*it2).second;
        }
    }
  //
  for(std::size_t i=0;i<sz1;i++)
    {
      int rowId=0;
      const std::vector< std::map<int,double> >& mat=_the_matrix_st[i];
      int curSrcId=_the_matrix_st_source_proc_id[i];
      std::vector<int>::iterator isItem1=std::find(_trg_proc_st2.begin(),_trg_proc_st2.end(),curSrcId);
      std::vector< std::map<int,double> >& denoM=_the_deno_st[i];
      denoM.resize(mat.size());
      if(isItem1==_trg_proc_st2.end() || curSrcId==myProcId)//item1 of step2 main algo. Simple, because rowId of mat are directly target ids.
        {
          int rowId=0;
          for(std::vector< std::map<int,double> >::const_iterator it1=mat.begin();it1!=mat.end();it1++,rowId++)
            for(std::map<int,double>::const_iterator it2=(*it1).begin();it2!=(*it1).end();it2++)
              denoM[rowId][(*it2).first]=deno[rowId];
        }
      else
        {
          std::vector<int>::iterator fnd=isItem1;
          int locId=std::distance(_trg_proc_st2.begin(),fnd);
          const DataArrayInt *trgIds=_trg_ids_st2[locId];
          const int *trgIds2=trgIds->getConstPointer();
          for(std::vector< std::map<int,double> >::const_iterator it1=mat.begin();it1!=mat.end();it1++,rowId++)
            for(std::map<int,double>::const_iterator it2=(*it1).begin();it2!=(*it1).end();it2++)
              denoM[rowId][(*it2).first]=deno[trgIds2[rowId]];
        }
    }
}

/*!
 * This method performs step #0/3 in serialization process.
 * \param count tells specifies nb of elems to send to corresponding proc id. size equal to _group.size().
 * \param offsets tells for a proc i where to start serialize#0 matrix. size equal to _group.size().
 * \param nbOfElemsSrc of size _group.size(). Comes from previous all2all call. tells how many srcIds per proc contains matrix for current proc.
 */
void OverlapMapping::serializeMatrixStep0ST(const int *nbOfElemsSrc, int *&bigArr, int *count, int *offsets,
                                            int *countForRecv, int *offsetsForRecv) const
{
  int grpSize=_group.size();
  std::fill<int *>(count,count+grpSize,0);
  int szz=0;
  int myProcId=_group.myRank();
  for(std::size_t i=0;i<_matrixes_st.size();i++)
    {
      if(_source_proc_id_st[i]==myProcId)// && _target_proc_id_st[i]!=myProcId
        {
          count[_target_proc_id_st[i]]=_matrixes_st[i].size()+1;
          szz+=_matrixes_st[i].size()+1;
        }
    }
  bigArr=new int[szz];
  offsets[0]=0;
  for(int i=1;i<grpSize;i++)
    offsets[i]=offsets[i-1]+count[i-1];
  for(std::size_t i=0;i<_matrixes_st.size();i++)
    {
      if(_source_proc_id_st[i]==myProcId)
        {
          int start=offsets[_target_proc_id_st[i]];
          int *work=bigArr+start;
          *work=0;
          const std::vector< std::map<int,double> >& mat=_matrixes_st[i];
          for(std::vector< std::map<int,double> >::const_iterator it=mat.begin();it!=mat.end();it++,work++)
            work[1]=work[0]+(*it).size();
        }
    }
  //
  offsetsForRecv[0]=0;
  for(int i=0;i<grpSize;i++)
    {
      if(nbOfElemsSrc[i]>0)
        countForRecv[i]=nbOfElemsSrc[i]+1;
      else
        countForRecv[i]=0;
      if(i>0)
        offsetsForRecv[i]=offsetsForRecv[i-1]+countForRecv[i-1];
    }
}

/*!
 * This method performs step#1 and step#2/3. It returns the size of expected array to get allToAllV.
 */
int OverlapMapping::serializeMatrixStep1ST(const int *nbOfElemsSrc, const int *recvStep0, const int *countStep0, const int *offsStep0,
                                           int *&bigArrI, double *&bigArrD, int *count, int *offsets,
                                           int *countForRecv, int *offsForRecv) const
{
  int grpSize=_group.size();
  int myProcId=_group.myRank();
  offsForRecv[0]=0;
  int szz=0;
  for(int i=0;i<grpSize;i++)
    {
      if(nbOfElemsSrc[i]!=0)
        countForRecv[i]=recvStep0[offsStep0[i]+nbOfElemsSrc[i]];
      else
        countForRecv[i]=0;
      szz+=countForRecv[i];
      if(i>0)
        offsForRecv[i]=offsForRecv[i-1]+countForRecv[i-1];
    }
  //
  std::fill(count,count+grpSize,0);
  offsets[0]=0;
  int fullLgth=0;
  for(std::size_t i=0;i<_matrixes_st.size();i++)
    {
      if(_source_proc_id_st[i]==myProcId)
        {
          const std::vector< std::map<int,double> >& mat=_matrixes_st[i];
          int lgthToSend=0;
          for(std::vector< std::map<int,double> >::const_iterator it=mat.begin();it!=mat.end();it++)
            lgthToSend+=(*it).size();
          count[_target_proc_id_st[i]]=lgthToSend;
          fullLgth+=lgthToSend;
        }
    }
  for(int i=1;i<grpSize;i++)
    offsets[i]=offsets[i-1]+count[i-1];
  //
  bigArrI=new int[fullLgth];
  bigArrD=new double[fullLgth];
  // feeding arrays
  fullLgth=0;
  for(std::size_t i=0;i<_matrixes_st.size();i++)
    {
      if(_source_proc_id_st[i]==myProcId)
        {
          const std::vector< std::map<int,double> >& mat=_matrixes_st[i];
          for(std::vector< std::map<int,double> >::const_iterator it1=mat.begin();it1!=mat.end();it1++)
            {
              int j=0;
              for(std::map<int,double>::const_iterator it2=(*it1).begin();it2!=(*it1).end();it2++,j++)
                {
                  bigArrI[fullLgth+j]=(*it2).first;
                  bigArrD[fullLgth+j]=(*it2).second;
                }
              fullLgth+=(*it1).size();
            }
        }
    }
  return szz;
}

/*!
 * This is the last step after all2Alls for matrix exchange.
 * _the_matrix_st is the final matrix : 
 *      - The first entry is srcId in current proc.
 *      - The second is the pseudo id of source proc (correspondance with true id is in attribute _the_matrix_st_source_proc_id and _the_matrix_st_source_ids)
 *      - the third is the srcId in the pseudo source proc
 */
void OverlapMapping::unserializationST(int nbOfTrgElems,
                                       const int *nbOfElemsSrcPerProc,//first all2all
                                       const int *bigArrRecv, const int *bigArrRecvCounts, const int *bigArrRecvOffs,//2nd all2all
                                       const int *bigArrRecv2, const double *bigArrDRecv2, const int *bigArrRecv2Count, const int *bigArrRecv2Offs)//3rd and 4th all2alls
{
  _the_matrix_st.clear();
  _the_matrix_st_source_proc_id.clear();
  //
  int grpSize=_group.size();
  for(int i=0;i<grpSize;i++)
    if(nbOfElemsSrcPerProc[i]!=0)
      _the_matrix_st_source_proc_id.push_back(i);
  int nbOfPseudoProcs=_the_matrix_st_source_proc_id.size();//_the_matrix_st_target_proc_id.size() contains number of matrix fetched remotely whose sourceProcId==myProcId
  _the_matrix_st.resize(nbOfPseudoProcs);
  //
  int j=0;
  for(int i=0;i<grpSize;i++)
    if(nbOfElemsSrcPerProc[i]!=0)
      {
        _the_matrix_st[j].resize(nbOfElemsSrcPerProc[i]);
        for(int k=0;k<nbOfElemsSrcPerProc[i];k++)
          {
            int offs=bigArrRecv[bigArrRecvOffs[i]+k];
            int lgthOfMap=bigArrRecv[bigArrRecvOffs[i]+k+1]-offs;
            for(int l=0;l<lgthOfMap;l++)
              _the_matrix_st[j][k][bigArrRecv2[bigArrRecv2Offs[i]+offs+l]]=bigArrDRecv2[bigArrRecv2Offs[i]+offs+l];
          }
        j++;
      }
}

/*!
 * This method should be called when all remote matrix with sourceProcId==thisProcId have been retrieved and are in 'this->_the_matrix_st' and 'this->_the_matrix_st_target_proc_id'
 * and 'this->_the_matrix_st_target_ids'.
 * This method finish the job of filling 'this->_the_matrix_st' and 'this->_the_matrix_st_target_proc_id' by putting candidates in 'this->_matrixes_st' into them.
 */
void OverlapMapping::finishToFillFinalMatrixST()
{
  int myProcId=_group.myRank();
  int sz=_matrixes_st.size();
  int nbOfEntryToAdd=0;
  for(int i=0;i<sz;i++)
    if(_source_proc_id_st[i]!=myProcId)
      nbOfEntryToAdd++;
  if(nbOfEntryToAdd==0)
    return ;
  int oldNbOfEntry=_the_matrix_st.size();
  int newNbOfEntry=oldNbOfEntry+nbOfEntryToAdd;
  _the_matrix_st.resize(newNbOfEntry);
  int j=oldNbOfEntry;
  for(int i=0;i<sz;i++)
    if(_source_proc_id_st[i]!=myProcId)
      {
        const std::vector<std::map<int,double> >& mat=_matrixes_st[i];
        _the_matrix_st[j]=mat;
        _the_matrix_st_source_proc_id.push_back(_source_proc_id_st[i]);
        j++;
      }
  _matrixes_st.clear();
}

/*!
 * This method performs the operation of target ids broadcasting.
 */
void OverlapMapping::prepareIdsToSendST()
{
  CommInterface commInterface=_group.getCommInterface();
  const MPIProcessorGroup *group=static_cast<const MPIProcessorGroup*>(&_group);
  const MPI_Comm *comm=group->getComm();
  int grpSize=_group.size();
  _source_ids_to_send_st.clear();
  _source_ids_to_send_st.resize(grpSize);
  INTERP_KERNEL::AutoPtr<int> nbsend=new int[grpSize];
  std::fill<int *>(nbsend,nbsend+grpSize,0);
  for(std::size_t i=0;i<_the_matrix_st_source_proc_id.size();i++)
    nbsend[_the_matrix_st_source_proc_id[i]]=_the_matrix_st_source_ids[i].size();
  INTERP_KERNEL::AutoPtr<int> nbrecv=new int[grpSize];
  commInterface.allToAll(nbsend,1,MPI_INT,nbrecv,1,MPI_INT,*comm);
  //
  INTERP_KERNEL::AutoPtr<int> nbsend2=new int[grpSize];
  std::copy((int *)nbsend,((int *)nbsend)+grpSize,(int *)nbsend2);
  INTERP_KERNEL::AutoPtr<int> nbsend3=new int[grpSize];
  nbsend3[0]=0;
  for(int i=1;i<grpSize;i++)
    nbsend3[i]=nbsend3[i-1]+nbsend2[i-1];
  int sendSz=nbsend3[grpSize-1]+nbsend2[grpSize-1];
  INTERP_KERNEL::AutoPtr<int> bigDataSend=new int[sendSz];
  for(std::size_t i=0;i<_the_matrix_st_source_proc_id.size();i++)
    {
      int offset=nbsend3[_the_matrix_st_source_proc_id[i]];
      std::copy(_the_matrix_st_source_ids[i].begin(),_the_matrix_st_source_ids[i].end(),((int *)nbsend3)+offset);
    }
  INTERP_KERNEL::AutoPtr<int> nbrecv2=new int[grpSize];
  INTERP_KERNEL::AutoPtr<int> nbrecv3=new int[grpSize];
  std::copy((int *)nbrecv,((int *)nbrecv)+grpSize,(int *)nbrecv2);
  nbrecv3[0]=0;
  for(int i=1;i<grpSize;i++)
    nbrecv3[i]=nbrecv3[i-1]+nbrecv2[i-1];
  int recvSz=nbrecv3[grpSize-1]+nbrecv2[grpSize-1];
  INTERP_KERNEL::AutoPtr<int> bigDataRecv=new int[recvSz];
  //
  commInterface.allToAllV(bigDataSend,nbsend2,nbsend3,MPI_INT,
                          bigDataRecv,nbrecv2,nbrecv3,MPI_INT,
                          *comm);
  for(int i=0;i<grpSize;i++)
    {
      if(nbrecv2[i]>0)
        {
          _source_ids_to_send_st[i].insert(_source_ids_to_send_st[i].end(),((int *)bigDataRecv)+nbrecv3[i],((int *)bigDataRecv)+nbrecv3[i]+nbrecv2[i]);
        }
    }
}

/*!
 * This method performs a transpose multiply of 'fieldInput' and put the result into 'fieldOutput'.
 * 'fieldInput' is expected to be the sourcefield and 'fieldOutput' the targetfield.
 */
void OverlapMapping::multiply(const MEDCouplingFieldDouble *fieldInput, MEDCouplingFieldDouble *fieldOutput) const
{
  int nbOfCompo=fieldInput->getNumberOfComponents();//to improve same number of components to test
  CommInterface commInterface=_group.getCommInterface();
  const MPIProcessorGroup *group=static_cast<const MPIProcessorGroup*>(&_group);
  const MPI_Comm *comm=group->getComm();
  int grpSize=_group.size();
  int myProcId=_group.myRank();
  //
  INTERP_KERNEL::AutoPtr<int> nbsend=new int[grpSize];
  INTERP_KERNEL::AutoPtr<int> nbsend2=new int[grpSize];
  INTERP_KERNEL::AutoPtr<int> nbrecv=new int[grpSize];
  INTERP_KERNEL::AutoPtr<int> nbrecv2=new int[grpSize];
  std::fill<int *>(nbsend,nbsend+grpSize,0);
  std::fill<int *>(nbrecv,nbrecv+grpSize,0);
  nbsend2[0]=0;
  nbrecv2[0]=0;
  std::vector<double> valsToSend;
  for(int i=0;i<grpSize;i++)
    {
      if(std::find(_proc_ids_to_send_vector_st.begin(),_proc_ids_to_send_vector_st.end(),i)!=_proc_ids_to_send_vector_st.end())
        {
          std::vector<int>::const_iterator isItem1=std::find(_src_proc_st2.begin(),_src_proc_st2.end(),i);
          MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> vals;
          if(isItem1!=_src_proc_st2.end())//item1 of step2 main algo
            {
              int id=std::distance(_src_proc_st2.begin(),isItem1);
              vals=fieldInput->getArray()->selectByTupleId(_src_ids_st2[id]->getConstPointer(),_src_ids_st2[id]->getConstPointer()+_src_ids_st2[id]->getNumberOfTuples());
            }
          else
            {//item0 of step2 main algo
              int id=std::distance(_src_ids_zip_proc_st2.begin(),std::find(_src_ids_zip_proc_st2.begin(),_src_ids_zip_proc_st2.end(),i));
              vals=fieldInput->getArray()->selectByTupleId(&(_src_ids_zip_st2[id])[0],&(_src_ids_zip_st2[id])[0]+_src_ids_zip_st2[id].size());
            }
          nbsend[i]=vals->getNbOfElems();
          valsToSend.insert(valsToSend.end(),vals->getConstPointer(),vals->getConstPointer()+nbsend[i]);
        }
      if(std::find(_proc_ids_to_recv_vector_st.begin(),_proc_ids_to_recv_vector_st.end(),i)!=_proc_ids_to_recv_vector_st.end())
        {
          std::vector<int>::const_iterator isItem0=std::find(_trg_proc_st2.begin(),_trg_proc_st2.end(),i);
          if(isItem0==_trg_proc_st2.end())//item1 of step2 main algo [ (0,1) computed on proc1 and Matrix-Vector on proc1 ]
            {
              std::vector<int>::const_iterator it1=std::find(_src_ids_proc_st2.begin(),_src_ids_proc_st2.end(),i);
              if(it1!=_src_ids_proc_st2.end())
                {
                  int id=std::distance(_src_ids_proc_st2.begin(),it1);
                  nbrecv[i]=_nb_of_src_ids_proc_st2[id]*nbOfCompo;
                }
              else if(i==myProcId)
                {
                  nbrecv[i]=fieldInput->getNumberOfTuplesExpected()*nbOfCompo;
                }
              else
                throw INTERP_KERNEL::Exception("Plouff ! send email to anthony.geay@cea.fr ! ");
            }
          else
            {//item0 of step2 main algo [ (2,1) computed on proc2 but Matrix-Vector on proc1 ] [(1,0) computed on proc1 but Matrix-Vector on proc0]
              int id=std::distance(_src_ids_zip_proc_st2.begin(),std::find(_src_ids_zip_proc_st2.begin(),_src_ids_zip_proc_st2.end(),i));
              nbrecv[i]=_src_ids_zip_st2[id].size()*nbOfCompo;
            }
        }
    }
  for(int i=1;i<grpSize;i++)
    {
      nbsend2[i]=nbsend2[i-1]+nbsend[i-1];
      nbrecv2[i]=nbrecv2[i-1]+nbrecv[i-1];
    }
  INTERP_KERNEL::AutoPtr<double> bigArr=new double[nbrecv2[grpSize-1]+nbrecv[grpSize-1]];
  commInterface.allToAllV(&valsToSend[0],nbsend,nbsend2,MPI_DOUBLE,
                          bigArr,nbrecv,nbrecv2,MPI_DOUBLE,*comm);
  fieldOutput->getArray()->fillWithZero();
  INTERP_KERNEL::AutoPtr<double> tmp=new double[nbOfCompo];
  for(int i=0;i<grpSize;i++)
    {
      if(nbrecv[i]>0)
        {
          double *pt=fieldOutput->getArray()->getPointer();
          std::vector<int>::const_iterator it=std::find(_the_matrix_st_source_proc_id.begin(),_the_matrix_st_source_proc_id.end(),i);
          if(it==_the_matrix_st_source_proc_id.end())
            throw INTERP_KERNEL::Exception("Big problem !");
          int id=std::distance(_the_matrix_st_source_proc_id.begin(),it);
          const std::vector< std::map<int,double> >& mat=_the_matrix_st[id];
          const std::vector< std::map<int,double> >& deno=_the_deno_st[id];
          std::vector<int>::const_iterator isItem0=std::find(_trg_proc_st2.begin(),_trg_proc_st2.end(),i);
          if(isItem0==_trg_proc_st2.end())//item1 of step2 main algo [ (0,1) computed on proc1 and Matrix-Vector on proc1 ]
            {
              int nbOfTrgTuples=mat.size();
              for(int j=0;j<nbOfTrgTuples;j++,pt+=nbOfCompo)
                {
                  const std::map<int,double>& mat1=mat[j];
                  const std::map<int,double>& deno1=deno[j];
                  std::map<int,double>::const_iterator it4=deno1.begin();
                  for(std::map<int,double>::const_iterator it3=mat1.begin();it3!=mat1.end();it3++,it4++)
                    {
                      std::transform(bigArr+nbrecv2[i]+((*it3).first)*nbOfCompo,bigArr+nbrecv2[i]+((*it3).first+1)*(nbOfCompo),(double *)tmp,std::bind2nd(std::multiplies<double>(),(*it3).second/(*it4).second));
                      std::transform((double *)tmp,(double *)tmp+nbOfCompo,pt,pt,std::plus<double>());
                    }
                }
            }
          else
            {//item0 of step2 main algo [ (2,1) computed on proc2 but Matrix-Vector on proc1 ]
              double *pt=fieldOutput->getArray()->getPointer();
              std::map<int,int> zipCor;
              int id=std::distance(_src_ids_zip_proc_st2.begin(),std::find(_src_ids_zip_proc_st2.begin(),_src_ids_zip_proc_st2.end(),i));
              const std::vector<int> zipIds=_src_ids_zip_st2[id];
              int newId=0;
              for(std::vector<int>::const_iterator it=zipIds.begin();it!=zipIds.end();it++,newId++)
                zipCor[*it]=newId;
              int id2=std::distance(_trg_proc_st2.begin(),std::find(_trg_proc_st2.begin(),_trg_proc_st2.end(),i));
              const DataArrayInt *tgrIds=_trg_ids_st2[id2];
              const int *tgrIds2=tgrIds->getConstPointer();
              int nbOfTrgTuples=mat.size();
              for(int j=0;j<nbOfTrgTuples;j++)
                {
                  const std::map<int,double>& mat1=mat[j];
                  const std::map<int,double>& deno1=deno[j];
                  std::map<int,double>::const_iterator it5=deno1.begin();
                  for(std::map<int,double>::const_iterator it3=mat1.begin();it3!=mat1.end();it3++,it5++)
                    {
                      std::map<int,int>::const_iterator it4=zipCor.find((*it3).first);
                      if(it4==zipCor.end())
                        throw INTERP_KERNEL::Exception("Hmmmmm send e mail to anthony.geay@cea.fr !");
                      std::transform(bigArr+nbrecv2[i]+((*it4).second)*nbOfCompo,bigArr+nbrecv2[i]+((*it4).second+1)*(nbOfCompo),(double *)tmp,std::bind2nd(std::multiplies<double>(),(*it3).second/(*it5).second));
                      std::transform((double *)tmp,(double *)tmp+nbOfCompo,pt+tgrIds2[j]*nbOfCompo,pt+tgrIds2[j]*nbOfCompo,std::plus<double>());
                    }
                }
            }
        }
    }
}

/*!
 * This method performs a transpose multiply of 'fieldInput' and put the result into 'fieldOutput'.
 * 'fieldInput' is expected to be the targetfield and 'fieldOutput' the sourcefield.
 */
void OverlapMapping::transposeMultiply(const MEDCouplingFieldDouble *fieldInput, MEDCouplingFieldDouble *fieldOutput)
{
}

/*!
 * This method should be called immediately after _the_matrix_st has been filled with remote computed matrix put in this proc for Matrix-Vector.
 * This method computes for these matrix the minimal set of source ids corresponding to the source proc id. 
 */
void OverlapMapping::updateZipSourceIdsForFuture()
{
  CommInterface commInterface=_group.getCommInterface();
  int myProcId=_group.myRank();
  int nbOfMatrixRecveived=_the_matrix_st_source_proc_id.size();
  for(int i=0;i<nbOfMatrixRecveived;i++)
    {
      int curSrcProcId=_the_matrix_st_source_proc_id[i];
      if(curSrcProcId!=myProcId)
        {
          const std::vector< std::map<int,double> >& mat=_the_matrix_st[i];
          _src_ids_zip_proc_st2.push_back(curSrcProcId);
          _src_ids_zip_st2.resize(_src_ids_zip_st2.size()+1);
          std::set<int> s;
          for(std::vector< std::map<int,double> >::const_iterator it1=mat.begin();it1!=mat.end();it1++)
            for(std::map<int,double>::const_iterator it2=(*it1).begin();it2!=(*it1).end();it2++)
              s.insert((*it2).first);
          _src_ids_zip_st2.back().insert(_src_ids_zip_st2.back().end(),s.begin(),s.end());
        }
    }
}

// #include <iostream>

// void OverlapMapping::printTheMatrix() const
// {
//   CommInterface commInterface=_group.getCommInterface();
//   const MPIProcessorGroup *group=static_cast<const MPIProcessorGroup*>(&_group);
//   const MPI_Comm *comm=group->getComm();
//   int grpSize=_group.size();
//   int myProcId=_group.myRank();
//   std::cerr << "I am proc #" << myProcId << std::endl;
//   int nbOfMat=_the_matrix_st.size();
//   std::cerr << "I do manage " << nbOfMat << "matrix : "<< std::endl;
//   for(int i=0;i<nbOfMat;i++)
//     {
//       std::cerr << "   - Matrix #" << i << " on source proc #" << _the_matrix_st_source_proc_id[i];
//       const std::vector< std::map<int,double> >& locMat=_the_matrix_st[i];
//       for(std::vector< std::map<int,double> >::const_iterator it1=locMat.begin();it1!=locMat.end();it1++)
//         {
//           for(std::map<int,double>::const_iterator it2=(*it1).begin();it2!=(*it1).end();it2++)
//             std::cerr << "(" << (*it2).first << "," << (*it2).second << "), ";
//           std::cerr << std::endl;
//         }
//     }
//   std::cerr << "*********" << std::endl;
// }
