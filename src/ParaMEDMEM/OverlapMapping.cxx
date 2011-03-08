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
 * This method stores from a matrix in format Source(rows)/Target(cols) for a source procId 'srcProcId' and for a target procId 'trgProcId'.
 * All ids (source and target) are in format of local ids. 
 */
void OverlapMapping::addContributionST(const std::vector< std::map<int,double> >& matrixST, const int *srcIds, const int *trgIds, int trgIdsLgth, int srcProcId, int trgProcId)
{
  int nbOfRows=matrixST.size();
  _matrixes_st.push_back(matrixST);
  _source_ids_st.resize(_source_ids_st.size()+1);
  _source_ids_st.back().insert(_source_ids_st.back().end(),srcIds,srcIds+nbOfRows);
  _source_proc_id_st.push_back(srcProcId);
  //
  _target_ids_st.resize(_target_ids_st.size()+1);
  _target_ids_st.back().insert(_target_ids_st.back().end(),trgIds,trgIds+trgIdsLgth);
  _target_proc_id_st.push_back(trgProcId);
}

/*!
 * 'procsInInteraction' gives the global view of interaction between procs.
 * In 'procsInInteraction' for a proc with id i, is in interaction with procs listed in procsInInteraction[i].
 *
 * This method is in charge to send matrixes in AlltoAll mode.
 * After the call of this method 'this' contains the matrixST for all source elements of the current proc and 
 * matrixTS for all target elements of current proc.
 */
void OverlapMapping::prepare(const std::vector< std::vector<int> >& procsInInteraction, int nbOfSrcElems)
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
  for(std::size_t i=0;i<_matrixes_st.size();i++)
    {
      if(_source_proc_id_st[i]!=myProcId)
        nbsend[_source_proc_id_st[i]]=_matrixes_st[i].size();
    }
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
                          *comm);// sending ids off sparse matrix (n+1 elems) + src ids (n elems)
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
  unserializationST(nbOfSrcElems,nbrecv,bigArrRecv,nbrecv1,nbrecv2,
                    bigArrRecv2,bigArrDRecv2,nbrecv3,nbrecv4);
  //finish to fill _the_matrix_st and _the_matrix_st_target_proc_id with already in place matrix in _matrixes_st
  finishToFillFinalMatrixST(nbOfSrcElems);
  //exchanging target ids for future sending
  prepareIdsToSendST();
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
      if(_source_proc_id_st[i]!=myProcId)
        {
          count[_source_proc_id_st[i]]=2*_matrixes_st[i].size()+1;
          szz+=2*_matrixes_st[i].size()+1;
        }
    }
  bigArr=new int[szz];
  offsets[0]=0;
  for(int i=1;i<grpSize;i++)
    offsets[i]=offsets[i-1]+count[i-1];
  for(std::size_t i=0;i<_matrixes_st.size();i++)
    {
      if(_source_proc_id_st[i]!=myProcId)
        {
          int start=offsets[_source_proc_id_st[i]];
          int *work=bigArr+start;
          *work=0;
          const std::vector< std::map<int,double> >& mat=_matrixes_st[i];
          for(std::vector< std::map<int,double> >::const_iterator it=mat.begin();it!=mat.end();it++,work++)
            work[1]=work[0]+(*it).size();
          std::copy(_source_ids_st[i].begin(),_source_ids_st[i].end(),work+1);
        }
    }
  //
  offsetsForRecv[0]=0;
  for(int i=0;i<grpSize;i++)
    {
      if(nbOfElemsSrc[i]>0)
        countForRecv[i]=2*nbOfElemsSrc[i]+1;
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
      if(_source_proc_id_st[i]!=myProcId)
        {
          const std::vector< std::map<int,double> >& mat=_matrixes_st[i];
          int lgthToSend=0;
          for(std::vector< std::map<int,double> >::const_iterator it=mat.begin();it!=mat.end();it++)
            lgthToSend+=(*it).size();
          count[_source_proc_id_st[i]]=lgthToSend;
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
      if(_source_proc_id_st[i]!=myProcId)
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
 *      - The second is the pseudo id of target proc (correspondance with true id is in attribute _the_matrix_st_target_proc_id and _the_matrix_st_target_ids)
 *      - the third is the trgId in the pseudo target proc
 */
void OverlapMapping::unserializationST(int nbOfSrcElems,
                                       const int *nbOfElemsSrcPerProc,//first all2all
                                       const int *bigArrRecv, const int *bigArrRecvCounts, const int *bigArrRecvOffs,//2nd all2all
                                       const int *bigArrRecv2, const double *bigArrDRecv2, const int *bigArrRecv2Count, const int *bigArrRecv2Offs)//3rd and 4th all2alls
{
  _the_matrix_st.clear();
  _the_matrix_st.resize(nbOfSrcElems);
  _the_matrix_st_target_proc_id.clear();
  //
  int grpSize=_group.size();
  for(int i=0;i<grpSize;i++)
    if(nbOfElemsSrcPerProc[i]!=0)
      _the_matrix_st_target_proc_id.push_back(i);
  int nbOfPseudoProcs=_the_matrix_st_target_proc_id.size();//_the_matrix_st_target_proc_id.size() contains number of matrix fetched remotely whose sourceProcId==myProcId
  _the_matrix_st_target_ids.resize(nbOfPseudoProcs);
  _the_matrix_st.resize(nbOfPseudoProcs);
  for(int i=0;i<nbOfPseudoProcs;i++)
    _the_matrix_st[i].resize(nbOfSrcElems);
  //
  int j=0;
  for(int i=0;i<grpSize;i++)
    if(nbOfElemsSrcPerProc[i]!=0)
      {
        std::set<int> targetIdsZip;// this zip is to reduce amount of data to send/rexcv on transposeMultiply target ids transfert
        for(int k=0;k<nbOfElemsSrcPerProc[i];k++)
          {
            int offs=bigArrRecv[bigArrRecvOffs[i]+k];
            int lgthOfMap=bigArrRecv[bigArrRecvOffs[i]+k+1]-offs;
            for(int l=0;l<lgthOfMap;l++)
              targetIdsZip.insert(bigArrRecv2[bigArrRecv2Offs[i]+offs+l]);
          }
        _the_matrix_st_target_ids[j].insert(_the_matrix_st_target_ids[j].end(),targetIdsZip.begin(),targetIdsZip.end());
        std::map<int,int> old2newTrgIds;
        int newNbTrg=0;
        for(std::set<int>::const_iterator it=targetIdsZip.begin();it!=targetIdsZip.end();it++,newNbTrg++)
          old2newTrgIds[*it]=newNbTrg;
        for(int k=0;k<nbOfElemsSrcPerProc[i];k++)
          {
            int srcId=bigArrRecv[bigArrRecvOffs[i]+nbOfElemsSrcPerProc[i]+1+k];
            int offs=bigArrRecv[bigArrRecvOffs[i]+k];
            int lgthOfMap=bigArrRecv[bigArrRecvOffs[i]+k+1]-offs;
            for(int l=0;l<lgthOfMap;l++)
              _the_matrix_st[j][srcId][old2newTrgIds[bigArrRecv2[bigArrRecv2Offs[i]+offs+l]]]=bigArrDRecv2[bigArrRecv2Offs[i]+offs+l];
          }
        j++;
      }
}

/*!
 * This method should be called when all remote matrix with sourceProcId==thisProcId have been retrieved and are in 'this->_the_matrix_st' and 'this->_the_matrix_st_target_proc_id'
 * and 'this->_the_matrix_st_target_ids'.
 * This method finish the job of filling 'this->_the_matrix_st' and 'this->_the_matrix_st_target_proc_id' by putting candidates in 'this->_matrixes_st' into them.
 */
void OverlapMapping::finishToFillFinalMatrixST(int nbOfSrcElems)
{
  int myProcId=_group.myRank();
  int sz=_matrixes_st.size();
  int nbOfEntryToAdd=0;
  for(int i=0;i<sz;i++)
    if(_source_proc_id_st[i]==myProcId)
      {
        nbOfEntryToAdd++;
        _the_matrix_st_target_proc_id.push_back(_target_proc_id_st[i]);
      }
  if(nbOfEntryToAdd==0)
    return ;
  int oldNbOfEntry=_the_matrix_st.size();
  int newNbOfEntry=oldNbOfEntry+nbOfEntryToAdd;
  _the_matrix_st.resize(newNbOfEntry);
  _the_matrix_st_target_ids.resize(newNbOfEntry);
  for(int i=oldNbOfEntry;i<newNbOfEntry;i++)
    _the_matrix_st[i].resize(nbOfSrcElems);
  int j=oldNbOfEntry;
  for(int i=0;i<sz;i++)
    if(_source_proc_id_st[i]==myProcId)
      {
        const std::vector<std::map<int,double> >& mat=_matrixes_st[i];
        const std::vector<int>& srcIds=_source_ids_st[i];
        int sz=srcIds.size();//assert srcIds.size()==mat.size()
        for(int k=0;k<sz;k++)
          {
            const std::map<int,double>& m2=mat[k];
            for(std::map<int,double>::const_iterator it=m2.begin();it!=m2.end();it++)
              _the_matrix_st[j][srcIds[k]][(*it).first]=(*it).second;
          }
        _the_matrix_st_target_ids[j].insert(_the_matrix_st_target_ids[j].end(),_target_ids_st[i].begin(),_target_ids_st[i].end());
        j++;
      }
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
  _target_ids_to_send_st.clear();
  _target_ids_to_send_st.resize(grpSize);
  INTERP_KERNEL::AutoPtr<int> nbsend=new int[grpSize];
  std::fill<int *>(nbsend,nbsend+grpSize,0);
  for(std::size_t i=0;i<_the_matrix_st_target_proc_id.size();i++)
    nbsend[_the_matrix_st_target_proc_id[i]]=_the_matrix_st_target_ids[i].size();
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
  for(std::size_t i=0;i<_the_matrix_st_target_proc_id.size();i++)
    {
      int offset=nbsend3[_the_matrix_st_target_proc_id[i]];
      std::copy(_the_matrix_st_target_ids[i].begin(),_the_matrix_st_target_ids[i].end(),((int *)nbsend3)+offset);
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
          _target_ids_to_send_st[i].insert(_target_ids_to_send_st[i].end(),((int *)bigDataRecv)+nbrecv3[i],((int *)bigDataRecv)+nbrecv3[i]+nbrecv2[i]);
        }
    }
}

/*!
 * This method performs a transpose multiply of 'fieldInput' and put the result into 'fieldOutput'.
 * 'fieldInput' is expected to be the targetfield and 'fieldOutput' the sourcefield.
 */
void OverlapMapping::transposeMultiply(const MEDCouplingFieldDouble *fieldInput, MEDCouplingFieldDouble *fieldOutput)
{
#if 0
  CommInterface commInterface=_group.getCommInterface();
  const MPIProcessorGroup *group=static_cast<const MPIProcessorGroup*>(&_group);
  const MPI_Comm *comm=group->getComm();
  int grpSize=_group.size();
  INTERP_KERNEL::AutoPtr<int> nbsend=new int[grpSize];
  std::fill<int *>(nbsend,nbsend+grpSize,0);
  for(std::size_t i=0;i<_the_matrix_st.size();i++)
    nbsend[_the_matrix_st_target_proc_id[i]]=_the_matrix_st_target_ids[i].size();
  INTERP_KERNEL::AutoPtr<int> nbrecv=new int[grpSize];
  commInterface.allToAll(nbsend,1,MPI_INT,nbrecv,1,MPI_INT,*comm);
  int nbOfCompo=fieldInput->getNumberOfComponents();//to improve same number of components
  std::transform((int *)nbsend,(int *)nbsend+grpSize,(int *)nbsend,std::bind2nd(std::multiplies<int>(),nbOfCompo));
  std::transform((int *)nbrecv,(int *)nbrecv+grpSize,(int *)nbrecv,std::bind2nd(std::multiplies<int>(),nbOfCompo));
  INTERP_KERNEL::AutoPtr<int> nbsend2=new int[grpSize];
  nbsend2[0]=0;
  for(int i=1;i<grpSize;i++)
    nbsend2[i]=nbsend2[i-1]+nbsend[i-1];
  int szToSend=nbsend2[grpSize-1]+nbsend[grpSize-1];
  INTERP_KERNEL::AutoPtr<double> nbsend3=new double[szToSend];
  double *work=nbsend3;
  for(std::size_t i=0;i<_the_matrix_st.size();i++)
    {
      MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ptr=fieldInput->getArray()->selectByTupleId(&_target_ids_st[i][0],&(_target_ids_st[i][0])+_target_ids_st[i].size());
      std::copy(ptr->getConstPointer(),ptr->getConstPointer()+nbsend[_target_proc_id_st[i]],work+nbsend2[_target_proc_id_st[i]]);
    }
  INTERP_KERNEL::AutoPtr<int> nbrecv3=new int[grpSize];
  nbrecv3[0]=0;
  for(int i=1;i<grpSize;i++)
    nbrecv3[i]=nbrecv3[i-1]+nbrecv[i-1];
  int szToFetch=nbsend2[grpSize-1]+nbrecv[grpSize-1];
  INTERP_KERNEL::AutoPtr<double> nbrecv2=new double[szToFetch];
#endif
  //
  int nbOfCompo=fieldInput->getNumberOfComponents();//to improve same number of components to test
  CommInterface commInterface=_group.getCommInterface();
  const MPIProcessorGroup *group=static_cast<const MPIProcessorGroup*>(&_group);
  const MPI_Comm *comm=group->getComm();
  int grpSize=_group.size();
  //
  INTERP_KERNEL::AutoPtr<int> nbsend=new int[grpSize];
  std::fill<int *>(nbsend,nbsend+grpSize,0);
  INTERP_KERNEL::AutoPtr<int> nbsend2=new int[grpSize];
  for(int i=0;i<grpSize;i++)
    nbsend[i]=nbOfCompo*_target_ids_to_send_st[i].size();
  nbsend2[0];
  for(int i=1;i<grpSize;i++)
    nbsend2[i]=nbsend2[i-1]+nbsend[i-1];
  int szToSend=nbsend2[grpSize-1]+nbsend[grpSize-1];
  INTERP_KERNEL::AutoPtr<double> nbsend3=new double[szToSend];
  INTERP_KERNEL::AutoPtr<int> nbrecv=new int[grpSize];
  INTERP_KERNEL::AutoPtr<int> nbrecv3=new int[grpSize];
  std::fill<int *>(nbrecv,nbrecv+grpSize,0);
  for(std::size_t i=0;i<_the_matrix_st_target_proc_id.size();i++)
    nbrecv[_the_matrix_st_target_proc_id[i]]=_the_matrix_st_target_ids[i].size()*nbOfCompo;
  nbrecv3[0]=0;
  for(int i=1;i<grpSize;i++)
    nbrecv3[i]=nbrecv3[i-1]+nbrecv[i-1];
  int szToRecv=nbrecv3[grpSize-1]+nbrecv[grpSize-1];
  double *work=nbsend3;
  for(int i=0;i<grpSize;i++)
    {
      MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ptr=fieldInput->getArray()->selectByTupleId(&(_target_ids_to_send_st[i][0]),
                                                                                                    &(_target_ids_to_send_st[i][0])+_target_ids_to_send_st[i].size());
      std::copy(ptr->getConstPointer(),ptr->getConstPointer()+nbsend[i],work+nbsend2[i]);
    }
  INTERP_KERNEL::AutoPtr<double> nbrecv2=new double[szToRecv];
  commInterface.allToAllV(nbsend3,nbsend,nbsend2,MPI_DOUBLE,
                          nbrecv2,nbrecv,nbrecv3,MPI_DOUBLE,
                          *comm);
  // deserialization
  fieldOutput->getArray()->fillWithZero();
  INTERP_KERNEL::AutoPtr<double> workZone=new double[nbOfCompo];
  for(std::size_t i=0;i<_the_matrix_st.size();i++)
    {
      double *res=fieldOutput->getArray()->getPointer();
      int targetProcId=_the_matrix_st_target_proc_id[i];
      const std::vector<std::map<int,double> >& m=_the_matrix_st[i];
      const std::vector<std::map<int,double> >& deno=_the_deno_st[i];
      const double *vecOnTargetProcId=((const double *)nbrecv2)+nbrecv3[targetProcId];
      std::size_t nbOfIds=m.size();
      for(std::size_t j=0;j<nbOfIds;j++,res+=nbOfCompo)
        {
          const std::map<int,double>& m2=m[j];
          const std::map<int,double>& deno2=deno[j];
          for(std::map<int,double>::const_iterator it=m2.begin();it!=m2.end();it++)
            {
              std::map<int,double>::const_iterator it2=deno2.find((*it).first);
              std::transform(vecOnTargetProcId+((*it).first*nbOfCompo),vecOnTargetProcId+((*it).first+1)*nbOfCompo,(double *)workZone,std::bind2nd(std::multiplies<double>(),(*it).second));
              std::transform((double *)workZone,(double *)workZone+nbOfCompo,(double *)workZone,std::bind2nd(std::multiplies<double>(),1./(*it2).second));
              std::transform((double *)workZone,(double *)workZone+nbOfCompo,res,res,std::plus<double>());
            }
        }
    }
}
