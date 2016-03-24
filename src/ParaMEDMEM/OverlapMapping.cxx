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
// Author : Anthony Geay (CEA/DEN)

#include "OverlapMapping.hxx"
#include "MPIProcessorGroup.hxx"

#include "MEDCouplingFieldDouble.hxx"
#include "MCAuto.hxx"

#include "InterpKernelAutoPtr.hxx"

#include <numeric>
#include <algorithm>

using namespace MEDCoupling;

OverlapMapping::OverlapMapping(const ProcessorGroup& group, const OverlapElementLocator & loc):
    _group(group),_locator(loc)
{
}

/*!
 * Keeps the link between a given a proc holding source mesh data, and the corresponding cell IDs.
 */
void OverlapMapping::keepTracksOfSourceIds(int procId, DataArrayInt *ids)
{
  ids->incrRef();
  _sent_src_ids[procId] = ids;
}

/*!
 * Same as keepTracksOfSourceIds() but for target mesh data.
*/
void OverlapMapping::keepTracksOfTargetIds(int procId, DataArrayInt *ids)
{
  ids->incrRef();
  _sent_trg_ids[procId] = ids;
}

/*!
 * This method stores in the local members the contribution coming from a matrix in format
 * Target(rows)/Source(cols) for a source procId 'srcProcId' and for a target procId 'trgProcId'.
 * All IDs received here (source and target) are in the format of local IDs.
 *
 * @param srcIds is null if the source mesh is on the local proc
 * @param trgIds is null if the source mesh is on the local proc
 *
 * One of the 2 is necessarily null (the two can be null together)
 */
void OverlapMapping::addContributionST(const std::vector< SparseDoubleVec >& matrixST, const DataArrayInt *srcIds, int srcProcId, const DataArrayInt *trgIds, int trgProcId)
{
  _matrixes_st.push_back(matrixST);
  _source_proc_id_st.push_back(srcProcId);
  _target_proc_id_st.push_back(trgProcId);
  if(srcIds)  // source mesh part is remote <=> srcProcId != myRank
      _nb_of_rcv_src_ids[srcProcId] = srcIds->getNumberOfTuples();
  else        // source mesh part is local
    {
      std::set<int> s;
      // For all source IDs (=col indices) in the sparse matrix:
      for(std::vector< SparseDoubleVec >::const_iterator it1=matrixST.begin();it1!=matrixST.end();it1++)
        for(SparseDoubleVec::const_iterator it2=(*it1).begin();it2!=(*it1).end();it2++)
          s.insert((*it2).first);
      vector<int> v(s.begin(), s.end());  // turn set into vector
      _src_ids_zip_comp[trgProcId] = v;
    }
}

/*!
 * This method is in charge to send matrices in AlltoAll mode.
 *
 * 'procsToSendField' gives the list of procs field data has to be sent to.
 * See OverlapElementLocator::computeBoundingBoxesAndTodoList()
 *
 * After the call of this method, 'this' contains the matrixST for all source cells of the current proc
 */
void OverlapMapping::prepare(const std::vector< int >& procsToSendField, int nbOfTrgElems)
{
#ifdef DEC_DEBUG
  printMatrixesST();
#endif

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

  //finish to fill _the_matrix_st with already in place matrix in _matrixes_st (local computation)
  finishToFillFinalMatrixST();

  //updating _src_ids_zip_st2 and _src_ids_zip_st2 with received matrix.
  fillSourceIdsZipReceivedForMultiply();
  // Prepare proc list for future field data exchange (mutliply()):
  _proc_ids_to_send_vector_st = procsToSendField;
  // Make some space on local proc:
  _matrixes_st.clear();

#ifdef DEC_DEBUG
  printTheMatrix();
#endif
}

///*!
// * Compute denominators for ExtensiveConservation interp.
// * TO BE REVISED: needs another communication since some bits are held non locally
// */
//void OverlapMapping::computeDenoGlobConstraint()
//{
//  _the_deno_st.clear();
//  std::size_t sz1=_the_matrix_st.size();
//  _the_deno_st.resize(sz1);
//  for(std::size_t i=0;i<sz1;i++)
//    {
//      std::size_t sz2=_the_matrix_st[i].size();
//      _the_deno_st[i].resize(sz2);
//      for(std::size_t j=0;j<sz2;j++)
//        {
//          double sum=0;
//          SparseDoubleVec& mToFill=_the_deno_st[i][j];
//          const SparseDoubleVec& m=_the_matrix_st[i][j];
//          for(SparseDoubleVec::const_iterator it=m.begin();it!=m.end();it++)
//            sum+=(*it).second;
//          for(SparseDoubleVec::const_iterator it=m.begin();it!=m.end();it++)
//            mToFill[(*it).first]=sum;
//        }
//    }
//    printDenoMatrix();
//}

///*! Compute integral denominators
// * TO BE REVISED: needs another communication since some source areas are held non locally
// */
//void OverlapMapping::computeDenoIntegral()
//{
//  _the_deno_st.clear();
//  std::size_t sz1=_the_matrix_st.size();
//  _the_deno_st.resize(sz1);
//  for(std::size_t i=0;i<sz1;i++)
//    {
//      std::size_t sz2=_the_matrix_st[i].size();
//      _the_deno_st[i].resize(sz2);
//      for(std::size_t j=0;j<sz2;j++)
//        {
//          SparseDoubleVec& mToFill=_the_deno_st[i][j];
//          for(SparseDoubleVec::const_iterator it=mToFill.begin();it!=mToFill.end();it++)
//            mToFill[(*it).first] = sourceAreas;
//        }
//    }
//    printDenoMatrix();
//}

/*! Compute rev integral denominators
  */
void OverlapMapping::computeDenoRevIntegral(const DataArrayDouble & targetAreas)
{
  _the_deno_st.clear();
  std::size_t sz1=_the_matrix_st.size();
  _the_deno_st.resize(sz1);
  const double * targetAreasP = targetAreas.getConstPointer();
  for(std::size_t i=0;i<sz1;i++)
    {
      std::size_t sz2=_the_matrix_st[i].size();
      _the_deno_st[i].resize(sz2);
      for(std::size_t j=0;j<sz2;j++)
        {
          SparseDoubleVec& mToFill=_the_deno_st[i][j];
          SparseDoubleVec& mToIterate=_the_matrix_st[i][j];
          for(SparseDoubleVec::const_iterator it=mToIterate.begin();it!=mToIterate.end();it++)
            mToFill[(*it).first] = targetAreasP[j];
        }
    }
//    printDenoMatrix();
}


/*!
 * Compute denominators for ConvervativeVolumic interp.
 */
void OverlapMapping::computeDenoConservativeVolumic(int nbOfTuplesTrg)
{
  int myProcId=_group.myRank();
  //
  _the_deno_st.clear();
  std::size_t sz1=_the_matrix_st.size();
  _the_deno_st.resize(sz1);
  std::vector<double> deno(nbOfTuplesTrg);
  // Fills in the vector indexed by target cell ID:
  for(std::size_t i=0;i<sz1;i++)
    {
      const std::vector< SparseDoubleVec >& mat=_the_matrix_st[i];
      int curSrcId=_the_matrix_st_source_proc_id[i];
      map < int, MCAuto<DataArrayInt> >::const_iterator isItem1 = _sent_trg_ids.find(curSrcId);
      int rowId=0;
      if(isItem1==_sent_trg_ids.end() || curSrcId==myProcId) // Local computation: simple, because rowId of mat are directly target cell ids.
        {
          for(std::vector< SparseDoubleVec >::const_iterator it1=mat.begin();it1!=mat.end();it1++,rowId++)
            for(SparseDoubleVec::const_iterator it2=(*it1).begin();it2!=(*it1).end();it2++)
              deno[rowId]+=(*it2).second;
        }
      else  // matrix was received, remote computation
        {
          const DataArrayInt *trgIds = (*isItem1).second;
          const int *trgIds2=trgIds->getConstPointer();
          for(std::vector< SparseDoubleVec >::const_iterator it1=mat.begin();it1!=mat.end();it1++,rowId++)
            for(SparseDoubleVec::const_iterator it2=(*it1).begin();it2!=(*it1).end();it2++)
              deno[trgIds2[rowId]]+=(*it2).second;
        }
    }
  // Broadcast the vector into a structure similar to the initial sparse matrix of numerators:
  for(std::size_t i=0;i<sz1;i++)
    {
      int rowId=0;
      const std::vector< SparseDoubleVec >& mat=_the_matrix_st[i];
      int curSrcId=_the_matrix_st_source_proc_id[i];
      map < int, MCAuto<DataArrayInt> >::const_iterator isItem1 = _sent_trg_ids.find(curSrcId);
      std::vector< SparseDoubleVec >& denoM=_the_deno_st[i];
      denoM.resize(mat.size());
      if(isItem1==_sent_trg_ids.end() || curSrcId==myProcId)//item1 of step2 main algo. Simple, because rowId of mat are directly target ids.
        {
          int rowId=0;
          for(std::vector< SparseDoubleVec >::const_iterator it1=mat.begin();it1!=mat.end();it1++,rowId++)
            for(SparseDoubleVec::const_iterator it2=(*it1).begin();it2!=(*it1).end();it2++)
              denoM[rowId][(*it2).first]=deno[rowId];
        }
      else
        {
          const DataArrayInt *trgIds = (*isItem1).second;
          const int *trgIds2=trgIds->getConstPointer();
          for(std::vector< SparseDoubleVec >::const_iterator it1=mat.begin();it1!=mat.end();it1++,rowId++)
            for(SparseDoubleVec::const_iterator it2=(*it1).begin();it2!=(*it1).end();it2++)
              denoM[rowId][(*it2).first]=deno[trgIds2[rowId]];
        }
    }
//  printDenoMatrix();
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
          const std::vector< SparseDoubleVec >& mat=_matrixes_st[i];
          for(std::vector< SparseDoubleVec >::const_iterator it=mat.begin();it!=mat.end();it++,work++)
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
 * It is where the locally computed matrices are serialized to be sent to adequate final proc.
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
          const std::vector< SparseDoubleVec >& mat=_matrixes_st[i];
          int lgthToSend=0;
          for(std::vector< SparseDoubleVec >::const_iterator it=mat.begin();it!=mat.end();it++)
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
          const std::vector< SparseDoubleVec >& mat=_matrixes_st[i];
          for(std::vector< SparseDoubleVec >::const_iterator it1=mat.begin();it1!=mat.end();it1++)
            {
              int j=0;
              for(SparseDoubleVec::const_iterator it2=(*it1).begin();it2!=(*it1).end();it2++,j++)
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
 * This method should be called when all remote matrix with sourceProcId==thisProcId have been retrieved and are
 * in 'this->_the_matrix_st' and 'this->_the_matrix_st_target_proc_id' and 'this->_the_matrix_st_target_ids'.
 * This method finish the job of filling 'this->_the_matrix_st' and 'this->_the_matrix_st_target_proc_id'
 * by putting candidates in 'this->_matrixes_st' into them (i.e. local computation result).
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
        const std::vector<SparseDoubleVec >& mat=_matrixes_st[i];
        _the_matrix_st[j]=mat;
        _the_matrix_st_source_proc_id.push_back(_source_proc_id_st[i]);
        j++;
      }
}


/*!
 * This method performs a transpose multiply of 'fieldInput' and put the result into 'fieldOutput'.
 * 'fieldInput' is expected to be the sourcefield and 'fieldOutput' the targetfield.
 */
void OverlapMapping::multiply(const MEDCouplingFieldDouble *fieldInput, MEDCouplingFieldDouble *fieldOutput, double default_val) const
{
  using namespace std;

  int nbOfCompo=fieldInput->getNumberOfComponents();//to improve same number of components to test
  CommInterface commInterface=_group.getCommInterface();
  const MPIProcessorGroup *group=static_cast<const MPIProcessorGroup*>(&_group);
  const MPI_Comm *comm=group->getComm();
  int grpSize=_group.size();
  int myProcID=_group.myRank();
  //
  INTERP_KERNEL::AutoPtr<int> nbsend=new int[grpSize];
  INTERP_KERNEL::AutoPtr<int> nbsend2=new int[grpSize];
  INTERP_KERNEL::AutoPtr<int> nbrecv=new int[grpSize];
  INTERP_KERNEL::AutoPtr<int> nbrecv2=new int[grpSize];
  fill<int *>(nbsend,nbsend+grpSize,0);
  fill<int *>(nbrecv,nbrecv+grpSize,0);
  nbsend2[0]=0;
  nbrecv2[0]=0;
  vector<double> valsToSend;

  /*
   * FIELD VALUE XCHGE:
   * We call the 'BB source IDs' (bounding box source IDs) the set of source cell IDs transmitted just based on the bounding box information.
   * This is potentially bigger than what is finally in the interp matrix and this is stored in _sent_src_ids.
   * We call 'interp source IDs' the set of source cell IDs with non null entries in the interp matrix. This is a sub-set of the above.
   */
  for(int procID=0;procID<grpSize;procID++)
    {
      /* SENDING part: compute field values to be SENT (and how many of them)
       *   - for all proc 'procID' in group
       *      * if procID == myProcID, send nothing
       *      * elif 'procID' in _proc_ids_to_send_vector_st (computed from the BB intersection)
       *        % if myProcID computed the job (myProcID, procID)
       *           => send only 'interp source IDs' field values (i.e. IDs stored in _src_ids_zip_comp)
       *        % else (=we just sent mesh data to procID, but have never seen the matrix, i.e. matrix was computed remotely by procID)
       *           => send 'BB source IDs' set of field values (i.e. IDs stored in _sent_src_ids)
       */
      if (procID == myProcID)
        nbsend[procID] = 0;
      else
        if(find(_proc_ids_to_send_vector_st.begin(),_proc_ids_to_send_vector_st.end(),procID)!=_proc_ids_to_send_vector_st.end())
          {
            MCAuto<DataArrayDouble> vals;
            if(_locator.isInMyTodoList(myProcID, procID))
              {
                map<int, vector<int> >::const_iterator isItem11 = _src_ids_zip_comp.find(procID);
                if (isItem11 == _src_ids_zip_comp.end())
                  throw INTERP_KERNEL::Exception("OverlapMapping::multiply(): internal error: SEND: unexpected end iterator in _src_ids_zip_comp!");
                const vector<int> & v = (*isItem11).second;
                int sz = v.size();
                vals=fieldInput->getArray()->selectByTupleId(&(v[0]),&(v[0])+sz);
              }
            else
              {
                map < int, MCAuto<DataArrayInt> >::const_iterator isItem11 = _sent_src_ids.find( procID );
                if (isItem11 == _sent_src_ids.end())
                  throw INTERP_KERNEL::Exception("OverlapMapping::multiply(): internal error: SEND: unexpected end iterator in _sent_src_ids!");
                vals=fieldInput->getArray()->selectByTupleId(*(*isItem11).second);
              }
            nbsend[procID] = vals->getNbOfElems();
            valsToSend.insert(valsToSend.end(),vals->getConstPointer(),vals->getConstPointer()+nbsend[procID]);
          }

      /* RECEIVE: compute number of field values to be RECEIVED
       *   - for all proc 'procID' in group
       *      * if procID == myProcID, rcv nothing
       *      * elif 'procID' in _proc_ids_to_recv_vector_st (computed from BB intersec)
       *        % if myProcID computed the job (procID, myProcID)
       *          => receive full set ('BB source IDs') of field data from proc #procID which has never seen the matrix
       *             i.e. prepare to receive the numb in _nb_of_rcv_src_ids
       *        % else (=we did NOT compute the job, hence procID has, and knows the matrix)
       *          => receive 'interp source IDs' set of field values
       */
      const std::vector< int > & _proc_ids_to_recv_vector_st = _the_matrix_st_source_proc_id;
      if (procID == myProcID)
        nbrecv[procID] = 0;
      else
        if(find(_proc_ids_to_recv_vector_st.begin(),_proc_ids_to_recv_vector_st.end(),procID)!=_proc_ids_to_recv_vector_st.end())
          {
            if(_locator.isInMyTodoList(procID, myProcID))
              {
                map <int,int>::const_iterator isItem11 = _nb_of_rcv_src_ids.find(procID);
                if (isItem11 == _nb_of_rcv_src_ids.end())
                  throw INTERP_KERNEL::Exception("OverlapMapping::multiply(): internal error: RCV: unexpected end iterator in _nb_of_rcv_src_ids!");
                nbrecv[procID] = (*isItem11).second;
              }
            else
              {
                map<int, vector<int> >::const_iterator isItem11 = _src_ids_zip_recv.find(procID);
                if (isItem11 == _src_ids_zip_recv.end())
                  throw INTERP_KERNEL::Exception("OverlapMapping::multiply(): internal error: RCV: unexpected end iterator in _src_ids_zip_recv!");
                nbrecv[procID] = (*isItem11).second.size()*nbOfCompo;
              }
          }
    }
  // Compute offsets in the sending/receiving array.
  for(int i=1;i<grpSize;i++)
    {
      nbsend2[i]=nbsend2[i-1]+nbsend[i-1];
      nbrecv2[i]=nbrecv2[i-1]+nbrecv[i-1];
    }
  INTERP_KERNEL::AutoPtr<double> bigArr=new double[nbrecv2[grpSize-1]+nbrecv[grpSize-1]];

#ifdef DEC_DEBUG
  stringstream scout;
  scout << "("  << myProcID << ") nbsend :" << nbsend[0] << "," << nbsend[1] << "," << nbsend[2] << "\n";
  scout << "("  << myProcID << ") nbrecv :" << nbrecv[0] << "," << nbrecv[1] << "," << nbrecv[2] << "\n";
  scout << "("  << myProcID << ") valsToSend: ";
  for (int iii=0; iii<valsToSend.size(); iii++)
    scout << ", " << valsToSend[iii];
  cout << scout.str() << "\n";
#endif

  /*
   * *********************** ALL-TO-ALL
   */
  commInterface.allToAllV(&valsToSend[0],nbsend,nbsend2,MPI_DOUBLE,
                          bigArr,nbrecv,nbrecv2,MPI_DOUBLE,*comm);
#ifdef DEC_DEBUG
  MPI_Barrier(MPI_COMM_WORLD);
  scout << "("  << myProcID << ") bigArray: ";
    for (int iii=0; iii<nbrecv2[grpSize-1]+nbrecv[grpSize-1]; iii++)
      scout << ", " << bigArr[iii];
  cout << scout.str() << "\n";
#endif

  /*
   * TARGET FIELD COMPUTATION (matrix-vec computation)
   */
  fieldOutput->getArray()->fillWithZero();
  INTERP_KERNEL::AutoPtr<double> tmp=new double[nbOfCompo];

  // By default field value set to default value - so mark which cells are hit
  INTERP_KERNEL::AutoPtr<bool> hit_cells = new bool[fieldOutput->getNumberOfTuples()];

  for(vector<int>::const_iterator itProc=_the_matrix_st_source_proc_id.begin(); itProc != _the_matrix_st_source_proc_id.end();itProc++)
  // For each source processor corresponding to a locally held matrix:
    {
      int srcProcID = *itProc;
      int id = distance(_the_matrix_st_source_proc_id.begin(),itProc);
      const vector< SparseDoubleVec >& mat =_the_matrix_st[id];
      const vector< SparseDoubleVec >& deno = _the_deno_st[id];

      /*   FINAL MULTIPLICATION
       *      * if srcProcID == myProcID, local multiplication without any mapping
       *         => for all target cell ID 'tgtCellID'
       *           => for all src cell ID 'srcCellID' in the sparse vector
       *             => tgtFieldLocal[tgtCellID] += srcFieldLocal[srcCellID] * matrix[tgtCellID][srcCellID] / deno[tgtCellID][srcCellID]
       */
      if (srcProcID == myProcID)
        {
          int nbOfTrgTuples=mat.size();
          double * targetBase = fieldOutput->getArray()->getPointer();
          for(int j=0; j<nbOfTrgTuples; j++)
            {
              const SparseDoubleVec& mat1=mat[j];
              const SparseDoubleVec& deno1=deno[j];
              SparseDoubleVec::const_iterator it5=deno1.begin();
              const double * localSrcField = fieldInput->getArray()->getConstPointer();
              double * targetPt = targetBase+j*nbOfCompo;
              for(SparseDoubleVec::const_iterator it3=mat1.begin();it3!=mat1.end();it3++,it5++)
                {
                  // Apply the multiplication for all components:
                  double ratio = (*it3).second/(*it5).second;
                  transform(localSrcField+((*it3).first)*nbOfCompo,
                            localSrcField+((*it3).first+1)*nbOfCompo,
                            (double *)tmp,
                            bind2nd(multiplies<double>(),ratio) );
                  // Accumulate with current value:
                  transform((double *)tmp,(double *)tmp+nbOfCompo,targetPt,targetPt,plus<double>());
                  hit_cells[j] = true;
                }
            }
        }

      if(nbrecv[srcProcID]<=0)  // also covers the preceding 'if'
        continue;

      /*      * if something was received
       *         %  if received matrix (=we didn't compute the job), this means that :
       *            1. we sent part of our targetIDs to srcProcID before, so that srcProcId can do the computation.
       *            2. srcProcID has sent us only the 'interp source IDs' field values
       *            => invert _src_ids_zip_recv -> 'revert_zip'
       *            => for all target cell ID 'tgtCellID'
       *              => mappedTgtID = _sent_trg_ids[srcProcID][tgtCellID]
       *              => for all src cell ID 'srcCellID' in the sparse vector
       *                 => idx = revert_zip[srcCellID]
       *                 => tgtFieldLocal[mappedTgtID] += rcvValue[srcProcID][idx] * matrix[tgtCellID][srcCellID] / deno[tgtCellID][srcCellID]
       */
      if(!_locator.isInMyTodoList(srcProcID, myProcID))
        {
          // invert _src_ids_zip_recv
          map<int,int> revert_zip;
          map<int, vector<int> >::const_iterator it11= _src_ids_zip_recv.find(srcProcID);
          if (it11 == _src_ids_zip_recv.end())
            throw INTERP_KERNEL::Exception("OverlapMapping::multiply(): internal error: MULTIPLY: unexpected end iterator in _src_ids_zip_recv!");
          const vector<int> & vec = (*it11).second;
          int newId=0;
          for(vector<int>::const_iterator it=vec.begin();it!=vec.end();it++,newId++)
            revert_zip[*it]=newId;
          map < int, MCAuto<DataArrayInt> >::const_iterator isItem24 = _sent_trg_ids.find(srcProcID);
          if (isItem24 == _sent_trg_ids.end())
            throw INTERP_KERNEL::Exception("OverlapMapping::multiply(): internal error: MULTIPLY: unexpected end iterator in _sent_trg_ids!");
          const DataArrayInt *tgrIdsDA = (*isItem24).second;
          const int *tgrIds = tgrIdsDA->getConstPointer();

          int nbOfTrgTuples=mat.size();
          double * targetBase = fieldOutput->getArray()->getPointer();
          for(int j=0;j<nbOfTrgTuples;j++)
            {
              const SparseDoubleVec& mat1=mat[j];
              const SparseDoubleVec& deno1=deno[j];
              SparseDoubleVec::const_iterator it5=deno1.begin();
              double * targetPt = targetBase+tgrIds[j]*nbOfCompo;
              for(SparseDoubleVec::const_iterator it3=mat1.begin();it3!=mat1.end();it3++,it5++)
                {
                  map<int,int>::const_iterator it4=revert_zip.find((*it3).first);
                  if(it4==revert_zip.end())
                    throw INTERP_KERNEL::Exception("OverlapMapping::multiply(): internal error: MULTIPLY: unexpected end iterator in revert_zip!");
                  double ratio = (*it3).second/(*it5).second;
                  transform(bigArr+nbrecv2[srcProcID]+((*it4).second)*nbOfCompo,
                            bigArr+nbrecv2[srcProcID]+((*it4).second+1)*nbOfCompo,
                            (double *)tmp,
                            bind2nd(multiplies<double>(),ratio) );
                  transform((double *)tmp,(double *)tmp+nbOfCompo,targetPt,targetPt,plus<double>());
                  hit_cells[tgrIds[j]] = true;
                }
            }
        }
      else
        /*         % else (=we computed the job and we received the 'BB source IDs' set of source field values)
         *            => for all target cell ID 'tgtCellID'
         *              => for all src cell ID 'srcCellID' in the sparse vector
         *                => tgtFieldLocal[tgtCellID] += rcvValue[srcProcID][srcCellID] * matrix[tgtCellID][srcCellID] / deno[tgtCellID][srcCellID]
         */
        {
          // Same loop as in the case srcProcID == myProcID, except that instead of working on local field data, we work on bigArr
          int nbOfTrgTuples=mat.size();
          double * targetBase = fieldOutput->getArray()->getPointer();
          for(int j=0;j<nbOfTrgTuples;j++)
            {
              const SparseDoubleVec& mat1=mat[j];
              const SparseDoubleVec& deno1=deno[j];
              SparseDoubleVec::const_iterator it5=deno1.begin();
              double * targetPt = targetBase+j*nbOfCompo;
              for(SparseDoubleVec::const_iterator it3=mat1.begin();it3!=mat1.end();it3++,it5++)
                {
                  // Apply the multiplication for all components:
                  double ratio = (*it3).second/(*it5).second;
                  transform(bigArr+nbrecv2[srcProcID]+((*it3).first)*nbOfCompo,
                            bigArr+nbrecv2[srcProcID]+((*it3).first+1)*nbOfCompo,
                            (double *)tmp,
                            bind2nd(multiplies<double>(),ratio));
                  // Accumulate with current value:
                  transform((double *)tmp,(double *)tmp+nbOfCompo,targetPt,targetPt,plus<double>());
                  hit_cells[j] = true;
                }
            }
        }
    }

  // Fill in default values for cells which haven't been hit:
  int i = 0;
  for(bool * hit_cells_ptr=hit_cells; i< fieldOutput->getNumberOfTuples(); hit_cells_ptr++,i++)
    if (!(*hit_cells_ptr))
      {
        double * targetPt=fieldOutput->getArray()->getPointer();
        fill(targetPt+i*nbOfCompo, targetPt+(i+1)*nbOfCompo, default_val);
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
 * This method should be called immediately after _the_matrix_st has been filled with remote computed matrix
 * put in this proc for Matrix-Vector.
 * It fills _src_ids_zip_recv (see member doc)
 */
void OverlapMapping::fillSourceIdsZipReceivedForMultiply()
{
  /* When it is called, only the bits received from other processors (i.e. the remotely executed jobs) are in the
    big matrix _the_matrix_st. */

  CommInterface commInterface=_group.getCommInterface();
  int myProcId=_group.myRank();
  int nbOfMatrixRecveived=_the_matrix_st_source_proc_id.size();
  for(int i=0;i<nbOfMatrixRecveived;i++)
    {
      int curSrcProcId=_the_matrix_st_source_proc_id[i];
      if(curSrcProcId!=myProcId)  // if =, data has been populated by addContributionST()
        {
          const std::vector< SparseDoubleVec >& mat=_the_matrix_st[i];
          std::set<int> s;
          for(std::vector< SparseDoubleVec >::const_iterator it1=mat.begin();it1!=mat.end();it1++)
            for(SparseDoubleVec::const_iterator it2=(*it1).begin();it2!=(*it1).end();it2++)
              s.insert((*it2).first);
          vector<int> vec(s.begin(),s.end());
          _src_ids_zip_recv[curSrcProcId] = vec;
        }
    }
}

#ifdef DEC_DEBUG
 void OverlapMapping::printTheMatrix() const
 {
   CommInterface commInterface=_group.getCommInterface();
   const MPIProcessorGroup *group=static_cast<const MPIProcessorGroup*>(&_group);
   const MPI_Comm *comm=group->getComm();
   int grpSize=_group.size();
   int myProcId=_group.myRank();
   std::stringstream oscerr;
   int nbOfMat=_the_matrix_st.size();
   oscerr << "(" <<  myProcId <<  ") I hold " << nbOfMat << " matrix(ces) : "<< std::endl;
   for(int i=0;i<nbOfMat;i++)
     {
       oscerr << "   - Matrix #" << i << " coming from source proc #" << _the_matrix_st_source_proc_id[i] << ":\n ";
       const std::vector< SparseDoubleVec >& locMat=_the_matrix_st[i];
       int j = 0;
       for(std::vector< SparseDoubleVec >::const_iterator it1=locMat.begin();it1!=locMat.end();it1++, j++)
         {
           oscerr << " Target Cell #" << j;
           for(SparseDoubleVec::const_iterator it2=(*it1).begin();it2!=(*it1).end();it2++)
             oscerr << " (" << (*it2).first << "," << (*it2).second << "), ";
           oscerr << std::endl;
         }
     }
   oscerr << "*********" << std::endl;

   // Hope this will be flushed in one go:
   std::cerr << oscerr.str() << std::endl;
//   if(myProcId != 0)
//     MPI_Barrier(MPI_COMM_WORLD);
 }

 void OverlapMapping::printMatrixesST() const
  {
    CommInterface commInterface=_group.getCommInterface();
    const MPIProcessorGroup *group=static_cast<const MPIProcessorGroup*>(&_group);
    const MPI_Comm *comm=group->getComm();
    int grpSize=_group.size();
    int myProcId=_group.myRank();
    std::stringstream oscerr;
    int nbOfMat=_matrixes_st.size();
    oscerr << "(" <<  myProcId <<  ") I hold " << nbOfMat << " LOCAL matrix(ces) : "<< std::endl;
    for(int i=0;i<nbOfMat;i++)
      {
        oscerr << "   - Matrix #" << i << ": (source proc #" << _source_proc_id_st[i] << " / tgt proc#" << _target_proc_id_st[i] << "): \n";
        const std::vector< SparseDoubleVec >& locMat=_matrixes_st[i];
        int j = 0;
        for(std::vector< SparseDoubleVec >::const_iterator it1=locMat.begin();it1!=locMat.end();it1++, j++)
          {
            oscerr << " Target Cell #" << j;
            for(SparseDoubleVec::const_iterator it2=(*it1).begin();it2!=(*it1).end();it2++)
              oscerr << " (" << (*it2).first << "," << (*it2).second << "), ";
            oscerr << std::endl;
          }
      }
    oscerr << "*********" << std::endl;

    // Hope this will be flushed in one go:
    std::cerr << oscerr.str() << std::endl;
  }

 void OverlapMapping::printDenoMatrix() const
   {
     CommInterface commInterface=_group.getCommInterface();
     const MPIProcessorGroup *group=static_cast<const MPIProcessorGroup*>(&_group);
     const MPI_Comm *comm=group->getComm();
     int grpSize=_group.size();
     int myProcId=_group.myRank();
     std::stringstream oscerr;
     int nbOfMat=_the_deno_st.size();
     oscerr << "(" <<  myProcId <<  ") I hold " << nbOfMat << " DENOMINATOR matrix(ces) : "<< std::endl;
     for(int i=0;i<nbOfMat;i++)
       {
         oscerr << "   - Matrix #" << i << " coming from source proc #" << _the_matrix_st_source_proc_id[i] << ": \n";
         const std::vector< SparseDoubleVec >& locMat=_the_deno_st[i];
         int j = 0;
         for(std::vector< SparseDoubleVec >::const_iterator it1=locMat.begin();it1!=locMat.end();it1++, j++)
           {
             oscerr << " Target Cell #" << j;
             for(SparseDoubleVec::const_iterator it2=(*it1).begin();it2!=(*it1).end();it2++)
               oscerr << " (" << (*it2).first << "," << (*it2).second << "), ";
             oscerr << std::endl;
           }
       }
     oscerr << "*********" << std::endl;

     // Hope this will be flushed in one go:
     std::cerr << oscerr.str() << std::endl;
   }
#endif
