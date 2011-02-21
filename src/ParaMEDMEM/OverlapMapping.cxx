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
#include "InterpKernelAutoPtr.hxx"

using namespace ParaMEDMEM;

OverlapMapping::OverlapMapping(const ProcessorGroup& group):_group(group)
{
}

/*!
 * This method stores from a matrix in format Source(rows)/Target(cols) for a source procId 'srcProcId' and for a target procId 'trgProcId'.
 * All ids (source and target) are in format of local ids. 
 */
void OverlapMapping::addContributionST(const std::vector< std::map<int,double> >& matrixST, const int *srcIds, int srcProcId, int trgProcId)
{
  int nbOfRows=matrixST.size();
  _matrixes_st.push_back(matrixST);
  _source_ids_st.resize(_source_ids_st.size()+1);
  _source_ids_st.back().insert(_source_ids_st.back().end(),srcIds,srcIds+nbOfRows);
  _source_proc_id_st.push_back(srcProcId);
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
  for(std::size_t i=0;i<_matrixes_st.size();i++)
    {
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
  commInterface.allToAllV(bigArr2,nbsend2,nbsend3,MPI_DOUBLE,
                          bigArrDRecv2,nbrecv3,nbrecv4,MPI_DOUBLE,
                          *comm);
  //finishing
  unserializationST(nbOfSrcElems,nbrecv,bigArrRecv,nbrecv1,nbrecv2,
                    bigArrRecv2,bigArrDRecv2,nbrecv3,nbrecv4);
}

/*!
 * Compute denominators.
 */
void OverlapMapping::computeDeno()
{
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
  for(std::size_t i=0;i<_matrixes_st.size();i++)
    {
      count[_source_proc_id_st[i]]=2*_matrixes_st[i].size()+1;
      szz+=2*_matrixes_st[i].size()+1;
    }
  bigArr=new int[szz];
  offsets[0]=0;
  for(int i=1;i<grpSize;i++)
    offsets[i]=offsets[i-1]+count[i-1];
  for(std::size_t i=0;i<_matrixes_st.size();i++)
    {
      int start=offsets[_source_proc_id_st[i]];
      int *work=bigArr+start;
      *work=0;
      const std::vector< std::map<int,double> >& mat=_matrixes_st[i];
      for(std::vector< std::map<int,double> >::const_iterator it=mat.begin();it!=mat.end();it++,work++)
        work[1]=work[0]+(*it).size();
      std::copy(_source_ids_st[i].begin(),_source_ids_st[i].end(),work);
    }
  //
  offsetsForRecv[0]=0;
  for(int i=0;i<grpSize;i++)
    {
      countForRecv[i]=2*nbOfElemsSrc[i]+1;
      if(i>0)
        offsetsForRecv[i]=offsetsForRecv[i-1]+2*nbOfElemsSrc[i]+1;
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
        offsForRecv[i]=offsForRecv[i-1]+countForRecv[i];
    }
  //
  std::fill(count,count+grpSize,0);
  offsets[0]=0;
  int fullLgth=0;
  for(std::size_t i=0;i<_matrixes_st.size();i++)
    {
      const std::vector< std::map<int,double> >& mat=_matrixes_st[i];
      int lgthToSend=0;
      for(std::vector< std::map<int,double> >::const_iterator it=mat.begin();it!=mat.end();it++)
        lgthToSend+=(*it).size();
      count[_source_proc_id_st[i]]=lgthToSend;
      fullLgth+=lgthToSend;
    }
  //
  bigArrI=new int[fullLgth];
  bigArrD=new double[fullLgth];
  // feeding arrays
  fullLgth=0;
  for(std::size_t i=0;i<_matrixes_st.size();i++)
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
  return szz;
}

/*!
 * This is the last step after all2Alls for matrix exchange.
 * _the_matrix_st is the final matrix : 
 *      - The first entry is srcId in current proc.
 *      - The second is the pseudo id of target proc (correspondance with true id is in attribute _the_matrix_st_target_proc_id)
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
      _the_matrix_st_target_proc_id.push_back(nbOfElemsSrcPerProc[i]);
  int nbOfPseudoProcs=_the_matrix_st_target_proc_id.size();
  for(int i=0;i<nbOfSrcElems;i++)
    _the_matrix_st.resize(nbOfPseudoProcs);
  //
  int j=0;
  for(int i=0;i<grpSize;i++)
    if(nbOfElemsSrcPerProc[i]!=0)
      {
        for(int k=0;k<nbOfElemsSrcPerProc[i];k++)
          {
            int srcId=bigArrRecv[bigArrRecvOffs[i]+nbOfElemsSrcPerProc[i]+1+k];
            int offs=bigArrRecv[bigArrRecvOffs[i]+k];
            int lgthOfMap=bigArrRecv[bigArrRecvOffs[i]+k+1]-offs;
            for(int l=0;l<lgthOfMap;l++)
              _the_matrix_st[srcId][j][bigArrRecv2[bigArrRecv2Offs[i]+offs+l]]=bigArrDRecv2[bigArrRecv2Offs[i]+offs+l];
          }
        j++;
      }
}
