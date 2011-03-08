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

#ifndef __OVERLAPMAPPING_HXX__
#define __OVERLAPMAPPING_HXX__

#include <vector>
#include <map>

namespace ParaMEDMEM
{
  class ProcessorGroup;
  class MEDCouplingFieldDouble;

  class OverlapMapping
  {
  public:
    OverlapMapping(const ProcessorGroup& group);
    void addContributionST(const std::vector< std::map<int,double> >& matrixST, const int *srcIds, const int *trgIds, int trgIdsLgth, int srcProcId, int trgProcId);
    void prepare(const std::vector< std::vector<int> >& procsInInteraction, int nbOfSrcElems);
    void computeDenoGlobConstraint();
    //
    void transposeMultiply(const MEDCouplingFieldDouble *fieldInput, MEDCouplingFieldDouble *fieldOutput);
  private:
    void serializeMatrixStep0ST(const int *nbOfElemsSrc, int *&bigArr, int *count, int *offsets,
                                int *countForRecv, int *offsetsForRecv) const;
    int serializeMatrixStep1ST(const int *nbOfElemsSrc, const int *recvStep0, const int *countStep0, const int *offsStep0,
                               int *&bigArrI, double *&bigArrD, int *count, int *offsets,
                               int *countForRecv, int *offsForRecv) const;
    void unserializationST(int nbOfSrcElems, const int *nbOfElemsSrcPerProc, const int *bigArrRecv, const int *bigArrRecvCounts, const int *bigArrRecvOffs,
                           const int *bigArrRecv2, const double *bigArrDRecv2, const int *bigArrRecv2Count, const int *bigArrRecv2Offs);
    void finishToFillFinalMatrixST(int nbOfSrcElems);
    void prepareIdsToSendST();
  private:
    const ProcessorGroup &_group;
    //! vector of matrixes the first entry correspond to source proc id in _source_ids_st
    std::vector< std::vector< std::map<int,double> > > _matrixes_st;
    std::vector< std::vector<int> > _source_ids_st;
    std::vector< int > _source_proc_id_st;
    std::vector< std::vector<int> > _target_ids_st;
    std::vector< int > _target_proc_id_st;
    //! the matrix for matrix-vector product. The first dimension the set of target procs that interacts with local source mesh. The second dimension correspond to nb of local source ids. 
    std::vector< std::vector< std::map<int,double> > > _the_matrix_st;
    std::vector< int > _the_matrix_st_target_proc_id;
    std::vector< std::vector<int> > _the_matrix_st_target_ids;
    std::vector< std::vector< std::map<int,double> > > _the_deno_st;
    //! this attribute is of size _group.size(); for each procId in _group _target_ids_to_send_st[procId] contains tupleId to send abroad
    std::vector< std::vector<int> > _target_ids_to_send_st;
  };
}

#endif
