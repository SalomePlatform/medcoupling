// Copyright (C) 2007-2015  CEA/DEN, EDF R&D
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

#ifndef __OVERLAPMAPPING_HXX__
#define __OVERLAPMAPPING_HXX__

#include "MEDCouplingAutoRefCountObjectPtr.hxx"

#include <vector>
#include <map>

namespace ParaMEDMEM
{
  class ProcessorGroup;
  class DataArrayInt;
  class MEDCouplingFieldDouble;

  class OverlapMapping
  {
  public:
    OverlapMapping(const ProcessorGroup& group);
    void keepTracksOfSourceIds(int procId, DataArrayInt *ids);
    void keepTracksOfTargetIds(int procId, DataArrayInt *ids);
    void addContributionST(const std::vector< std::map<int,double> >& matrixST, const DataArrayInt *srcIds, int srcProcId, const DataArrayInt *trgIds, int trgProcId);
    void prepare(const std::vector< std::vector<int> >& procsInInteraction, int nbOfTrgElems);
    void computeDenoConservativeVolumic(int nbOfTuplesTrg);
    void computeDenoGlobConstraint();
    //
    void multiply(const MEDCouplingFieldDouble *fieldInput, MEDCouplingFieldDouble *fieldOutput) const;
    void transposeMultiply(const MEDCouplingFieldDouble *fieldInput, MEDCouplingFieldDouble *fieldOutput);
  private:
    void serializeMatrixStep0ST(const int *nbOfElemsSrc, int *&bigArr, int *count, int *offsets,
                                int *countForRecv, int *offsetsForRecv) const;
    int serializeMatrixStep1ST(const int *nbOfElemsSrc, const int *recvStep0, const int *countStep0, const int *offsStep0,
                               int *&bigArrI, double *&bigArrD, int *count, int *offsets,
                               int *countForRecv, int *offsForRecv) const;
    void unserializationST(int nbOfTrgElems, const int *nbOfElemsSrcPerProc, const int *bigArrRecv, const int *bigArrRecvCounts, const int *bigArrRecvOffs,
                           const int *bigArrRecv2, const double *bigArrDRecv2, const int *bigArrRecv2Count, const int *bigArrRecv2Offs);
    void finishToFillFinalMatrixST();
    void prepareIdsToSendST();
    void updateZipSourceIdsForFuture();
    //void printTheMatrix() const;
  private:
    const ProcessorGroup &_group;
    //! vector of ids
    std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > _src_ids_st2;//item #1
    std::vector< int > _src_proc_st2;//item #1
    std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > _trg_ids_st2;//item #0
    std::vector< int > _trg_proc_st2;//item #0
    std::vector< int > _nb_of_src_ids_proc_st2;//item #1
    std::vector< int > _src_ids_proc_st2;//item #1
    std::vector< std::vector<int> > _src_ids_zip_st2;//same size as _src_ids_zip_proc_st2. Sorted. specifies for each id the corresponding ids to send. This is for item0 of Step2 of main algorithm
    std::vector< int > _src_ids_zip_proc_st2;
    //! vector of matrixes the first entry correspond to source proc id in _source_ids_st
    std::vector< std::vector< std::map<int,double> > > _matrixes_st;
    std::vector< std::vector<int> > _source_ids_st;
    std::vector< int > _source_proc_id_st;
    std::vector< std::vector<int> > _target_ids_st;
    std::vector< int > _target_proc_id_st;
    //! the matrix for matrix-vector product. The first dimension the set of target procs that interacts with local source mesh. The second dimension correspond to nb of local source ids. 
    std::vector< std::vector< std::map<int,double> > > _the_matrix_st;
    std::vector< int > _the_matrix_st_source_proc_id;
    std::vector< std::vector<int> > _the_matrix_st_source_ids;
    std::vector< std::vector< std::map<int,double> > > _the_deno_st;
    //! this attribute stores the proc ids that wait for data from this proc ids for matrix-vector computation
    std::vector< int > _proc_ids_to_send_vector_st;
    std::vector< int > _proc_ids_to_recv_vector_st;
    //! this attribute is of size _group.size(); for each procId in _group _source_ids_to_send_st[procId] contains tupleId to send abroad
    std::vector< std::vector<int> > _source_ids_to_send_st;
  };
}

#endif
