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

  typedef std::map<int,double> SparseDoubleVec;

  /*!
   * Internal class, not part of the public API.
   *
   * Used by the impl of OverlapInterpolationMatrix, plays an equivalent role than what the NxM_Mapping
   * does for the InterpolationMatrix.
   *
   */
  class OverlapMapping
  {
  public:

    OverlapMapping(const ProcessorGroup& group);
    void keepTracksOfSourceIds(int procId, DataArrayInt *ids);
    void keepTracksOfTargetIds(int procId, DataArrayInt *ids);
    void addContributionST(const std::vector< SparseDoubleVec >& matrixST, const DataArrayInt *srcIds, int srcProcId, const DataArrayInt *trgIds, int trgProcId);
    void prepare(const std::vector< int >& procsToSendField, int nbOfTrgElems);
    void computeDenoConservativeVolumic(int nbOfTuplesTrg);
    void computeDenoGlobConstraint();
    //
    void multiply(const MEDCouplingFieldDouble *fieldInput, MEDCouplingFieldDouble *fieldOutput) const;
    void transposeMultiply(const MEDCouplingFieldDouble *fieldInput, MEDCouplingFieldDouble *fieldOutput);
  private:
    void fillProcToSendRcvForMultiply(const std::vector< int >& procsToSendField);
    void serializeMatrixStep0ST(const int *nbOfElemsSrc, int *&bigArr, int *count, int *offsets,
                                int *countForRecv, int *offsetsForRecv) const;
    int serializeMatrixStep1ST(const int *nbOfElemsSrc, const int *recvStep0, const int *countStep0, const int *offsStep0,
                               int *&bigArrI, double *&bigArrD, int *count, int *offsets,
                               int *countForRecv, int *offsForRecv) const;
    void unserializationST(int nbOfTrgElems, const int *nbOfElemsSrcPerProc, const int *bigArrRecv, const int *bigArrRecvCounts, const int *bigArrRecvOffs,
                           const int *bigArrRecv2, const double *bigArrDRecv2, const int *bigArrRecv2Count, const int *bigArrRecv2Offs);
    void finishToFillFinalMatrixST();
    void updateZipSourceIdsForFuture();

    // Debug
//    void printMatrixesST() const;
//    void printTheMatrix() const;
  private:
    const ProcessorGroup &_group;
    /**! Vector of DAInt of cell identifiers. The 2 following class members work in pair. For a proc ID i,
     * first member gives an old2new map for the local part of the source mesh that has been sent.
     * Second member gives proc ID.  */
    std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > _sent_src_ids_st2;
    //! see above _sent_src_ids_st2
    std::vector< int > _sent_src_proc_st2;

    //! See _src_ids_st2 and _sent_src_proc_st2. Same for target mesh.
    std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > _sent_trg_ids_st2;
    //! See _src_ids_st2 and _sent_src_proc_st2. Same for target mesh.
    std::vector< int > _sent_trg_proc_st2;


    /**! Vector of matrixes (partial interpolation ratios), result of the local interpolator run.
     * Indexing shared with _source_proc_id_st, and _target_proc_id_st.   */
    std::vector< std::vector< SparseDoubleVec > > _matrixes_st;
    //! See _matrixes_st - vec of source proc IDs
    std::vector< int > _source_proc_id_st;
    //! See _matrixes_st - vec of target proc IDs
    std::vector< int > _target_proc_id_st;

    //! Vector of remote remote proc IDs for source mesh. Indexing shared with _nb_of_src_ids_proc_st2
    std::vector< int > _src_ids_proc_st2;
    //! Number of cells in the mesh/mapping received from the remote proc i for source mesh. See _src_ids_proc_st2 above
    std::vector< int > _nb_of_src_ids_proc_st2;

    /**! Specifies for each remote proc ID (given in _src_ids_zip_proc_st2 below) the corresponding local
     * source cell IDs to use/send. Same indexing as _src_ids_zip_proc_st2. Sorted.
     * On a given proc, those two members contain exactly the same set of cell identifiers as what is given
     * in the locally held interpolation matrices.   */
    std::vector< std::vector<int> > _src_ids_zip_st2;
    //! Vector of remote proc ID to which the local source mapping above corresponds. See _src_ids_zip_st2 above.
    std::vector< int > _src_ids_zip_proc_st2;

    /**! THE matrix for matrix-vector product. The first dimension is indexed in the set of target procs
    * that interacts with local source mesh. The second dim is the pseudo id of source proc.
    * Same indexing as _the_matrix_st_source_proc_id  */
    std::vector< std::vector< SparseDoubleVec > > _the_matrix_st;
    //! See _the_matrix_st above. List of source proc IDs contributing to _the_matrix_st
    std::vector< int > _the_matrix_st_source_proc_id;

    //! Proc IDs to which data will be sent (originating this current proc) for matrix-vector computation
    std::vector< int > _proc_ids_to_send_vector_st;
    //! Proc IDs from which data will be received (on this current proc) for matrix-vector computation
    std::vector< int > _proc_ids_to_recv_vector_st;

    // Denominators (computed from the numerator matrix)
    std::vector< std::vector< SparseDoubleVec > > _the_deno_st;

//    std::vector< std::vector<int> > _the_matrix_st_source_ids;
//    //! this attribute is of size _group.size(); for each procId in _group _source_ids_to_send_st[procId] contains tupleId to send abroad
//    std::vector< std::vector<int> > _source_ids_to_send_st;
  };
}

#endif
