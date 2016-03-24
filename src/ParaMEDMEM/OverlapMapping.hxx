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

#ifndef __OVERLAPMAPPING_HXX__
#define __OVERLAPMAPPING_HXX__

#include "MCAuto.hxx"
#include "OverlapElementLocator.hxx"

#include <vector>
#include <map>
//#define DEC_DEBUG

namespace MEDCoupling
{
  class ProcessorGroup;
  class DataArrayInt;
  class MEDCouplingFieldDouble;

  using namespace std;
  typedef map<int,double> SparseDoubleVec;

  /*!
   * Internal class, not part of the public API.
   *
   * Used by the impl of OverlapInterpolationMatrix, plays an equivalent role than what the NxM_Mapping
   * does for the InterpolationMatrix.
   */
  class OverlapMapping
  {
  public:

    OverlapMapping(const ProcessorGroup& group, const OverlapElementLocator& locator);
    void keepTracksOfSourceIds(int procId, DataArrayInt *ids);
    void keepTracksOfTargetIds(int procId, DataArrayInt *ids);
    void addContributionST(const vector< SparseDoubleVec >& matrixST, const DataArrayInt *srcIds, int srcProcId, const DataArrayInt *trgIds, int trgProcId);
    void prepare(const vector< int >& procsToSendField, int nbOfTrgElems);
    void computeDenoConservativeVolumic(int nbOfTuplesTrg);
//    void computeDenoIntegralGlobConstraint();
//    void computeDenoIntegral();
    void computeDenoRevIntegral(const DataArrayDouble & targetAreas);
    //
    void multiply(const MEDCouplingFieldDouble *fieldInput, MEDCouplingFieldDouble *fieldOutput, double default_val) const;
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
    void fillSourceIdsZipReceivedForMultiply();

#ifdef DEC_DEBUG
    void printMatrixesST() const;
    void printTheMatrix() const;
    void printDenoMatrix() const;
#endif
  private:
    const ProcessorGroup &_group;
    const OverlapElementLocator& _locator;

    /**! Map of DAInt of cell identifiers. For a proc ID i,
     * gives an old2new map for the local part of the source mesh that has been sent to proc#i, just based on the
     * bounding box computation (this is potentially a larger set than what is finally in the interp matrix).
     * Second member gives proc ID.  */
    map < int, MCAuto<DataArrayInt> > _sent_src_ids;

    //! See _sent_src_ids. Same for target mesh.
    map < int, MCAuto<DataArrayInt> > _sent_trg_ids;

    /**! Vector of matrixes (partial interpolation ratios), result of the LOCAL interpolator run.
     * Indexing shared with _source_proc_id_st, and _target_proc_id_st.   */
    vector< vector< SparseDoubleVec > > _matrixes_st;
    //! See _matrixes_st - vec of source proc IDs
    vector< int > _source_proc_id_st;
    //! See _matrixes_st - vec of target proc IDs
    vector< int > _target_proc_id_st;

    /**! Number of received source mesh IDs at mesh data exchange.
     Counting the number of IDs suffices, as we just need this to prepare the receive side, when doing the final vector matrix multiplication.
     First dimension is the remote proc ID from which we received. */
    map <int, int > _nb_of_rcv_src_ids;

    /**! Specifies for each (target) remote proc ID (first dim of the map) the corresponding
     * source cell IDs to use.
     * This information is stored from the *locally* COMPuted matrices, and corresponds hence to field value that will need to
     * sent later on, if this matrix bit itself is sent aways.  */
    map<int, vector<int> > _src_ids_zip_comp;

    /**! Same idea as _src_ids_zip_comp above, but for RECEIVED matrix. */
    map<int, vector<int> > _src_ids_zip_recv;

    /**! THE matrix for matrix-vector product. The first dimension is indexed in the set of target procs
    * that interacts with local source mesh. The second dim is the target cell ID.
    * Same indexing as _the_matrix_st_source_proc_id and _the_deno_st.
    * We don't use a map here to be more efficient in the final matrix-vector computation which requires the joint
    * taversal of _the_matrix_st and _the_deno_st.
    * This matrix is filled after receival from other procs, contrary to _matrixes_st which contains local computations.*/
    vector< vector< SparseDoubleVec > > _the_matrix_st;
    //! See _the_matrix_st above. List of source proc IDs contributing to _the_matrix_st
    vector< int > _the_matrix_st_source_proc_id;
    // Denominators (computed from the numerator matrix). As for _the_matrix_st it is paired with _the_matrix_st_source_proc_id
    vector< vector< SparseDoubleVec > > _the_deno_st;

    //! Proc IDs to which data will be sent (originating this current proc) for matrix-vector computation
    vector< int > _proc_ids_to_send_vector_st;
  };
}

#endif
