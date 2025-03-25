// Copyright (C) 2007-2025  CEA, EDF
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
// Author : Anthony Geay (EDF R&D)

#include "MEDCouplingFieldDiscretizationOnNodesFE.hxx"
#include "InterpKernelDenseMatrix.hxx"
#include "InterpKernelGaussCoords.hxx"
#include "InterpKernelRootsMultiDim.hxx"
#include "InterpolationHelper.txx"
#include "MEDCouplingNormalizedUnstructuredMesh.txx"
#include "MEDCouplingUMesh.hxx"
#include <array>
#include <sstream>

using namespace MEDCoupling;

const char MEDCouplingFieldDiscretizationOnNodesFE::REPR[] = "FE";

std::string
MEDCouplingFieldDiscretizationOnNodesFE::getStringRepr() const {
  return std::string(REPR);
}

const char *
MEDCouplingFieldDiscretizationOnNodesFE::getRepr() const {
  return MEDCouplingFieldDiscretizationOnNodesFE::REPR;
}

void
MEDCouplingFieldDiscretizationOnNodesFE::reprQuickOverview(std::ostream &stream) const {
  stream << "NodeFE spatial discretization.";
}

MCAuto<MEDCouplingFieldDiscretization>
MEDCouplingFieldDiscretizationOnNodesFE::aggregate(std::vector<const MEDCouplingFieldDiscretization *> &fds) const {
  return EasyAggregate<MEDCouplingFieldDiscretizationOnNodesFE>(fds);
}

bool
MEDCouplingFieldDiscretizationOnNodesFE::isEqualIfNotWhy(const MEDCouplingFieldDiscretization *other, double eps, std::string &reason) const {
  if (!other) {
    reason = "other spatial discretization is NULL, and this spatial discretization (Node FE) is defined.";
    return false;
  }
  const MEDCouplingFieldDiscretizationOnNodesFE *otherC = dynamic_cast<const MEDCouplingFieldDiscretizationOnNodesFE *>(other);
  bool ret = otherC != 0;
  if (!ret)
    reason = "Spatial discrtization of this is ON_NODES_FE, which is not the case of other.";
  return ret;
}

/*!
 * This method is simply called by MEDCouplingFieldDiscretization::deepCopy. It performs the deep copy of \a this.
 *
 * \sa MEDCouplingFieldDiscretization::deepCopy.
 */
MEDCouplingFieldDiscretization *
MEDCouplingFieldDiscretizationOnNodesFE::clone() const {
  return new MEDCouplingFieldDiscretizationOnNodesFE;
}

void
MEDCouplingFieldDiscretizationOnNodesFE::checkCompatibilityWithNature(NatureOfField nat) const {
  if (nat != IntensiveMaximum)
    throw INTERP_KERNEL::Exception("Invalid nature for NodeFE field : expected IntensiveMaximum !");
}

MEDCouplingFieldDouble *
MEDCouplingFieldDiscretizationOnNodesFE::getMeasureField(const MEDCouplingMesh *mesh, bool isAbs) const {
  if (!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationOnNodesFE::getMeasureField : mesh instance specified is NULL !");
  throw INTERP_KERNEL::Exception("getMeasureField on MEDCouplingFieldDiscretizationOnNodesFE : not implemented yet !");
}

template <int SPACEDIM>
bool
IsInside(const MEDCouplingGaussLocalization &gl, const std::vector<double> &ptsInCell, const double locInReal[SPACEDIM], double locInRef[SPACEDIM]) {
  constexpr double EPS_IN_OUT = 1e-12;
  std::size_t nbPtsInCell(ptsInCell.size() / SPACEDIM);
  bool ret(false);
  std::vector<double> vini;
  const double *refCoo(gl.getRefCoords().data());
  INTERP_KERNEL::NormalizedCellType ct(gl.getType());

  class Functor {
  private:
    const MEDCouplingGaussLocalization *_gl;
    std::size_t _nb_pts_in_cell;
    const double *_pts_in_cell;
    const double *_point;

  public:
    Functor(const MEDCouplingGaussLocalization &gl, std::size_t nbPtsInCell, const double *ptsInCell, const double point[SPACEDIM])
        : _gl(&gl), _nb_pts_in_cell(nbPtsInCell), _pts_in_cell(ptsInCell), _point(point) {}
    std::vector<double>
    operator()(const std::vector<double> &x) {
      MEDCouplingGaussLocalization gl(_gl->getType(), _gl->getRefCoords(), x, {1.0});
      MCAuto<DataArrayDouble> shapeFunc = gl.getShapeFunctionValues();
      const double *shapeFuncPtr(shapeFunc->begin());
      std::vector<double> ret(SPACEDIM, 0);
      for (std::size_t iPt = 0; iPt < _nb_pts_in_cell; ++iPt) {
        for (short iDim = 0; iDim < SPACEDIM; ++iDim)
          ret[iDim] += shapeFuncPtr[iPt] * _pts_in_cell[SPACEDIM * iPt + iDim];
      }
      for (std::size_t ii = 0; ii < SPACEDIM; ii++)
        ret[ii] -= _point[ii];
      return ret;
    }
  };

  Functor func(gl, nbPtsInCell, ptsInCell.data(), locInReal);

  auto myJacobian = [&gl, nbPtsInCell, ptsInCell](const std::vector<double> &x, const std::vector<double> &, INTERP_KERNEL::DenseMatrix &jacobian) {
    MEDCouplingGaussLocalization mygl(gl.getType(), gl.getRefCoords(), x, {1.0});
    MCAuto<DataArrayDouble> shapeFunc = mygl.getDerivativeOfShapeFunctionValues();
    for (std::size_t i = 0; i < SPACEDIM; ++i)
      for (std::size_t j = 0; j < SPACEDIM; ++j) {
        double res = 0.0;
        for (std::size_t k = 0; k < nbPtsInCell; ++k)
          res += ptsInCell[k * SPACEDIM + i] * shapeFunc->getIJ(0, SPACEDIM * k + j);
        jacobian[i][j] = res;
      }
  };

  // loop on refcoords as initialization point for Newton algo. vini is the initialization vector of Newton.
  for (std::size_t attemptId = 0; attemptId < nbPtsInCell + 1; ++attemptId) {
    // TODO: optim, only use barycenter.
    // Moreover for robustesse begin with linear element or use split in simplex
    // to have a good initial guess.
    // Maybe need more test.
    if (attemptId == 0) {
      vini = INTERP_KERNEL::GaussInfo::GetReferenceCoordinatesOfBarycenterOf(ct, refCoo);
    } else {
      vini = std::vector<double>(refCoo + (attemptId - 1) * SPACEDIM, refCoo + (attemptId)*SPACEDIM);
    }

    try {
      bool check(true);
      // if not converge with Newton after 50 iterations, go to the next pt
      INTERP_KERNEL::SolveWithNewtonWithJacobian(vini, check, func, myJacobian, EPS_IN_OUT, 50);
      ret =
          (check == false); // looks strange but OK regarding newt (SolveWithNewton) at page 387 of numerical recipes for semantic of check parameter
    } catch (INTERP_KERNEL::Exception &ex) {
      ret = false;
    } // Something get wrong during Newton process
    if (ret) { // Newton has converged. Now check if it converged to a point inside cell
      // TODO: add in a class method ?
      if (INTERP_KERNEL::GaussInfo::IsInOrOutForReference(ct, refCoo, vini.data(), EPS_IN_OUT)) {
        // converged but locInReal has been detected outside of cell
        for (int i = 0; i < SPACEDIM; i++) {
          locInRef[i] = vini[i];
        }
        return true;
      }
    }
  }
  std::fill(locInRef, locInRef + SPACEDIM, std::numeric_limits<double>::max());
  return false;
}

/* Search closest point on the surface unsing orthogonal projection.
 *  See mmnewt.F90 routine of code_aster
 */
template <int SPACEDIM>
bool
ClosestPoint(const MEDCouplingGaussLocalization &gl, const std::vector<double> &ptsInCell, const double locInReal[SPACEDIM],
             double locInRef[SPACEDIM], double &dist) {
  constexpr double EPS_IN_OUT = 1e-12;
  std::size_t nbPtsInCell(ptsInCell.size() / SPACEDIM);
  bool ret(false);
  std::vector<double> vini;
  const double *refCoo(gl.getRefCoords().data());
  INTERP_KERNEL::NormalizedCellType ct(gl.getType());

  auto tangent = [nbPtsInCell, ptsInCell](const MCAuto<DataArrayDouble> &dshapeFunc) {
    std::array<std::array<double, SPACEDIM>, SPACEDIM - 1> tang;

    for (std::size_t j = 0; j < SPACEDIM - 1; ++j) {
      for (std::size_t i = 0; i < SPACEDIM; ++i) {
        tang[j][i] = 0.0;
      }
    }

    for (std::size_t k = 0; k < nbPtsInCell; ++k) {
      for (std::size_t j = 0; j < SPACEDIM - 1; ++j) {
        for (std::size_t i = 0; i < SPACEDIM; ++i) {
          tang[j][i] += ptsInCell[k * SPACEDIM + i] * dshapeFunc->getIJ(0, (SPACEDIM - 1) * k + j);
        }
      }
    }

    return tang;
  };

  auto dist_pt = [nbPtsInCell, ptsInCell, locInReal](const MCAuto<DataArrayDouble> &shapeFunc) {
    const double *shapeFuncPtr(shapeFunc->begin());
    std::vector<double> posi(SPACEDIM, 0);
    for (std::size_t iPt = 0; iPt < nbPtsInCell; ++iPt) {
      for (short iDim = 0; iDim < SPACEDIM; ++iDim)
        posi[iDim] += shapeFuncPtr[iPt] * ptsInCell[SPACEDIM * iPt + iDim];
    }
    for (std::size_t ii = 0; ii < SPACEDIM; ii++) {
      posi[ii] -= locInReal[ii];
    }

    return posi;
  };

  auto distance = [&gl, &dist_pt](const std::vector<double> &x) {
    MEDCouplingGaussLocalization mygl(gl.getType(), gl.getRefCoords(), x, {1.0});
    MCAuto<DataArrayDouble> shapeFunc = mygl.getShapeFunctionValues();

    const auto _pt = dist_pt(shapeFunc);

    double _dist = 0.;
    for (int i = 0; i < SPACEDIM; i++) {
      _dist += _pt[i] * _pt[i];
    }
    return std::sqrt(_dist);
  };

  auto myResidual = [&gl, nbPtsInCell, ptsInCell, locInReal, &tangent, &dist_pt](const std::vector<double> &x) {
    MEDCouplingGaussLocalization mygl(gl.getType(), gl.getRefCoords(), x, {1.0});
    MCAuto<DataArrayDouble> shapeFunc = mygl.getShapeFunctionValues();
    MCAuto<DataArrayDouble> dshapeFunc = mygl.getDerivativeOfShapeFunctionValues();

    const auto tang = tangent(dshapeFunc);

    std::vector<double> posi = dist_pt(shapeFunc);

    std::vector<double> ret(SPACEDIM - 1, 0);
    for (std::size_t jj = 0; jj < SPACEDIM - 1; jj++) {
      ret[jj] = 0.0;
      for (std::size_t ii = 0; ii < SPACEDIM; ii++) {
        ret[jj] += posi[ii] * tang[jj][ii];
      }
    }

    return ret;
  };

  auto myJacobian = [&gl, &tangent](const std::vector<double> &x, const std::vector<double> &, INTERP_KERNEL::DenseMatrix &jacobian) {
    MEDCouplingGaussLocalization mygl(gl.getType(), gl.getRefCoords(), x, {1.0});
    MCAuto<DataArrayDouble> dshapeFunc = mygl.getDerivativeOfShapeFunctionValues();

    const auto tang = tangent(dshapeFunc);

    for (std::size_t i = 0; i < SPACEDIM - 1; ++i) {
      for (std::size_t j = i; j < SPACEDIM - 1; ++j) {
        double res = 0.0;
        // No local curvature for the moment
        for (std::size_t k = 0; k < SPACEDIM; ++k)
          res += tang[i][k] * tang[j][k];
        jacobian[i][j] = res;
        jacobian[j][i] = res;
      }
    }
  };

  bool find_one = false;
  std::vector<double> vtmp;
  // loop on refcoords as initialization point for Newton algo. vini is the initialization vector of Newton.
  for (std::size_t attemptId = 0; attemptId < nbPtsInCell + 1; ++attemptId) {
    if (attemptId == 0) {
      vini = INTERP_KERNEL::GaussInfo::GetReferenceCoordinatesOfBarycenterOf(ct, refCoo);
    } else {
      vini = std::vector<double>(refCoo + (attemptId - 1) * (SPACEDIM - 1), refCoo + (attemptId) * (SPACEDIM - 1));
    }

    try {
      bool check(true);
      // if not converge with Newton after 50 iterations, go to the next pt
      INTERP_KERNEL::SolveWithNewtonWithJacobian(vini, check, myResidual, myJacobian, EPS_IN_OUT, 50);
      ret = (check == false); // looks strange but OK regarding newt (SolveWithNewton) at page 387 of numerical recipes for semantic of check
                              // parameter
    } catch (INTERP_KERNEL::Exception &ex) {
      ret = false;
    } // Something get wrong during Newton process
    if (ret) { // Newton has converged. Now check if it converged to a point inside cell
      if (INTERP_KERNEL::GaussInfo::IsInOrOutForReference(ct, refCoo, vini.data(), EPS_IN_OUT)) {
        // converged but locInReal has been detected outside of cell
        for (int i = 0; i < (SPACEDIM - 1); i++) {
          locInRef[i] = vini[i];
        }
        dist = distance(vini);
        return true;
      }
      if (!find_one) {
        find_one = true;
        vtmp = vini;
      }
    }
  }
  // No point inside but at least one on the border
  if (find_one) {
    INTERP_KERNEL::GaussInfo::AdapatCoorForReference(ct, refCoo, vtmp.data());
    for (int i = 0; i < (SPACEDIM - 1); i++) {
      locInRef[i] = vtmp[i];
    }
    dist = distance(vtmp);
    return true;
  }
  std::fill(locInRef, locInRef + (SPACEDIM - 1), std::numeric_limits<double>::max());
  dist = std::numeric_limits<double>::max();
  return false;
}

void
MEDCouplingFieldDiscretizationOnNodesFE::getValueOn(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *ptsCoo,
                                                    double *res) const {
  MCAuto<DataArrayDouble> res2(this->getValueOnMulti(arr, mesh, ptsCoo, 1));
  std::copy(res2->begin(), res2->end(), res);
}

void
MEDCouplingFieldDiscretizationOnNodesFE::GetRefCoordOfListOfPts(
    const MEDCouplingUMesh *umesh, const double *ptsCoo, mcIdType nbOfPts,
    std::function<void(const MEDCouplingGaussLocalization &, const std::vector<mcIdType> &)> customFunc) {
  switch (umesh->getSpaceDimension()) {
  case 1:
    MEDCouplingFieldDiscretizationOnNodesFE::GetRefCoordOfListOfNDPtsInND<1>(umesh, ptsCoo, nbOfPts, customFunc);
    break;

  case 2:
    MEDCouplingFieldDiscretizationOnNodesFE::GetRefCoordOfListOfNDPtsInND<2>(umesh, ptsCoo, nbOfPts, customFunc);
    break;

  case 3:
    MEDCouplingFieldDiscretizationOnNodesFE::GetRefCoordOfListOfNDPtsInND<3>(umesh, ptsCoo, nbOfPts, customFunc);
    break;
  default:
    THROW_IK_EXCEPTION("GetRefCoordOfListOfPts : invalid dimension !")
    break;
  }
}

void
MEDCouplingFieldDiscretizationOnNodesFE::GetClosestRefCoordOfListOfPts(
    const MEDCouplingUMesh *umesh, const double *ptsCoo, mcIdType nbOfPts,
    std::function<void(const MEDCouplingGaussLocalization &, const std::vector<mcIdType> &)> customFunc) {
  switch (umesh->getSpaceDimension()) {
  case 2:
    MEDCouplingFieldDiscretizationOnNodesFE::GetClosestRefCoordOfListOfNDPtsInND<2>(umesh, ptsCoo, nbOfPts, customFunc);
    break;

  case 3:
    MEDCouplingFieldDiscretizationOnNodesFE::GetClosestRefCoordOfListOfNDPtsInND<3>(umesh, ptsCoo, nbOfPts, customFunc);
    break;
  default:
    THROW_IK_EXCEPTION("GetClosestRefCoordOfListOfPts : invalid dimension !")
    break;
  }
}

template <int SPACEDIM>
void
MEDCouplingFieldDiscretizationOnNodesFE::GetRefCoordOfListOfNDPtsInND(
    const MEDCouplingUMesh *umesh, const double *ptsCoo, mcIdType nbOfPts,
    std::function<void(const MEDCouplingGaussLocalization &, const std::vector<mcIdType> &)> customFunc) {
  if (SPACEDIM != umesh->getSpaceDimension()) {
    THROW_IK_EXCEPTION("GetRefCoordOfListOfNDPtsInND : mesh and parameter have a different dimension !")
  }
  MEDCouplingNormalizedUnstructuredMesh<SPACEDIM, SPACEDIM> mesh_wrapper(umesh);
  auto tree = INTERP_KERNEL::BuildBBTree<MEDCouplingNormalizedUnstructuredMesh<SPACEDIM, SPACEDIM>, SPACEDIM>(mesh_wrapper);
  for (mcIdType iPt = 0; iPt < nbOfPts; ++iPt) {
    const std::vector<double> gsCoo(ptsCoo + iPt * SPACEDIM, ptsCoo + (iPt + 1) * SPACEDIM);
    GetRefCoordOfListOf1PtInND(umesh, tree, gsCoo, customFunc);
  }
}

template <typename BBTree>
void
MEDCouplingFieldDiscretizationOnNodesFE::GetRefCoordOfListOf1PtInND(
    const MEDCouplingUMesh *umesh, const BBTree &tree, const std::vector<double> &ptCoor,
    std::function<void(const MEDCouplingGaussLocalization &, const std::vector<mcIdType> &)> customFunc) {

  constexpr int SPACEDIM = BBTree::dimension;
  const double *coordsOfMesh(umesh->getCoords()->begin());
  std::vector<mcIdType> elems;
  tree.getElementsAroundPoint(ptCoor.data(), elems);
  for( mcIdType cellId : elems )
  {
    INTERP_KERNEL::NormalizedCellType gt(umesh->getTypeOfCell(cellId));
    std::vector<mcIdType> conn;
    umesh->getNodeIdsOfCell(cellId, conn);
    MCAuto<DataArrayDouble> refCoo(MEDCouplingGaussLocalization::GetDefaultReferenceCoordinatesOf(gt));
    std::vector<double> refCooCpp(refCoo->begin(), refCoo->end());
    MEDCouplingGaussLocalization gl(gt, refCooCpp, std::vector<double>(SPACEDIM, 0.), {1.});
    std::vector<double> ptsInCell;
    ptsInCell.reserve(conn.size() * gl.getDimension());
    std::for_each(conn.cbegin(), conn.cend(), [SPACEDIM, coordsOfMesh, &ptsInCell](mcIdType c) {
      ptsInCell.insert(ptsInCell.end(), coordsOfMesh + c * SPACEDIM, coordsOfMesh + (c + 1) * SPACEDIM);
    });
    std::vector<double> locInRef(SPACEDIM);
    if (IsInside<SPACEDIM>(gl, ptsInCell, ptCoor.data(), locInRef.data())) {
      gl.setGaussCoords(locInRef);
      customFunc(gl, conn);
      // By continuity of the field, we can stop here.
      return;
    }
  }

  std::array<std::string, 3> dir = {"X", "Y", "Z"};
  std::string mess;
  for (mcIdType i = 0; i < SPACEDIM; i++) {
    mess += " " + dir[i] + "=" + std::to_string(ptCoor[i]);
  }
  THROW_IK_EXCEPTION("GetRefCoordOfListOf1PtInND on MEDCouplingFieldDiscretizationOnNodesFE : fail to locate point " << mess << " !");
}

template <int SPACEDIM>
void
MEDCouplingFieldDiscretizationOnNodesFE::GetClosestRefCoordOfListOfNDPtsInND(
    const MEDCouplingUMesh *umesh, const double *ptsCoo, mcIdType nbOfPts,
    std::function<void(const MEDCouplingGaussLocalization &, const std::vector<mcIdType> &)> customFunc) {
  if (SPACEDIM != umesh->getSpaceDimension()) {
    THROW_IK_EXCEPTION("GetClosestRefCoordOfListOfNDPtsInND : mesh and parameter have a different dimension !")
  }
  MEDCouplingNormalizedUnstructuredMesh<SPACEDIM, SPACEDIM - 1> mesh_wrapper(umesh);
  auto tree = INTERP_KERNEL::BuildBBTree<MEDCouplingNormalizedUnstructuredMesh<SPACEDIM, SPACEDIM - 1>, SPACEDIM>(mesh_wrapper);
  for (mcIdType iPt = 0; iPt < nbOfPts; ++iPt) {
    const std::vector<double> gsCoo(ptsCoo + iPt * SPACEDIM, ptsCoo + (iPt + 1) * SPACEDIM);
    GetClosestRefCoordOfListOf1PtInND(umesh, tree, gsCoo, customFunc);
  }
}

template <typename BBTree>
void
MEDCouplingFieldDiscretizationOnNodesFE::GetClosestRefCoordOfListOf1PtInND(
    const MEDCouplingUMesh *umesh, const BBTree &tree, const std::vector<double> &ptCoor,
    std::function<void(const MEDCouplingGaussLocalization &, const std::vector<mcIdType> &)> customFunc, const double dist_max) {

  constexpr int SPACEDIM = BBTree::dimension;
  const double *coordsOfMesh(umesh->getCoords()->begin());
  // EDF31461 : TODO : An algorithm could be implemented to automatically reduce number of candidates.
  //std::vector<mcIdType> elems;
  //tree.getElementsAroundPoint(ptCoor.data(), elems);
  bool found = false;
  double dist_min = std::numeric_limits<double>::max(), dist_loc;
  MEDCouplingGaussLocalization gl_min(INTERP_KERNEL::NORM_SEG2);
  mcIdType cell_id;

  //for( mcIdType cellId : elems ) // 
  for (mcIdType cellId = 0; cellId < umesh->getNumberOfCells(); ++cellId)
  {
    INTERP_KERNEL::NormalizedCellType gt(umesh->getTypeOfCell(cellId));
    std::vector<mcIdType> conn;
    umesh->getNodeIdsOfCell(cellId, conn);
    MCAuto<DataArrayDouble> refCoo(MEDCouplingGaussLocalization::GetDefaultReferenceCoordinatesOf(gt));
    std::vector<double> refCooCpp(refCoo->begin(), refCoo->end());
    MEDCouplingGaussLocalization gl(gt, refCooCpp, std::vector<double>(SPACEDIM - 1, 0.), {1.});
    std::vector<double> ptsInCell;
    ptsInCell.reserve(conn.size() * SPACEDIM);
    std::for_each(conn.cbegin(), conn.cend(), [SPACEDIM, coordsOfMesh, &ptsInCell](mcIdType c) {
      ptsInCell.insert(ptsInCell.end(), coordsOfMesh + c * SPACEDIM, coordsOfMesh + (c + 1) * SPACEDIM);
    });
    std::vector<double> locInRef(SPACEDIM - 1);
    if (ClosestPoint<SPACEDIM>(gl, ptsInCell, ptCoor.data(), locInRef.data(), dist_loc)) {
      if (dist_loc <= dist_max && dist_loc < dist_min) {
        found = true;
        cell_id = cellId;
        gl_min = gl;
        gl_min.setGaussCoords(locInRef);
        dist_min = dist_loc;
      }
    }
  }
  if (found) {
    std::vector<mcIdType> conn;
    umesh->getNodeIdsOfCell(cell_id, conn);
    customFunc(gl_min, conn);
  } else {
    std::array<std::string, 3> dir = {"X", "Y", "Z"};
    std::string mess;
    for (mcIdType i = 0; i < SPACEDIM; i++) {
      mess += " " + dir[i] + "=" + std::to_string(ptCoor[i]);
    }
    THROW_IK_EXCEPTION("GetClosestRefCoordOfListOf1PtInND on MEDCouplingFieldDiscretizationOnNodesFE : fail to locate point " << mess << " !");
  }
}

const MEDCouplingUMesh *
MEDCouplingFieldDiscretizationOnNodesFE::checkConfig(const MEDCouplingMesh *mesh) const {
  const MEDCouplingUMesh *umesh(dynamic_cast<const MEDCouplingUMesh *>(mesh));
  if (!umesh)
    THROW_IK_EXCEPTION("checkConfig : not implemented yet for type != MEDCouplingUMesh !");
  if (umesh->getSpaceDimension() != umesh->getMeshDimension())
    THROW_IK_EXCEPTION("checkConfig : implemented only for meshes with spacedim == meshdim !");
  umesh->checkConsistency();
  return umesh;
}

const MEDCouplingUMesh *
MEDCouplingFieldDiscretizationOnNodesFE::checkConfigSurf(const MEDCouplingMesh *mesh) const {
  const MEDCouplingUMesh *umesh(dynamic_cast<const MEDCouplingUMesh *>(mesh));
  if (!umesh)
    THROW_IK_EXCEPTION("checkConfigSurf : not implemented yet for type != MEDCouplingUMesh !");
  if (umesh->getSpaceDimension() - 1 != umesh->getMeshDimension())
    THROW_IK_EXCEPTION("checkConfigSurf : implemented only for meshes with spacedim - 1 = meshdim !");
  umesh->checkConsistency();
  return umesh;
}

DataArrayDouble *
MEDCouplingFieldDiscretizationOnNodesFE::getValueOnMulti(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc,
                                                         mcIdType nbOfTargetPoints) const {
  if (!arr || !arr->isAllocated())
    throw INTERP_KERNEL::Exception("getValueOnMulti : input array is null or not allocated !");
  mcIdType nbOfRows = getNumberOfMeshPlaces(mesh);
  const MEDCouplingUMesh *umesh = checkConfig(mesh);
  std::size_t nbCompo(arr->getNumberOfComponents());
  MCAuto<DataArrayDouble> ret(DataArrayDouble::New());
  ret->alloc(nbOfTargetPoints, nbCompo);
  double *res(ret->getPointer());

  auto arrayFeeder = [arr, &res, nbCompo](const MEDCouplingGaussLocalization &gl, const std::vector<mcIdType> &conn) {
    MCAuto<DataArrayDouble> resVector(gl.getShapeFunctionValues());
    {
      std::for_each(res, res + nbCompo, [](double &v) { v = 0.0; });
      for (std::size_t iComp = 0; iComp < nbCompo; ++iComp)
        for (int iPt = 0; iPt < gl.getNumberOfPtsInRefCell(); ++iPt) {
          {
            res[iComp] += resVector->getIJ(0, iPt) * arr->getIJ(conn[iPt], iComp);
          }
        }
      res += nbCompo;
    }
  };

  GetRefCoordOfListOfPts(umesh, loc, nbOfTargetPoints, arrayFeeder);
  return ret.retn();
}

/*!
 * Returns for each \a nbOfPoints point in \a loc its coordinate in reference element.
 */
MCAuto<DataArrayDouble>
MEDCouplingFieldDiscretizationOnNodesFE::getCooInRefElement(const MEDCouplingMesh *mesh, const double *loc, mcIdType nbOfPoints) const {
  const MEDCouplingUMesh *umesh = checkConfig(mesh);
  MCAuto<DataArrayDouble> ret(DataArrayDouble::New());
  const mcIdType dim = mesh->getSpaceDimension();
  ret->alloc(nbOfPoints, dim);
  double *retPtr(ret->getPointer());

  auto arrayFeeder = [&retPtr, &dim](const MEDCouplingGaussLocalization &gl, const std::vector<mcIdType> &conn) {
    std::vector<double> resVector(gl.getGaussCoords());
    {
      std::copy(resVector.begin(), resVector.end(), retPtr);
      retPtr += dim;
    }
  };

  GetRefCoordOfListOfPts(umesh, loc, nbOfPoints, arrayFeeder);
  return ret;
}

/*!
 * Returns for each \a nbOfPoints point in \a loc its coordinate in reference element.
 */
MCAuto<DataArrayDouble>
MEDCouplingFieldDiscretizationOnNodesFE::getClosestCooInRefElement(const MEDCouplingMesh *mesh, const double *loc, mcIdType nbOfPoints) const {
  const MEDCouplingUMesh *umesh = checkConfigSurf(mesh);
  MCAuto<DataArrayDouble> ret(DataArrayDouble::New());
  const mcIdType dim = mesh->getSpaceDimension();
  ret->alloc(nbOfPoints, dim - 1);
  double *retPtr(ret->getPointer());

  auto arrayFeeder = [&retPtr, &dim](const MEDCouplingGaussLocalization &gl, const std::vector<mcIdType> &conn) {
    std::vector<double> resVector(gl.getGaussCoords());
    {
      std::copy(resVector.begin(), resVector.end(), retPtr);
      retPtr += dim - 1;
    }
  };

  GetClosestRefCoordOfListOfPts(umesh, loc, nbOfPoints, arrayFeeder);
  return ret;
}

void
MEDCouplingFieldDiscretizationOnNodesFE::computeCrudeMatrix(const MEDCouplingUMesh *srcMesh, const double *ptsCoo, const mcIdType nbOfPts,
                                                            std::vector<std::map<mcIdType, double>> &matrix, const INTERP_KERNEL::FEInterpolationOptions& option)
{
  switch (srcMesh->getSpaceDimension()) {
  case 1:
    MEDCouplingFieldDiscretizationOnNodesFE::computeCrudeMatrixNd<1>(srcMesh, ptsCoo, nbOfPts, matrix, option);
    break;

  case 2:
    MEDCouplingFieldDiscretizationOnNodesFE::computeCrudeMatrixNd<2>(srcMesh, ptsCoo, nbOfPts, matrix, option);
    break;

  case 3:
    MEDCouplingFieldDiscretizationOnNodesFE::computeCrudeMatrixNd<3>(srcMesh, ptsCoo, nbOfPts, matrix, option);
    break;
  default:
    THROW_IK_EXCEPTION("computeCrudeMatrix : invalid dimension !")
    break;
  }
}

template <int SPACEDIM>
void
MEDCouplingFieldDiscretizationOnNodesFE::computeCrudeMatrixNd(const MEDCouplingUMesh *srcMesh, const double *ptsCoo, const mcIdType nbOfPts,
                                                              std::vector<std::map<mcIdType, double>> &matrix, const INTERP_KERNEL::FEInterpolationOptions& option)
{
  // This code is largely duplicated
  if (SPACEDIM != srcMesh->getSpaceDimension()) {
    THROW_IK_EXCEPTION("computeCrudeMatrixNd : mesh and parameter have a different dimension !")
  }

  // Prepare matrix
  matrix.clear();
  matrix.resize(nbOfPts);
  mcIdType rowId = 0;

  auto matrixFeeder = [&matrix, &rowId](const MEDCouplingGaussLocalization &gl, const std::vector<mcIdType> &conn) {
    auto &row = matrix[rowId++];
    const MCAuto<DataArrayDouble> resVector(gl.getShapeFunctionValues());
    for (int iPt = 0; iPt < gl.getNumberOfPtsInRefCell(); ++iPt) {
      row[conn[iPt]] = resVector->getIJ(0, iPt);
    }
  };
  // This projection is named COLLOCATION in PROJ_CHAMP in code_aster.

  const double dist_max = option.getProjectionMaxDistance();
  const bool projOnSurf( option.getProjectionOnSurfStatus() ), use_dist_max( option.getMaxDistanceStatus() );

  if (srcMesh->getMeshDimension() == srcMesh->getSpaceDimension()) {
    // Compute bounding box
    MEDCouplingNormalizedUnstructuredMesh<SPACEDIM, SPACEDIM> mesh_wrapper(srcMesh);
    auto tree = INTERP_KERNEL::BuildBBTree<MEDCouplingNormalizedUnstructuredMesh<SPACEDIM, SPACEDIM>, SPACEDIM>(mesh_wrapper);

    MEDCouplingUMesh *srcSkinMesh = nullptr;
    if (projOnSurf) {
      srcSkinMesh = srcMesh->computeSkin();
    }

    for (mcIdType iPt = 0; iPt < nbOfPts; ++iPt) {
      const std::vector<double> gsCoo(ptsCoo + iPt * SPACEDIM, ptsCoo + (iPt + 1) * SPACEDIM);
      try {
        GetRefCoordOfListOf1PtInND(srcMesh, tree, gsCoo, matrixFeeder);
      } catch (INTERP_KERNEL::Exception & /*e*/) {
        if (projOnSurf) {
          try {
            GetClosestRefCoordOfListOf1PtInND(srcSkinMesh, tree, gsCoo, matrixFeeder, dist_max);
          } catch (INTERP_KERNEL::Exception & /*e*/) {
            // Nothing to do since the point is to far - go to next line.
            if (use_dist_max) {
              rowId++;
            } else {
              throw;
            }
          }
        } else {
          // Nothing to do since point is not inside volume - go to next line.
          rowId++;
        }
      }
    }
  } else if (srcMesh->getMeshDimension() == srcMesh->getSpaceDimension() - 1) {
    // Compute bounding box
    MEDCouplingNormalizedUnstructuredMesh<SPACEDIM, SPACEDIM - 1> mesh_wrapper(srcMesh);
    auto tree = INTERP_KERNEL::BuildBBTree<MEDCouplingNormalizedUnstructuredMesh<SPACEDIM, SPACEDIM - 1>, SPACEDIM>(mesh_wrapper);
    for (mcIdType iPt = 0; iPt < nbOfPts; ++iPt) {
      const std::vector<double> gsCoo(ptsCoo + iPt * SPACEDIM, ptsCoo + (iPt + 1) * SPACEDIM);
      try {
        GetClosestRefCoordOfListOf1PtInND(srcMesh, tree, gsCoo, matrixFeeder, dist_max);
      } catch (INTERP_KERNEL::Exception & /*e*/) {
        // Nothing to do since the point is to far - go to next line.
        if (use_dist_max) {
          rowId++;
        } else {
          throw;
        }
      }
    }
  } else {
    THROW_IK_EXCEPTION("computeCrudeMatrixNd : incompatible dimension between space and mesh !")
  }
}
