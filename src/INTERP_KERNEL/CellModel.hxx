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

#ifndef __CELLMODEL_INTERP_KERNEL_HXX__
#define __CELLMODEL_INTERP_KERNEL_HXX__

#include "INTERPKERNELDefines.hxx"

#include "NormalizedUnstructuredMesh.hxx"

#include <map>

namespace INTERP_KERNEL
{
  class DiameterCalculator;
  class OrientationInverter;
  
  /*!
   * This class describes all static elements (different from polygons and polyhedron) 3D, 2D and 1D.
   */
  class CellModel
  {
  public:
    static const unsigned MAX_NB_OF_SONS=8;
    static const unsigned MAX_NB_OF_NODES_PER_ELEM=30;
    static const unsigned MAX_NB_OF_LITTLE_SONS=12;
  private:
    CellModel(NormalizedCellType type);
    static void buildUniqueInstance();
  public:
    INTERPKERNEL_EXPORT static const CellModel& GetCellModel(NormalizedCellType type);
    INTERPKERNEL_EXPORT NormalizedCellType getEnum() const { return _type; }
    INTERPKERNEL_EXPORT const char *getRepr() const;
    INTERPKERNEL_EXPORT bool isExtruded() const { return _is_extruded; }
    INTERPKERNEL_EXPORT bool isDynamic() const { return _dyn; }
    INTERPKERNEL_EXPORT bool isQuadratic() const { return _quadratic; }
    INTERPKERNEL_EXPORT unsigned getDimension() const { return _dim; }
    INTERPKERNEL_EXPORT bool isCompatibleWith(NormalizedCellType type) const;
    INTERPKERNEL_EXPORT bool isSimplex() const { return _is_simplex; }
    //! sonId is in C format.
    INTERPKERNEL_EXPORT const unsigned *getNodesConstituentTheSon(unsigned sonId) const { return _sons_con[sonId]; }
    INTERPKERNEL_EXPORT const unsigned *getNodesConstituentTheLittleSon(unsigned littleSonId) const { return _little_sons_con[littleSonId]; }
    INTERPKERNEL_EXPORT bool getOrientationStatus(unsigned lgth, const int *conn1, const int *conn2) const;
    INTERPKERNEL_EXPORT unsigned getNumberOfNodes() const { return _nb_of_pts; }
    INTERPKERNEL_EXPORT unsigned getNumberOfSons() const { return _nb_of_sons; }
    INTERPKERNEL_EXPORT unsigned getNumberOfSons2(const int *conn, int lgth) const;
    INTERPKERNEL_EXPORT unsigned getNumberOfEdgesIn3D(const int *conn, int lgth) const;
    INTERPKERNEL_EXPORT unsigned getNumberOfMicroEdges() const;
    INTERPKERNEL_EXPORT unsigned getNumberOfNodesConstituentTheSon(unsigned sonId) const { return _nb_of_sons_con[sonId]; }
    INTERPKERNEL_EXPORT unsigned getNumberOfNodesConstituentTheSon2(unsigned sonId, const int *nodalConn, int lgth) const;
    INTERPKERNEL_EXPORT NormalizedCellType getExtrudedType() const { return _extruded_type; }
    INTERPKERNEL_EXPORT NormalizedCellType getCorrespondingPolyType() const;
    INTERPKERNEL_EXPORT NormalizedCellType getReverseExtrudedType() const { return _reverse_extruded_type; }
    INTERPKERNEL_EXPORT NormalizedCellType getLinearType() const { return _linear_type; }
    INTERPKERNEL_EXPORT NormalizedCellType getQuadraticType() const { return _quadratic_type; }
    INTERPKERNEL_EXPORT NormalizedCellType getQuadraticType2() const { return _quadratic_type2; }
    INTERPKERNEL_EXPORT NormalizedCellType getSonType(unsigned sonId) const { return _sons_type[sonId]; }
    INTERPKERNEL_EXPORT NormalizedCellType getSonType2(unsigned sonId) const;
    INTERPKERNEL_EXPORT unsigned fillSonCellNodalConnectivity(int sonId, const int *nodalConn, int *sonNodalConn) const;
    INTERPKERNEL_EXPORT unsigned fillSonCellNodalConnectivity2(int sonId, const int *nodalConn, int lgth, int *sonNodalConn, NormalizedCellType& typeOfSon) const;
    INTERPKERNEL_EXPORT unsigned fillSonCellNodalConnectivity4(int sonId, const int *nodalConn, int lgth, int *sonNodalConn, NormalizedCellType& typeOfSon) const;
    INTERPKERNEL_EXPORT unsigned fillSonEdgesNodalConnectivity3D(int sonId, const int *nodalConn, int lgth, int *sonNodalConn, NormalizedCellType& typeOfSon) const;
    INTERPKERNEL_EXPORT unsigned fillMicroEdgeNodalConnectivity(int sonId, const int *nodalConn, int *sonNodalConn, NormalizedCellType& typeOfSon) const;
    INTERPKERNEL_EXPORT void changeOrientationOf2D(int *nodalConn, unsigned int sz) const;
    INTERPKERNEL_EXPORT void changeOrientationOf1D(int *nodalConn, unsigned int sz) const;
    INTERPKERNEL_EXPORT DiameterCalculator *buildInstanceOfDiameterCalulator(int spaceDim) const;
    INTERPKERNEL_EXPORT OrientationInverter *buildOrientationInverter() const;
  private:
    bool _dyn;
    bool _quadratic;
    bool _is_simplex;
    bool _is_extruded;
    unsigned _dim;
    unsigned _nb_of_pts;
    unsigned _nb_of_sons;
    unsigned _nb_of_little_sons;
    NormalizedCellType _type;
    NormalizedCellType _extruded_type;
    NormalizedCellType _reverse_extruded_type;
    NormalizedCellType _linear_type;
    NormalizedCellType _quadratic_type;
    NormalizedCellType _quadratic_type2;
    unsigned _sons_con[MAX_NB_OF_SONS][MAX_NB_OF_NODES_PER_ELEM];
    unsigned _little_sons_con[MAX_NB_OF_LITTLE_SONS][3];
    unsigned _nb_of_sons_con[MAX_NB_OF_SONS];
    NormalizedCellType _sons_type[MAX_NB_OF_SONS];
    static std::map<NormalizedCellType,CellModel> _map_of_unique_instance;
    static const char *CELL_TYPES_REPR[];
  };
}

#endif
