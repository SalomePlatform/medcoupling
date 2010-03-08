//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
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
#ifndef __CELLMODEL_INTERP_KERNEL_HXX__
#define __CELLMODEL_INTERP_KERNEL_HXX__

#include "INTERPKERNELDefines.hxx"

#include "NormalizedUnstructuredMesh.hxx"

#include <map>

namespace INTERP_KERNEL
{
  /*!
   * This class descibes all static elements (different from polygons and polyhedron) 3D, 2D and 1D.
   */
  class INTERPKERNEL_EXPORT CellModel
  {
  public:
    static const unsigned MAX_NB_OF_SONS=6;
    static const unsigned MAX_NB_OF_NODES_PER_ELEM=30;
  private:
    CellModel(NormalizedCellType type);
    static void buildUniqueInstance();
  public:
    static const CellModel& getCellModel(NormalizedCellType type);
    bool isDynamic() const { return _dyn; }
    bool isQuadratic() const { return _quadratic; }
    unsigned getDimension() const { return _dim; }
    bool isCompatibleWith(NormalizedCellType type) const;
    //! sonId is in C format.
    const unsigned *getNodesConstituentTheSon(unsigned sonId) const { return _sons_con[sonId]; }
    unsigned getNumberOfNodes() const { return _nb_of_pts; }
    unsigned getNumberOfSons() const { return _nb_of_sons; }
    unsigned getNumberOfSons2(const int *conn, int lgth) const;
    unsigned getNumberOfNodesConstituentTheSon(unsigned sonId) const { return _nb_of_sons_con[sonId]; }
    unsigned getNumberOfNodesConstituentTheSon2(unsigned sonId, const int *nodalConn, int lgth) const;
    NormalizedCellType getExtrudedType() const { return _extruded_type; }
    NormalizedCellType getSonType(unsigned sonId) const { return _sons_type[sonId]; }
    NormalizedCellType getSonType2(unsigned sonId) const;
    unsigned fillSonCellNodalConnectivity(int sonId, const int *nodalConn, int *sonNodalConn) const;
    unsigned fillSonCellNodalConnectivity2(int sonId, const int *nodalConn, int lgth, int *sonNodalConn, NormalizedCellType& typeOfSon) const;
  private:
    bool _dyn;
    bool _quadratic;
    unsigned _dim;
    unsigned _nb_of_pts;
    unsigned _nb_of_sons;
    NormalizedCellType _type;
    NormalizedCellType _extruded_type;
    unsigned _sons_con[MAX_NB_OF_SONS][MAX_NB_OF_NODES_PER_ELEM];
    unsigned _nb_of_sons_con[MAX_NB_OF_SONS];
    NormalizedCellType _sons_type[MAX_NB_OF_SONS];
    static std::map<NormalizedCellType,CellModel> _map_of_unique_instance;
  };
}

#endif
