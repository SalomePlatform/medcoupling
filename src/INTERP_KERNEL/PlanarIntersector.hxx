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
#ifndef _PLANAR_INTERSECTOR_HXX_
#define _PLANAR_INTERSECTOR_HXX_

#include <vector>

namespace INTERP_KERNEL
{
  class TranslationRotationMatrix;
  
  template<class MyMeshType>
  class INTERPKERNEL_EXPORT PlanarIntersector
  {
  public:
    static const int SPACEDIM=MyMeshType::MY_SPACEDIM;
    static const int MESHDIM=MyMeshType::MY_MESHDIM;
    typedef typename MyMeshType::MyConnType ConnType;
    static const NumberingPolicy numPol=MyMeshType::My_numPol;
  public:
    //! \addtogroup InterpKerGrpIntPlan @{
    PlanarIntersector(double dimCaracteristic, double precision, double medianPlane, bool doRotate, int printLevel);
    //! @}
    virtual ~PlanarIntersector() { }
  
    /*!
     * \addtogroup InterpKerGrpIntPlan
     * @{
     */
    /*!
     * Tool for cell intersection, result is always positive.
     */
    virtual double intersectCells(int icell_A, int icell_B, int nb_NodesA, int nb_NodesB) = 0;
    //! @}
    //Tool for cell filtering
    void createBoundingBoxes(const MyMeshType& mesh, std::vector<double>& bbox);
    void adjustBoundingBoxes(std::vector<double>& bbox, double Surf3DAdjustmentEps);
    inline void getElemBB(double* bb, const MyMeshType& mesh, ConnType iP, ConnType nb_nodes);
  protected :
    static int projection(std::vector<double>& Coords_A, std::vector<double>& Coords_B, 
                           int nb_NodesA, int nb_NodesB, double epsilon, double median_plane, bool do_rotate);
    static void rotate3DTriangle( double* PP1, double*PP2, double*PP3,
                                  TranslationRotationMatrix& rotation_matrix);
  protected:
    double _dimCaracteristic;
    double _precision;
    double _medianPlane;
    bool _doRotate;
    int _printLevel;
  };
}

#endif
