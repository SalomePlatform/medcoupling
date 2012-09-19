// Copyright (C) 2007-2012  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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

#ifndef __INTERPOLATIONOPTIONS_HXX__
#define __INTERPOLATIONOPTIONS_HXX__

#include "INTERPKERNELDefines.hxx"

#include <string>

namespace INTERP_KERNEL
{
  typedef enum { Triangulation, Convex, Geometric2D, PointLocator } IntersectionType;
  /// Type describing the different ways in which the hexahedron can be split into tetrahedra.
  /// The PLANAR_* policies persume that each face is to be considered planar, while the general
  /// policies make no such hypothesis. The integer at the end gives the number of tetrahedra
  /// that result from the split.
  typedef enum  { PLANAR_FACE_5 = 5, PLANAR_FACE_6 = 6, GENERAL_24 = 24, GENERAL_48 = 48 } SplittingPolicy;
  
  /*!
   * \class InterpolationOptions
   * Class defining the options for all interpolation algorithms.
   * 
   * List of options, possible values and default values can be found on this page:
   * \ref InterpKerIntersectors
   */
  class INTERPKERNEL_EXPORT InterpolationOptions
  {
  private:
    int _print_level ;
    IntersectionType _intersection_type;
    double _precision;
    double _median_plane ;
    bool _do_rotate ;
    //! this measure is relative to the caracteristic dimension
    double _bounding_box_adjustment ;
    //! this measure is absolute \b not relative to the cell size
    double _bounding_box_adjustment_abs ;
    double _max_distance_for_3Dsurf_intersect;
    int _orientation ;
    bool _measure_abs;
    SplittingPolicy _splitting_policy ;
    bool _P1P0_bary_method; // issue 0020440

  public:
    InterpolationOptions() { init(); }
    int getPrintLevel() const { return _print_level; }
    void setPrintLevel(int pl) { _print_level=pl; }

    IntersectionType getIntersectionType() const { return _intersection_type; }
    void setIntersectionType(IntersectionType it) { _intersection_type=it; }
    std::string getIntersectionTypeRepr() const;

    double getPrecision() const { return _precision; }
    void setPrecision(double p) { _precision=p; }

    double getMedianPlane() const { return _median_plane; }
    void setMedianPlane(double mp) { _median_plane=mp; }
    
    bool getDoRotate() const { return _do_rotate; }
    void setDoRotate( bool dr) { _do_rotate = dr; }
    
    double getBoundingBoxAdjustment() const { return _bounding_box_adjustment; }
    void setBoundingBoxAdjustment(double bba) { _bounding_box_adjustment=bba; }

    double getBoundingBoxAdjustmentAbs() const { return _bounding_box_adjustment_abs; }
    void setBoundingBoxAdjustmentAbs(double bba) { _bounding_box_adjustment_abs=bba; }
    
    double getMaxDistance3DSurfIntersect() const { return _max_distance_for_3Dsurf_intersect; }
    void setMaxDistance3DSurfIntersect(double bba) { _max_distance_for_3Dsurf_intersect=bba; }

    int getOrientation() const { return _orientation; }
    void setOrientation(int o) { _orientation=o; }

    bool getMeasureAbsStatus() const { return _measure_abs; }
    void setMeasureAbsStatus(bool newStatus) { _measure_abs=newStatus; }
    
    SplittingPolicy getSplittingPolicy() const { return _splitting_policy; }
    void setSplittingPolicy(SplittingPolicy sp) { _splitting_policy=sp; }
    std::string getSplittingPolicyRepr() const;

    void setP1P0BaryMethod(bool isP1P0) { _P1P0_bary_method=isP1P0; }
    bool getP1P0BaryMethod() const { return _P1P0_bary_method; }

    std::string filterInterpolationMethod(const std::string& meth) const;

    void init()
    {  
      _print_level=0;
      _intersection_type=Triangulation;
      _precision=1e-12;
      _median_plane=DFT_MEDIAN_PLANE;
      _do_rotate=true;
      _bounding_box_adjustment=DFT_SURF3D_ADJ_EPS;
      _bounding_box_adjustment_abs=0.;
      _max_distance_for_3Dsurf_intersect=DFT_MAX_DIST_3DSURF_INTERSECT;
      _orientation=0;
      _measure_abs=true;
      _splitting_policy=GENERAL_48;
      _P1P0_bary_method=false;
    }
    bool setInterpolationOptions(long print_level,
                                 std::string intersection_type,
                                 double precision,
                                 double median_plane,
                                 bool do_rotate,
                                 double bounding_box_adjustment,
                                 double bounding_box_adjustment_abs,
                                 double max_distance_for_3Dsurf_intersect,
                                 long orientation,
                                 bool measure_abs,
                                 std::string splitting_policy,
                                 bool P1P0_bary_method );
    void copyOptions(const InterpolationOptions & other) { *this = other; }
    bool setOptionDouble(const std::string& key, double value);
    bool setOptionInt(const std::string& key, int value);
    bool setOptionString(const std::string& key, const std::string& value);
    std::string printOptions() const;
  private:
    static const double DFT_MEDIAN_PLANE;
    static const double DFT_SURF3D_ADJ_EPS;
    static const double DFT_MAX_DIST_3DSURF_INTERSECT;
  public:
    static const char PRECISION_STR[];
    static const char MEDIANE_PLANE_STR[];
    static const char BOUNDING_BOX_ADJ_STR[];
    static const char BOUNDING_BOX_ADJ_ABS_STR[];
    static const char MAX_DISTANCE_3DSURF_INSECT_STR[];
    static const char PRINT_LEV_STR[];
    static const char DO_ROTATE_STR[];
    static const char ORIENTATION_STR[];
    static const char MEASURE_ABS_STR[];
    static const char INTERSEC_TYPE_STR[];
    static const char SPLITTING_POLICY_STR[];
    static const char TRIANGULATION_INTERSECT2D_STR[];
    static const char CONVEX_INTERSECT2D_STR[];
    static const char GEOMETRIC_INTERSECT2D_STR[];
    static const char POINTLOCATOR_INTERSECT_STR[];
    static const char PLANAR_SPLIT_FACE_5_STR[];
    static const char PLANAR_SPLIT_FACE_6_STR[];
    static const char GENERAL_SPLIT_24_STR[];
    static const char GENERAL_SPLIT_48_STR[];
  };

}
#endif
