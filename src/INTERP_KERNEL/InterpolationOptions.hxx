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
#ifndef __INTERPOLATIONOPTIONS_HXX__
#define __INTERPOLATIONOPTIONS_HXX__


namespace INTERP_KERNEL {
  typedef enum { Triangulation, Convex, Geometric2D } IntersectionType;
  /// Type describing the different ways in which the hexahedron can be split into tetrahedra.
  /// The PLANAR_* policies persume that each face is to be considered planar, while the general
  /// policies make no such hypothesis. The integer at the end gives the number of tetrahedra
  /// that result from the split.
  typedef enum  { PLANAR_FACE_5 = 5, PLANAR_FACE_6 = 6, GENERAL_24 = 24, GENERAL_48 = 48 } SplittingPolicy;

  
  class InterpolationOptions{
  private :
    int _print_level ;
    IntersectionType _intersection_type;
    double _precision;
    double _median_plane ;
    bool _do_rotate ;
    double _bounding_box_adjustment ;
    int _orientation ;
    SplittingPolicy _splitting_policy ;

  public:
    InterpolationOptions() { init(); }
    int getPrintLevel() const { return _print_level; }
    void setPrintLevel(int pl) { _print_level=pl; }

    IntersectionType getIntersectionType() const { return InterpolationOptions::_intersection_type; }
    void setIntersectionType(IntersectionType it) { InterpolationOptions::_intersection_type=it; }

    double getPrecision() const { return InterpolationOptions::_precision; }
    void setPrecision(double p) { InterpolationOptions::_precision=p; }

    double getMedianPlane() { return InterpolationOptions::_median_plane; }
    void setMedianPlane(double mp) { InterpolationOptions::_median_plane=mp; }
    
    bool getDoRotate() { return InterpolationOptions::_do_rotate; }
    void setDoRotate( bool dr) { InterpolationOptions::_do_rotate = dr; }
    
    double getBoundingBoxAdjustment() { return InterpolationOptions::_bounding_box_adjustment; }
    void setBoundingBoxAdjustment(double bba) { InterpolationOptions::_bounding_box_adjustment=bba; }
    
    int getOrientation() { return InterpolationOptions::_orientation; }
    void setOrientation(int o) { InterpolationOptions::_orientation=o; }
    
    SplittingPolicy getSplittingPolicy() { return _splitting_policy; }
    void setSplittingPolicy(SplittingPolicy sp) { _splitting_policy=sp; }
    void init()
    {  
      _print_level=0;
      _intersection_type=Triangulation;
      _precision=1e-12;;
      _median_plane=0.5;
      _do_rotate=true;
      _bounding_box_adjustment=0.1;
      _orientation=0;
      _splitting_policy=GENERAL_48;
    }
  };

}
#endif
