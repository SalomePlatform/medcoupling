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

#include "InterpolationOptions.hxx"

#include <sstream>

const double INTERP_KERNEL::InterpolationOptions::DFT_MEDIAN_PLANE=0.5;

const double INTERP_KERNEL::InterpolationOptions::DFT_SURF3D_ADJ_EPS=1.e-4;

const double INTERP_KERNEL::InterpolationOptions::DFT_MAX_DIST_3DSURF_INTERSECT=-1.;

const char INTERP_KERNEL::InterpolationOptions::PRECISION_STR[]="Precision";

const char INTERP_KERNEL::InterpolationOptions::MEDIANE_PLANE_STR[]="MedianPlane";

const char INTERP_KERNEL::InterpolationOptions::BOUNDING_BOX_ADJ_STR[]="BoundingBoxAdjustment";

const char INTERP_KERNEL::InterpolationOptions::BOUNDING_BOX_ADJ_ABS_STR[]="BoundingBoxAdjustmentAbs";

const char INTERP_KERNEL::InterpolationOptions::MAX_DISTANCE_3DSURF_INSECT_STR[]="MaxDistance3DSurfIntersect";

const char INTERP_KERNEL::InterpolationOptions::PRINT_LEV_STR[]="PrintLevel";

const char INTERP_KERNEL::InterpolationOptions::DO_ROTATE_STR[]="DoRotate";

const char INTERP_KERNEL::InterpolationOptions::ORIENTATION_STR[]="Orientation";

const char INTERP_KERNEL::InterpolationOptions::MEASURE_ABS_STR[]="MeasureAbs";

const char INTERP_KERNEL::InterpolationOptions::INTERSEC_TYPE_STR[]="IntersectionType";

const char INTERP_KERNEL::InterpolationOptions::SPLITTING_POLICY_STR[]="SplittingPolicy";

const char INTERP_KERNEL::InterpolationOptions::TRIANGULATION_INTERSECT2D_STR[]="Triangulation";

const char INTERP_KERNEL::InterpolationOptions::CONVEX_INTERSECT2D_STR[]="Convex";

const char INTERP_KERNEL::InterpolationOptions::GEOMETRIC_INTERSECT2D_STR[]="Geometric2D";

const char INTERP_KERNEL::InterpolationOptions::POINTLOCATOR_INTERSECT_STR[]="PointLocator";

const char INTERP_KERNEL::InterpolationOptions::PLANAR_SPLIT_FACE_5_STR[]="PLANAR_FACE_5";

const char INTERP_KERNEL::InterpolationOptions::PLANAR_SPLIT_FACE_6_STR[]="PLANAR_FACE_6";

const char INTERP_KERNEL::InterpolationOptions::GENERAL_SPLIT_24_STR[]="GENERAL_24";

const char INTERP_KERNEL::InterpolationOptions::GENERAL_SPLIT_48_STR[]="GENERAL_48";

std::string INTERP_KERNEL::InterpolationOptions::getIntersectionTypeRepr() const
{
  if(_intersection_type==INTERP_KERNEL::Triangulation)
    return std::string(TRIANGULATION_INTERSECT2D_STR);
  else if(_intersection_type==INTERP_KERNEL::Convex)
    return std::string(CONVEX_INTERSECT2D_STR);
  else if(_intersection_type==INTERP_KERNEL::Geometric2D)
    return std::string(GEOMETRIC_INTERSECT2D_STR);
  else if(_intersection_type==INTERP_KERNEL::PointLocator)
    return std::string(POINTLOCATOR_INTERSECT_STR);
  else
    return std::string("UNKNOWN_INTERSECT_TYPE");
}

bool INTERP_KERNEL::InterpolationOptions::setOptionDouble(const std::string& key, double value)
{
  if(key==PRECISION_STR) 
    {
      setPrecision(value);
      return true;
    }
  else if(key==MEDIANE_PLANE_STR) 
    {
      setMedianPlane(value);
      return true;
    }
  else if(key==BOUNDING_BOX_ADJ_STR) 
    {
      setBoundingBoxAdjustment(value);
      return true;
    }
  else if(key==BOUNDING_BOX_ADJ_ABS_STR) 
    {
      setBoundingBoxAdjustmentAbs(value);
      return true;
    }
  else if(key==MAX_DISTANCE_3DSURF_INSECT_STR) 
    {
      setMaxDistance3DSurfIntersect(value);
      return true;
    }
  else
    return false;
}

bool INTERP_KERNEL::InterpolationOptions::setOptionInt(const std::string& key, int value)
{
  if(key==PRINT_LEV_STR) 
    {
      setPrintLevel(value);
      return true;
    }
    else if(key==DO_ROTATE_STR) 
      {
        setDoRotate(value != 0);
        return true;
      }
    else if(key==ORIENTATION_STR) 
      {
        setOrientation(value);
        return true;
      }
    else if(key==MEASURE_ABS_STR)
      {
        bool valBool=(value!=0);
        setMeasureAbsStatus(valBool);
        return true;
      }
    else
      return false;
}

bool INTERP_KERNEL::InterpolationOptions::setOptionString(const std::string& key, const std::string& value)
{
  if(key==INTERSEC_TYPE_STR) 
    {
      if(value==TRIANGULATION_INTERSECT2D_STR)
        {
          setIntersectionType(INTERP_KERNEL::Triangulation);
          return true;
        }
      else if(value==CONVEX_INTERSECT2D_STR)
        {
          setIntersectionType(INTERP_KERNEL::Convex);
          return true;
        }
      else if(value==GEOMETRIC_INTERSECT2D_STR)
        {
          setIntersectionType(INTERP_KERNEL::Geometric2D);
          return true;
        }
      else if(value==POINTLOCATOR_INTERSECT_STR)
        {
          setIntersectionType(INTERP_KERNEL::PointLocator);
          return true;
        }
    }
  else if(key==SPLITTING_POLICY_STR) 
    {
      if(value==PLANAR_SPLIT_FACE_5_STR)
        {
          setSplittingPolicy(INTERP_KERNEL::PLANAR_FACE_5);
          return true;
        }
      else if(value==PLANAR_SPLIT_FACE_6_STR)
        {
          setSplittingPolicy(INTERP_KERNEL::PLANAR_FACE_6);
          return true;
        }
      else if(value==GENERAL_SPLIT_24_STR)
        {
          setSplittingPolicy(INTERP_KERNEL::GENERAL_24);
          return true;
        }
      else if(value==GENERAL_SPLIT_48_STR)
        {
          setSplittingPolicy(INTERP_KERNEL::GENERAL_48);
          return true;
        }
      else
        return false;
    }
  return false;
}

std::string INTERP_KERNEL::InterpolationOptions::getSplittingPolicyRepr() const
{
  if(_splitting_policy==INTERP_KERNEL::PLANAR_FACE_5)
    return std::string(PLANAR_SPLIT_FACE_5_STR);
  else if(_splitting_policy==INTERP_KERNEL::PLANAR_FACE_6)
    return std::string(PLANAR_SPLIT_FACE_6_STR);
  else if(_splitting_policy==INTERP_KERNEL::GENERAL_24)
    return std::string(GENERAL_SPLIT_24_STR);
  else if(_splitting_policy==INTERP_KERNEL::GENERAL_48)
    return std::string(GENERAL_SPLIT_48_STR);
  else
    return std::string("UNKNOWN_SPLITTING_POLICY");
}

std::string INTERP_KERNEL::InterpolationOptions::filterInterpolationMethod(const std::string& meth) const
{
  if ( _P1P0_bary_method && meth == "P1P0" )
    return "P1P0Bary";
  return meth;
}

bool INTERP_KERNEL::InterpolationOptions::setInterpolationOptions(long print_level,
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
                                                                  bool P1P0_bary_method )
{
  _print_level=print_level;
  _precision=precision;
  _median_plane=median_plane;
  _do_rotate=do_rotate;
  _bounding_box_adjustment=bounding_box_adjustment;
  _bounding_box_adjustment_abs=bounding_box_adjustment_abs;
  _max_distance_for_3Dsurf_intersect=max_distance_for_3Dsurf_intersect;
  _orientation=orientation;
  _measure_abs=measure_abs;
  _P1P0_bary_method=P1P0_bary_method;
  return(setOptionString(INTERSEC_TYPE_STR,intersection_type) && setOptionString(SPLITTING_POLICY_STR,splitting_policy));
}

std::string INTERP_KERNEL::InterpolationOptions::printOptions() const
{
  std::ostringstream oss; oss.precision(15); oss << "Interpolation Options ******" << std::endl;
  oss << "Print level : " << _print_level << std::endl;
  oss << "Intersection type : " << getIntersectionTypeRepr() << std::endl;
  oss << "Precision : " << _precision << std::endl;
  oss << "Median plane : " << _median_plane << std::endl;
  oss << "Do Rotate status : " << std::boolalpha << _do_rotate << std::endl;
  oss << "Bounding box adj : " << _bounding_box_adjustment << std::endl;
  oss << "Bounding box adj abs : " << _bounding_box_adjustment_abs << std::endl;
  oss << "Max distance for 3DSurf intersect : " << _max_distance_for_3Dsurf_intersect << std::endl;
  oss << "Orientation : " << _orientation << std::endl;
  oss << "Measure abs : " << _measure_abs << std::endl;
  oss << "Splitting policy : " << getSplittingPolicyRepr() << std::endl;
  oss << "P1P0 Barycentric method : " << _P1P0_bary_method << std::endl;
  oss << "****************************" << std::endl;
  return oss.str();
}
