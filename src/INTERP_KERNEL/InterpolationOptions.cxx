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
#include "InterpolationOptions.hxx"

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

const char INTERP_KERNEL::InterpolationOptions::POINTLOCATOR_INTERSECT2D_STR[]="PointLocator2D";

const char INTERP_KERNEL::InterpolationOptions::PLANAR_SPLIT_FACE_5_STR[]="PLANAR_FACE_5";

const char INTERP_KERNEL::InterpolationOptions::PLANAR_SPLIT_FACE_6_STR[]="PLANAR_FACE_6";

const char INTERP_KERNEL::InterpolationOptions::GENERAL_SPLIT_24_STR[]="GENERAL_24";

const char INTERP_KERNEL::InterpolationOptions::GENERAL_SPLIT_48_STR[]="GENERAL_48";

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
        setDoRotate(value);
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
      }
    else
      return false;
}

bool INTERP_KERNEL::InterpolationOptions::setOptionString(const std::string& key, std::string& value)
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
      else if(value==POINTLOCATOR_INTERSECT2D_STR)
        {
          setIntersectionType(INTERP_KERNEL::PointLocator2D);
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

std::string INTERP_KERNEL::InterpolationOptions::filterInterpolationMethod(const std::string& meth) const
{
  if ( _P1P0_bary_method && meth == "P1P0" )
    return "P1P0Bary";
  return meth;
}
