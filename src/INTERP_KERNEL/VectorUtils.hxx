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

#ifndef __VECTORUTILS_HXX__
#define __VECTORUTILS_HXX__

#include <sstream>
#include <numeric>
#include <string>
#include <cmath>
#include <map>

namespace INTERP_KERNEL
{
  /// Precision used for tests of 3D part of INTERP_KERNEL
  const double VOL_PREC = 1.0e-6;
  
  /// Default relative tolerance in epsilonEqualRelative
  const double DEFAULT_REL_TOL = 1.0e-6;
  
  /// Default absolute tolerance in epsilonEqual and epsilonEqualRelative
  const double DEFAULT_ABS_TOL = 5.0e-12;

  /**
   * @param a first point. Should point on a array of size at least equal to SPACEDIM.
   * @param b second point. Should point on a array of size at least equal to SPACEDIM.
   */
  template<int SPACEDIM>
  inline double getDistanceBtw2Pts(const double *a, const double *b)
  {
    double ret2=0.;
    for(int i=0;i<SPACEDIM;i++)
      ret2+=(a[i]-b[i])*(a[i]-b[i]);
    return sqrt(ret2);
  }

  // -------------------------------------------------------------------
  // Math operations for vectors represented by double[3] - arrays  
  // -------------------------------------------------------------------
  
  /**
   * Copies a double[3] vector from src to dest
   *
   * @param src   source vector
   * @param dest  destination vector
   *
   */
  inline void copyVector3(const double* src, double* dest)
  {
    for(int i = 0 ; i < 3 ; ++i)
      dest[i] = src[i];
  }
  
  /**
   * Creates a string representation of a double[3] vector
   *
   * @param  pt  a 3-vector
   * @return a string of the form [x, y, z]
   */
  inline const std::string vToStr(const double* pt)
  {
    std::stringstream ss(std::ios::out);
    ss << "[" << pt[0] << ", " << pt[1] << ", " << pt[2] << "]";
    return ss.str();
  }

  /**
   * Adds a double[3] - vector to another one.
   *
   * @param v     vector v
   * @param res   vector in which to store the result res + v.
   */
  inline void add(const double* v, double* res)
  {
    res[0] += v[0];
    res[1] += v[1];
    res[2] += v[2];
  }

  /**
   * Calculates the cross product of two double[3] - vectors.
   *
   * @param v1    vector v1
   * @param v2    vector v2
   * @param res   vector in which to store the result v1 x v2. It should not be one of v1 and v2.
   */
  inline void cross(const double* v1, const double* v2,double* res)
  {
    res[0] = v1[1]*v2[2] - v1[2]*v2[1];
    res[1] = v1[2]*v2[0] - v1[0]*v2[2];
    res[2] = v1[0]*v2[1] - v1[1]*v2[0];
  }

  /**
   * Calculates the dot product of two double[3] - vectors
   *
   * @param v1   vector v1
   * @param v2   vector v2
   * @return   dot (scalar) product v1.v2
   */
  inline double dot(const double* v1, const double* v2)
  {
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
  }

  /**
   * Calculates norm of a double[3] vector
   *
   * @param v  a vector v
   * @return euclidean norm of v
   */
  inline double norm(const double* v)
  {
    return sqrt(dot(v,v));
  }

  /**
   * Compares doubles using an absolute tolerance
   * This is suitable mainly for comparisons with 0.0
   * 
   * @param x         first value
   * @param y         second value
   * @param errTol    maximum allowed absolute difference that is to be treated as equality
   * @return  true if |x - y| < errTol, false otherwise
   */
  inline bool epsilonEqual(const double x, const double y, const double errTol = DEFAULT_ABS_TOL)
  {
    return y < x ? x - y < errTol : y - x < errTol;
    //    return std::fabs(x - y) < errTol;
  }


  /**
   * Test whether two 3D vectors are colinear. The two vectors are expected to be of unit norm (not checked)
   * Implemented by checking that the norm of the cross product is null.
   */
  inline bool isColinear3D(const double *v1, const double *v2, const double eps = DEFAULT_ABS_TOL)
  {
    double cros[3];
    cross(v1, v2, cros);
    return epsilonEqual(dot(cros, cros), 0.0, eps);
  }


  /**
   * Compares doubles using a relative tolerance
   * This is suitable mainly for comparing larger values to each other. Before performing the relative test,
   * an absolute test is performed to guard from problems when comparing to 0.0
   * 
   * @param x         first value
   * @param y         second value
   * @param relTol    maximum allowed relative difference that is to be treated as equality
   * @param absTol    maximum allowed absolute difference that is to be treated as equality
   * @return  true if |x - y| <= absTol or |x - y|/max(|x|,|y|) <= relTol, false otherwise
   */
  inline bool epsilonEqualRelative(const double x, const double y, const double relTol = DEFAULT_REL_TOL, const double absTol = DEFAULT_ABS_TOL)
  {
    // necessary for comparing values close to zero
    // in order to avoid division by very small numbers
    if(std::fabs(x - y) < absTol)
      {
        return true;
      }

    const double relError = std::fabs((x - y) / std::max(std::fabs(x), std::fabs(y)));

    return relError < relTol;
  }

}

#endif
