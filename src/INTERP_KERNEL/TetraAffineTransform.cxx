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

#include "TetraAffineTransform.hxx"
#include "VectorUtils.hxx"

#include <cmath>
#include <cstring>
#include <iostream>

#include "Log.hxx"

namespace INTERP_KERNEL
{
  /////////////////////////////////////////////////////////////////////////////////////////
  /// PUBLIC INTERFACE METHODS                                               //////////////
  /////////////////////////////////////////////////////////////////////////////////////////

  /**
   * Constructor
   * Create the TetraAffineTransform object from the tetrahedron 
   * with corners specified in pts. If the tetrahedron is degenerate or almost degenerate, 
   * construction succeeds, but the determinant of the transform is set to 0.
   *
   * @param pts  a 4x3 matrix containing 4 points (P0X,P0Y,P0Z,P1X,P1Y,P1Z,...) of 3 coordinates each
   */
  TetraAffineTransform::TetraAffineTransform(const double *pts)
  {
    
    LOG(2,"Creating transform from tetraeder : ");
    
    // three last points -> linear transform
    for(int i = 0; i < 3 ; ++i)
      {
        for(int j = 0 ; j < 3 ; ++j)
          {
            // NB we insert columns, not rows
            _linear_transform[3*j + i] = pts[(i+1)*3+j] - pts[j];
          }
      }

    // remember _linear_transform for the reverse transformation
    memcpy( _back_linear_transform, _linear_transform, 9*sizeof(double));
    memcpy( _back_translation,      pts,          3*sizeof(double));
    
    calculateDeterminant();
    
    LOG(3, "determinant before inverse = " << _determinant);
    
    // check that tetra is non-planar -> determinant is not zero
    // otherwise set _determinant to zero to signal caller that transformation did not work
    if(epsilonEqual(_determinant, 0.0))
      {
        _determinant = 0.0;
        return;
      }

    // we need the inverse transform
    invertLinearTransform();
    
    // first point -> translation
    // calculate here because translation takes place in "transformed space",
    // or in other words b = -A*O where A is the linear transform
    // and O is the position vector of the point that is mapped onto the origin
    for(int i = 0 ; i < 3 ; ++i)
      {
        _translation[i] = -(_linear_transform[3*i]*pts[0] + _linear_transform[3*i+1]*pts[1] + _linear_transform[3*i+2]*pts[2]) ;
      }
    
    // precalculate determinant (again after inversion of transform)
    calculateDeterminant();

#ifdef INVERSION_SELF_CHECK
    // debugging : check that applying the inversed transform to the original points
    // gives us the unit tetrahedron
    LOG(4, "transform determinant is " << _determinant);
    LOG(4, "*Self-check : Applying transformation to original points ... ");
    for(int i = 0; i < 4 ; ++i)
      {
        double v[3];
        apply(v, pts+3*i);
        LOG(4, vToStr(v))
          for(int j = 0; j < 3; ++j)
            {
              assert(epsilonEqual(v[j], (3*i+j == 3 || 3*i+j == 7 || 3*i+j == 11 ) ? 1.0 : 0.0));
            }
      }
    
    LOG(4, " ok");
#endif
  }

  /**
   * Calculates the transform of point srcPt and stores the result in destPt.
   * If destPt == srcPt, then srcPt is overwritten safely.
   *  
   *
   * @param destPt  double[3] in which to store the transformed point
   * @param srcPt   double[3] containing coordinates of points to be transformed
   *
   */
  void TetraAffineTransform::apply(double* destPt, const double* srcPt) const
  {
    double* dest = destPt;
    
    // are we self-allocating ?
    const bool selfAllocation = (destPt == srcPt);
    
    if(selfAllocation)
      {
        // alloc temporary memory
        dest = new double[3];
       
        LOG(6, "Info : Self-affectation in TetraAffineTransform::apply");
      }
    
    for(int i = 0 ; i < 3 ; ++i)
      {
        // matrix - vector multiplication
        dest[i] = _linear_transform[3*i] * srcPt[0] + _linear_transform[3*i + 1] * srcPt[1] + _linear_transform[3*i + 2] * srcPt[2];
       
        // translation
        dest[i] += _translation[i];
      }

    if(selfAllocation)
      {
        // copy result back to destPt
        for(int i = 0 ; i < 3 ; ++i)
          {
            destPt[i] = dest[i];
          }
        delete[] dest;
      }
  }

  /**
   * Calculates the reverse transform of point srcPt and stores the result in destPt.
   * If destPt == srcPt, then srcPt is overwritten safely.
   *
   * @param destPt  double[3] in which to store the transformed point
   * @param srcPt   double[3] containing coordinates of points to be transformed
   */
  void TetraAffineTransform::reverseApply(double* destPt, const double* srcPt) const
  {
    double* dest = destPt;
    
    // are we self-allocating ?
    const bool selfAllocation = (destPt == srcPt);
    
    if(selfAllocation)
      {
        // alloc temporary memory
        dest = new double[3];
       
        LOG(6, "Info : Self-affectation in TetraAffineTransform::reverseApply");
      }
    
    for(int i = 0 ; i < 3 ; ++i)
      {
        // matrix - vector multiplication
        dest[i] = _back_linear_transform[3*i] * srcPt[0] + _back_linear_transform[3*i + 1] * srcPt[1] + _back_linear_transform[3*i + 2] * srcPt[2];

        // translation
        dest[i] += _back_translation[i];
      }

    if(selfAllocation)
      {
        // copy result back to destPt
        for(int i = 0 ; i < 3 ; ++i)
          {
            destPt[i] = dest[i];
          }
        delete[] dest;
      }
  }

  /**
   * Returns the determinant of the linear part A of the transform.
   *
   * @return determinant of the transform
   *
   */
  double TetraAffineTransform::determinant() const
  {
    return _determinant;
  }

  /**
   * Outputs  to std::cout the matrix A and the vector b
   * of the transform Ax + b
   *
   */
  void TetraAffineTransform::dump() const
  {
    std::cout << "A = " << std::endl << "[";
    for(int i = 0; i < 3; ++i)
      {
        std::cout << _linear_transform[3*i] << ", " << _linear_transform[3*i + 1] << ", " << _linear_transform[3*i + 2];
        if(i != 2 ) std::cout << std::endl;
      }
    std::cout << "]" << std::endl;
    
    std::cout << "b = " << "[" << _translation[0] << ", " << _translation[1] << ", " << _translation[2] << "]" << std::endl;
  }

  /////////////////////////////////////////////////////////////////////////////////////////
  /// PRIVATE  METHODS                                                       //////////////
  /////////////////////////////////////////////////////////////////////////////////////////

  /**
   * Calculates the inverse of the matrix A, stored in _linear_transform
   * by LU-factorization and substitution
   *
   */
  void TetraAffineTransform::invertLinearTransform()
  {
    //{ we copy the matrix for the lu-factorization
    // maybe inefficient    
    double lu[9];
    for(int i = 0 ; i < 9; ++i)
      {
        lu[i] = _linear_transform[i];
      }
    
    // calculate LU factorization
    int idx[3];
    factorizeLU(lu, idx);
    
    // calculate inverse by forward and backward substitution
    // store in _linear_transform
    // NB _linear_transform cannot be overwritten with lu in the loop
    for(int i = 0 ; i < 3 ; ++i)
      {
        // form standard base vector i
        const double b[3] = 
          {
            double ( int(i == 0) ),
            double ( int(i == 1) ),
            double ( int(i == 2) ),
          };

        LOG(6,  "b = [" << b[0] << ", " << b[1] << ", " << b[2] << "]");
       
        double y[3];
        forwardSubstitution(y, lu, b, idx);
       
        double x[3];
        backwardSubstitution(x, lu, y, idx);
       
        // copy to _linear_transform matrix
        // NB : this is a column operation, so we cannot 
        // do this directly when we calculate x
        for(int j = 0 ; j < 3 ; j++)
          {
            _linear_transform[3*j + i] = x[idx[j]];
          }
      }
  }
  
  /**
   * Updates the member _determinant of the matrix A of the transformation.
   *
   */
  void TetraAffineTransform::calculateDeterminant()
  {
    const double subDet[3] = 
      {
        _linear_transform[4] * _linear_transform[8] - _linear_transform[5] * _linear_transform[7],
        _linear_transform[3] * _linear_transform[8] - _linear_transform[5] * _linear_transform[6],
        _linear_transform[3] * _linear_transform[7] - _linear_transform[4] * _linear_transform[6]
      };

    _determinant = _linear_transform[0] * subDet[0] - _linear_transform[1] * subDet[1] + _linear_transform[2] * subDet[2]; 
  }

  
  /////////////////////////////////////////////////
  /// Auxiliary methods for inverse calculation ///
  /////////////////////////////////////////////////


  /**
   * Calculates the LU-factorization of the matrix A (_linear_transform)
   * and stores it in lu. Since partial pivoting is used, there are 
   * row swaps. This is represented by the index permutation vector idx : to access element
   * (i,j) of lu, use lu[3*idx[i] + j]
   *
   * @param lu double[9] in which to store LU-factorization
   * @param idx int[3] in which to store row permutation vector
   */
  void TetraAffineTransform::factorizeLU(double* lu, int* idx) const
  {
    // 3 x 3 LU factorization
    // initialise idx
    for(int i = 0 ; i < 3 ; ++i)
      {
        idx[i] = i;
      }
            
    for(int k = 0; k < 2 ; ++k)
      {
         
        // find pivot
        int i = k;
        double max = std::fabs(lu[3*idx[k] + k]);
        int row = i;
        while(i < 3)
          {
            if(std::fabs(lu[3*idx[i] + k]) > max)
              {
                max = fabs(lu[3*idx[i] + k]);
                row = i;
              }
            ++i;
          }
             
        // swap rows in index vector
        int tmp = idx[k];
        idx[k] = idx[row];
        idx[row] = tmp;
      
        // calculate row
        for(int j = k + 1 ; j < 3 ; ++j)
          {
            // l_jk = u_jk / u_kk
            lu[3*idx[j] + k] /= lu[3*idx[k] + k];
            for(int s = k + 1; s < 3 ; ++s)
              {
                // case s = k will always become zero, and then be replaced by
                // the l-value
                // there's no need to calculate this explicitly

                // u_js = u_js - l_jk * u_ks
                lu[3*idx[j] + s] -= lu[3*idx[j] + k] * lu[3*idx[k] + s] ;
              }
          }
      }
  }

  /**
   * Solves the system Lx = b, where L is lower unit-triangular (ones on the diagonal)
   * 
   * @param x   double[3] in which the solution is stored
   * @param lu  double[9] containing the LU-factorization
   * @param b   double[3] containing the right-hand side
   * @param idx int[3] containing the permutation vector associated with lu
   *
   */
  void TetraAffineTransform::forwardSubstitution(double* x, const double* lu, const double* b, const int* idx) const
  {
    // divisions are not carried out explicitly because lu does not store
    // the diagonal values of L, which are all 1
    // Compare with backwardSubstitution()
    x[idx[0]] = ( b[idx[0]] ); // / lu[3*idx[0]]
    x[idx[1]] = ( b[idx[1]] - lu[3*idx[1]] * x[idx[0]] ); // / lu[3*idx[1] + 1];
    x[idx[2]] = ( b[idx[2]] - lu[3*idx[2]] * x[idx[0]] - lu[3*idx[2] + 1] * x[idx[1]] ); // / lu[3*idx[2] + 2];
  }

  /**
   * Solves the system Ux = b, where U is upper-triangular 
   * 
   * @param x   double[3] in which the solution is stored
   * @param lu  double[9] containing the LU-factorization
   * @param b   double[3] containing the right-hand side
   * @param idx int[3] containing the permutation vector associated with lu
   *
   */
  void TetraAffineTransform::backwardSubstitution(double* x, const double* lu, const double* b, const int* idx) const
  {
    x[idx[2]] = ( b[idx[2]] ) / lu[3*idx[2] + 2];
    x[idx[1]] = ( b[idx[1]] - lu[3*idx[1] + 2] * x[idx[2]] ) / lu[3*idx[1] + 1];
    x[idx[0]] = ( b[idx[0]] - lu[3*idx[0] + 1] * x[idx[1]] - lu[3*idx[0] + 2] * x[idx[2]] ) / lu[3*idx[0]];
  }
  
}
