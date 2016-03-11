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

#ifndef __TETRA_AFFINE_TRANSFORM_HXX__
#define __TETRA_AFFINE_TRANSFORM_HXX__

#include "INTERPKERNELDefines.hxx"

#undef INVERSION_SELF_CHECK // debugging : check that calculated inverse is correct

namespace INTERP_KERNEL
{
  /**
   * \brief Class representing an affine transformation x -> Ax + b that transforms a given tetrahedron
   * into the unit tetrahedron.
   *
   */
  class INTERPKERNEL_EXPORT TetraAffineTransform
  {

  public:
    TetraAffineTransform(const double *pts);

    void apply(double* destPt, const double* srcPt) const;

    void reverseApply(double* destPt, const double* srcPt) const;

    double determinant() const;

    void dump() const;

  private:

    void invertLinearTransform();

    void calculateDeterminant();

    void factorizeLU(double* lu, int* idx) const;
      
    void forwardSubstitution(double* x, const double* lu, const double* b, const int* idx) const;

    void backwardSubstitution(double* x, const double* lu, const double* b, const int* idx) const;

    // The affine transformation Ax + b is represented with _linear_transformation containing the elements of
    // A in row-first ordering and _translation containing the elements of b

    /// 3x3 matrix A in affine transform x -> Ax + b
    double _linear_transform[9];

    /// 3x1 vector b in affine transform x -> Ax + b
    double _translation[3];

    /// The determinant of the matrix A is calculated at the construction of the object and cached.
    double _determinant;

    /// 3x3 matrix AT is transposed A matrix used for y -> ATy - c transformation
    double _back_linear_transform[9];

    /// 3x1 vector c in affine transform y -> ATy - c
    double _back_translation[3];

  };
}

#endif
