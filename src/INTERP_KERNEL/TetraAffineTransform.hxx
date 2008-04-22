#ifndef __TETRA_AFFINE_TRANSFORM_HXX__
#define __TETRA_AFFINE_TRANSFORM_HXX__

#undef INVERSION_SELF_CHECK // debugging : check that calculated inverse is correct

namespace INTERP_KERNEL
{

  /**
   * \brief Class representing an affine transformation x -> Ax + b that transforms a given tetrahedron
   * into the unit tetrahedron.
   *
   */
  class TetraAffineTransform
  {

  public:
    TetraAffineTransform(const double** pts);

    void apply(double* destPt, const double* srcPt) const;

    double determinant() const;

    void dump() const;

  private:

    void invertLinearTransform();

    void calculateDeterminant();

    void factorizeLU(double* lu, int* idx) const;
      
    void forwardSubstitution(double* x, const double* lu, const double* b, const int* idx) const;

    void backwardSubstitution(double* x, const double* lu, const double* b, const int* idx) const;

    // The affine transformation Ax + b is represented with _linearTransformation containing the elements of
    // A in row-first ordering and _translation containing the elements of b

    /// 3x3 matrix A in affine transform x -> Ax + b
    double _linearTransform[9];

    /// 3x1 vector b in affine transform x -> Ax + b
    double _translation[3];

    /// The determinant of the matrix A is calculated at the construction of the object and cached.
    double _determinant;

  };


};


#endif
