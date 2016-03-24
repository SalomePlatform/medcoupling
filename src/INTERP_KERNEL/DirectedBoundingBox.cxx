// Copyright (C) 2009-2016  OPEN CASCADE
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
// File      : DirectedBoundingBox.cxx
// Created   : Mon Apr 12 14:41:22 2010
// Author    : Edward AGAPOV (eap)

#include "DirectedBoundingBox.hxx"

#include "InterpolationUtils.hxx"

#define __TENSOR(i,j) tensor[(i)*_dim+(j)]
#define __AXIS(i)     (&_axes[(i)*_dim])
#define __MIN(i)      _minmax[i*2]
#define __MAX(i)      _minmax[i*2+1]
#define __MYID        (long(this)%10000)
#define __DMP(msg) \
  //  cout << msg << endl

using namespace std;

namespace
{
  //================================================================================
  /*!
   * \brief Add point coordinates to inertia tensor in 3D space
   */
  //================================================================================

  inline void addPointToInertiaTensor3D(const double*   coord,
                                        const double*   gc,
                                        vector<double>& tensor)
  {
    // we fill the upper triangle of tensor only
    const int _dim = 3;
    double x = coord[0] - gc[0], y = coord[1] - gc[1], z = coord[2] - gc[2];
    __TENSOR(0,0) += y*y + z*z;
    __TENSOR(1,1) += x*x + z*z;
    __TENSOR(2,2) += x*x + y*y;
    __TENSOR(0,1) -= x*y;
    __TENSOR(0,2) -= x*z;
    __TENSOR(1,2) -= y*z;
  }
  //================================================================================
  /*!
   * \brief Add point coordinates to inertia tensor in 2D space
   */
  //================================================================================

  inline void addPointToInertiaTensor2D(const double*   coord,
                                        const double*   gc,
                                        vector<double>& tensor)
  {
    // we fill the upper triangle of tensor only
    const int _dim = 2;
    double x = coord[0] - gc[0], y = coord[1] - gc[1];
    __TENSOR(0,0) += y*y;
    __TENSOR(1,1) += x*x;
    __TENSOR(0,1) -= x*y;
  }

  //================================================================================
  /*!
   * \brief Find eigenvectors of tensor using Jacobi's method
   */
  //================================================================================

  bool JacobiEigenvectorsSearch( const int _dim, vector<double>& tensor, vector<double>& _axes)
  {
    if ( _dim == 3 )
      {
        __DMP( "Tensor : {"
               << "{ "<<__TENSOR(0,0) << ", "<<__TENSOR(0,1) << ", "<<__TENSOR(0,2) << "} "
               << "{ "<<__TENSOR(1,0) << ", "<<__TENSOR(1,1) << ", "<<__TENSOR(1,2) << "} "
               << "{ "<<__TENSOR(2,0) << ", "<<__TENSOR(2,1) << ", "<<__TENSOR(2,2) << "}} ");
      }
    else
      {
        __DMP( "Tensor : {"
          << "{ "<<__TENSOR(0,0) << ", "<<__TENSOR(0,1) << "} "
          << "{ "<<__TENSOR(1,0) << ", "<<__TENSOR(1,1) << "}} ");
      }

    const int maxRot = 5*_dim*_dim; // limit on number of rotations
    const double tol = 1e-9;

    // set _axes to identity
    int i,j;
    for ( i = 0; i < _dim; ++i )
      for ( j = 0; j < _dim; ++j )
        __AXIS(i)[j] = ( i==j ? 1. : 0 );

    bool solved = false;
    for ( int iRot = 0; iRot < maxRot; ++ iRot )
      {
        // find max off-diagonal element of the tensor
        int k = 0, l = 0;
        double max = 0.;
        for ( i = 0; i < _dim-1; ++i )
          for ( j = i+1; j < _dim; ++j )
            if ( fabs( __TENSOR(i,j)) > max )
              max = fabs( __TENSOR(i,j) ), k = i, l = j;
        solved = ( max < tol );
        if ( solved )
          break;

        // Rotate to make __TENSOR(k,l) == 0

        double diff = __TENSOR(l,l) - __TENSOR(k,k);
        double t; // tangent of rotation angle
        if ( fabs(__TENSOR(k,l)) < abs(diff)*1.0e-36)
          {
            t = __TENSOR(k,l)/diff;
          }
        else
          {
            double phi = diff/(2.0*__TENSOR(k,l));
            t = 1.0/(abs(phi) + sqrt(phi*phi + 1.0));
            if ( phi < 0.0) t = -t;
          }
        double c = 1.0/sqrt(t*t + 1.0); // cosine of rotation angle
        double s = t*c; // sine of rotation angle
        double tau = s/(1.0 + c);
        __TENSOR(k,k) -= t*__TENSOR(k,l);
        __TENSOR(l,l) += t*__TENSOR(k,l);
        __TENSOR(k,l) = 0.0;

#define __ROTATE(T,r1,c1,r2,c2) \
{ \
int i1 = r1*_dim+c1, i2 = r2*_dim+c2; \
double t1 = T[i1], t2 = T[i2]; \
T[i1] -= s * ( t2 + tau * t1);\
T[i2] += s * ( t1 - tau * t2);\
}
        for ( i = 0; i < k; ++i )       // Case of i < k
          __ROTATE(tensor, i,k,i,l);

        for ( i = k+1; i < l; ++i )     // Case of k < i < l
          __ROTATE(tensor, k,i,i,l);

        for ( i = l + 1; i < _dim; ++i )   // Case of i > l
          __ROTATE(tensor, k,i,l,i);

        for ( i = 0; i < _dim; ++i )       // Update transformation matrix
          __ROTATE(_axes, i,k,i,l);
      }

    __DMP( "Solved = " << solved );
    if ( _dim == 3 ) {
      __DMP( " Eigen " << __TENSOR(0,0)<<", "<<__TENSOR(1,1)<<", "<<__TENSOR(2,2) );
      for ( int ii=0; ii <3; ++ii )
        __DMP( ii << ": " << __AXIS(ii)[0] << ", " << __AXIS(ii)[1] << ", " << __AXIS(ii)[2] );
    }
    else {
      __DMP( " Eigen " << __TENSOR(0,0) << ", " << __TENSOR(1,1) );
      for ( int ii=0; ii <2; ++ii )
        __DMP( ii << ": " << __AXIS(ii)[0] << ", " << __AXIS(ii)[1] );
    }

    return solved;
  }

  //================================================================================
  /*!
   * \brief Return true if two minmaxes do not intersect
   */
  //================================================================================

  inline bool isMinMaxOut(const double* minmax1,
                          const double* minmax2,
                          int           dim)
  {
    for ( int i = 0; i < dim; ++i )
      {
        if ( minmax1[i*2] > minmax2[i*2+1] ||
             minmax1[i*2+1] < minmax2[i*2] )
          return true;
      }
    return false;
  }

} // noname namespace

namespace INTERP_KERNEL
{

  //================================================================================
  /*!
   * \brief Creates empty box intended to further initalization via setData()
   */
  //================================================================================

  DirectedBoundingBox::DirectedBoundingBox():_dim(0)
  {
  }

  //================================================================================
  /*!
   * \brief Creates bounding box of a mesh
   *  \param pts - coordinates of points in full interlace
   *  \param numPts - number of points in the mesh
   *  \param dim - space dimension
   */
  //================================================================================

  DirectedBoundingBox::DirectedBoundingBox(const double*  pts,
                                           const unsigned numPts,
                                           const unsigned dim)
    : _dim(dim), _axes(dim*dim), _minmax(2*dim)
  {
    // init box extremities
    for ( unsigned i = 0; i < _dim; ++i )
      _minmax[1+i*2] = -numeric_limits<double>::max(),
        _minmax[i*2] =  numeric_limits<double>::max();

    if ( numPts < 1 ) return;

    __DMP( "DirectedBoundingBox " << __MYID );

    const double* coord = pts;
    const double* coordEnd = coord + numPts * dim;

    // compute gravity center of points
    double gc[3] = {0,0,0};
    if ( dim > 1 )
      {
        for ( coord = pts; coord < coordEnd; )
          for ( int i = 0; i < (int)dim; ++i )
            gc[i] += *coord++;
        for ( int j = 0; j < (int)dim; ++j )
          gc[j] /= numPts;

      }

    // compute axes and box extremities
    vector<double> tensor( dim * dim, 0.);
    switch ( dim )
      {
      case 3:
        for ( coord = pts; coord < coordEnd; coord += dim )
          addPointToInertiaTensor3D( coord, gc, tensor );

        //computeAxes3D(tensor);
        JacobiEigenvectorsSearch(_dim, tensor, _axes);

        for ( coord = pts; coord < coordEnd; coord += dim )
          addPointToBox( coord );

        break;

      case 2:
        for ( coord = pts; coord < coordEnd; coord += dim )
          addPointToInertiaTensor2D( coord, gc, tensor );

        //computeAxes2D(tensor);
        JacobiEigenvectorsSearch(_dim, tensor, _axes);

        for ( coord = pts; coord < coordEnd; coord += dim )
          addPointToBox( coord );

        break;

      default:
        for ( coord = pts; coord < coordEnd; coord += dim )
          {
            if ( *coord < _minmax[0] ) _minmax[0] = *coord;
            if ( *coord > _minmax[1] ) _minmax[1] = *coord;
          }
      }
  }

  //================================================================================
  /*!
   * \brief Creates bounding box of an element
   *  \param pts - coordinates of points of element
   *  \param numPts - number of points in the element
   *  \param dim - space dimension
   */
  //================================================================================

  DirectedBoundingBox::DirectedBoundingBox(const double** pts,
                                           const unsigned numPts,
                                           const unsigned dim)
    : _dim(dim), _axes(dim*dim), _minmax(2*dim)
  {
    // init box extremities
    for ( unsigned i = 0; i < _dim; ++i )
      _minmax[1+i*2] = -numeric_limits<double>::max(),
        _minmax[i*2] =  numeric_limits<double>::max();

    if ( numPts < 1 ) return;

    __DMP( "DirectedBoundingBox " << __MYID );

    // compute gravity center of points
    double gc[3] = {0,0,0};
    if ( dim > 1 )
      {
        for ( unsigned i = 0; i < numPts; ++i )
          for ( int j = 0; j < (int)dim; ++j )
            gc[j] += pts[i][j];
        for ( int j = 0; j < (int)dim; ++j )
          gc[j] /= numPts;
      }

    // compute axes and box extremities
    vector<double> tensor( dim * dim, 0.);
    switch ( dim )
      {
      case 3:
        for ( unsigned i = 0; i < numPts; ++i )
          addPointToInertiaTensor3D( pts[i], gc, tensor );

        //computeAxes3D(tensor);
        JacobiEigenvectorsSearch(_dim, tensor, _axes);

        for ( unsigned i = 0; i < numPts; ++i )
          addPointToBox( pts[i] );

        break;
      case 2:
        for ( unsigned i = 0; i < numPts; ++i )
          addPointToInertiaTensor2D( pts[i], gc, tensor );

        //computeAxes2D(tensor);
        JacobiEigenvectorsSearch(_dim, tensor, _axes);

        for ( unsigned i = 0; i < numPts; ++i )
          addPointToBox( pts[i] );

        break;
      default:
        for ( unsigned i = 0; i < numPts; ++i )
          {
            if ( pts[i][0] < _minmax[0] ) _minmax[0] = pts[i][0];
            if ( pts[i][0] > _minmax[1] ) _minmax[1] = pts[i][0];
          }
        _axes[0] = 1.0;
      }
  }

  //================================================================================
  /*!
   * \brief Compute eigenvectors of inertia tensor
   */
  //================================================================================

  // void DirectedBoundingBox::computeAxes3D(const std::vector<double>& tensor)
//   {
//     // compute principal moments of inertia which are eigenvalues of the tensor
//     double eig[3];
//     {
//       // coefficients of polynomial equation det(tensor-eig*I) = 0
//       double a = -1;
//       double b = __TENSOR(0,0)+__TENSOR(1,1)+__TENSOR(2,2);
//       double c =
//         __TENSOR(0,1)*__TENSOR(0,1) +
//         __TENSOR(0,2)*__TENSOR(0,2) +
//         __TENSOR(1,2)*__TENSOR(1,2) -
//         __TENSOR(0,0)*__TENSOR(1,1) -
//         __TENSOR(0,0)*__TENSOR(2,2) -
//         __TENSOR(1,1)*__TENSOR(2,2);
//       double d =
//         __TENSOR(0,0)*__TENSOR(1,1)*__TENSOR(2,2) -
//         __TENSOR(0,0)*__TENSOR(1,2)*__TENSOR(1,2) -
//         __TENSOR(1,1)*__TENSOR(0,2)*__TENSOR(0,2) -
//         __TENSOR(2,2)*__TENSOR(0,1)*__TENSOR(0,1) +
//         __TENSOR(0,1)*__TENSOR(0,2)*__TENSOR(1,2)*2;

//       // find eigenvalues which are roots of characteristic polynomial
//       double x = (3*c/a - b*b/(a*a))/3;
//       double y = (2*b*b*b/(a*a*a) - 9*b*c/(a*a) + 27*d/a)/27;
//       double z = y*y/4 + x*x*x/27;

//       double i = sqrt(y*y/4 - z) + 1e-300;
//       double j = -pow(i,1/3.);
//       double y2 = -y/(2*i);
//       if ( y2 > 1.0) y2 = 1.; else if ( y2 < -1.0) y2 = -1.;
//       double k = acos(y2);
//       double m = cos(k/3);
//       double n = sqrt(3)*sin(k/3);
//       double p = -b/(3*a);

//       eig[0] = -2*j*m + p;
//       eig[1] = j *(m + n) + p;
//       eig[2] = j *(m - n) + p;
//     }
//     // compute eigenvector of the tensor at each eigenvalue
//     // by solving system [tensor-eig*I]*[axis] = 0
//     bool ok = true;
//     __DMP( "Tensor : {"
//          << "{ "<<__TENSOR(0,0) << ", "<<__TENSOR(0,1) << ", "<<__TENSOR(0,2) << "} "
//          << "{ "<<__TENSOR(1,0) << ", "<<__TENSOR(1,1) << ", "<<__TENSOR(1,2) << "} "
//          << "{ "<<__TENSOR(2,0) << ", "<<__TENSOR(2,1) << ", "<<__TENSOR(2,2) << "}} ");
//     for ( int i = 0; i < 3 && ok; ++i ) // loop on 3 eigenvalues
//       {
//         // [tensor-eig*I]
//         double T[3][3]=
//           {{ __TENSOR(0,0)-eig[i],__TENSOR(0,1),       __TENSOR(0,2),      },
//            { __TENSOR(0,1),       __TENSOR(1,1)-eig[i],__TENSOR(1,2),      },
//            { __TENSOR(0,2),       __TENSOR(1,2),       __TENSOR(2,2)-eig[i]}};
//         // The determinant of T is zero, so that the equations are not linearly independent.
//         // Therefore, we assign an arbitrary value (1.) to i-th component of eigenvector
//         // and use two of the equations to compute the other two components
//         double M[2][3], sol[2];
//         for ( int j = 0, c = 0; j < 3; ++j )
//           if ( i == j )
//             M[0][2] = -T[0][j], M[1][2] = -T[1][j];
//           else
//             M[0][c] = T[0][j], M[1][c] = T[1][j], c++;

//         ok = solveSystemOfEquations<2>( M, sol );

//         double* eigenVec = __AXIS(i);
//         for ( int j = 0, c = 0; j < 3; ++j )
//           eigenVec[j] = ( i == j ) ? 1. : sol[c++];

//         // normilize
//         double size = sqrt(eigenVec[0]*eigenVec[0] +
//                            eigenVec[1]*eigenVec[1] +
//                            eigenVec[2]*eigenVec[2] );
//         if ((ok = (size > numeric_limits<double>::min() )))
//           {
//             eigenVec[0] /= size;
//             eigenVec[1] /= size;
//             eigenVec[2] /= size;
//           }
//       }
//     if ( !ok )
//       {
//         __DMP( " solve3EquationSystem() - KO " );
//         _axes = vector<double>( _dim*_dim, 0);
//         __AXIS(0)[0] = __AXIS(1)[1] = __AXIS(2)[2] = 1.;
//       }
//     __DMP( " Eigen " << eig[0] << ", " << eig[1] << ", " << eig[2] );
//     for ( int i=0; i <3; ++i )
//       __DMP( i << ": " << __AXIS(i)[0] << ", " << __AXIS(i)[1] << ", " << __AXIS(i)[2] );

//     double* a0 = __AXIS(0), *a1 = __AXIS(1);
//     double cross[3] = { a0[1]*a1[2]-a1[1]*a0[2],
//                         a0[2]*a1[0]-a1[2]*a0[0],
//                         a0[0]*a1[1]-a1[0]*a0[1] };
//     __DMP( " Cross a1^a2 " << cross[0] << ", " << cross[1] << ", " << cross[2] );
//   }

  //================================================================================
  /*!
   * \brief Compute eigenvectors of inertia tensor
   */
  //================================================================================

  // void DirectedBoundingBox::computeAxes2D(const std::vector<double>& tensor)
//   {
//     // compute principal moments of inertia which are eigenvalues of the tensor
//     // by solving square equation det(tensor-eig*I)
//     double X = (__TENSOR(0,0)+__TENSOR(1,1))/2;
//     double Y = sqrt(4*__TENSOR(0,1)*__TENSOR(0,1) +
//                     (__TENSOR(0,0)-__TENSOR(1,1)) * (__TENSOR(0,0)-__TENSOR(1,1)))/2;
//     double eig[2] =
//       {
//         X + Y,
//         X - Y
//       };
//     // compute eigenvector of the tensor at each eigenvalue
//     // by solving system [tensor-eig*I]*[axis] = 0
//     bool ok = true;
//     for ( int i = 0; i < 2 && ok; ++i )
//       {
//         // [tensor-eig*I]
//         double T[2][2]=
//           {{ __TENSOR(0,0)-eig[i],__TENSOR(0,1)        },
//            { __TENSOR(0,1),       __TENSOR(1,1)-eig[i] }};

//         // The determinant of T is zero, so that the equations are not linearly independent.
//         // Therefore, we assign an arbitrary value (1.) to i-th component of eigenvector
//         // and use one equation to compute the other component
//         double* eigenVec = __AXIS(i);
//         eigenVec[i] = 1.;
//         int j = 1-i;
//         if ((ok = ( fabs( T[j][j] ) > numeric_limits<double>::min() )))
//           eigenVec[j] = -T[j][i] / T[j][j];
//       }
//     if ( !ok )
//       {
//         _axes = vector<double>( _dim*_dim, 0);
//         __AXIS(0)[0] = __AXIS(1)[1] = 1.;
//       }
//   }

  //================================================================================
  /*!
   * \brief Convert point coordinates into local coordinate system of the box
   */
  //================================================================================

  void DirectedBoundingBox::toLocalCS(const double* p, double* pLoc) const
  {
    switch ( _dim )
      {
      case 3:
        pLoc[0] = dotprod<3>( p, __AXIS(0));
        pLoc[1] = dotprod<3>( p, __AXIS(1));
        pLoc[2] = dotprod<3>( p, __AXIS(2));
        break;
      case 2:
        pLoc[0] = dotprod<2>( p, __AXIS(0));
        pLoc[1] = dotprod<2>( p, __AXIS(1));
        break;
      default:
        pLoc[0] = p[0];
      }
  }

  //================================================================================
  /*!
   * \brief Convert point coordinates from local coordinate system of the box to global CS
   */
  //================================================================================

  void DirectedBoundingBox::fromLocalCS(const double* p, double* pGlob) const
  {
    switch ( _dim )
      {
      case 3:
        pGlob[0] = p[0] * __AXIS(0)[0] + p[1] * __AXIS(1)[0] + p[2] * __AXIS(2)[0];
        pGlob[1] = p[0] * __AXIS(0)[1] + p[1] * __AXIS(1)[1] + p[2] * __AXIS(2)[1];
        pGlob[2] = p[0] * __AXIS(0)[2] + p[1] * __AXIS(1)[2] + p[2] * __AXIS(2)[2];
        break;
      case 2:
        pGlob[0] = p[0] * __AXIS(0)[0] + p[1] * __AXIS(1)[0];
        pGlob[1] = p[0] * __AXIS(0)[1] + p[1] * __AXIS(1)[1];
        break;
      default:
        pGlob[0] = p[0];
      }
  }

  //================================================================================
  /*!
   * \brief Enlarge box size by given value
   */
  //================================================================================

  void DirectedBoundingBox::enlarge(const double tol)
  {
    for ( unsigned i = 0; i < _dim; ++i )
      __MIN(i) -= tol, __MAX(i) += tol;
  }

  //================================================================================
  /*!
   * \brief Return coordinates of corners of bounding box
   */
  //================================================================================

  void DirectedBoundingBox::getCorners(std::vector<double>& corners,
                                       const double*        minmax) const
  {
    int iC, nbCorners = 1;
    for ( int i=0;i<(int)_dim;++i ) nbCorners *= 2;
    corners.resize( nbCorners * _dim );
    // each coordinate is filled with either min or max, nbSwap is number of corners
    // after which min and max swap
    int nbSwap = nbCorners/2;
    for ( unsigned i = 0; i < _dim; ++i )
      {
        iC = 0;
        while ( iC < nbCorners )
          {
            for (int j = 0; j < nbSwap; ++j, ++iC ) corners[iC*_dim+i] = minmax[i*2];
            for (int j = 0; j < nbSwap; ++j, ++iC ) corners[iC*_dim+i] = minmax[i*2+1];
          }
        nbSwap /= 2;
      }
  }

  //================================================================================
  /*!
   * \brief Test if this box intersects with the other
   *  \retval bool - true if there is no intersection
   */
  //================================================================================

  bool DirectedBoundingBox::isDisjointWith(const DirectedBoundingBox& box) const
  {
    if ( _dim < 1 || box._dim < 1 ) return false; // empty box includes all
    if ( _dim == 1 )
      return isMinMaxOut( &box._minmax[0], &this->_minmax[0], _dim );

    // boxes are disjoined if their minmaxes in local CS of either of boxes do not intersect
    for ( int isThisCS = 0; isThisCS < 2; ++isThisCS )
      {
        const DirectedBoundingBox* axisBox   = isThisCS ? this : &box;
        const DirectedBoundingBox* cornerBox = isThisCS ? &box : this;

        // find minmax of cornerBox in the CS of axisBox

        DirectedBoundingBox mmBox((double*)0,0,_dim); //!< empty box with CS == axisBox->_axes
        mmBox._axes = axisBox->_axes;

        vector<double> corners;
        getCorners( corners, &cornerBox->_minmax[0] );

        double globCorner[3];
        for ( int iC = 0, nC = corners.size()/_dim; iC < nC; ++iC)
          {
            cornerBox->fromLocalCS( &corners[iC*_dim], globCorner );
            mmBox.addPointToBox( globCorner );
          }
        if ( isMinMaxOut( &mmBox._minmax[0], &axisBox->_minmax[0], _dim ))
          return true;
      }
    return false;
  }

  //================================================================================
  /*!
   * \brief Test if this box intersects with an non-directed box
   *  \retval bool - true if there is no intersection
   */
  //================================================================================

  bool DirectedBoundingBox::isDisjointWith(const double* box) const
  {
    if ( _dim < 1 ) return false; // empty box includes all
    if ( _dim == 1 )
      return isMinMaxOut( &_minmax[0], box, _dim );

    // boxes are disjoined if their minmaxes in local CS of either of boxes do not intersect

    // compare minmaxes in locals CS of this directed box
    {
      vector<double> cornersOther;
      getCorners( cornersOther, box );
      DirectedBoundingBox mmBox((double*)0,0,_dim); //!< empty box with CS == this->_axes
      mmBox._axes = this->_axes;
      for ( int iC = 0, nC = cornersOther.size()/_dim; iC < nC; ++iC)
        mmBox.addPointToBox( &cornersOther[iC*_dim] );

      if ( isMinMaxOut( &mmBox._minmax[0], &this->_minmax[0], _dim ))
        return true;
    }

    // compare minmaxes in global CS
    {
      vector<double> cornersThis;
      getCorners( cornersThis, &_minmax[0] );
      DirectedBoundingBox mmBox((double*)0,0,_dim); //!< initailized _minmax
      double globCorner[3];
      for ( int iC = 0, nC = cornersThis.size()/_dim; iC < nC; ++iC)
        {
          fromLocalCS( &cornersThis[iC*_dim], globCorner );
          for ( int i = 0; i < (int)_dim; ++i )
            {
              if ( globCorner[i] < mmBox._minmax[i*2] )   mmBox._minmax[i*2] = globCorner[i];
              if ( globCorner[i] > mmBox._minmax[i*2+1] ) mmBox._minmax[i*2+1] = globCorner[i];
            }
        }
      if ( isMinMaxOut( &mmBox._minmax[0], box, _dim ))
        return true;
    }
    return false;
  }

  //================================================================================
  /*!
   * \brief Return true if given point is out of this box
   */
  //================================================================================

  bool DirectedBoundingBox::isOut(const double* point) const
  {
    if ( _dim < 1 ) return false; // empty box includes all

    double pLoc[3];
    toLocalCS( point, pLoc );
    bool out = isLocalOut( pLoc );
#ifdef _DEBUG_
    switch (_dim)
      {
      case 3:
        __DMP(__MYID<<": "<<point[0]<<", "<<point[1]<<", "<<point[2]<<" "<<(out?"OUT":"IN"));break;
      case 2:
        __DMP(__MYID<<": "<<point[0]<<", "<<point[1]<<" "<<(out?"OUT":"IN"));break;
      case 1:
        __DMP(__MYID<<": "<<point[0]<<" "<<(out?"OUT":"IN"));break;
      }
#endif
    return out;
  }

  //================================================================================
  /*!
   * \brief Return array of internal data
   */
  //================================================================================

  vector<double> DirectedBoundingBox::getData() const
  {
    vector<double> data(1, _dim);
    if ( _dim > 0 )
    {
      data.insert( data.end(), &_axes[0], &_axes[0] + _axes.size());
      data.insert( data.end(), &_minmax[0], &_minmax[0] + _minmax.size());
    }
    if ( data.size() < (unsigned)dataSize( _dim ))
      data.resize( dataSize( _dim ), 0 );
    return data;
  }

  //================================================================================
  /*!
   * \brief Initializes self with data retrieved via getData()
   */
  //================================================================================

  void DirectedBoundingBox::setData(const double* data)
  {
    _dim = unsigned( *data++ );
    if ( _dim > 0 )
      {
        _axes.assign( data, data+_dim*_dim ); data += _dim*_dim;
        _minmax.assign( data, data+2*_dim );
      }
    else
      {
        _axes.clear();
        _minmax.clear();
      }
  }

  //================================================================================
  /*!
   * \brief Return size of internal data returned by getData() depending on space dim
   */
  //================================================================================

  int DirectedBoundingBox::dataSize(int dim)
  {
    return 1 + dim*dim + 2*dim; // : _dim + _axes + _minmax
  }
}
