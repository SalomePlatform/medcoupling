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
// Author : Anthony Geay (CEA/DEN)

#include "InterpKernelMatrixTools.hxx"
#include "InterpKernelException.hxx"
#include "InterpKernelAutoPtr.hxx"

#include <sstream>
#include <algorithm>

namespace INTERP_KERNEL
{
  /*
   *  Computes the dot product of two vectors.
   *  This routine uses unrolled loops for increments equal to one.
   *
   *  Reference:
   *
   *    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
   *    LINPACK User's Guide,
   *    SIAM, 1979,
   *    ISBN13: 978-0-898711-72-1,
   *    LC: QA214.L56.
   *
   *    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
   *    Basic Linear Algebra Subprograms for Fortran Usage,
   *    Algorithm 539, 
   *    ACM Transactions on Mathematical Software, 
   *    Volume 5, Number 3, September 1979, pages 308-323.
   *
   *  \param [in] n the number of entries in the vectors.
   *  \param [in] dx the first vector.
   *  \param [in] incx the increment between successive entries in \a dx.
   *  \param [in] dy the second vector.
   *  \param [in] incy the increment between successive entries in \a dy.
   *  \return the sum of the product of the corresponding entries of \a dx and \a dy.
   */
  double ddot(int n, const double *dx, int incx, const double *dy, int incy)
  {
    double dtemp=0.0;
    int i,ix,iy,m;
    if(n<=0)
      return dtemp;
    // Code for unequal increments or equal increments not equal to 1.
    if(incx!=1 || incy!=1)
      {
        if (incx>=0)
          ix=0;
        else
          ix=(-n+1)*incx;

        if(incy>=0)
          iy=0;
        else
          iy=(-n+1)*incy;
        for(i=0;i<n;i++,ix+=incx,iy+=incy)
          dtemp+=dx[ix]*dy[iy];
      }
    //Code for both increments equal to 1.
    else
      {
        m=n%5;
        for(i=0;i<m;i++)
          dtemp+=dx[i]*dy[i];
        for (i=m;i<n;i+=5)
          dtemp+=dx[i]*dy[i]+dx[i+1]*dy[i+1]+dx[i+2]*dy[i+2]+dx[i+3]*dy[i+3]+dx[i+4]*dy[i+4];
      }
    return dtemp;
  }


  void dscal(int n, double sa, double *x, int incx)
  {
    int i,ix,m;

    if(n<=0) { }
    else if(incx==1)
      {
        m=n%5;
        for(i=0;i<m;i++ )
          x[i]=sa*x[i];

        for(i=m;i<n;i+=5)
          {
            x[i]  =sa*x[i];
            x[i+1]=sa*x[i+1];
            x[i+2]=sa*x[i+2];
            x[i+3]=sa*x[i+3];
            x[i+4]=sa*x[i+4];
          }
      }
    else
      {
        if(0<=incx)
          ix=0;
        else
          ix=(-n+1)*incx;

        for(i=0;i<n;i++,ix+=incx)
          x[ix]=sa*x[ix];
      }
  }

  void daxpy(int n, double da, const double *dx, int incx, double *dy, int incy)
  {
    int i,ix,iy,m;
    if (n<=0)
      return;
    if (da==0.0)
      return;
    // Code for unequal increments or equal increments not equal to 1.
    if(incx!=1 || incy!=1)
      {
        if(0<=incx)
          ix=0;
        else
          ix=(-n+1)*incx;

        if(0<=incy)
          iy=0;
        else
          iy=(-n+1)*incy;

        for(i=0;i<n;i++,ix+=incx,iy+=incy)
          dy[iy]=dy[iy]+da*dx[ix];
      }
    // Code for both increments equal to 1.
    else
      {
        m=n%4;
        for (i=0;i<m;i++)
          dy[i]=dy[i]+da*dx[i];
        for(i=m;i<n;i+=4)
          {
            dy[i  ]=dy[i  ]+da*dx[i  ];
            dy[i+1]=dy[i+1]+da*dx[i+1];
            dy[i+2]=dy[i+2]+da*dx[i+2];
            dy[i+3]=dy[i+3]+da*dx[i+3];
          }
      }
  }

  double r8_abs(double x)
  {
    if(x>=0.0)
      return x;
    else
      return -x;
  }

  void dswap(int n, double *x, int incx, double *y, int incy)
  {
    int i,ix,iy,m;
    double temp;

    if(n<=0) { }
    else if(incx==1 && incy==1)
      {
        m=n%3;
        for(i=0;i<m;i++)
          { temp=x[i]; x[i]=y[i]; y[i]=temp; }
        for(i=m;i<n;i+=3)
          { 
            temp=x[i]; x[i]=y[i]; y[i]=temp;
            temp=x[i+1]; x[i+1]=y[i+1]; y[i+1]=temp;
            temp=x[i+2]; x[i+2]=y[i+2]; y[i+2]=temp;
          }
      }
    else
      {
        if(0<=incx)
          ix=0;
        else
          ix=(-n+1)*incx;
        if(0<=incy)
          iy=0;
        else
          iy=(-n+1)*incy;
        for(i=0;i<n;i++,ix+=incx,iy+=incy)
          {
            temp=x[ix]; x[ix]=y[iy];
          }
      }
  }

  /*
   *  Finds the index of the vector element of maximum absolute value.
   *  \warning This index is a 1-based index, not a 0-based index !
   *
   *  Reference:
   *    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
   *    LINPACK User's Guide,
   *    SIAM, 1979,
   *    ISBN13: 978-0-898711-72-1,
   *    LC: QA214.L56.
   *
   *    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
   *    Basic Linear Algebra Subprograms for Fortran Usage,
   *    Algorithm 539,
   *    ACM Transactions on Mathematical Software,
   *    Volume 5, Number 3, September 1979, pages 308-323.
   *
   *    \param [in] n the number of entries in the vector.
   *    \param [in] dx the vector to be examined.
   *    \param [in] incx the increment between successive entries of SX.
   *    \return the index of the element of maximum absolute value (in C convention).
   */
  int idamax(int n, const double *dx, int incx)
  {
    double dmax;
    int i,ix,value;
    value=-1;
    if ( n < 1 || incx <= 0 )
      return value;
    value=0;
    if(n==1)
      return value;
    if(incx==1)
      {
        dmax=r8_abs(dx[0]);
        for(i=1;i<n;i++)
          {
            if(dmax<r8_abs(dx[i]))
              {
                value=i;
                dmax=r8_abs(dx[i]);
              }
          }
      }
    else
      {
        ix=0;
        dmax=r8_abs(dx[0]);
        ix+=incx;
        for(i=1;i<n;i++)
          {
            if(dmax<r8_abs(dx[ix]))
              {
                value=i;
                dmax=r8_abs(dx[ix]);
              }
            ix+=incx;
          }
      }
    return value;
  }


  /*
   *  Purpose:
   *
   *     factors a real general matrix.
   *
   *  Reference:
   *
   *    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
   *    LINPACK User's Guide,
   *    SIAM, (Society for Industrial and Applied Mathematics),
   *    3600 University City Science Center,
   *    Philadelphia, PA, 19104-2688.
   *    ISBN 0-89871-172-X
   *
   *  Parameters:
   *
   *  \param [in,out] a a matrix of size LDA*N. On intput, the matrix to be factored.
   *    On output, an upper triangular matrix and the multipliers used to obtain
   *    it.  The factorization can be written A=L*U, where L is a product of
   *    permutation and unit lower triangular matrices, and U is upper triangular.
   *
   *  \param [in] lda the leading dimension of \a a.
   *  \param [in] n the order of the matrix \a a.
   *  \param [out] ipvt the pivot indices of size \a n.
   *  \return the singularity indicator.
   *  - 0, normal value.
   *  - K, if U(K-1,K-1) == 0.  This is not an error condition for this subroutine,
   *    but it does indicate that DGESL or DGEDI will divide by zero if called.
   */
  int dgefa(double *a, int lda, int n, int *ipvt)
  {
    int info=0;
    int l;
    double t;
    // Gaussian elimination with partial pivoting.
    for(int k=0;k<n-1;k++)
      {
        //  Find L=pivot index.
        l=idamax(n-k,a+k+k*lda,1)+k;
        ipvt[k]=l;
        // Zero pivot implies this column already triangularized.
        if(a[l+k*lda]==0.0)
          {
            info=k;
            continue;
          }
        //Interchange if necessary.
        if(l!=k)
          {
            t=a[l+k*lda];
            a[l+k*lda]=a[k+k*lda];
            a[k+k*lda]=t;
          }
        // Compute multipliers.
        t=-1.0/a[k+k*lda];
        dscal(n-k-1,t,a+k+1+k*lda,1);
        // Row elimination with column indexing.
        for(int j=k+1;j<n;j++)
          {
            t=a[l+j*lda];
            if(l!=k)
              {
                a[l+j*lda]=a[k+j*lda];
                a[k+j*lda]=t;
              }
            daxpy(n-k-1,t,a+k+1+k*lda,1,a+k+1+j*lda,1);
          }
      }
    ipvt[n-1]=n-1;
    if(a[n-1+(n-1)*lda]==0.0)
      info=n;
    return info;
  }

  /*
   *  - Purpose: computes inverse (job=1) of a matrix factored by dgefa.
   *
   *  - Discussion: A division by zero will occur if the input factor contains
   *    a zero on the diagonal. It will not occur if the subroutines are called correctly
   *    and dgfa has set info to 0.
   *
   *  Reference:
   *
   *    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
   *    LINPACK User's Guide,
   *    SIAM, (Society for Industrial and Applied Mathematics),
   *    3600 University City Science Center,
   *    Philadelphia, PA, 19104-2688.
   *    ISBN 0-89871-172-X
   *
   *  Parameters:
   *
   *    \param [in,out] a matrix of size lda*n, on input, the LU factor information, as output of dgfa. On output, the inverse matrix.
   *    \param [in] lda, the leading dimension of the array \a a.
   *    \param [in] n, the order of the matrix \a a.
   *    \param [in] ipvt, the pivot vector from dgfa.
   *    \param [in,out] work a work array of size at least equal to \a n.
   */
  void dgedi(double *a, int lda, int n, const int *ipvt, double *work)
  {
    double t;
    for(int k=0;k<n;k++)
      {
        a[k+k*lda]=1.0/a[k+k*lda];
        t=-a[k+k*lda];
        dscal(k,t,a+0+k*lda,1);
        for(int j=k+1;j<n;j++)
          {
            t=a[k+j*lda];
            a[k+j*lda]=0.0;
            daxpy(k+1,t,a+0+k*lda,1,a+0+j*lda,1);
          }
      }
    // Form inverse(U) * inverse(L).
    for(int k=n-2;k>=0;k--)
      {
        for(int i=k+1;i<n;i++)
          {
            work[i]=a[i+k*lda];
            a[i+k*lda]=0.0;
          }

        for(int j=k+1;j<n;j++)
          {
            t=work[j];
            daxpy(n,t,a+0+j*lda,1,a+0+k*lda,1);
          }
        int l=ipvt[k];
        if(l!=k-1)
          dswap(n,a+0+k*lda,1,a+0+l*lda,1);
      }
  }

  void matrixProduct(const double *A, int n1, int p1, const double *B, int n2, int p2, double *C)
  {
    if(p1!=n2)
      {
        std::ostringstream oss; oss << "matrixProduct : the size of input matrix are not coherent the nb of cols of input matrix #0 is " <<  p1 << " whereas the number of rows of input matrix #1 is " << n2 << " !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    for(int i=0;i<n1;i++)
      {
        for(int j=0;j<p2;j++)
          {
            C[i*p2+j] = 0.;
            for(int k=0;k<p1;k++)
              C[i*p2+j]+=A[i*p1+k]*B[k*p2+j];
          }
      }
  }

  void inverseMatrix(const double *A, int n, double *iA)
  {
    INTERP_KERNEL::AutoPtr<int> ipvt=new int[n];
    INTERP_KERNEL::AutoPtr<double> work=new double[n*n];
    std::copy(A,A+n*n,iA);
    dgefa(iA,n,n,ipvt);
    dgedi(iA,n,n,ipvt,work);
  }
}
