// Copyright (C) 2022-2024  CEA, EDF
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

// Implementation coming from Numerical Recipes in C of 1994 (version 2.04) p 384
// Root finding and non linear sets of equation using line searches and backtracking

#include "InterpKernelDenseMatrix.hxx"
#include "InterpKernelLUDecomp.hxx"
#include "InterpKernelQRDecomp.hxx"
#include "InterpKernelException.hxx"
#include "MCIdType.hxx"

#include <vector>
#include <limits>
#include <cmath>

template<class T>
inline T sqr(const T a) { return a*a; }

namespace INTERP_KERNEL
{
  template <class T>
  void LineSearches(const std::vector<double>& xold, const double fold, const std::vector<double>& g, std::vector<double> &p,
  std::vector<double> &x, double &f, const double stpmax, bool &check, T &func)
  {
    const double ALF=1.0e-4, TOLX=std::numeric_limits<double>::epsilon();
    double a,alam,alam2=0.0,alamin,b,disc,f2=0.0;
    double rhs1,rhs2,slope=0.0,sum=0.0,temp,test,tmplam;
    mcIdType i,n=xold.size();
    check=false;
    for (i=0;i<n;i++) sum += p[i]*p[i];
    sum=std::sqrt(sum);
    if (sum > stpmax)
      for (i=0;i<n;i++)
        p[i] *= stpmax/sum;
    for (i=0;i<n;i++)
      slope += g[i]*p[i];
    if (slope >= 0.0) THROW_IK_EXCEPTION("Roundoff problem in LineSearches.");
    test=0.0;
    for (i=0;i<n;i++)
    {
      temp=std::abs(p[i])/std::fmax(std::abs(xold[i]),1.0);
      if (temp > test) test=temp;
    }
    alamin=TOLX/test;
    alam=1.0;
    for (;;)
    {
      for (i=0;i<n;i++) x[i]=xold[i]+alam*p[i];
      f=func(x);
      if (alam < alamin)
      {
        for (i=0;i<n;i++) x[i]=xold[i];
        check=true;
        return;
      } else if (f <= fold+ALF*alam*slope) return;
      else
      {
        if (alam == 1.0)
          tmplam = -slope/(2.0*(f-fold-slope));
        else
        {
          rhs1=f-fold-alam*slope;
          rhs2=f2-fold-alam2*slope;
          a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
          b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
          if (a == 0.0) tmplam = -slope/(2.0*b);
          else {
              disc=b*b-3.0*a*slope;
              if (disc < 0.0) tmplam=0.5*alam;
              else if (b <= 0.0) tmplam=(-b+std::sqrt(disc))/(3.0*a);
              else tmplam=-slope/(b+std::sqrt(disc));
          }
          if (tmplam>0.5*alam)
              tmplam=0.5*alam;
        }
      }
      alam2=alam;
      f2 = f;
      alam=std::fmax(tmplam,0.1*alam);
    }
  }
  template <class T>
  class JacobianCalculator
  {
  private:
    const double EPS;
    T &func;
  public:
    JacobianCalculator(T &funcc) : EPS(1.0e-8),func(funcc) {}
    INTERP_KERNEL::DenseMatrix operator() (const std::vector<double>& x, const std::vector<double>& fvec)
    {
      mcIdType n=x.size();
      INTERP_KERNEL::DenseMatrix df(n,n);
      std::vector<double> xh=x;
      for (mcIdType j=0;j<n;j++)
      {
        double temp=xh[j];
        double h=EPS*std::abs(temp);
        if (h == 0.0) h=EPS;
        xh[j]=temp+h;
        h=xh[j]-temp;
        std::vector<double> f=func(xh);
        xh[j]=temp;
        for (mcIdType i=0;i<n;i++)
          df[i][j]=(f[i]-fvec[i])/h;
      }
      return df;
    }
  };

  template <class T>
  class FMin
  {
  private:
    std::vector<double> fvec;
    T &func;
    mcIdType n;
  public:
    FMin(T &funcc) : func(funcc){}
    double operator() (const std::vector<double>& x)
    {
      n=x.size();
      double sum=0;
      fvec=func(x);
      for (mcIdType i=0;i<n;i++) sum += sqr(fvec[i]);
      return 0.5*sum;
    }
    std::vector<double>& getVector() { return fvec; }
  };

  /*!
   * check is false on normal return.
   * check is true if the routine has converged to a local minimum.
   */
  template <class T>
  void SolveWithNewtonWithJacobian(std::vector<double> &x, bool &check, T &vecfunc,
    std::function<void(const std::vector<double>&, const std::vector<double>&, INTERP_KERNEL::DenseMatrix&)> jacobian)
  {
    const mcIdType MAXITS=200;
    const double TOLF=1.0e-8,TOLMIN=1.0e-12,STPMX=100.0;
    const double TOLX=std::numeric_limits<double>::epsilon();
    mcIdType i,j,its,n=x.size();
    double den,f,fold,stpmax,sum,temp,test;
    std::vector<double> g(n),p(n),xold(n);
    INTERP_KERNEL::DenseMatrix fjac(n,n);
    FMin<T> fmin(vecfunc);
    std::vector<double> &fvec=fmin.getVector();
    f=fmin(x);
    test=0.0;
    for (i=0;i<n;i++)
      if (std::abs(fvec[i]) > test) test=std::abs(fvec[i]);
    if (test < 0.01*TOLF)
    {
      check=false;
      return;
    }
    sum=0.0;
    for (i=0;i<n;i++) sum += sqr(x[i]);
    stpmax=STPMX*std::fmax(std::sqrt(sum),double(n));
    for (its=0;its<MAXITS;its++)
    {
      jacobian(x,fvec,fjac);//fjac=fdjac(x,fvec);
      for (i=0;i<n;i++)
      {
        sum=0.0;
        for (j=0;j<n;j++) sum += fjac[j][i]*fvec[j];
        g[i]=sum;
      }
      for (i=0;i<n;i++) xold[i]=x[i];
      fold=f;
      for (i=0;i<n;i++) p[i] = -fvec[i];
      INTERP_KERNEL::LUDecomp alu(fjac);
      alu.solve(p,p);
      LineSearches(xold,fold,g,p,x,f,stpmax,check,fmin);
      test=0.0;
      for (i=0;i<n;i++)
          if (std::abs(fvec[i]) > test) test=std::abs(fvec[i]);
      if (test < TOLF)
      {
        check=false;
        return;
      }
      if (check)
      {
        test=0.0;
        den=std::fmax(f,0.5*double(n));
        for (i=0;i<n;i++) {
            temp=std::abs(g[i])*std::fmax(std::abs(x[i]),1.0)/den;
            if (temp > test) test=temp;
        }
        check=(test < TOLMIN);
        if( check )
          return;
        else
          continue;
      }
      test=0.0;
      for (i=0;i<n;i++)
      {
        temp=(std::abs(x[i]-xold[i]))/std::fmax(std::abs(x[i]),1.0);
        if (temp > test) test=temp;
      }
      if (test < TOLX)
          return;
    }
    THROW_IK_EXCEPTION("MAXITS exceeded in SolveWithNewtonWithJacobian");
  }

  /*!
   * check is false on normal return.
   * check is true if the routine has converged to a local minimum.
   */
  template <class T>
  void SolveWithNewton(std::vector<double> &x, bool &check, T &vecfunc)
  {
    JacobianCalculator<T> fdjac(vecfunc);
    auto myJacobian = [&fdjac,vecfunc](const std::vector<double>& x, const std::vector<double>& fvec, INTERP_KERNEL::DenseMatrix& fjac)
    {
      fjac = fdjac(x,fvec);
    };
    SolveWithNewtonWithJacobian(x,check,vecfunc,myJacobian);
  }
}
