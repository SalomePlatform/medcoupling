// Copyright (C) 2022-2025  CEA, EDF
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
void LineSearches(const std::vector<double> &xold, std::vector<double> &dx, std::vector<double> &x,
                  T &func, const double ALF = 0.1, const mcIdType MAXIT = 10) {
  const double TOLX = std::numeric_limits<double>::epsilon();
  const mcIdType n = xold.size();
  std::vector<double> &fvec = func.getVector();
  double p0, p1, f1, f0, rho_0, rho_1;
  double rho, rho_neg, rho_pos, rho_opt, rho_new, rho_cur;
  double f, f_opt, f_cur;
  bool b_pos;

  // fixed paramters - from code_aster
  const double rho_min = 1e-2, rho_max = 10., rho_excl = 0.9e-2;
  const double parmul = 3.0;

  // Doc:
  // https://codeaster.pages.pleiade.edf.fr/doc/docaster/manuals/man_r/r5/r5.03.01/Recherche_lin_aire.html
  // METHODE="MIXTE"

  // Compute residual.dot(increment) (and update solution)
  auto _f = [&fvec, &func, &dx, &xold, &n, &x](const double &rho) {
    for (mcIdType j = 0; j < n; j++)
      x[j] = xold[j] + rho * dx[j];
    const auto norm = func(x);

    double f = 0.0;
    for (mcIdType j = 0; j < n; j++)
      f += fvec[j] * dx[j];
    return f;
  };

  // project bound on admissible interval
  auto _proj = [rho_min, rho_max, rho_excl](double &rho) {
    const double rho_tmp = rho;
    if (rho_tmp < rho_min) {
      rho = rho_min;
    }
    if (rho_tmp > rho_max) {
      rho = rho_max;
    }
    if (rho_tmp < 0.0 && rho_tmp >= -rho_excl) {
      rho = -rho_excl;
    }
    if (rho_tmp >= 0 && rho_tmp <= rho_excl) {
      rho = rho_excl;
    }
  };

  // initial values
  const double f_old = _f(0.0);
  const double f_cvg = ALF * std::abs(f_old);
  const double sens = (f_old <= 0.0) ? 1.0 : -1.0;

  rho_opt = 1.0, rho = sens * 1.0;
  rho_neg = 0.0, rho_pos = std::numeric_limits<double>::signaling_NaN();
  f_opt = 10e100;

  rho_0 = 0.0, rho_1 = rho_0;
  f0 = sens * f_old, f1 = f0;

  b_pos = false;

  for (mcIdType its = 0; its < MAXIT; its++) {
    // Compute new residual
    f = _f(rho);

    rho_cur = sens * rho;
    f_cur = sens * f;

    // Store value
    rho_0 = rho_1, f0 = f1;
    rho_1 = rho_cur, f1 = f_cur;

    // Update bounds
    if (f_cur < 0.0) {
      rho_neg = rho_cur;
    } else {
      b_pos = true;
      rho_pos = rho_cur;
    }

    // Optimal solution until now ?
    if (std::abs(f_cur) < std::abs(f_opt)) {
      rho_opt = rho_cur;
      f_opt = f_cur;
      _proj(rho_opt);
    }

    // Search maximal bound
    if (b_pos) {
      if (std::abs(f1) >= std::abs(f0)) {
        // f is not decreased - use dichotomie
        rho_new = 0.5 * (rho_neg + rho_pos);
      } else {
        // linear interpolation
        if (std::abs(rho_1 - rho_0) > TOLX) {
          p1 = (f1 - f0) / (rho_1 - rho_0);
          p0 = f0 - p1 * rho_0;

          if (std::abs(p1) <= std::abs(f0) / (rho_pos + rho_0)) {
            rho_new = 0.5 * (rho_neg + rho_pos);
          } else {
            rho_new = -p0 / p1;
          }
        } else {
          // failed
          break;
        }
      }
    } else {
      rho_new = parmul * rho_cur;
    }

    // minimal bound
    if (rho_new < rho_neg) {
      if (b_pos) {
        rho_new = 0.5 * (rho_neg + rho_pos);
      } else {
        // failed
        break;
      }
    }

    // maximal bound
    if (b_pos && rho_new > rho_pos) {
      rho_new = 0.5 * (rho_neg + rho_pos);
    }

    // project bound
    _proj(rho_new);

    // update
    rho = sens * rho_new;

    // Test convergence ?
    if (std::abs(f_opt) <= f_cvg) {
      break;
    }
  }

  /* Return optimal value */
  f = _f(rho_opt);
}
template <class T> class JacobianCalculator {
private:
  const double EPS;
  T &func;

public:
  JacobianCalculator(T &funcc) : EPS(1.0e-8), func(funcc) {}
  INTERP_KERNEL::DenseMatrix operator()(const std::vector<double> &x,
                                        const std::vector<double> &fvec) {
    mcIdType n = x.size();
    INTERP_KERNEL::DenseMatrix df(n, n);
    std::vector<double> xh = x;
    for (mcIdType j = 0; j < n; j++) {
      double temp = xh[j];
      double h = EPS * std::abs(temp);
      if (h == 0.0)
        h = EPS;
      xh[j] = temp + h;
      h = xh[j] - temp;
      std::vector<double> f = func(xh);
      xh[j] = temp;
      for (mcIdType i = 0; i < n; i++)
        df[i][j] = (f[i] - fvec[i]) / h;
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
  void SolveWithNewtonWithJacobian(
      std::vector<double> &x, bool &check, T &vecfunc,
      std::function<void(const std::vector<double> &, const std::vector<double> &, INTERP_KERNEL::DenseMatrix &)> jacobian,
      const double TOLF = 1.0e-8,
      const mcIdType MAXITS = 200) {
    const double TOLX=std::numeric_limits<double>::epsilon();
    mcIdType i,j,its,n=x.size();
    double f, temp, test;
    std::vector<double> p(n), xold(n);
    INTERP_KERNEL::DenseMatrix fjac(n,n);
    FMin<T> fmin(vecfunc);
    std::vector<double> &fvec=fmin.getVector();
    f=fmin(x);
    test = 0.0;
    for (i=0;i<n;i++)
      if (std::abs(fvec[i]) > test) test=std::abs(fvec[i]);
    if (test < 0.01*TOLF)
    {
      check=false;
      return;
    }
    for (its=0;its<MAXITS;its++)
    {
      jacobian(x, fvec, fjac);
      for (i = 0; i < n; i++) xold[i] = x[i];
      for (i=0;i<n;i++) p[i] = -fvec[i];
      INTERP_KERNEL::LUDecomp alu(fjac);
      alu.solve(p,p);
      LineSearches(xold, p, x, fmin);
      test = 0.0;
      for (i = 0; i < n; i++) {
        if (std::abs(fvec[i]) > test)
          test = std::abs(fvec[i]);
      }
      if (test < TOLF) {
        check=false;
        return;
      }
      test=0.0;
      for (i=0;i<n;i++)
      {
        temp=(std::abs(x[i]-xold[i]))/std::fmax(std::abs(x[i]),1.0);
        if (temp > test) test=temp;
      }
      if (test < TOLX) {
        check = false;
        return;
      }
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
