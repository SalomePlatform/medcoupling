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

// Implementation coming from Numerical Recipes in C of 1994 (version 2.04)

#include "InterpKernelQRDecomp.hxx"
#include "InterpKernelException.hxx"

#include <cmath>

using namespace INTERP_KERNEL;

template<class T>
inline T sqr(const T a) { return a*a; }

template<class T>
inline T sign(const T &a, const T &b) { return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a); }

QRDecomp::QRDecomp(const INTERP_KERNEL::DenseMatrix& a): n(a.nrows()), qt(n,n), r(a), sing(false)
{
	mcIdType i,j,k;
	std::vector<double> c(n), d(n);
	double scale,sigma,sum,tau;
	for (k=0;k<n-1;k++) {
		scale=0.0;
		for (i=k;i<n;i++) scale=std::fmax(scale,std::abs(r[i][k]));
		if (scale == 0.0) {
			sing=true;
			c[k]=d[k]=0.0;
		} else {
			for (i=k;i<n;i++) r[i][k] /= scale;
			for (sum=0.0,i=k;i<n;i++) sum += sqr(r[i][k]);
			sigma=sign(std::sqrt(sum),r[k][k]);
			r[k][k] += sigma;
			c[k]=sigma*r[k][k];
			d[k] = -scale*sigma;
			for (j=k+1;j<n;j++) {
				for (sum=0.0,i=k;i<n;i++) sum += r[i][k]*r[i][j];
				tau=sum/c[k];
				for (i=k;i<n;i++) r[i][j] -= tau*r[i][k];
			}
		}
	}
	d[n-1]=r[n-1][n-1];
	if (d[n-1] == 0.0) sing=true;
	for (i=0;i<n;i++) {
		for (j=0;j<n;j++) qt[i][j]=0.0;
		qt[i][i]=1.0;
	}
	for (k=0;k<n-1;k++) {
		if (c[k] != 0.0) {
			for (j=0;j<n;j++) {
				sum=0.0;
				for (i=k;i<n;i++)
					sum += r[i][k]*qt[i][j];
				sum /= c[k];
				for (i=k;i<n;i++)
					qt[i][j] -= sum*r[i][k];
			}
		}
	}
	for (i=0;i<n;i++) {
		r[i][i]=d[i];
		for (j=0;j<i;j++) r[i][j]=0.0;
	}
}

void QRDecomp::solve(const std::vector<double>& b, std::vector<double> &x)
{
	qtmult(b,x);
	rsolve(x,x);
}

void QRDecomp::qtmult(const std::vector<double>& b, std::vector<double> &x)
{
	mcIdType i,j;
	double sum;
	for (i=0;i<n;i++)
  {
		sum = 0.;
		for (j=0;j<n;j++) sum += qt[i][j]*b[j];
		x[i] = sum;
	}
}

void QRDecomp::rsolve(const std::vector<double>& b, std::vector<double> &x)
{
	mcIdType i,j;
	double sum;
	if (sing) THROW_IK_EXCEPTION("attempting solve in a singular QR");
	for (i=n-1;i>=0;i--) {
		sum=b[i];
		for (j=i+1;j<n;j++) sum -= r[i][j]*x[j];
		x[i]=sum/r[i][i];
	}
}
void QRDecomp::update(const std::vector<double>& u, const std::vector<double>& v)
{
	mcIdType i,k;
	std::vector<double> w(u);
	for (k=n-1;k>=0;k--)
		if (w[k] != 0.0) break;
	if (k < 0) k=0;
	for (i=k-1;i>=0;i--) {
		rotate(i,w[i],-w[i+1]);
		if (w[i] == 0.0)
			w[i]=std::abs(w[i+1]);
		else if (std::abs(w[i]) > std::abs(w[i+1]))
			w[i]=std::abs(w[i])*std::sqrt(1.0+sqr(w[i+1]/w[i]));
		else w[i]=std::abs(w[i+1])*std::sqrt(1.0+sqr(w[i]/w[i+1]));
	}
	for (i=0;i<n;i++) r[0][i] += w[0]*v[i];
	for (i=0;i<k;i++)
		rotate(i,r[i][i],-r[i+1][i]);
	for (i=0;i<n;i++)
		if (r[i][i] == 0.0) sing=true;
}

void QRDecomp::rotate(const mcIdType i, const double a, const double b)
{
	mcIdType j;
	double c,fact,s,w,y;
	if (a == 0.0)
  {
		c=0.0;
		s=(b >= 0.0 ? 1.0 : -1.0);
	}
  else if (std::abs(a) > std::abs(b))
  {
		fact=b/a;
		c=sign(1.0/std::sqrt(1.0+(fact*fact)),a);
		s=fact*c;
	}
  else
  {
		fact=a/b;
		s=sign(1.0/std::sqrt(1.0+(fact*fact)),b);
		c=fact*s;
	}
	for (j=i;j<n;j++)
  {
		y=r[i][j];
		w=r[i+1][j];
		r[i][j]=c*y-s*w;
		r[i+1][j]=s*y+c*w;
	}
	for (j=0;j<n;j++)
  {
		y=qt[i][j];
		w=qt[i+1][j];
		qt[i][j]=c*y-s*w;
		qt[i+1][j]=s*y+c*w;
	}
}
