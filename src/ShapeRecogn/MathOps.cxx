// Copyright (C) 2024  CEA, EDF
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

#include "MathOps.hxx"
#include "MCIdType.hxx"

#include <algorithm>
#if defined(_MSC_VER)
  // see https://github.com/xianyi/OpenBLAS/issues/3661
  #define _CRT_USE_C_COMPLEX_H
  #include <complex.h>
  #define LAPACK_COMPLEX_CUSTOM
  #define lapack_complex_float _Fcomplex
  #define lapack_complex_double _Dcomplex
  #include <openblas/lapacke.h>
  #include <openblas/cblas.h>
#else
  #include <lapacke.h>
  #include <cblas.h>
#endif
#include <iostream>
#include <cfloat>
#include <cmath>
#include <cstdlib>

using namespace MEDCoupling;

std::vector<double>
MathOps::lstsq(std::vector<double> &a, const std::vector<double> &b)
{
    int m = (int)b.size();
    int n = 3;
    int nrhs = 1;
    return lstsq(a, b, m, n, nrhs);
}

std::vector<double>
MathOps::lstsq(std::vector<double> &a, const std::vector<double> &b, int m, int n, int nrhs)
{
    int ldb = std::max<int>(m, n);
    int lds = std::min<int>(m, n);
    std::vector<double> x(ldb, 0.0);
    for (size_t i = 0; i < b.size(); ++i) x[i] = b[i];
    std::vector<double> s(lds, 0.0);
    double rcond = DBL_EPSILON * (double)ldb;  // same value as numpy.linalg.lstsq
    int rank = 0;
    int info = LAPACKE_dgelsd(LAPACK_COL_MAJOR, m, n, nrhs, a.data(), m, x.data(), ldb, s.data(), rcond, &rank);
    return x;
}

std::vector<double>
MathOps::lstsqRow(std::vector<double> &a, const std::vector<double> &b)
{
    auto m = b.size();
    std::size_t n = 3;
    std::size_t nrhs = 1;
    std::size_t ldb = std::max<std::size_t>(m, n);
    std::size_t lds = std::min<std::size_t>(m, n);
    std::vector<double> x(ldb, 0.0);
    for (size_t i = 0; i < b.size(); ++i) x[i] = b[i];
    std::vector<double> s(lds, 0.0);
    double rcond = DBL_EPSILON * (double)ldb;  // same value as numpy.linalg.lstsq
    int rank = 0;
    int info = LAPACKE_dgelsd(
        LAPACK_ROW_MAJOR,
        FromIdType<int>(m),
        FromIdType<int>(n),
        FromIdType<int>(nrhs),
        a.data(),
        FromIdType<int>(n),
        x.data(),
        FromIdType<int>(nrhs),
        s.data(),
        rcond,
        &rank
    );
    return x;
}
std::array<double, 3>
MathOps::cross(const std::array<double, 3> &a, const std::array<double, 3> &b)
{
    std::array<double, 3> c{0.0, 0.0, 0.0};
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
    return c;
}

std::array<double, 3>
MathOps::normalize(const std::array<double, 3> &a)
{
    std::array<double, 3> an{0.0, 0.0, 0.0};
    double n = computeNorm(a);
    for (size_t i = 0; i < 3; i++) an[i] = a[i] / n;
    return an;
}

double
MathOps::computeNorm(const std::array<double, 3> &a)
{
    double n = 0;
    for (size_t i = 0; i < 3; i++) n += a[i] * a[i];
    return sqrt(n);
}

std::vector<double>
MathOps::computeNorm(const std::vector<double> &a)
{
    size_t s = a.size() / 3;
    std::vector<double> n(s, 0.0);
    for (size_t i = 0; i < s; i++)
    {
        for (size_t j = 0; j < 3; j++) n[i] += a[3 * i + j] * a[3 * i + j];
        n[i] = sqrt(n[i]);
    }
    return n;
}

double
MathOps::dot(const std::array<double, 3> &a, const std::array<double, 3> &b)
{
    double d = 0.0;
    for (size_t i = 0; i < 3; i++) d += a[i] * b[i];
    return d;
}

std::vector<double>
MathOps::dot(const std::vector<double> &a, const std::array<double, 3> &b)
{
    size_t nbNodes = a.size() / 3;
    std::vector<double> d(nbNodes, 0.0);
    cblas_dgemv(
        CBLAS_LAYOUT::CblasRowMajor,
        CBLAS_TRANSPOSE::CblasNoTrans,
        (int)nbNodes,
        3,
        1.0,
        a.data(),
        3,
        b.data(),
        1,
        0.0,
        d.data(),
        1
    );
    return d;
}

double
MathOps::mean(const std::vector<double> &values)
{
    double mean = 0.0;
    for (double value : values) mean += value;
    return mean / double(values.size());
}

std::array<double, 3>
MathOps::meanCoordinates(const std::vector<double> &coordinates)
{
    std::array<double, 3> coordsMean{0.0, 0.0, 0.0};
    size_t nbNodes = coordinates.size() / 3;
    for (size_t nodeId = 0; nodeId < nbNodes; ++nodeId)
    {
        coordsMean[0] += coordinates[3 * nodeId];
        coordsMean[1] += coordinates[3 * nodeId + 1];
        coordsMean[2] += coordinates[3 * nodeId + 2];
    }
    coordsMean[0] /= double(nbNodes);
    coordsMean[1] /= double(nbNodes);
    coordsMean[2] /= double(nbNodes);
    return coordsMean;
}

std::array<double, 9>
MathOps::computeCov(const std::vector<double> &coordinates)
{
    std::array<double, 9> covMatrix;
    covMatrix.fill(0);
    size_t nbNodes = coordinates.size() / 3;
    // Center the coordinates
    std::vector<double> coordsCentered(coordinates);
    std::array<double, 3> coordsMean = meanCoordinates(coordinates);
    for (size_t nodeId = 0; nodeId < nbNodes; ++nodeId)
    {
        coordsCentered[3 * nodeId] -= coordsMean[0];
        coordsCentered[3 * nodeId + 1] -= coordsMean[1];
        coordsCentered[3 * nodeId + 2] -= coordsMean[2];
    }
    cblas_dgemm(
        CBLAS_LAYOUT::CblasColMajor,
        CBLAS_TRANSPOSE::CblasNoTrans,
        CBLAS_TRANSPOSE::CblasTrans,
        3,
        3,
        (int)nbNodes,
        1.0,
        coordsCentered.data(),
        3,
        coordsCentered.data(),
        3,
        0.0,
        covMatrix.data(),
        3
    );
    for (size_t i = 0; i < 9; ++i)
    {
        covMatrix[i] /= (double)nbNodes - 1;
    }
    return covMatrix;
}

std::array<double, 9>
MathOps::computePCA(const std::vector<double> &coordinates)
{
    std::array<double, 9> covMatrix = computeCov(coordinates);
    std::array<double, 3> eigenValues{0.0, 0.0, 0.0};
    LAPACKE_dsyevd(LAPACK_COL_MAJOR, 'V', 'U', 3, covMatrix.data(), 3, eigenValues.data());
    return covMatrix;
}

std::array<double, 3>
MathOps::computePCAFirstAxis(const std::vector<double> &coordinates)
{
    std::array<double, 9> pca = computePCA(coordinates);
    // The eignvalues are in ascending order so the first axis correspond to the last vector
    return {pca[6], pca[7], pca[8]};
}

std::array<double, 3>
MathOps::computePCAThirdAxis(const std::vector<double> &coordinates)
{
    std::array<double, 9> pca = computePCA(coordinates);
    // The eignvalues are in ascending order so the third axis correspond to the first vector
    return {pca[0], pca[1], pca[2]};
}

double
MathOps::computeQuantile(const std::vector<double> &values, double q)
{
    std::vector<double> sortedValues(values);
    std::sort(sortedValues.begin(), sortedValues.end());
    double pos = q * double(sortedValues.size() - 1);
    size_t index = static_cast<size_t>(pos);
    if (std::abs(pos - double(index)) < 1e-12)
        return sortedValues[index];
    else
    {
        double frac = pos - double(index);
        return sortedValues[index] * (1 - frac) + sortedValues[index + 1] * frac;
    }
}

double
MathOps::computeAngle(const std::array<double, 3> &direction, std::array<double, 3> axis)
{
    double angle = dot(direction, axis);
    double normAxis = computeNorm(axis);
    double normDirection = computeNorm(direction);
    angle /= normDirection * normAxis;
    if (fabs(angle) >= 1)
        return 0;
    else
        return acos(angle);
}

std::vector<double>
MathOps::computeAngles(const std::vector<double> &directions, std::array<double, 3> axis)
{
    size_t nbDirections = directions.size() / 3;
    std::vector<double> angles(nbDirections, 0.0);
    cblas_dgemv(
        CBLAS_LAYOUT::CblasRowMajor,
        CBLAS_TRANSPOSE::CblasNoTrans,
        int(nbDirections),
        3,
        1.0,
        directions.data(),
        3,
        axis.data(),
        1,
        0.0,
        angles.data(),
        1
    );
    double normAxis = computeNorm(axis);
    std::vector<double> normDirections = computeNorm(directions);
    for (size_t i = 0; i < nbDirections; ++i)
    {
        angles[i] /= normAxis * normDirections[i];
        if (fabs(angles[i]) >= 1.0)
            angles[i] = 0.0;
        angles[i] = acos(angles[i]);
    }
    return angles;
}

double
MathOps::computeOrientedAngle(
    const std::array<double, 3> &normal, const std::array<double, 3> &vector1, const std::array<double, 3> &vector2
)
{
    double angle = computeAngle(vector1, vector2);
    if (dot(cross(vector1, vector2), normal) >= 0.0)
        return angle;
    else
        return -angle;
}

double
MathOps::computeVariance(std::vector<double> values)
{
    size_t n = values.size();
    double m = mean(values);
    double d2 = 0;
    for (size_t i = 0; i < n; ++i) d2 += pow(values[i] - m, 2);
    return d2 / (double)n;
}

std::array<double, 6>
MathOps::computeBaseFromNormal(std::array<double, 3> normal)
{
    std::array<double, 3> n_normal = normalize(normal);
    std::array<double, 3> s;
    std::array<double, 1> u;
    std::array<double, 9> v;
    LAPACKE_dgesdd(LAPACK_COL_MAJOR, 'A', 1, 3, n_normal.data(), 1, s.data(), u.data(), 1, v.data(), 3);
    std::array<double, 3> v1 = {v[3], v[4], v[5]};
    std::array<double, 3> v2 = {v[6], v[7], v[8]};
    std::array<double, 3> u1 = normalize(v1);
    std::array<double, 3> u2;
    u2.fill(0.0);
    double innerProd = dot(u1, v2);
    for (size_t i = 0; i < 3; ++i) u2[i] = v2[i] - innerProd * u1[i];
    u2 = normalize(u2);
    double sign = dot(cross(u1, u2), normal);
    if (sign < 0)
    {
        for (size_t i = 0; i < 3; ++i) u1[i] *= -1;
    }
    return {u1[0], u1[1], u1[2], u2[0], u2[1], u2[2]};
}
