// Copyright (C) 2024-2026  CEA, EDF
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

#pragma once

#include <vector>
#include <array>
#include "ShapeRecognDefines.hxx"
namespace MEDCoupling
{
class SHAPE_RECOGNITION_EXPORT MathOps
{
   public:
    static std::vector<double> lstsq(std::vector<double> &a, const std::vector<double> &b);
    static std::vector<double> lstsq(
        std::vector<double> &a, const std::vector<double> &b, int m, int n = 3, int nrhs = 1
    );
    static std::vector<double> lstsqRow(std::vector<double> &a, const std::vector<double> &b);
    static std::array<double, 3> cross(const std::array<double, 3> &a, const std::array<double, 3> &b);
    static std::array<double, 3> normalize(const std::array<double, 3> &a);
    static double computeNorm(const std::array<double, 3> &a);
    static std::vector<double> computeNorm(const std::vector<double> &a);
    static double dot(const std::array<double, 3> &a, const std::array<double, 3> &b);
    static std::vector<double> dot(const std::vector<double> &a, const std::array<double, 3> &b);
    static double mean(const std::vector<double> &values);
    static std::array<double, 3> meanCoordinates(const std::vector<double> &coordinates);
    static std::array<double, 9> computeCov(const std::vector<double> &coordinates);
    static std::array<double, 9> computePCA(const std::vector<double> &coordinates);
    static std::array<double, 3> computePCAFirstAxis(const std::vector<double> &coordinates);
    static std::array<double, 3> computePCAThirdAxis(const std::vector<double> &coordinates);
    static double computeQuantile(const std::vector<double> &values, double q);
    static double computeAngle(const std::array<double, 3> &direction, std::array<double, 3> axis);
    static std::vector<double> computeAngles(const std::vector<double> &directions, std::array<double, 3> axis);
    static double computeOrientedAngle(
        const std::array<double, 3> &normal, const std::array<double, 3> &vector1, const std::array<double, 3> &vector2
    );
    static double computeVariance(std::vector<double> values);
    static std::array<double, 6> computeBaseFromNormal(std::array<double, 3> normal);
};
}  // namespace MEDCoupling
