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

namespace MEDCoupling
{
// Nodes
constexpr double EPSILON_PRIMITIVE = 0.005;
// Areas
// - Match
constexpr double TOL_MATCH_CYLINDER = 0.05;
constexpr double TOL_MATCH_SPHERE = 0.05;
// - Relative distance
constexpr double DELTA_PLANE = 0.05;
constexpr double DELTA_SPHERE = 0.05;
constexpr double DELTA_CYLINDER = 0.05;
constexpr double DELTA_CONE = 0.05;
// - Invalid Zones
constexpr double THETA_MAX_CYLINDER = 0.05;
// - Thresholds
constexpr int THRESHOLD_MIN_NB_NODES = 5;
constexpr int THRESHOLD_MAX_NB_AREAS = 500;
}  // namespace MEDCoupling
