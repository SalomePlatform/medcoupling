# Copyright (C) 2015  CEA/DEN, EDF R&D
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
# See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#

ADD_TEST(TestMEDCoupling TestMEDCoupling)
SET_TESTS_PROPERTIES(TestMEDCoupling PROPERTIES LABELS "${COMPONENT_NAME}")

ADD_TEST(TestMEDCouplingRemapper TestMEDCouplingRemapper)
SET_TESTS_PROPERTIES(TestMEDCouplingRemapper PROPERTIES LABELS "${COMPONENT_NAME}")

ADD_TEST(TestMEDCouplingExamples TestMEDCouplingExamples)
SET_TESTS_PROPERTIES(TestMEDCouplingExamples PROPERTIES LABELS "${COMPONENT_NAME}")
