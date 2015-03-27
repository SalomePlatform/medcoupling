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

ADD_TEST(MEDCouplingBasicsTest python MEDCouplingBasicsTest.py)
SET_TESTS_PROPERTIES(MEDCouplingBasicsTest PROPERTIES LABELS "${COMPONENT_NAME}")

ADD_TEST(MEDCouplingExamplesTest python MEDCouplingExamplesTest.py)
SET_TESTS_PROPERTIES(MEDCouplingExamplesTest PROPERTIES LABELS "${COMPONENT_NAME}")

ADD_TEST(MEDCouplingRemapperTest python MEDCouplingRemapperTest.py)
SET_TESTS_PROPERTIES(MEDCouplingRemapperTest PROPERTIES LABELS "${COMPONENT_NAME}")

# if numpy is used
ADD_TEST(MEDCouplingNumPyTest python MEDCouplingNumPyTest.py)
SET_TESTS_PROPERTIES(MEDCouplingNumPyTest PROPERTIES LABELS "${COMPONENT_NAME}")

ADD_TEST(MEDCouplingPickleTest python MEDCouplingPickleTest.py)
SET_TESTS_PROPERTIES(MEDCouplingPickleTest PROPERTIES LABELS "${COMPONENT_NAME}")
