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

ADD_TEST(MEDLoaderTest python MEDLoaderTest.py)
SET_TESTS_PROPERTIES(MEDLoaderTest PROPERTIES LABELS "${COMPONENT_NAME}")

ADD_TEST(MEDLoaderTest2 python MEDLoaderTest2.py)
SET_TESTS_PROPERTIES(MEDLoaderTest2 PROPERTIES LABELS "${COMPONENT_NAME}")

ADD_TEST(MEDLoaderTest3 python MEDLoaderTest3.py)
SET_TESTS_PROPERTIES(MEDLoaderTest3 PROPERTIES LABELS "${COMPONENT_NAME}")

ADD_TEST(MEDLoaderTest4 python MEDLoaderTest4.py)
SET_TESTS_PROPERTIES(MEDLoaderTest4 PROPERTIES LABELS "${COMPONENT_NAME}")

ADD_TEST(MEDLoaderExamplesTest python MEDLoaderExamplesTest.py)
SET_TESTS_PROPERTIES(MEDLoaderExamplesTest PROPERTIES LABELS "${COMPONENT_NAME}")

ADD_TEST(SauvLoaderTest python SauvLoaderTest.py)
SET_TESTS_PROPERTIES(SauvLoaderTest PROPERTIES LABELS "${COMPONENT_NAME}")

# if numpy is used
ADD_TEST(MEDLoaderCouplingTrainingSession python MEDLoaderCouplingTrainingSession.py)
SET_TESTS_PROPERTIES(MEDLoaderCouplingTrainingSession PROPERTIES LABELS "${COMPONENT_NAME}")
