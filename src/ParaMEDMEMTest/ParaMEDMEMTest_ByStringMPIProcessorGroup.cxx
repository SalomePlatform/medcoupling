// Copyright (C) 2007-2023  CEA, EDF
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

#include "ParaMEDMEMTest.hxx"
#include <cppunit/TestAssert.h>
#include "CommInterface.hxx"
#include "ProcessorGroup.hxx"
#include "ByStringMPIProcessorGroup.hxx"

#include <string>

// use this define to enable lines, execution of which leads to Segmentation Fault
#define ENABLE_FAULTS

// use this define to enable CPPUNIT asserts and fails, showing bugs
#define ENABLE_FORCED_FAILURES

using namespace std;
using namespace MEDCoupling;

/*
 * Check methods defined in MPPIProcessorGroup.hxx
 *
 (+) ByStringMPIProcessorGroup(const CommInterface& interface);
 (+) ByStringMPIProcessorGroup(const CommInterface& interface, std::string& codeTag, const MPI_Comm& world_comm );
 (+) ByStringMPIProcessorGroup(ByStringMPIProcessorGroup& other );
*/

void
ParaMEDMEMTest::testByStringMPIProcessorGroup_constructor()
{
    CommInterface comm_interface;
    ByStringMPIProcessorGroup *group = new ByStringMPIProcessorGroup(comm_interface);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    CPPUNIT_ASSERT_EQUAL(size, group->size());
    int size2;
    const MPI_Comm *communicator = group->getComm();
    MPI_Comm_size(*communicator, &size2);
    CPPUNIT_ASSERT_EQUAL(size, size2);
    delete group;
}

void
ParaMEDMEMTest::testByStringMPIProcessorGroup_stringconstructor()
{
    int size, rankId;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankId);

    if (size != 3)
        return;

    std::string myTag;
    if (rankId == 0 || rankId == 2)
        myTag = "group0";
    else
        myTag = "gr1";

    CommInterface comm_interface;
    ByStringMPIProcessorGroup *group = new ByStringMPIProcessorGroup(comm_interface, myTag, MPI_COMM_WORLD);
    ByStringMPIProcessorGroup *copygroup = new ByStringMPIProcessorGroup(*group);
    CPPUNIT_ASSERT(group);
    CPPUNIT_ASSERT(copygroup);

    std::set<int> ranksInGroup = group->getProcIDs();
    std::set<int> ranksInCopiedGroup = group->getProcIDs();
    if (rankId == 0 || rankId == 2)
    {
        CPPUNIT_ASSERT_EQUAL((int)ranksInGroup.size(), 2);
        CPPUNIT_ASSERT_EQUAL((int)ranksInCopiedGroup.size(), 2);
    }
    else
    {
        CPPUNIT_ASSERT_EQUAL((int)ranksInGroup.size(), 1);
        CPPUNIT_ASSERT_EQUAL((int)ranksInCopiedGroup.size(), 1);
    }
    CPPUNIT_ASSERT(group->contains(rankId));
    CPPUNIT_ASSERT(copygroup->contains(rankId));
    delete group;
    delete copygroup;
}
