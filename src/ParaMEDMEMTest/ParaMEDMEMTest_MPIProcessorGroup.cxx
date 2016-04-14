// Copyright (C) 2007-2016  CEA/DEN, EDF R&D
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
#include "MPIProcessorGroup.hxx"
#include "InterpolationUtils.hxx"

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
 (+) MPIProcessorGroup(const CommInterface& interface);
 (+) MPIProcessorGroup(const CommInterface& interface, set<int> proc_ids);
 (u) MPIProcessorGroup (const ProcessorGroup& proc_group, set<int> proc_ids);
 (+) MPIProcessorGroup(const CommInterface& interface,int pstart, int pend);
 (+) virtual ~MPIProcessorGroup();
 (+) virtual ProcessorGroup* fuse (const ProcessorGroup&) const;
 (u) void intersect (ProcessorGroup&){};
 (+) int myRank() const {int rank; MPI_Comm_rank(_comm,&rank); return rank;}
 (+) bool containsMyRank() const { int rank; MPI_Group_rank(_group, &rank); return (rank!=MPI_UNDEFINED);}
 (+) int translateRank(const ProcessorGroup* group, int rank) const;
 (+) const MPI_Comm* getComm() const {return &_comm;}
 (+) ProcessorGroup* createComplementProcGroup() const;
 (o) ProcessorGroup* createProcGroup() const;
   
*/
 
void ParaMEDMEMTest::testMPIProcessorGroup_constructor()
{
  CommInterface comm_interface;
  MPIProcessorGroup* group = new MPIProcessorGroup(comm_interface);;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  CPPUNIT_ASSERT_EQUAL(size,group->size());
  int size2;
  const MPI_Comm* communicator=group->getComm();
  MPI_Comm_size(*communicator, &size2);
  CPPUNIT_ASSERT_EQUAL(size,size2);
  delete group;

  set <int> procs;

  procs.insert(0);
  procs.insert(1);
  if (size==1)
    CPPUNIT_ASSERT_THROW(group=new MPIProcessorGroup(comm_interface,procs),INTERP_KERNEL::Exception);
  else
    {
      CPPUNIT_ASSERT_NO_THROW(  group=new MPIProcessorGroup(comm_interface,procs));
      CPPUNIT_ASSERT_EQUAL (group->size(),2);
      delete group;
    }

  //throws because plast<pfirst
  CPPUNIT_ASSERT_THROW(group=new MPIProcessorGroup(comm_interface,1,0),INTERP_KERNEL::Exception);
  //throws because plast is beyond size-1
  CPPUNIT_ASSERT_THROW(group=new MPIProcessorGroup(comm_interface,0,size),INTERP_KERNEL::Exception);
  if (size>1)
    {
      group=new MPIProcessorGroup(comm_interface,0,size-2);
      CPPUNIT_ASSERT_EQUAL(group->size(),size-1);
      delete group;
    }

}
 
void ParaMEDMEMTest::testMPIProcessorGroup_boolean()
{
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  CommInterface comm_interface;
  MPIProcessorGroup group(comm_interface,0,0);
  MPIProcessorGroup group2(comm_interface,size-1,size-1);
  ProcessorGroup* group_fuse=group.fuse(group2);
  int group_fuse_size=(size==1)?1:2;
  CPPUNIT_ASSERT_EQUAL(group_fuse_size,group_fuse->size());
 
  ProcessorGroup* group_complement=((MPIProcessorGroup*)group_fuse)->createComplementProcGroup();
  CPPUNIT_ASSERT_EQUAL(group_complement->size(),size-group_fuse_size);
  
  delete group_fuse;
  delete group_complement;

  //intersect not implemented yet
  //   if (size>1)
  //   {
  //     MPIProcessorGroup group3(comm_interface,0,size-2);
  //     MPIProcessorGroup group4(comm_interface,1,size-1);
  //     group3.intersect(group4);
  //     CPPUNIT_ASSERT_EQUAL(group3.size(),size-2);
  //   }
}

void ParaMEDMEMTest::testMPIProcessorGroup_rank()
{
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  CommInterface comm_interface;
  MPIProcessorGroup group(comm_interface,0,0);
  MPIProcessorGroup group2(comm_interface,size-1,size-1);
  ProcessorGroup* group_fuse=group2.fuse(group);
  
  if (group.containsMyRank())
    CPPUNIT_ASSERT_EQUAL (group.myRank(), rank);

  if (group2.containsMyRank())
    {
      int trank=group_fuse->translateRank(&group2,0);
      if (size==1)
        CPPUNIT_ASSERT_EQUAL(trank,0);
      else  
        CPPUNIT_ASSERT_EQUAL(trank,1);
    }
  delete group_fuse;
}
