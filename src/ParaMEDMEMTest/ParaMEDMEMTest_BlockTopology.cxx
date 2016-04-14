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

#include "InterpolationUtils.hxx"
#include "CommInterface.hxx"
#include "ProcessorGroup.hxx"
#include "MPIProcessorGroup.hxx"
#include "Topology.hxx"
#include "BlockTopology.hxx"

#include <string>

// use this define to enable lines, execution of which leads to Segmentation Fault
#define ENABLE_FAULTS

// use this define to enable CPPUNIT asserts and fails, showing bugs
#define ENABLE_FORCED_FAILURES


using namespace std;
using namespace MEDCoupling;
 
/*
 * Check methods defined in BlockTopology.hxx
 *
  BlockTopology(){};
  BlockTopology(const ProcessorGroup& group, const MEDMEM::GRID& grid); 
  BlockTopology(const BlockTopology& geom_topo, const ComponentTopology& comp_topo);
  (+) BlockTopology(const ProcessorGroup& group, int nb_elem);
  virtual ~BlockTopology();
  (+) inline int getNbElements()const;
  (+) inline int getNbLocalElements() const;
  const ProcessorGroup* getProcGroup()const {return _proc_group;};
  (+) inline std::pair<int,int> globalToLocal (const int) const ;
  (+) inline int localToGlobal (const std::pair<int,int>) const;
  (+) std::vector<std::pair<int,int> > getLocalArrayMinMax() const ;
  (+) int getDimension() const {return _dimension;};
  (+) void serialize(int* & serializer, int& size) const ;
  (+) void unserialize(const int* serializer, const CommInterface& comm_interface);
  
 */
 
void ParaMEDMEMTest::testBlockTopology_constructor()
{
  //test constructor
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  CommInterface interface;
  MPIProcessorGroup group(interface);
  BlockTopology blocktopo(group,1);
  CPPUNIT_ASSERT_EQUAL(1,blocktopo.getNbLocalElements());
  CPPUNIT_ASSERT_EQUAL(size,blocktopo.getNbElements());
  CPPUNIT_ASSERT_EQUAL(1,blocktopo.getDimension());
  
  //checking access methods
  BlockTopology blocktopo2(group,2);
  std::pair<int,int> local= blocktopo2.globalToLocal(0);
  CPPUNIT_ASSERT_EQUAL(local.first,0);
  CPPUNIT_ASSERT_EQUAL(local.second,0);
  int global=blocktopo2.localToGlobal(local);
  CPPUNIT_ASSERT_EQUAL(global,0);
  
  local = blocktopo2.globalToLocal(1);
  CPPUNIT_ASSERT_EQUAL(local.first,0);
  CPPUNIT_ASSERT_EQUAL(local.second,1);
  global=blocktopo2.localToGlobal(local);
  CPPUNIT_ASSERT_EQUAL(global,1);
  
  local = blocktopo2.globalToLocal(2*size-1);
  CPPUNIT_ASSERT_EQUAL(local.first,size-1);
  CPPUNIT_ASSERT_EQUAL(local.second,1);
  global=blocktopo2.localToGlobal(local);
  CPPUNIT_ASSERT_EQUAL(global,2*size-1);

  std::vector<std::pair<int,int> > bounds = blocktopo2.getLocalArrayMinMax();
  int vecsize = bounds.size();
  CPPUNIT_ASSERT_EQUAL(1,vecsize);
  CPPUNIT_ASSERT_EQUAL(2*rank, (bounds[0]).first);
  CPPUNIT_ASSERT_EQUAL(2*rank+2, (bounds[0]).second);
 }
 
void ParaMEDMEMTest::testBlockTopology_serialize()
{

  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  CommInterface interface;
  MPIProcessorGroup group(interface);
  BlockTopology blocktopo(group,3);

//testing the serialization process that is used to transfer a
//block topology via a MPI_Send/Recv comm  
  BlockTopology blocktopo_recv;
  int* serializer;
  int sersize;
  blocktopo.serialize(serializer,sersize);
  blocktopo_recv.unserialize(serializer,interface);
  CPPUNIT_ASSERT_EQUAL(blocktopo.getNbElements(),blocktopo_recv.getNbElements());
  delete [] serializer;
}
