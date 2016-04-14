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

#ifndef __PROCESSORGROUP_HXX__
#define __PROCESSORGROUP_HXX__

#include "CommInterface.hxx"

#include <set>

namespace MEDCoupling
{
  /*!
   * Abstract class defining a group of processors (computation nodes) in a parallel run of a code.
   *
   * See the non-abstract child \ref MPIProcessorGroup-det "MPIProcessorGroup"
   */
  class ProcessorGroup
  {
  public:
  
    ProcessorGroup(const CommInterface& interface):_comm_interface(interface) { }
    ProcessorGroup(const CommInterface& interface, std::set<int> proc_ids):
      _comm_interface(interface),_proc_ids(proc_ids) { }
    ProcessorGroup (const ProcessorGroup& proc_group, std::set<int> proc_ids):
      _comm_interface(proc_group.getCommInterface()),_proc_ids(proc_ids) { }
    ProcessorGroup (const ProcessorGroup& other):
      _comm_interface(other.getCommInterface()),_proc_ids(other._proc_ids) { }
    ProcessorGroup (const CommInterface& interface, int start, int end);
    virtual ~ProcessorGroup() { }
    virtual ProcessorGroup *deepCopy() const = 0;
    virtual ProcessorGroup* fuse (const ProcessorGroup&) const = 0;
    virtual void intersect (ProcessorGroup&) = 0;
    bool contains(int rank) const { return _proc_ids.find(rank)!=_proc_ids.end(); }
    virtual bool containsMyRank() const = 0;
    int size() const  { return _proc_ids.size(); }
    const CommInterface& getCommInterface()const { return _comm_interface; }
    virtual int myRank() const = 0;
    virtual int translateRank(const ProcessorGroup*, int) const = 0;
    virtual ProcessorGroup* createComplementProcGroup() const = 0;
    virtual ProcessorGroup* createProcGroup() const = 0;
    virtual const std::set<int>& getProcIDs()const  { return _proc_ids; } 
  protected:
    const CommInterface _comm_interface;
    std::set<int> _proc_ids;
  };
}

#endif
