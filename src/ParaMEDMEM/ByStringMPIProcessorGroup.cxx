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

#include "ByStringMPIProcessorGroup.hxx"

#include <cstddef>
#include <string>
#include "MPIProcessorGroup.hxx"
#include "mpi.h"

using namespace std;


namespace MEDCoupling
{
  /*!
   \class ByStringMPIProcessorGroup

   The ByStringMPIProcessorGroup implements a derived version of MPIProcessorGroup. 

   Groups are formed from MPI ranks with the same simCodeTag. Two trivial cases: 
    - All simCodeTag are equal, then one group is formed from all mpi ranks in the communicator
    - All simCodeTag are different, then n-MPI groups are formed.
  */


  /*! Build and return a map between different simCodeTag and the set of MPI ranks ids (based on the passed communicator) grouped by the same identifier.
    \param interface CommInterface object giving access to the MPI communication layer
    \param simCodeTag the string identifiying the tag for the group.
    \param world_comm mpi communicator
    \return Map relating unique simCodeTag and the group of MPI ranks ids belowing to that group
  */
  static std::map<std::string,std::set<int>> DefineSetIdByStringName( const CommInterface& interface, const std::string& simCodeTag, const MPI_Comm& world_comm )
  {
    int size_world;
    int rank_world;
    interface.commSize(world_comm,&size_world);
    interface.commRank(world_comm,&rank_world);

    std::map<std::string,std::set<int>> myRanksSet;

    std::vector<int> displacement(size_world, 0 );
    std::vector<int> words_size(size_world);

    int stringSize = (int) simCodeTag.size();
    interface.allGather( &stringSize, 1, MPI_INT, words_size.data(), 1, MPI_INT, world_comm );

    for (size_t rank = 1; rank < words_size.size(); rank++)
      displacement[ rank ] = words_size[ rank - 1 ] + displacement[ rank - 1 ];

    char globalnames[displacement[size_world-1]];

    interface.allGatherV( simCodeTag.c_str(), stringSize, MPI_CHAR, &globalnames, 
                          words_size.data(), displacement.data(), MPI_CHAR, world_comm );

    for (size_t rank = 0; rank < size_world; rank++)
    {
      std::string strByRank( &globalnames[displacement[ rank ]], words_size[ rank ] );
      myRanksSet[ strByRank ].insert( (int)rank );
    }
    return myRanksSet;
  }

  /*! 
   * Creates a processor group that is based on all the
   processors of MPI_COMM_WORLD .This routine must be called by all processors in MPI_COMM_WORLD.
   \param interface CommInterface object giving access to the MPI
   communication layer
  */
  ByStringMPIProcessorGroup::ByStringMPIProcessorGroup(const CommInterface& interface):
    MPIProcessorGroup(interface)
  {
  }

  /*! Creates a processor group based in the simCodeTag passed.

    \param interface CommInterface object giving access to the MPI
    communication layer
    \param simCodeTag the string identifiying the tag for the group.
    \param world_comm mpi communicator
  */
  ByStringMPIProcessorGroup::ByStringMPIProcessorGroup(const CommInterface& interface, const std::string& simCodeTag, const MPI_Comm& world_comm ):
    MPIProcessorGroup(interface, DefineSetIdByStringName( interface, simCodeTag, world_comm ), simCodeTag, world_comm )
  {
  } 

  ByStringMPIProcessorGroup::ByStringMPIProcessorGroup(const ByStringMPIProcessorGroup& other):
    MPIProcessorGroup(other)
  {
  }

  ByStringMPIProcessorGroup::~ByStringMPIProcessorGroup()
  = default;

  ByStringMPIProcessorGroup *ByStringMPIProcessorGroup::deepCopy() const
  {
    return new ByStringMPIProcessorGroup(*this);
  }

}
