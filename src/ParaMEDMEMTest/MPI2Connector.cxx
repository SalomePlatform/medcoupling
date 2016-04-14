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

#include "MPI2Connector.hxx"

#include <iostream>
#include <cstring>

#ifndef WIN32
#include <unistd.h>
#endif

MPI2Connector::MPI2Connector()
{
  MPI_Comm_size( MPI_COMM_WORLD, &_nb_proc );
  MPI_Comm_rank( MPI_COMM_WORLD, &_num_proc );
}

MPI2Connector::~MPI2Connector()
{
}

MPI_Comm MPI2Connector::remoteMPI2Connect(const std::string& service)
{
  int i;
  char port_name[MPI_MAX_PORT_NAME];
  char port_name_clt[MPI_MAX_PORT_NAME];
  std::ostringstream msg;
  MPI_Comm icom;

  if( service.size() == 0 )
    {
      msg << "[" << _num_proc << "] You have to give a service name !";
      std::cerr << msg.str().c_str() << std::endl;
      throw std::exception();
    }

  _srv = false;

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
  if( _num_proc == 0 )
    { 
      /* rank 0 try to be a server. If service is already published, try to be a cient */
      MPI_Open_port(MPI_INFO_NULL, port_name); 
      if ( MPI_Lookup_name((char*)service.c_str(), MPI_INFO_NULL, port_name_clt) == MPI_SUCCESS )
        {
          std::cerr << "[" << _num_proc << "] I get the connection with " << service << " at " << port_name_clt << std::endl;
          MPI_Close_port( port_name );
        }
      else if ( MPI_Publish_name((char*)service.c_str(), MPI_INFO_NULL, port_name) == MPI_SUCCESS )
        {
          _srv = true;
          _port_name = port_name;
          std::cerr << "[" << _num_proc << "] service " << service << " available at " << port_name << std::endl;
        }      
      else if ( MPI_Lookup_name((char*)service.c_str(), MPI_INFO_NULL, port_name_clt) == MPI_SUCCESS )
        {
          std::cerr << "[" << _num_proc << "] I get the connection with " << service << " at " << port_name_clt << std::endl;
          MPI_Close_port( port_name );
        }
      else
        {
          msg << "[" << _num_proc << "] Error on connection with " << service << " at " << port_name_clt;
          std::cerr << msg.str().c_str() << std::endl;
          throw std::exception();
        }
    }
  else
    {
      i=0;
      /* Waiting rank 0 publish name and try to be a client */
      while ( i != TIMEOUT  ) 
        {
          sleep(1);
          if ( MPI_Lookup_name((char*)service.c_str(), MPI_INFO_NULL, port_name_clt) == MPI_SUCCESS )
            {
              std::cerr << "[" << _num_proc << "] I get the connection with " << service << " at " << port_name_clt << std::endl;
              break;
            }
          i++;
        }
      if(i==TIMEOUT)
        {
          msg << "[" << _num_proc << "] Error on connection with " << service << " at " << port_name_clt;
          std::cerr << msg.str().c_str() << std::endl;
          throw std::exception();
        }
    }
  MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
  
  /* If rank 0 is server, all processes call MPI_Comm_accept */
  /* If rank 0 is not server, all processes call MPI_Comm_connect */
  int srv = (int)_srv;
  MPI_Bcast(&srv,1,MPI_INT,0,MPI_COMM_WORLD);
  _srv = (bool)srv;
  if ( _srv )
    MPI_Comm_accept( port_name, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &icom );
  else
    MPI_Comm_connect(port_name_clt, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &icom );

  /* create global communicator: servers have low index in global communicator*/
  MPI_Intercomm_merge(icom,!_srv,&_gcom);

  /* only rank 0 can be server for unpublish name */
  if(_num_proc != 0) _srv = false;

  return _gcom;

}

void MPI2Connector::remoteMPI2Disconnect(const std::string& service)
{
  std::ostringstream msg;

  if( service.size() == 0 )
    {
      msg << "[" << _num_proc << "] You have to give a service name !";
      std::cerr << msg.str().c_str() << std::endl;
      throw std::exception();
    }

  MPI_Comm_disconnect( &_gcom ); 
  if ( _srv )
    {

      char port_name[MPI_MAX_PORT_NAME];
      strcpy(port_name,_port_name.c_str());

      MPI_Unpublish_name((char*)service.c_str(), MPI_INFO_NULL, port_name); 
      std::cerr << "[" << _num_proc << "] " << service << ": close port " << _port_name << std::endl;
      MPI_Close_port( port_name ); 
    }
  
}

