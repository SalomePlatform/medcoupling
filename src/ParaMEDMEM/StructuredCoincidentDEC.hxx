//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
#ifndef STRUCTUREDCOINCIDENTDEC_HXX_
#define STRUCTUREDCOINCIDENTDEC_HXX_

#include "DEC.hxx"
#include "BlockTopology.hxx"


namespace ParaMEDMEM
{
class DEC;
class BlockTopology;
class StructuredCoincidentDEC: public DEC
{
public:
	StructuredCoincidentDEC();
	StructuredCoincidentDEC( ProcessorGroup& source, ProcessorGroup& target);
	virtual ~StructuredCoincidentDEC();
	void synchronize();
	void recvData();
	void sendData();
	void prepareSourceDE();
	void prepareTargetDE();

private :
	void synchronizeTopology();
	void broadcastTopology(BlockTopology*&, int tag);

	BlockTopology* _toposource;
	BlockTopology* _topotarget;
	int* _sendcounts;
	int* _recvcounts;
	int* _senddispls;
	int* _recvdispls;
	double* _recvbuffer;
	double* _sendbuffer;
};

}

#endif /*STRUCTUREDCOINCIDENTDEC_HXX_*/
	
