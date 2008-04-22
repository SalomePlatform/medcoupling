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
	
