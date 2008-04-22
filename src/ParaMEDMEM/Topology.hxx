#ifndef TOPOLOGY_HXX_
#define TOPOLOGY_HXX_

#include <utility>

using namespace std;
namespace ParaMEDMEM
{
  class ProcessorGroup;

class Topology
{
public:
	Topology(){}
	virtual ~Topology(){}
//	virtual std::pair<int,int> globalToLocal (const int) const =0;
//	virtual int localToGlobal (const std::pair<int,int>) const =0;
	virtual int getNbElements() const=0;
	virtual int getNbLocalElements() const =0;
	virtual const ProcessorGroup* getProcGroup()const =0;
};

}

#endif /*TOPOLOGY_HXX_*/
