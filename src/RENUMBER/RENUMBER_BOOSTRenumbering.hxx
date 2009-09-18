#ifndef BOOSTRENUMBERING_HXX_
#define BOOSTRENUMBERING_HXX_

#include "RENUMBER_Renumbering.hxx"

class BOOSTRenumbering:public Renumbering
{
public:
  virtual void renumber(const int* graph,const int* index_graph,int nb_cell,std::vector<int>& iperm,std::vector<int>& perm);
};

#endif /*BOOSTRENUMBERING_HXX_*/
