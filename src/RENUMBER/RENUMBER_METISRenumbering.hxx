#ifndef METISRENUMBERING_HXX_
#define METISRENUMBERING_HXX_

#include "RENUMBER_Renumbering.hxx"

class METISRenumbering:public Renumbering
{
public:
  virtual void renumber(const int* graph,const int* index_graph,int nb_cell,std::vector<int>& iperm,std::vector<int>& perm);
};

#endif /*METISRENUMBERING_HXX_*/
