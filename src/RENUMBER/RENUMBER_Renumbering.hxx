#ifndef RENUMBERING_HXX_
#define RENUMBERING_HXX_
#include <vector>

class Renumbering
{
public:
  virtual void renumber(const int* graphe,const int* index_graphe,int nb_cell,std::vector<int>& iperm,std::vector<int>& perm)=0;
}; 

#endif /*RENUMBERING_HXX_*/
