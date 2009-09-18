extern "C"
{
#include "metis.h"
}

#include "RENUMBER_METISRenumbering.hxx"

void METISRenumbering::renumber(const int* graph,const int* index_graph,int nb_cell,std::vector<int>& iperm,std::vector<int>& perm)
{
  iperm.resize(nb_cell,0);
  perm.resize(nb_cell,0);
  int num_flag=1;
  int options=0;
  METIS_NodeND(&nb_cell,(int*)index_graph,(int*)graph,&num_flag,&options,&iperm[0],&perm[0]);
}
