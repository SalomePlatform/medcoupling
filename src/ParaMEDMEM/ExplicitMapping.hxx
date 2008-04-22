#ifndef EXPLICIT_MAPPING_HXX_
#define EXPLICIT_MAPPING_HXX_

#include <vector>
#include <map>
#include <set>

namespace ParaMEDMEM
{
  class ExplicitMapping
  {
  public:

    ExplicitMapping():_numbers(0), _domains(0), _commbuffer(0) {};

    ~ExplicitMapping()
    {
      if (_domains!=0) delete[] _domains;
      if (_numbers!=0) delete[] _numbers;
      if (_commbuffer!=0) delete[] _commbuffer;
    }
    void pushBackElem(pair<int,int> idistant) {
      _mapping.push_back(idistant);
    }
    void  setDistantElem(int ilocal, pair<int,int> idistant)
    {
      _mapping[ilocal]=idistant;
    }

    int nbDistantDomains()
    {
      if (_distant_domains.empty())
	{
	 for (vector <pair<int,int> >::const_iterator iter= _mapping.begin();
	       iter!=_mapping.end();
	       iter++)
	    _distant_domains.insert(iter->first);
	}
      return _distant_domains.size();
    }
    
    pair <int,int> getDistantNumbering(int ielem)const
    {
      return _mapping[ielem];
    }
    
   
    int getDistantDomain(int i)
    {
      if (_domains==0)
	computeNumbers();

      return _domains[i];
    }

    int getNbDistantElems(int i)
    {
      if (_numbers==0)
	computeNumbers();
      return _numbers[i];	  
    }

    int* serialize(int idproc)
    {
      _commbuffer=new int[_mapping.size()*2];
      vector<int> offsets(_distant_domains.size());
      offsets[0]=0;
      for (int i=1; i<_distant_domains.size();i++)
	offsets[i]=offsets[i-1]+_numbers[i-1];
      
      for (int i=0; i< _mapping.size(); i++)
	{
	  int offset= offsets[_mapping[i].first];
	  _commbuffer[offset*2]=idproc;
	  _commbuffer[offset*2+1]=_mapping[i].second;
	  offsets[_mapping[i].first]++;
	}
      return _commbuffer;
    }

    void unserialize(int nbprocs, int* sizes,int nbtarget, int* targetrank, int* commbuffer)
    {
      int total_size=0;
      for (int i=0; i< nbprocs; i++)
	total_size+=sizes[i];
      
      _mapping.resize(total_size);
      _buffer_index=new int[total_size];
      int indmap=0;
      for (int i=0; i<nbprocs; i++)
	for (int ielem=0; ielem<sizes[i]; ielem++)
	{
	  _mapping[indmap].first=i;
	  _mapping[indmap].second=commbuffer[indmap*2+1];
	  _buffer_index[indmap]=commbuffer[indmap*2+1];
	  indmap++;
	}	
      _numbers=new int [nbtarget];
      _domains=new int [nbtarget];
      
      int index=0;      
      for (int i=0; i<nbtarget; i++)
	{
	  if (sizes[targetrank[i]]>0)
	    {
	      _numbers[index]=sizes[targetrank[i]];
	      _domains[index]=i;
	      index++;
	    }
	}
      _send_counts=new int[nbprocs];
      for (int i=0; i<nbprocs; i++)
	_send_counts[i]=sizes[i];
    }

    int* getBufferIndex() const { return _buffer_index;}
    int* getCounts() const{return _send_counts;}
  private:
    vector <pair<int,int> > _mapping;
    set<int> _distant_domains;
    int* _numbers;
    int* _domains;
    int* _commbuffer;
    int* _buffer_index;
    int* _send_counts;

    void computeNumbers()
    {
      map <int,int> counts;
      if (_numbers==0)
	{
	  _numbers=new int[nbDistantDomains()];
	  _domains=new int[nbDistantDomains()];
	  for (int i=0; i< _mapping.size(); i++)
	    {
	      if ( counts.find(_mapping[i].first) == counts.end())
		counts.insert(make_pair(_mapping[i].first,1));
	      else
		(counts[_mapping[i].first])++;
	    }
	  int counter=0;
	  for (map<int,int>::const_iterator iter=counts.begin(); 
	       iter!=counts.end(); 
	       iter++)
	    {
	      _numbers[counter]=iter->second;
	      _domains[counter]=iter->first;
	      counter++;
	    }
	}
    }
  };

}
#endif
