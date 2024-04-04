// Copyright (C) 2007-2024  CEA, EDF
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

#include "ExplicitMapping.hxx"
#include <utility>
#include "MCIdType.hxx"
#include <vector>
#include <map>

namespace MEDCoupling
{

  ExplicitMapping::ExplicitMapping():
        _mapping(), _distant_domains(),
        _numbers(nullptr), _domains(nullptr), _comm_buffer(nullptr),
        _buffer_index(nullptr), _send_counts(nullptr)
  { }

  ExplicitMapping::~ExplicitMapping()
  {
    if (_domains!=nullptr) delete[] _domains;
    if (_numbers!=nullptr) delete[] _numbers;
    if (_comm_buffer!=nullptr) delete[] _comm_buffer;
  }

  void ExplicitMapping::pushBackElem(std::pair<int,mcIdType> idistant)
  {
    _mapping.push_back(idistant);
  }

  void  ExplicitMapping::setDistantElem(mcIdType ilocal, std::pair<int,mcIdType> idistant)
  {
    _mapping[ilocal]=idistant;
  }

  int ExplicitMapping::nbDistantDomains()
  {
    if (_distant_domains.empty())
      {
        for (const auto & iter : _mapping)
          _distant_domains.insert(iter.first);
      }
    return (int)_distant_domains.size();
  }

  std::pair <int,mcIdType> ExplicitMapping::getDistantNumbering(mcIdType ielem)const
  {
    return _mapping[ielem];
  }

  int ExplicitMapping::getDistantDomain(int i)
  {
    if (_domains==nullptr)
      computeNumbers();

    return _domains[i];
  }

  int ExplicitMapping::getNbDistantElems(int i)
  {
    if (_numbers==nullptr)
      computeNumbers();
    return _numbers[i];
  }

  mcIdType* ExplicitMapping::serialize(int idproc)
  {
    _comm_buffer=new mcIdType[_mapping.size()*2];
    std::vector<mcIdType> offsets(_distant_domains.size());
    offsets[0]=0;
    for (int i=1; i<(int)_distant_domains.size();i++)
      offsets[i]=offsets[i-1]+_numbers[i-1];

    for (auto & i : _mapping)
      {
        mcIdType const offset= offsets[i.first];
        _comm_buffer[offset*2]=idproc;
        _comm_buffer[offset*2+1]=i.second;
        offsets[i.first]++;
      }
    return _comm_buffer;
  }

  void ExplicitMapping::unserialize(int nbprocs, int* sizes,int nbtarget, int* targetrank, mcIdType* commbuffer)
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
          _buffer_index[indmap]=(int)commbuffer[indmap*2+1];
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

  void ExplicitMapping::computeNumbers()
  {
    std::map <int,int> counts;
    if (_numbers==nullptr)
      {
        _numbers=new int[nbDistantDomains()];
        _domains=new int[nbDistantDomains()];
        for (auto & i : _mapping)
          {
            if ( counts.find(i.first) == counts.end())
              counts.insert(std::make_pair(i.first,1));
            else
              (counts[i.first])++;
          }
        int counter=0;
        for (std::map<int,int>::const_iterator iter=counts.begin();
            iter!=counts.end();
            iter++)
          {
            _numbers[counter]=iter->second;
            _domains[counter]=iter->first;
            counter++;
          }
      }
  }

}
