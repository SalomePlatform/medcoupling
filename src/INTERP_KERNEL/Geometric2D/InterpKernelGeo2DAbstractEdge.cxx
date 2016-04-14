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
// Author : Anthony Geay (CEA/DEN)

#include "InterpKernelGeo2DAbstractEdge.hxx"
#include "InterpKernelGeo2DComposedEdge.hxx"
#include "InterpKernelGeo2DElementaryEdge.hxx"

using namespace INTERP_KERNEL;

IteratorOnComposedEdge::IteratorOnComposedEdge():_list_handle(0)
{
}

IteratorOnComposedEdge::IteratorOnComposedEdge(ComposedEdge *compEdges):_list_handle(compEdges->getListBehind()) 
{
  first(); 
}

void IteratorOnComposedEdge::operator=(const IteratorOnComposedEdge& other)
{
  _deep_it=other._deep_it;
  _list_handle=other._list_handle;
}

void IteratorOnComposedEdge::last()
{
  _deep_it=_list_handle->end();
  _deep_it--;
}

void IteratorOnComposedEdge::nextLoop()
{
  _deep_it++;
  if(_deep_it==_list_handle->end())
    first();
}

void IteratorOnComposedEdge::previousLoop()
{
  if(_deep_it!=_list_handle->begin())
    _deep_it--;
  else
    last();
}

bool IteratorOnComposedEdge::goToNextInOn(bool direction, int& i, int nbMax)
{
  TypeOfEdgeLocInPolygon loc=current()->getLoc();
  if(direction)
    {
      while(loc==FULL_OUT_1 && i<nbMax)
        {
          nextLoop(); i++;
          loc=current()->getLoc();
        }
      if(i==nbMax)
        return false;
      return true;
    }
  else
    {
      while(loc==FULL_OUT_1 && i<nbMax)
        {
          previousLoop(); i++;
          loc=current()->getLoc();
        }
      if(i==nbMax)
        return false;
      while(loc!=FULL_OUT_1 && i<nbMax)
        {
          previousLoop(); i++;
          loc=current()->getLoc();
        }
      nextLoop(); i--;
      return true;
    }
}

void IteratorOnComposedEdge::assignMySelfToAllElems(ComposedEdge *elems)
{
  std::list<ElementaryEdge *> *myList=elems->getListBehind();
  for(std::list<ElementaryEdge *>::iterator iter=myList->begin();iter!=myList->end();iter++)
    (*iter)->getIterator()=(*this);
}

void IteratorOnComposedEdge::insertElemEdges(ComposedEdge *elems, bool changeMySelf)
{
  std::list<ElementaryEdge *> *myListToInsert=elems->getListBehind();
  std::list<ElementaryEdge *>::iterator iter=myListToInsert->begin();
  *_deep_it=*iter;
  _deep_it++;
  iter++;
  int sizeOfMyList=myListToInsert->size();
  _list_handle->insert(_deep_it,iter,myListToInsert->end());
  if(!changeMySelf)
    {
      for(int i=0;i<sizeOfMyList;i++)
        _deep_it--;
    }
}

