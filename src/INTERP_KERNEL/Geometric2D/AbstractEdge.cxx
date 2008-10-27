#include "AbstractEdge.hxx"
#include "ComposedEdge.hxx"
#include "ElementaryEdge.hxx"

using namespace INTERP_KERNEL;

IteratorOnComposedEdge::IteratorOnComposedEdge():_listHandle(0)
{
}

IteratorOnComposedEdge::IteratorOnComposedEdge(ComposedEdge *compEdges):_listHandle(compEdges->getListBehind()) 
{
  first(); 
}

void IteratorOnComposedEdge::operator=(const IteratorOnComposedEdge& other)
{
  _deepIt=other._deepIt;
  _listHandle=other._listHandle;
}

void IteratorOnComposedEdge::last()
{
  _deepIt=_listHandle->end();
  _deepIt--;
}

void IteratorOnComposedEdge::nextLoop()
{
  _deepIt++;
  if(_deepIt==_listHandle->end())
    first();
}

void IteratorOnComposedEdge::previousLoop()
{
  if(_deepIt!=_listHandle->begin())
    _deepIt--;
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
  *_deepIt=*iter;
  _deepIt++;
  iter++;
  int sizeOfMyList=myListToInsert->size();
  _listHandle->insert(_deepIt,iter,myListToInsert->end());
  if(!changeMySelf)
    {
      for(int i=0;i<sizeOfMyList;i++)
        _deepIt--;
    }
}
