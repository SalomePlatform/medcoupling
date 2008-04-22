#include "AbstractEdge.hxx"
#include "ComposedEdge.hxx"
#include "ElementaryEdge.hxx"

using namespace INTERP_KERNEL;

IteratorOnComposedEdge::IteratorOnComposedEdge(ComposedEdge *cont):_container(cont) 
{
  first(); 
}

void IteratorOnComposedEdge::operator=(const IteratorOnComposedEdge& other)
{
  _container=other._container;
  for(ItOnFixdLev it=0;it<MAX_INTERSCT_DEPH;it++)
    _current[it]=other._current[it];
}

void IteratorOnComposedEdge::last()
{
  ItOnFixdLev delta=1;
  _container->getLastElementary(delta);
  _current[0]=delta;
  AbstractEdge *cur=_container;
  for(ItOnFixdLev it=1;it<=delta;it++)
    {
      _current[it]=cur->size()-1;
      cur=(AbstractEdge *)(*cur)[cur->size()-1];
    }
}

void IteratorOnComposedEdge::first()
{
  _current[0]=1;
  _current[1]=0;
}

void IteratorOnComposedEdge::next()
{
  updateNumbering();
  bool levelToIncr=false;
  do
    {
      if(getLowestDealing()->size()-1>_current[_current[0]])
        {
          _current[_current[0]]++;
          levelToIncr=true;
        }
      else
	_current[0]--;
    }
  while(_current[0]!=0 && !levelToIncr);
  if(levelToIncr)
    updateNumbering();
}

void IteratorOnComposedEdge::nextLoop()
{
  updateNumbering();
  bool levelToIncr=false;
  do
    {
      if(getLowestDealing()->size()-1>_current[_current[0]])
        {
          _current[_current[0]]++;
          levelToIncr=true;
        }
      else
	_current[0]--;
    }
  while(_current[0]!=0 && !levelToIncr);
  if(levelToIncr)
    updateNumbering();
  else
    first();
}

void IteratorOnComposedEdge::previousLoop()
{
  bool levelToIncr=false;
  do
    {
      if(_current[_current[0]]>0)
        {
          _current[_current[0]]--;
          levelToIncr=true;
	  ItOnFixdLev delta=0;
	  AbstractEdge *curLevel=getLowestDealing();
          AbstractEdge *coarseElem=(*curLevel)[_current[_current[0]]];
          if(dynamic_cast<ComposedEdge *>(coarseElem))
            {
              ((ComposedEdge *)coarseElem)->getLastElementary(++delta);
              for(ItOnFixdLev it=1;it<=delta;it++)
                {
                  _current[_current[0]+it]=coarseElem->size()-1;
                  curLevel=(AbstractEdge *)(*coarseElem)[coarseElem->size()-1];
                }
              _current[0]+=delta;
            }
        }
      else
	_current[0]--;
    }
  while(_current[0]!=0 && !levelToIncr);
  if(!levelToIncr)
    last();
}

bool IteratorOnComposedEdge::finished() const
{
  return _current[0]==0;
}

AbstractEdge *IteratorOnComposedEdge::currentDirect() const
{
  return (AbstractEdge *)(*getLowestDealing())[_current[_current[0]]];
}

ElementaryEdge* &IteratorOnComposedEdge::updateNumbering()
{
  ItOnFixdLev delta=0;
  AbstractEdge *& valToTest=(*getLowestDealing())[_current[_current[0]]];
  ElementaryEdge **ret;
  if(dynamic_cast<ElementaryEdge *>(valToTest))
    ret=(ElementaryEdge **)&valToTest;
  else
    ret=&(valToTest->getFirstElementary(++delta));
  if(delta==0)
    return *ret;
  else
    {
      for(ItOnFixdLev it=1;it<=delta;it++)
        _current[_current[0]+it]=0;
      _current[0]+=delta;
    }
  return *ret;
}

AbstractEdge *IteratorOnComposedEdge::getLowestDealing() const
{
  if(_current[0]==0)
    return 0;
  AbstractEdge *ret=(AbstractEdge *)_container;
  for(ItOnFixdLev iter=1;iter<_current[0];iter++)
    ret=(AbstractEdge *)(*ret)[_current[iter]];
  return ret;
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
