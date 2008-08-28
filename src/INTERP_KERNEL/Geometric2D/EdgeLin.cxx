#include "EdgeLin.hxx"
#include "Node.hxx"
#include "InterpolationUtils.hxx"

using namespace std;
using namespace INTERP_KERNEL;

namespace INTERP_KERNEL
{
  extern const unsigned MAX_SIZE_OF_LINE_XFIG_FILE=1024;
}

SegSegIntersector::SegSegIntersector(const EdgeLin& e1, const EdgeLin& e2):SameTypeIntersector(e1,e2)
{
  _matrix[0]=(*(e2.getStartNode()))[0]-(*(e2.getEndNode()))[0];
  _matrix[1]=(*(e1.getEndNode()))[0]-(*(e1.getStartNode()))[0];
  _matrix[2]=(*(e2.getStartNode()))[1]-(*(e2.getEndNode()))[1];
  _matrix[3]=(*(e1.getEndNode()))[1]-(*(e1.getStartNode()))[1];
  _col[0]=_matrix[3]*(*(e1.getStartNode()))[0]-_matrix[1]*(*(e1.getStartNode()))[1];
  _col[1]=-_matrix[2]*(*(e2.getStartNode()))[0]+_matrix[0]*(*(e2.getStartNode()))[1];
  //Little trick to avoid problems if 'e1' and 'e2' are colinears and along Ox or Oy axes.
  if(fabs(_matrix[3])>fabs(_matrix[1]))
    _ind=0;
  else
    _ind=1;
}

/*!
 * Must be called when 'this' and 'other' have been detected to be at least colinear. Typically they are overlapped.
 * Must be called after call of areOverlappedOrOnlyColinears.
 */
bool SegSegIntersector::haveTheySameDirection() const
{
  return (_matrix[_ind?1:0]>0. && _matrix[_ind?3:2]>0.) || (_matrix[_ind?1:0]<0. && _matrix[_ind?3:2]<0.);
}

/*!
 * Precondition start and end must be so that there predecessor was in the same direction than 'e1'
 */
void SegSegIntersector::getPlacements(Node *start, Node *end, TypeOfLocInEdge& whereStart, TypeOfLocInEdge& whereEnd, MergePoints& commonNode) const
{
  getCurveAbscisse(start,whereStart,commonNode);
  getCurveAbscisse(end,whereEnd,commonNode);
}

void SegSegIntersector::getCurveAbscisse(Node *node, TypeOfLocInEdge& where, MergePoints& commonNode) const
{
  bool obvious;
  obviousCaseForCurvAbscisse(node,where,commonNode,obvious);
  if(obvious)
    return ;
  int index = !_ind;
  double ret=((*node)[index]-(*_e1.getStartNode())[index])/((*_e1.getEndNode())[index]-(*_e1.getStartNode())[index]);
  if(ret>0. && ret <1.)
    where=INSIDE;
  else if(ret<0.)
    where=OUT_BEFORE;
  else
    where=OUT_AFTER;
}

/*!
 * areColinears method should be called before with a returned colinearity equal to false to avoid bad news.
 */
std::list< IntersectElement > SegSegIntersector::getIntersectionsCharacteristicVal() const
{
  std::list< IntersectElement > ret;
  double x=_matrix[0]*_col[0]+_matrix[1]*_col[1];
  double y=_matrix[2]*_col[0]+_matrix[3]*_col[1];
  //Only one intersect point possible
  Node *node=new Node(x,y);
  node->declareOn();
  bool i_1S=_e1.getStartNode()->isEqual(*node);
  bool i_1E=_e1.getEndNode()->isEqual(*node);
  bool i_2S=_e2.getStartNode()->isEqual(*node);
  bool i_2E=_e2.getEndNode()->isEqual(*node);
  ret.push_back(IntersectElement(_e1.getCharactValue(*node),
                                 _e2.getCharactValue(*node),
                                 i_1S,i_1E,i_2S,i_2E,node,_e1,_e2,keepOrder()));
  return ret;
}

/*!
 * retrieves if segs are colinears.
 * WARNING !!! Contrary to areOverlappedOrOnlyColinears method, this method use an
 * another precision to detect colinearity !
 */
bool SegSegIntersector::areColinears() const
{
  double determinant=_matrix[0]*_matrix[3]-_matrix[1]*_matrix[2];
  return fabs(determinant)<QUADRATIC_PLANAR::_arcDetectionPrecision;
}

/*!
 * Should be called \b once ! non const method.
 * \param whereToFind specifies the box where final seek should be done. Essentially it is used for caracteristic reason.
 * \param colinearity returns if regarding QUADRATIC_PLANAR::_precision ; e1 and e2 are colinears
 *                    If true 'this' is modified ! So this method be called once above all if true is returned for this parameter.
 * \param areOverlapped if colinearity if true, this parameter looks if e1 and e2 are overlapped.
 */
void SegSegIntersector::areOverlappedOrOnlyColinears(const Bounds *whereToFind, bool& colinearity, bool& areOverlapped)
{
  double determinant=_matrix[0]*_matrix[3]-_matrix[1]*_matrix[2];
  if(fabs(determinant)>2.*QUADRATIC_PLANAR::_precision)//2*_precision due to max of offset on _start and _end
    {
      colinearity=false; areOverlapped=false;
      _matrix[0]/=determinant; _matrix[1]/=determinant; _matrix[2]/=determinant; _matrix[3]/=determinant;
    }
  else
    {
      colinearity=true;
      //retrieving initial matrix
      double tmp=_matrix[0]; _matrix[0]=_matrix[3]; _matrix[3]=tmp;
      _matrix[1]=-_matrix[1]; _matrix[2]=-_matrix[2];
      //
      double deno=sqrt(_matrix[0]*_matrix[0]+_matrix[1]*_matrix[1]);
      double x=(*(_e1.getStartNode()))[0]-(*(_e2.getStartNode()))[0];
      double y=(*(_e1.getStartNode()))[1]-(*(_e2.getStartNode()))[1];
      areOverlapped=fabs((_matrix[1]*y+_matrix[0]*x)/deno)<QUADRATIC_PLANAR::_precision;
    }
}

EdgeLin::EdgeLin(std::istream& lineInXfig)
{
  char currentLine[MAX_SIZE_OF_LINE_XFIG_FILE];
  lineInXfig.getline(currentLine,MAX_SIZE_OF_LINE_XFIG_FILE);
  _start=new Node(lineInXfig);
  _end=new Node(lineInXfig);
  updateBounds();
}

EdgeLin::EdgeLin(Node *start, Node *end, bool direction):Edge(start,end,direction)
{
  updateBounds();
}

EdgeLin::EdgeLin(double sX, double sY, double eX, double eY):Edge(sX,sY,eX,eY)
{
  updateBounds();
}

EdgeLin::~EdgeLin()
{
}

/*!
 * Characteristic for edges is relative position btw 0.;1.
 */
bool EdgeLin::isIn(double characterVal) const
{
  return characterVal>0. && characterVal<1.;
}

Node *EdgeLin::buildRepresentantOfMySelf() const
{
  return new Node(((*(_start))[0]+(*(_end))[0])/2.,((*(_start))[1]+(*(_end))[1])/2.);
}

double EdgeLin::getCharactValue(const Node& node) const
{
  double car1_1x=node[0]-(*(_start))[0]; double car1_2x=(*(_end))[0]-(*(_start))[0];
  double car1_1y=node[1]-(*(_start))[1]; double car1_2y=(*(_end))[1]-(*(_start))[1];
  return (car1_1x*car1_2x+car1_1y*car1_2y)/(car1_2x*car1_2x+car1_2y*car1_2y);
}

void EdgeLin::dumpInXfigFile(std::ostream& stream, bool direction, int resolution, const Bounds& box) const
{
  stream << "2 1 0 1 ";
  fillXfigStreamForLoc(stream);
  stream << " 7 50 -1 -1 0.000 0 0 -1 0 0 2" << endl;
  direction?_start->dumpInXfigFile(stream,resolution,box):_end->dumpInXfigFile(stream,resolution,box);
  direction?_end->dumpInXfigFile(stream,resolution,box):_start->dumpInXfigFile(stream,resolution,box);
  stream << endl;
}

void EdgeLin::update(Node *m)
{
  updateBounds();
}

double EdgeLin::getNormSq() const
{
  return _start->distanceWithSq(*_end);
}

double EdgeLin::getAreaOfZone() const
{
  return ((*_start)[0]-(*_end)[0])*((*_start)[1]+(*_end)[1])/2.;
}

void EdgeLin::getBarycenter(double *bary) const
{
  bary[0]=((*_start)[0]+(*_end)[0])/2.;
  bary[1]=((*_start)[1]+(*_end)[1])/2.;
}

double EdgeLin::getCurveLength() const
{
  double x=(*_start)[0]-(*_end)[0];
  double y=(*_start)[1]-(*_end)[1];
  return sqrt(x*x+y*y);
}

Edge *EdgeLin::buildEdgeLyingOnMe(Node *start, Node *end, bool direction) const
{
  return new EdgeLin(start,end,direction);
}

/*!
 * No precision should be introduced here. Just think as if precision was perfect.
 */
void EdgeLin::updateBounds()
{
  _bounds.setValues(fmin((*_start)[0],(*_end)[0]),fmax((*_start)[0],(*_end)[0]),fmin((*_start)[1],(*_end)[1]),fmax((*_start)[1],(*_end)[1]));
}
