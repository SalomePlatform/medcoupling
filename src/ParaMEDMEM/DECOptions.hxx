// Copyright (C) 2007-2014  CEA/DEN, EDF R&D
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

#ifndef __DECOPTIONS_HXX__
#define __DECOPTIONS_HXX__

#include <string>

namespace ParaMEDMEM
{
  //Enum describing the allToAll method used in the communication pattern
  typedef enum { Native, PointToPoint } AllToAllMethod;
  typedef enum { WithoutTimeInterp, LinearTimeInterp } TimeInterpolationMethod;

  class DECOptions
  {
  protected:
    std::string _method;
    bool _asynchronous;
    TimeInterpolationMethod _timeInterpolationMethod;
    AllToAllMethod _allToAllMethod;
    bool _forcedRenormalization;
  public:
    DECOptions():_method("P0"),
                 _asynchronous(false),
                 _timeInterpolationMethod(WithoutTimeInterp),
                 _allToAllMethod(Native),
                 _forcedRenormalization(false)
    {
    }
    
    DECOptions(const DECOptions& deco)
    {
      _method=deco._method;
      _timeInterpolationMethod=deco._timeInterpolationMethod;
      _asynchronous=deco._asynchronous;
      _forcedRenormalization=deco._forcedRenormalization;
      _allToAllMethod=deco._allToAllMethod;
    }
    
    const std::string& getMethod() const { return _method; }
    void setMethod(const char *m) { _method=m; }

    TimeInterpolationMethod getTimeInterpolationMethod() const { return DECOptions::_timeInterpolationMethod; }
    void setTimeInterpolationMethod(TimeInterpolationMethod it) { DECOptions::_timeInterpolationMethod=it; }

    bool getForcedRenormalization() const { return DECOptions::_forcedRenormalization; }
    void setForcedRenormalization( bool dr) { DECOptions::_forcedRenormalization = dr; }

    bool getAsynchronous() const { return DECOptions::_asynchronous; }
    void setAsynchronous( bool dr) { DECOptions::_asynchronous = dr; }
     
    AllToAllMethod getAllToAllMethod() const { return _allToAllMethod; }
    void setAllToAllMethod(AllToAllMethod sp) { _allToAllMethod=sp; }
  };
}

#endif
