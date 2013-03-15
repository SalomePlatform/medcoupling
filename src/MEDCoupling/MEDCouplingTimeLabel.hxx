// Copyright (C) 2007-2012  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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

#ifndef __PARAMEDMEM_TIMELABEL_HXX__
#define __PARAMEDMEM_TIMELABEL_HXX__

#include "MEDCoupling.hxx"

namespace ParaMEDMEM
{
  /*!
   * Class representing a label of time of the lastely modified part of this.
   * More _time is high more the object has been modified recently.
   */
  class MEDCOUPLING_EXPORT TimeLabel
  {
  public:
    TimeLabel& operator=(const TimeLabel& other);
    //! This method should be called when write access has been done on this.
    void declareAsNew() const;
    //! This method should be called on high level classes as Field or Mesh to take into acount modifications done in aggregates objects.
    virtual void updateTime() const = 0;
    unsigned int getTimeOfThis() const { return _time; }
  protected:
    TimeLabel();
    virtual ~TimeLabel();
    void updateTimeWith(const TimeLabel& other) const;
    void forceTimeOfThis(const TimeLabel& other) const;
  private:
    static unsigned int GLOBAL_TIME;
    mutable unsigned int _time;
  };
}

#endif
