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

#include "MEDCouplingTimeLabel.hxx"

using namespace ParaMEDMEM;

unsigned int TimeLabel::GLOBAL_TIME=0;

TimeLabel::TimeLabel():_time(GLOBAL_TIME++)
{
}

TimeLabel::~TimeLabel()
{
}

TimeLabel& TimeLabel::operator=(const TimeLabel& other)
{
  _time=GLOBAL_TIME++;
  return *this;
}

void TimeLabel::declareAsNew() const
{
  _time=GLOBAL_TIME++;
}

void TimeLabel::updateTimeWith(const TimeLabel& other) const
{
  if(_time<other._time)
    _time=other._time;
}
