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

#ifndef __DECOPTIONS_HXX__
#define __DECOPTIONS_HXX__

#include <string>

namespace MEDCoupling
{
  //! Enum describing the allToAll method used in the communication pattern
  typedef enum { Native, PointToPoint } AllToAllMethod;
  //! Enum describing the time interpolation method
  typedef enum { WithoutTimeInterp, LinearTimeInterp } TimeInterpolationMethod;

  /*!
   This class groups the various options accepted by all \ref para-dec "DECs" (which all inherit from %DECOptions).

   The following code excerpt shows how to set options on a %DEC :

   \code
   InterpKernelDEC dec(source_group,target_group);
   dec.setForcedRenormalization(true);
   dec.attachLocalField(field);
   dec.synchronize();
   if (source_group.containsMyRank())
     dec.sendData();
   else
     dec.recvData();
   \endcode
   *
   *
   */
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
    

    /*!
     * \sa setMethod()
     */
    const std::string& getMethod() const { return _method; }
    /*!
     * Set interpolation method. Defaults to "P0".
     */
    void setMethod(const char *m) { _method=m; }

    /*!
     * \sa setTimeInterpolationMethod()
     */
    TimeInterpolationMethod getTimeInterpolationMethod() const { return DECOptions::_timeInterpolationMethod; }
    /*!
     * Set time interpolation method. Default to WithoutTimeInterp.
     */
    void setTimeInterpolationMethod(TimeInterpolationMethod it) { DECOptions::_timeInterpolationMethod=it; }

    /*!
     * \sa setForcedRenormalization()
     */
    bool getForcedRenormalization() const { return DECOptions::_forcedRenormalization; }

    /*!
     * Force renormalization of the field after it has been received so that the total sum
     * of the field values are the same on both the sending and the receiving side. Defaults to
     * false.
     */
    void setForcedRenormalization( bool dr) { DECOptions::_forcedRenormalization = dr; }


    /*!
     * \sa setAsynchronous()
     */
    bool getAsynchronous() const { return DECOptions::_asynchronous; }

    /*!
     * Switch to asynchronous data transfer mode. Default is false.
     */
    void setAsynchronous( bool dr) { DECOptions::_asynchronous = dr; }
     
    /*!
     * \sa setAllToAllMethod()
     */
    AllToAllMethod getAllToAllMethod() const { return _allToAllMethod; }
    /*!
     * Set the broadcast method for synchronisation processes. Default to Native.
     */
    void setAllToAllMethod(AllToAllMethod sp) { _allToAllMethod=sp; }
  };
}

#endif
