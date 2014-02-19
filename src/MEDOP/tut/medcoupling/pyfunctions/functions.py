#!/usr/bin/env python
# Copyright (C) 2012-2014  CEA/DEN, EDF R&D
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
# See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#


class Function:
    def __init__(self, **kwargs):
        self.kwargs = kwargs

    def function(self, x, **kwargs):
        # This should be implemented in a derived class
        raise Runtime("function is not implemented yet")

    def __call__(self, x):
        # The argument can be a scalar or a list, we have to check
        # that first.
        if isIterable(x):
            y = map(self,x)
        else:
            y = self.function(x, **self.kwargs)
        return y

def isIterable(x):
    """
    This returns True if the parameter is an iterable, a list or an
    array, or any collection object.
    """
    try:
        len(x)
        return True
    except TypeError, e:
        return False

#
# =====================================================
# Implementation of standard functions. All function are normalized
# function: the xrange is [0,1], the yrange is [0,1]
import numpy
from scipy.constants import pi

# Note that in this situation, the create another constructor because
# the parameters can be deduced from one single parameter xlimit. The
# constructor must create a kwargs dictionary that map the arguments
# of the method "function"
class FuncConique(Function):
    def __init__(self,xlimit):
        a = 1./(xlimit*xlimit-2*xlimit+1)
        b = -2.*a
        c = a
        d = 1/(xlimit*xlimit)
        # We call the super constructor to redefine the kwarg
        # attribute, so that it fits with the arguments of the method
        # "function":
        Function.__init__(self,xlimit=xlimit,a=a,b=b,c=c,d=d)
        # NOTE: Instead of calling the super constructor, we could
        # redefine directly the kwargs attribute:
        #self.kwargs = {"xlimit":xlimit,
        #               "a":a, "b":b,
        #               "c":c, "d":d}
    
    def function(self,x,xlimit,a,b,c,d):
        if x<xlimit:
            y=d*x*x
        else:
            y=a*x*x+b*x+c
        return y


class FuncChapeau(Function):
    def function(self,x,xlimit):
        if x<xlimit:
            y=x/xlimit
        else:
            y=(x-1)/(xlimit-1)
        return y

class FuncStiffExp(Function):
    """
    xlimit : the x position of the top of the function
    stiffness : the higher it is, the stiffer the function is
    """
    def function(self,x,xlimit,stiffness):
        if x<xlimit:
            y=numpy.exp(stiffness*(x-xlimit))
        else:
            y=numpy.exp(-stiffness*(x-xlimit))
        return y

class FuncCosinus(Function):
    def __init__(self,nbPeriods):
        # The pulsation w must be choosen so that w*xmax=n*2pi where
        # xmax=1 and n is an integer that corresponds to the number of
        # oscilations on the xrange [0,xmax].
        w=nbPeriods*2*pi
        Function.__init__(self,w=w)
    
    def function(self,x,w):
        y=numpy.cos(w*x)
        return y

class FuncStiffPulse(Function):
    def __init__(self,xlimit, stiffness, nbPeriods):
        self.stiffexp=FuncStiffExp(xlimit=xlimit,stiffness=stiffness)
        self.cosinus=FuncCosinus(nbPeriods=nbPeriods)
        Function.__init__(self)

    def function(self,x):
        y=self.stiffexp(x)*numpy.abs(self.cosinus(x))
        return y

class FuncHeaviside(Function):
    def function(self,x,xlimit):
        if x<xlimit:
            y=0
        else:
            y=1
        return y

class FuncPorte(Function):
    def function(self,x,xinf,xsup):
        if x<xinf or x>xsup:
            y=0
        else:
            y=1
        return y
        
import lagrange
class FuncLagrange(Function):
    def __init__(self,points):
        """
        @points : a dictionary whose keys are x values and values are
        y values to be considered as fixed points for interpolation.  
        """
        Function.__init__(self)
        self.polynom = lagrange.lagrange(points)

    def function(self,x):
        return self.polynom(x)

#
# =====================================================
# Unit tests functions
# =====================================================
#
class MyFunction(Function):
    def function(self,x,a,b):
        y=a*x+b
        return y

def TEST_Function():
    # The parameters of the constructor of MyFunction must be
    # consistent with the kwargs parameters of the method function of
    # the class MyFunction (it must map exactly).
    f=MyFunction(a=3.,b=7.)

    x=2
    y_ref = 3.*x+7.
    y_res = f(x)
    print y_ref
    print y_res
    if y_ref != y_res:
        print "ERR"
    else:
        print "OK"

def TEST_Function_withIterable():
    f=MyFunction(a=3.,b=1.)
    
    arrX = [0., 1., 2., 3.]
    arrY = f(arrX)

    arrY_ref = [1., 4., 7., 10.]
    print "arrY res =%s"%arrY
    print "arrY ref =%s"%arrY_ref
    
def TEST_FuncConique():
    f=FuncConique(xlimit=0.3)
    from plotter import plot
    plot(f)

def TEST_FuncChapeau():
    f=FuncChapeau(xlimit=0.3)
    from plotter import plot
    plot(f)

def TEST_FuncStiffExp():
    f=FuncStiffExp(xlimit=0.3,stiffness=20.)
    from plotter import plot
    plot(f)

def TEST_FuncCosinus():
    f=FuncCosinus(nbPeriods=20)
    from plotter import plot
    plot(f, step=0.001)

def TEST_FuncStiffPulse():
    f=FuncStiffPulse(xlimit=0.3,stiffness=50,nbPeriods=15)
    from plotter import plot
    plot(f, step=0.001)

def TEST_FuncHeaviside():
    f=FuncHeaviside(xlimit=0.3)
    from plotter import plot
    plot(f)

def TEST_FuncPorte():
    f=FuncPorte(xinf=0.3,xsup=0.4)
    from plotter import plot
    plot(f)

def TEST_customize_01():
    f=FuncStiffPulse(xlimit=0.3,stiffness=40,nbPeriods=20)

    # One can customize the final function as follow (in this example,
    # a linear transform)
    def myfunc(x):
        y=5*f(x)+2
        return y
    
    from plotter import plot
    plot(myfunc, step=0.001)

def TEST_customize_02():
    f=FuncHeaviside(xlimit=0.3)

    # One can customize the final function as follow (in this example,
    # reverse of heaviside)
    def myfunc(x):
        y=1-f(x)
        return y
    
    from plotter import plot
    plot(myfunc)

def TEST_FuncLagrange():
    points = {0.:5, 0.2:10, 0.9:10, 0.6:21, 1:8} 
    f=FuncLagrange(points)
    from plotter import plot
    plot(f)

if __name__ == "__main__":
    TEST_Function()
    TEST_Function_withIterable()
    #TEST_FuncConique()
    #TEST_FuncChapeau()
    #TEST_FuncStiffExp()
    #TEST_FuncCosinus()
    #TEST_FuncStiffPulse()
    #TEST_FuncHeaviside()
    #TEST_FuncPorte()
    #TEST_customize_01()
    #TEST_customize_02()
    #TEST_FuncLagrange()
