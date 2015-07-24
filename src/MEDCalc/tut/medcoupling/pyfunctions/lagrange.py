#!/usr/bin/env python
# Copyright (C) 2012-2015  CEA/DEN, EDF R&D
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

import scipy
def lagrange(points):
    '''
    This returns a polynom that fits the points specified in the input
    dictionary (lagrange interpolation). In this dictionary, the keys
    are x value and the values are y corresponding values
    (i.e. y=polynom(x)). The polynom is a scipy polynom and then is a
    callable (can be used as a function).
    '''
    tmp = scipy.poly1d([0])
    result=scipy.poly1d([0])
    
    for i in points.keys():
        numerator=scipy.poly1d([1])
        denom = 1.0
        for j in points.keys():
            if (i != j):
                tmp = scipy.poly1d([1,-j])
                numerator = numerator * tmp
                denom = denom * (i - j)
        tmp = (numerator/denom) * points.get(i)
        result = result + tmp

    return result

def points_usingfunction(arrX,function):
    points={}
    for x in arrX:
        points[x] = function(x)
    return points

def points_usingarray(arrX,arrY):
    points={}
    for i in range(len(arrX)):
        x=arrX[i]
        y=arrY[i]
        points[x] = y
    return points

def sortdict(points):
    # Sort this dictionary by keys and returns 2 lists, the list of X
    # and the list of Y, the whole ordered by X
    keys = points.keys()
    keys.sort()
    return keys, [points[key] for key in keys]

import pylab
import numpy
def plot(function, start=0., stop=1., step=0.01):
    """
    The parameter function must be a callable.
    """
    arrX=numpy.arange(start, stop, step, dtype='float64')
    # function is a callable
    arrY=map(function,arrX)
    pylab.plot(arrX, arrY)
    pylab.show()


# ---
# The points does not need to be ordered by x values
def TEST_lagrange_01():
    input = {0.:5, 0.2:10, 0.9:10, 0.6:21, 1:8} 
    polynom = lagrange(input)
    print polynom 
    plot(function=polynom, start=0., stop=1., step=0.001)

def TEST_lagrange_02():
    input = {0.:0., 0.5:1., 1.:0.} 
    polynom = lagrange(input)
    print polynom 
    plot(function=polynom, start=0., stop=1., step=0.001)

# ---
# One can create the input dictionary  from arrays
def TEST_lagrange_usingarrays_01():
    arrX = [0., 0.2, 0.9, 0.6, 1] 
    arrY = [5, 10, 10, 21, 8]
    input = points_usingarray(arrX,arrY)
    polynom = lagrange(input)
    print polynom 
    plot(function=polynom, start=0., stop=1., step=0.001)

# Another example using numpy
def TEST_lagrange_usingarrays_02():
    arrX=numpy.arange(start=0., stop=1., step=0.1, dtype='float64')
    arrY=numpy.zeros(len(arrX), dtype='float64')
    arrY[3]=2
    input = points_usingarray(arrX,arrY)
    polynom = lagrange(input)
    print polynom 
    plot(function=polynom, start=0., stop=1., step=0.001)

# ---
# One can create the input dictionary  from a function applied to an
# array of X values

# simple method for mathematical functions
def TEST_lagrange_usingfunction_01():
    arrX=numpy.arange(start=0., stop=1., step=0.1, dtype='float64')
    arrY=numpy.cos(10*arrX)
    input = points_usingarray(arrX,arrY)
    polynom = lagrange(input)
    print polynom
    plot(function=polynom, start=0., stop=1., step=0.001)

# General method
xlimit=0.8
def chapeau(x):
    if x<xlimit:
        y=x
    else:
        y=2*xlimit-x
    return y

def TEST_lagrange_usingfunction_01():
    arrX=numpy.arange(start=0., stop=1., step=0.1, dtype='float64')
    input = points_usingfunction(arrX,chapeau)
    polynom = lagrange(input)
    print polynom
    plot(function=polynom, start=0., stop=1., step=0.001)


if __name__ == "__main__":
    #TEST_lagrange_01()
    TEST_lagrange_02()
    #TEST_lagrange_usingarrays_01()
    #TEST_lagrange_usingarrays_02()
    #TEST_lagrange_usingfunction_01()
    #TEST_lagrange_usingfunction_01()
