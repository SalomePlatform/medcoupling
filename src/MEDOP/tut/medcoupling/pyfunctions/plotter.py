#!/usr/bin/env python

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
