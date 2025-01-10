#! /usr/bin/env python3
#  -*- coding: utf-8 -*-
# Copyright (C) 2021-2025  CEA, EDF
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
# Author : Anthony GEAY (EDF R&D)

from iterative_stats.iterative_mean import IterativeMean
from iterative_stats.iterative_variance import IterativeVariance
from iterative_stats.iterative_covariance import IterativeCovariance
from iterative_stats.sensitivity.sensitivity_saltelli import IterativeSensitivitySaltelli

import numpy as np
try:
    import openturns as ot
    have_ot = True
except ImportError:
    have_ot = False


class IterativeFieldMoments:
    """
    Compute statistics over fields.
    
    Parameters
    ----------
    enable_mean : bool, optional
        Whether to aggregate mean
        Default is True
    enable_variance : bool, optional
        Whether to aggregate variance
        Default is True
    enable_covariance : bool, optional
        Whether to aggregate covariance
        Default is True
    """

    def __init__(self, enable_mean=True, enable_variance=True, enable_covariance=True):
        self._dim = None
        self._size = None
        self._enable_mean = enable_mean
        self._enable_variance = enable_variance
        self._enable_covariance = enable_covariance

    def increment(self, field):
        """
        Increment state.

        Parameters
        ----------
        field : medcoupling.MEDCouplingFieldDouble
            Field to compute statistics from
        """
        values = field.getArray()

        if self._dim is None:
            self._dim = values.getNumberOfComponents()
            self._size = values.getNumberOfTuples()
            if have_ot:
                self._agg_mean = [ot.IterativeMoments(1, self._dim) for k in range(self._size)] if self._enable_mean else None
                self._agg_variance = [ot.IterativeMoments(2, self._dim) for k in range(self._size)] if self._enable_variance else None
            else:
                self._agg_mean = [IterativeMean(dim=self._dim) for k in range(self._size)] if self._enable_mean else None
                self._agg_variance = [IterativeVariance(dim=self._dim) for k in range(self._size)] if self._enable_variance else None

            if self._enable_covariance:
                n_tri = (self._dim * (self._dim + 1) // 2)
                self._agg_covariance = [[IterativeCovariance(dim=1) for i in range(n_tri)] for k in range(self._size)]
            else:
                self._agg_covariance = None

        if values.getNumberOfComponents() != self._dim:
            raise ValueError(f"Incorrect number of components {values.getNumberOfComponents()} expected {self._dim}")
        if values.getNumberOfTuples() != self._size:
            raise ValueError(f"Incorrect number of tuples {values.getNumberOfTuples()} expected {self._size}")
        for k in range(self._size):
            tk = values.getTuple(k)
            if self._enable_mean:
                self._agg_mean[k].increment(tk)
            if self._enable_variance:
                self._agg_variance[k].increment(tk)
            if self._enable_covariance:
                for i in range(self._dim):
                    for j in range(i + 1):
                        self._agg_covariance[k][i * (i + 1) // 2 + j].increment(tk[i], tk[j])

    def mean(self):
        """
        Mean.

        Returns
        -------
        mean : numpy.array of shape (n, d)
            Mean field
        """
        if self._agg_mean is None:
            raise ValueError(f"No data aggregated")

        mean = np.zeros((self._size, self._dim))
        for i in range(self._size):
            if have_ot:
                mean[i] = self._agg_mean[i].getMean()
            else:
                mean[i] = self._agg_mean[i].get_stats()
        return mean

    def variance(self):
        """
        Variance.

        Returns
        -------
        variance : numpy.array
            Variance per component
        """
        if self._agg_variance is None:
            raise ValueError(f"No data aggregated")

        variance = np.zeros((self._size, self._dim))
        for i in range(self._size):
            if have_ot:
                variance[i] = self._agg_variance[i].getVariance()
            else:
                variance[i] = self._agg_variance[i].get_stats()
        return variance

    def stddev(self):
        """
        Standard deviation.

        Returns
        -------
        stddev : numpy.array of shape (n, d)
            Standard deviation per component
        """
        return np.sqrt(self.variance())

    def covariance(self):
        """
        Compute the variance-covariance matrix.

        Returns
        -------
        cov : numpy.array of shape (n, d, d)
            Variance-covariance matrix
        """
        if self._agg_covariance is None:
            raise ValueError(f"No data aggregated")

        covariance = np.zeros((self._size, self._dim, self._dim))
        for k in range(self._size):
            for i in range(self._dim):
                for j in range(i + 1):
                    covij = self._agg_covariance[k][i * (i + 1) // 2 + j].get_stats()[0]
                    covariance[k, i, j] = covij
                    if i != j:
                        covariance[k, j, i] = covij
        return covariance


class IterativeFieldSobol:
    """
    Iterative Sobol indices.
    
    Parameters
    ----------
    nb_params : int
        Number of parameters of the field
    """
    def __init__(self, nb_parms: int):
        if nb_parms < 2:
            raise ValueError(f"Got {nb_parms} parameters, expected at least 2")
        self._nb_parms = nb_parms
        self._dim = None
        self._size = None
        self._agg_sobol = None

    def increment(self, fields):
        """
        Update data.

        Parameters
        ----------
        fields : List[medcoupling.MEDCouplingFieldDouble]
            List of fields for each pick-freeze combination of the parameters
        """

        if len(fields) != self._nb_parms + 2:
            raise ValueError(f"Got {len(fields)} fields, expected {self._nb_parms + 2}")
        if self._dim == None:
            self._dim = fields[0].getArray().getNumberOfComponents()
            self._size = fields[0].getArray().getNumberOfTuples()
            self._agg_sobol = [IterativeSensitivitySaltelli(self._nb_parms, dim=self._dim) for k in range(self._size)]

        for k in range(self._size):
            tks = np.array([field.getArray().getTuple(k) for field in fields])
            self._agg_sobol[k].increment(tks.reshape((self._nb_parms + 2,)))

    def indices(self, k):
        """
        Sobol indices accessor.

        Parameters
        ----------
        k : int
            Tuple index in the field

        Returns
        -------
        first, total : np.array
            First and total order Sobol' indices of shape (nb_parms, dim)
        """
        if self._agg_sobol is None:
            raise ValueError("No data aggregated")
        result = self._agg_sobol[k].getFirstOrderIndices(), self._agg_sobol[k].getTotalOrderIndices()
        return result
