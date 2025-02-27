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

import medcoupling as mc
from MEDCouplingIterativeStatistics import IterativeFieldMoments, IterativeFieldSobol
from iterative_stats.experimental_design.experiment import AbstractExperiment
from numpy.testing import assert_allclose
import numpy as np
import unittest


class MEDCouplingIterativeStatisticsTest(unittest.TestCase):
    def test1(self):
        #! [UG_MEDCouplingIterativeStatistics_1]
        # check moments on simple linear 3-d field defined on an 1-d mesh
        import medcoupling as mc
        from MEDCouplingIterativeStatistics import IterativeFieldMoments

        # 1-d mesh
        size = int(1e1)
        coords = [i / size for i in range(size)]
        arrX = mc.DataArrayDouble(coords, size, 1)
        mesh = mc.MEDCouplingCMesh()
        mesh.setCoords(arrX)

        # compute statistics
        istats = IterativeFieldMoments()

        # aggregate 3-d fields
        sampling_size = 10
        dim = 3
        for i in range(sampling_size):
            coeff = 1.0 + i / sampling_size
            field = mesh.fillFromAnalytic(mc.ON_NODES, dim, f"2*x*IVec + 3*x*{coeff}*JVec + 5*KVec")
            istats.increment(field)

        # check moments on last node
        mean = istats.mean()
        variance = istats.variance()
        stddev = istats.stddev()
        covariance = istats.covariance()
        #! [UG_MEDCouplingIterativeStatistics_1]

        print(f"mean={mean[-1]}")
        assert mean.shape == (size, dim), f"wrong shape: {mean.shape}"
        assert_allclose(mean[-1], [1.8, 3.915, 5.0])

        print(f"variance={variance[-1]}")
        assert variance.shape == (size, dim), f"wrong shape: {variance.shape}"
        assert_allclose(variance[-1], [0., 0.66825, 0.], atol=1e-9)

        print(f"stddev={stddev[-1]}")
        assert stddev.shape == (size, dim), f"wrong shape: {stddev.shape}"
        assert_allclose(stddev[-1], [0., 0.8174656, 0.], atol=1e-9)

        print(f"covariance={covariance[-1]}")
        assert covariance.shape == (size, dim, dim), f"wrong shape: {covariance.shape}"
        assert_allclose(
            covariance[-1].flatten(), [0.0, 0.0, 0.0, 0.0, 0.66825, 0.0, 0.0, 0.0, 0.0]
        )

    def test2(self):
        # check ref values with OT-only
        try:
            import openturns as ot
        except ImportError:
            return

        size = int(1e1)
        mesh = ot.RegularGrid(0.0, 1.0 / size, size)

        sampling_size = 10
        sample = ot.Sample(sampling_size, 3)
        for i in range(sampling_size):
            coeff = 1.0 + i / sampling_size
            f = ot.SymbolicFunction(["x"], ["2*x", f"3*x*{coeff}", "5"])
            values = f(mesh.getVertices())
            sample[i] = values[-1]

        # check moments on last node
        mean = sample.computeMean()
        print(f"mean={mean}")
        assert_allclose(np.asarray(mean), [1.8, 3.915, 5.0], atol=1e-9)

        variance = sample.computeVariance()
        print(f"variance={variance}")
        assert_allclose(np.asarray(variance), [0., 0.668250, 0.], atol=1e-9)

        covariance = sample.computeCovariance()
        print(f"covariance={covariance}")
        assert_allclose(
            np.asarray(covariance).flatten(), [0.0, 0.0, 0.0, 0.0, 0.66825, 0.0, 0.0, 0.0, 0.0],
            atol=1e-9
        )

    def test3(self):
        # test on bigger mesh

        # 1-d mesh
        size = int(1e4)
        coords = [i / size for i in range(size)]
        arrX = mc.DataArrayDouble(coords, size, 1)
        mesh = mc.MEDCouplingCMesh()
        mesh.setCoords(arrX)

        # compute mean only
        istats = IterativeFieldMoments(enable_variance=False, enable_covariance=False)

        # aggregate 2-d fields
        sampling_size = 10
        dim = 2
        for i in range(sampling_size):
            coeff = 1.0 + i / sampling_size
            field = mesh.fillFromAnalytic(mc.ON_NODES, dim, f"2*x*IVec + 3*x*{coeff}*JVec")
            istats.increment(field)

        mean = istats.mean()
        print(f"mean={mean[-1]}")
        assert mean.shape == (size, dim), f"wrong shape: {mean.shape}"
        assert_allclose(mean[-1], [1.9998 , 4.349565])

    def test4(self):
        #! [UG_MEDCouplingIterativeStatistics_2]

        # Sobol' test over an 1-d mesh
        import medcoupling as mc
        from MEDCouplingIterativeStatistics import IterativeFieldSobol
        from iterative_stats.experimental_design.experiment import AbstractExperiment

        # 1-d mesh [-1; -1]
        size = int(1e2)
        coords = [-1 + 2*i / size for i in range(size)]
        arrX = mc.DataArrayDouble(coords, size, 1)
        mesh = mc.MEDCouplingCMesh()
        mesh.setCoords(arrX)

        # class to sample 2 parameters: in [0,1]^2 for nb_sim repetitions in a pick-freeze manner
        nb_parms = 2
        nb_sim = 300
        class ABSampler(AbstractExperiment):
            def __init__(self):
                super(ABSampler, self).__init__(nb_parms, nb_sim)
                np.random.seed(42)
            def draw(self):
                return np.random.rand(1, self.nb_parms)

        # aggregate the nb_sim simulations
        isobol = IterativeFieldSobol(nb_parms)
        for pf_sample in ABSampler().generator():
            # build the fields of the time-series f(a, b)_x=a+bx^2 for each (a,b) tuple in the pick-freeze sample
            fields = [mesh.fillFromAnalytic(mc.ON_NODES, 1, f"({a}+{b}*x*x)*IVec") for a,b in pf_sample]
            # note that for each simulation we need to instanciate nb_parms+2 fields
            isobol.increment(fields)

        # retrieve the Sobol' indices of f(a, b)_x in the beginning (x=-1) and middle (x=0) tuples of the time grid
        first_beg, total_beg = isobol.indices(0)
        print(f"at 0  : first={first_beg} total={total_beg}")

        first_mid, total_mid = isobol.indices(size // 2)
        print(f"at mid: first={first_mid} total={total_mid}")

        #! [UG_MEDCouplingIterativeStatistics_2]

        assert_allclose(first_beg[0], [0.59680837, 0.56303733])
        assert_allclose(total_beg[0], [0.52188301, 0.47895543])

        assert_allclose(first_mid[0], [0.94759752, -0.00683533], rtol=1e-6)
        assert_allclose(total_mid[0], [1.00739829e+00, -9.41396650e-07])


if __name__ == "__main__":
    import logging
    from iterative_stats.utils.logger import logger
    logger.level = logging.INFO
    unittest.main()
