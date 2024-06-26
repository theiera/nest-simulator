# -*- coding: utf-8 -*-
#
# test_threads.py
#
# This file is part of NEST.
#
# Copyright (C) 2004 The NEST Initiative
#
# NEST is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# NEST is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with NEST.  If not, see <http://www.gnu.org/licenses/>.

"""
UnitTests for multithreaded pynest
"""

import unittest
import nest


@nest.ll_api.check_stack
class ThreadTestCase(unittest.TestCase):
    """Tests for multi-threading"""

    def nest_multithreaded(self):
        """Return True, if we have a thread-enabled NEST, False otherwise"""

        return nest.ll_api.sli_func("statusdict/threading :: (no) eq not")

    def test_Threads(self):
        """Multiple threads"""

        if not self.nest_multithreaded():
            self.skipTest("NEST was compiled without multi-threading")

        nest.ResetKernel()
        self.assertEqual(nest.local_num_threads, 1)

        nest.local_num_threads = 8
        n = nest.Create("iaf_psc_alpha", 8)
        st = list(nest.GetStatus(n, "vp"))
        st.sort()
        self.assertEqual(st, [0, 1, 2, 3, 4, 5, 6, 7])

    def test_ThreadsGetConnections(self):
        """GetConnections with threads"""

        if not self.nest_multithreaded():
            self.skipTest("NEST was compiled without multi-threading")

        nest.ResetKernel()
        nest.local_num_threads = 8
        pre = nest.Create("iaf_psc_alpha")
        post = nest.Create("iaf_psc_alpha", 6)

        nest.Connect(pre, post)

        conn = nest.GetConnections(pre)
        # Because of threading, targets may be in a different order than
        # in post, so we sort the vector.
        targets = list(conn.get("target"))
        targets.sort()

        self.assertEqual(targets, post.tolist())

    def test_ThreadsGetEvents(self):
        """Gathering events across threads"""

        if not self.nest_multithreaded():
            self.skipTest("NEST was compiled without multi-threading")

        threads = (1, 2, 4, 8)

        n_events_sr = []
        n_events_vm = []

        N = 128
        Simtime = 1000.0

        for t in threads:
            nest.ResetKernel()
            nest.local_num_threads = t

            # force a lot of spike events
            n = nest.Create("iaf_psc_alpha", N, {"I_e": 2000.0})
            sr = nest.Create("spike_recorder")
            vm = nest.Create("voltmeter")

            nest.Connect(n, sr)
            nest.Connect(vm, n)

            nest.Simulate(Simtime)

            n_events_sr.append(nest.GetStatus(sr, "n_events")[0])
            n_events_vm.append(nest.GetStatus(vm, "n_events")[0])

        ref_vm = N * (Simtime - 1)
        ref_sr = n_events_sr[0]

        # could be done more elegantly with any(), ravel(),
        # but we dont want to be dependent on numpy et al
        [self.assertEqual(x, ref_vm) for x in n_events_vm]
        [self.assertEqual(x, ref_sr) for x in n_events_sr]


def suite():
    suite = unittest.makeSuite(ThreadTestCase, "test")
    return suite


def run():
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())


if __name__ == "__main__":
    run()
