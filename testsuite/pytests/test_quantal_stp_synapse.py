# -*- coding: utf-8 -*-
#
# test_quantal_stp_synapse.py
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

# This script compares the two variants of the Tsodyks/Markram synapse in NEST.

import nest
import numpy
import unittest


@nest.ll_api.check_stack
class QuantalSTPSynapseTestCase(unittest.TestCase):
    """Compare quantal_stp_synapse with its deterministic equivalent."""

    def test_QuantalSTPSynapse(self):
        """Compare quantal_stp_synapse with its deterministic equivalent"""
        nest.ResetKernel()
        nest.rng_seed = 1
        nest.set_verbosity(100)
        n_syn = 12  # number of synapses in a connection
        n_trials = 100  # number of measurement trials

        # parameter set for facilitation
        fac_params = {"U": 0.03, "u": 0.03, "tau_fac": 500.0, "tau_rec": 200.0, "weight": 1.0}

        # Here we assign the parameter set to the synapse models
        t1_params = fac_params  # for tsodyks2_synapse
        t2_params = t1_params.copy()  # for furhmann_synapse

        t1_params["synapse_model"] = "tsodyks2_synapse"

        t2_params["n"] = n_syn
        t2_params["weight"] = 1.0 / n_syn
        t2_params["synapse_model"] = "quantal_stp_synapse"

        source = nest.Create("spike_generator")
        source.spike_times = [
            30.0,
            60.0,
            90.0,
            120.0,
            150.0,
            180.0,
            210.0,
            240.0,
            270.0,
            300.0,
            330.0,
            360.0,
            390.0,
            900.0,
        ]

        parrot = nest.Create("parrot_neuron")
        neuron = nest.Create("iaf_psc_exp", 2, params={"tau_syn_ex": 3.0, "tau_m": 70.0})

        # We must send spikes via parrot because devices cannot
        # connect through plastic synapses
        # See #478.
        nest.Connect(source, parrot)
        nest.Connect(parrot, neuron[:1], syn_spec=t1_params)
        nest.Connect(parrot, neuron[1:], syn_spec=t2_params)

        voltmeter = nest.Create("voltmeter", 2)

        t_tot = 1500.0

        # the following is a dry run trial so that the synapse dynamics is
        # idential in all subsequent trials.

        nest.Simulate(t_tot)

        # Now we connect the voltmeters
        nest.Connect(voltmeter[:1], neuron[:1])
        nest.Connect(voltmeter[1:], neuron[1:])

        for t in range(n_trials):
            t_net = nest.biological_time
            nest.SetStatus(source, {"origin": t_net})
            nest.Simulate(t_tot)

        nest.Simulate(0.1)  # flush the last voltmeter events from the queue

        vm = numpy.array(voltmeter[1].events["V_m"])
        vm_reference = numpy.array(voltmeter[0].events["V_m"])

        assert len(vm) % n_trials == 0
        n_steps = int(len(vm) / n_trials)
        vm.shape = (n_trials, n_steps)
        vm_reference.shape = (n_trials, n_steps)

        vm_mean = numpy.mean(vm, axis=0)
        vm_ref_mean = numpy.mean(vm_reference, axis=0)

        error = numpy.sqrt((vm_ref_mean - vm_mean) ** 2)
        self.assertLess(numpy.max(error), 4.0e-4)


def suite():
    suite = unittest.makeSuite(QuantalSTPSynapseTestCase, "test")
    return suite


def run():
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())


if __name__ == "__main__":
    run()
