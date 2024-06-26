# -*- coding: utf-8 -*-
#
# connex.py
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
Circular mask and flat probability
----------------------------------

Create two populations on a 30x30 grid of iaf_psc_alpha neurons,
connect with circular mask, flat probability,
visualize.

BCCN Tutorial @ CNS*09
Hans Ekkehard Plesser, UMB
"""

import nest
import matplotlib.pyplot as plt
import numpy as np

nest.ResetKernel()

pos = nest.spatial.grid(shape=[30, 30], extent=[3.0, 3.0])

#######################################################
# create and connect two populations
a = nest.Create("iaf_psc_alpha", positions=pos)
b = nest.Create("iaf_psc_alpha", positions=pos)

cdict = {"rule": "pairwise_bernoulli", "p": 0.5, "mask": {"circular": {"radius": 0.5}}}

nest.Connect(a, b, conn_spec=cdict, syn_spec={"weight": nest.random.uniform(0.5, 2.0)})

#################################################################
# plot targets of neurons in different grid locations

# first, clear existing figure, get current figure
plt.clf()
fig = plt.gcf()

# plot targets of two source neurons into same figure, with mask
for src_index in [30 * 15 + 15, 0]:
    # obtain node id for center
    src = a[src_index : src_index + 1]
    nest.PlotTargets(src, b, mask=cdict["mask"], fig=fig)

# beautify
plt.axes().set_xticks(np.arange(-1.5, 1.55, 0.5))
plt.axes().set_yticks(np.arange(-1.5, 1.55, 0.5))
plt.grid(True)
plt.axis([-2.0, 2.0, -2.0, 2.0])
plt.axes().set_aspect("equal", "box")
plt.title("Connection targets")

plt.show()

# plt.savefig('connex.pdf')
