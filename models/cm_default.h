/*
 *  cm_default.h
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef CM_DEFAULT_H
#define CM_DEFAULT_H

// Includes from nestkernel:
#include "archiving_node.h"
#include "event.h"
#include "nest_types.h"
#include "universal_data_logger.h"

#include "cm_tree.h"
#include "cm_compartmentcurrents.h"

namespace nest
{

/* BeginUserDocs: neuron, compartmental model

Short description
+++++++++++++++++

A neuron model with user-defined dendrite structure.
Currently, AMPA, GABA or AMPA+NMDA receptors.

Description
+++++++++++

`cm_default` is an implementation of a compartmental model. The structure of the
neuron -- soma, dendrites, axon -- is user-defined at runtime by adding
compartments through `nest.SetStatus()`. Each compartment can be assigned
receptors, also through `nest.SetStatus()`.

The default model is passive, but sodium and potassium currents can be added
by passing non-zero conductances 'g_Na' and 'g_K' with the parameter dictionary
when adding compartments. Receptors can be AMPA and/or NMDA (excitatory), and
GABA (inhibitory). Ion channel and receptor currents to the compartments can be
customized through NESTML

Usage
+++++
The structure of the dendrite is user defined. Thus after creation of the neuron
in the standard manner

.. code-block:: Python
    cm =  nest.Create('cm_default')

compartments can using the `nest.SetStatus()` function

.. code-block:: Python
    nest.SetStatus(cm, {"compartments": [{"parent_idx": -1,
                                          "params": {"e_L": -65.}},
                                          {"parent_idx": 0,
                                          "params": {"e_L": -60., "g_C": 0.02}})

Each compartment is assigned an index, corresponding to the order in which they
were added. Subsequently, compartment indices are used to specify parent
compartments in the tree, or are used to assign receptors to the compartments.
By convention, the first compartment is the root (soma), which has no parent.
In this case, `parent_index` is ``-1``.

To add receptors to the model, also use `nest.SetStatus()`

.. code-block:: Python
    nest.SetStatus(cm, {"receptors": {"comp_idx": 1,
                                      "receptor_type": "AMPA",
                                      "params": {"e_AMPA": 0., "tau_AMPA": 3.})

Similar to compartments, each receptor is assigned an index, starting at 0 and
corresponding to the order in which they are added. This index is used
subsequently to connect synapses to the receptor

.. code-block:: Python
    nest.Connect(pre, cm_model, syn_spec={
        'synapse_model': 'static_synapse', 'weight': 5., 'delay': 0.5,
        'receptor_type': 2})

Note:
  In the `nest.SetStatus()` call, the 'receptor_type' entry is string that
  specifies the type of receptor. In the `nest.Connect()` call, the
  'receptor_type' entry is and integer that specifies the receptor index.

Note:
  The "compartments" resp. "receptors" entries can be a dict or a list
  of dicts containing compartment resp. receptor details. When a dict is provided,
  a single compartment resp. receptor is added to the model. When a list of dicts
  is provided, multiple compartments resp. receptors are added with a single
  `nest.SetStatus()` call.

Compartment voltages can be recorded. To do so, create a multimeter in the
standard manner but specify the to be recorded voltages as
'v_comp{compartment_index}'. Ion channels state variables can be recorded as well,
using the syntax '{state_variable_name}{compartment_index}'. For receptor state
variables, use the receptor index '{state_variable_name}{receptor_index}' i.e.

.. code-block:: Python
    mm = nest.Create('multimeter', 1, {'record_from': ['v_comp0'}, ...})

Current generators can be connected to the model. In this case, the receptor
type is the compartment index, i.e.

.. code-block:: Python
    dc = nest.Create('dc_generator', {...})
    nest.Connect(dc, cm, syn_spec={..., 'receptor_type': 0}

Parameters
++++++++++

The following parameters can be set in the status dictionary.

=========== ======= ===========================================================
 V_th       mV      Spike threshold (default: -55.0 mV)
=========== ======= ===========================================================

The following parameters can be set using the `AddCompartment` function

=========== ======= ===============================================================
 C_m        uF      Capacitance of compartment (default: 1 uF)
 g_C        uS      Coupling conductance with parent compartment (default: 0.01 uS)
 g_L        uS      Leak conductance of the compartment (default: 0.1 uS)
 e_L        mV      Leak reversal of the compartment (default: -70. mV)
=========== ======= ===============================================================

Ion channels and receptor types for the default model are hardcoded.
For ion channels, there is a Na-channel and a K-channel. Parameters can be set
by specifying the following entries in the `SetStatus` dictionary argument:

=========== ======= ===========================================================
 gbar_Na    uS      Maximal conductance Na channel (default: 0 uS)
 e_Na       mV      Reversal Na channel default (default: 50 mV)
 gbar_K     uS      Maximal conductance K channel (default: 0 uS)
 e_K        mV      Reversal K channel (default: 50 mV)
=========== ======= ===========================================================

For receptors, the choice is from 'AMPA', 'GABA' or 'NMDA' or 'AMPA_NMDA'.
Ion channels and receptor types can be customized with NESTML.

If `receptor_type` is AMPA

=========== ======= ===========================================================
 e_AMPA     mV      AMPA reversal (default 0 mV)
 tau_r_AMPA ms      AMPA rise time (default .2 ms)
 tau_d_AMPA ms      AMPA decay time (default 3. ms)
=========== ======= ===========================================================

If `receptor_type` is GABA

=========== ======= ===========================================================
 e_GABA     mV      GABA reversal (default 0 mV)
 tau_r_GABA ms      GABA rise time (default .2 ms)
 tau_d_GABA ms      GABA decay time (default 10. ms)
=========== ======= ===========================================================

If `receptor_type` is NMDA

=========== ======= ===========================================================
 e_NMDA     mV      NMDA reversal (default 0 mV)
 tau_r_NMDA ms      NMDA rise time (default .2 ms)
 tau_d_NMDA ms      NMDA decay time (default 43. ms)
=========== ======= ===========================================================

If `receptor_type` is AMPA_NMDA

============ ======= ===========================================================
 e_AMPA_NMDA mV      NMDA reversal (default 0 mV)
 tau_r_AMPA  ms      AMPA rise time (default .2 ms)
 tau_d_AMPA  ms      AMPA decay time (default 3. ms)
 tau_r_NMDA  ms      NMDA rise time (default .2 ms)
 tau_d_NMDA  ms      NMDA decay time (default 43. ms)
 NMDA_ratio  (1)     Ratio of NMDA versus AMPA channels
============ ======= ===========================================================

Sends
+++++

SpikeEvent

Receives
++++++++

SpikeEvent, CurrentEvent, DataLoggingRequest

References
++++++++++

Data-driven reduction of dendritic morphologies with preserved dendro-somatic responses
WAM Wybo, J Jordan, B Ellenberger, UM Mengual, T Nevian, W Senn
Elife 10, `e60936 <https://elifesciences.org/articles/60936>`_

See also
++++++++

NEURON simulator ;-D

EndUserDocs*/

class cm_default : public ArchivingNode
{

public:
  cm_default();
  cm_default( const cm_default& );

  using Node::handle;
  using Node::handles_test_event;

  port send_test_event( Node&, rport, synindex, bool );

  void handle( SpikeEvent& );
  void handle( CurrentEvent& );
  void handle( DataLoggingRequest& );

  port handles_test_event( SpikeEvent&, rport );
  port handles_test_event( CurrentEvent&, rport );
  port handles_test_event( DataLoggingRequest&, rport );

  void get_status( DictionaryDatum& ) const;
  void set_status( const DictionaryDatum& );

private:
  void add_compartment_( DictionaryDatum& dd );
  void add_receptor_( DictionaryDatum& dd );

  void init_recordables_pointers_();
  void calibrate();

  void update( Time const&, const long, const long );

  CompTree c_tree_;
  std::vector< RingBuffer > syn_buffers_;

  // To record variables with DataAccessFunctor
  double
  get_state_element( size_t elem )
  {
    return *recordables_values[ elem ];
  };

  // The next classes need to be friends to access the State_ class/member
  friend class DataAccessFunctor< cm_default >;
  friend class DynamicRecordablesMap< cm_default >;
  friend class DynamicUniversalDataLogger< cm_default >;

  /*
  internal ordering of all recordables in a vector
  the vector 'recordables_values' stores pointers to all state variables
  present in the model
  */
  std::vector< Name > recordables_names;
  std::vector< double* > recordables_values;

  //! Mapping of recordables names to access functions
  DynamicRecordablesMap< cm_default > recordablesMap_;
  //! Logger for all analog data
  DynamicUniversalDataLogger< cm_default > logger_;

  double V_th_;
};


inline port
nest::cm_default::send_test_event( Node& target, rport receptor_type, synindex, bool )
{
  SpikeEvent e;
  e.set_sender( *this );
  return target.handles_test_event( e, receptor_type );
}

inline port
cm_default::handles_test_event( SpikeEvent&, rport receptor_type )
{
  if ( ( receptor_type < 0 ) or ( receptor_type >= static_cast< port >( syn_buffers_.size() ) ) )
  {
    throw IncompatibleReceptorType( receptor_type, get_name(), "SpikeEvent" );
  }
  return receptor_type;
}

inline port
cm_default::handles_test_event( CurrentEvent&, rport receptor_type )
{
  // if get_compartment returns nullptr, raise the error
  if ( not c_tree_.get_compartment( long( receptor_type ), c_tree_.get_root(), 0 ) )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return receptor_type;
}

inline port
cm_default::handles_test_event( DataLoggingRequest& dlr, rport receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return logger_.connect_logging_device( dlr, recordablesMap_ );
}

} // namespace

#endif /* #ifndef CM_DEFAULT_H */
