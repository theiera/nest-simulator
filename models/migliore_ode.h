/*
 *  migliore_ode.h
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

#ifndef MIGLIORE_ODE_H
#define MIGLIORE_ODE_H
#include "config.h"

#ifdef HAVE_GSL
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

// Includes from nestkernel:
#include "archiving_node.h"
#include "connection.h"
#include "event.h"
#include "nest_types.h"
#include "ring_buffer.h"
#include "universal_data_logger.h"

namespace nest
{

/* BeginUserDocs: neuron, integrate-and-fire

Short description
+++++++++++++++++

Migliore_Ode neuron model

Description
+++++++++++

Implementation of the spiking neuron model introduced by Migliore_Ode
[1]_. The dynamics are given by:

.. math::

   dV_m/dt &= 0.04 V_m^2 + 5 V_m + 140 - u + I
   du/dt &= a (b V_m - u)


.. math::

   &\text{if}\;\;\; V_m \geq V_{th}:\\
   &\;\;\;\; V_m \text{ is set to } c\\
   &\;\;\;\; u \text{ is incremented by } d\\
   & \, \\
   &v \text{ jumps on each spike arrival by the weight of the spike}

As published in [1]_, the numerics differs from the standard forward Euler
technique in two ways:

1) the new value of :math:`u` is calculated based on the new value of
   :math:`V_m`, rather than the previous value
2) the variable :math:`V_m` is updated using a time step half the size of that
   used to update variable :math:`u`.

This model offers both forms of integration, they can be selected using the
boolean parameter ``consistent_integration``. To reproduce some results
published on the basis of this model, it is necessary to use the published form
of the dynamics. In this case, ``consistent_integration`` must be set to false.
For all other purposes, it is recommended to use the standard technique for
forward Euler integration. In this case, ``consistent_integration`` must be set
to true (default).

For a detailed analysis and discussion of the numerical issues in the original publication, see [2]_.

Parameters
++++++++++

The following parameters can be set in the status dictionary.

======================= =======  ==============================================
 V_m                    mV       Membrane potential
 U_m                    mV       Membrane potential recovery variable
 V_th                   mV       Spike threshold
 I_e                    pA       Constant input current (R=1)
 V_min                  mV       Absolute lower value for the membrane potential
 a                      real     Describes time scale of recovery variable
 b                      real     Sensitivity of recovery variable
 c                      mV       After-spike reset value of V_m
 d                      mV       After-spike reset value of U_m
 consistent_integration boolean  Use standard integration technique
======================= =======  ==============================================

References
++++++++++

.. [1] Izhikevich EM. (2003). Simple model of spiking neurons. IEEE Transactions
       on Neural Networks, 14:1569-1572. DOI: https://doi.org/10.1109/TNN.2003.820440

.. [2] Pauli R, Weidel P, Kunkel S, Morrison A (2018). Reproducing polychronization: A guide to maximizing
       the reproducibility of spiking network models. Frontiers in Neuroinformatics, 12.
       DOI: https://www.frontiersin.org/article/10.3389/fninf.2018.00046

Sends
+++++

SpikeEvent

Receives
++++++++

SpikeEvent, CurrentEvent, DataLoggingRequest

See also
++++++++

iaf_psc_delta, mat2_psc_exp

EndUserDocs */

// forwards the declaration of the function
extern "C" int migliore_ode_dynamics(double, const double*, double*, void* );

class migliore_ode : public ArchivingNode
{

public:
  migliore_ode();
  migliore_ode( const migliore_ode& );

  /**
   * Import sets of overloaded virtual functions.
   * @see Technical Issues / Virtual Functions: Overriding, Overloading, and
   * Hiding
   */
  using Node::handle;
  using Node::handles_test_event;

  size_t send_test_event( Node&, size_t, synindex, bool ) override;

  void handle( SpikeEvent& ) override;
  void handle( CurrentEvent& ) override;
  void handle( DataLoggingRequest& ) override;

  size_t handles_test_event( SpikeEvent&, size_t ) override;
  size_t handles_test_event( CurrentEvent&, size_t ) override;
  size_t handles_test_event( DataLoggingRequest&, size_t ) override;

  void get_status( DictionaryDatum& ) const override;
  void set_status( const DictionaryDatum& ) override;

private:
  // friend class RecordablesMap< migliore_ode >;
  // friend class UniversalDataLogger< migliore_ode >;

  void init_buffers_() override;
  void pre_run_hook() override;

  double tagliorette(double);
  // double migliV(double, double, double, 
  // 		double, double, double, 
  // 		double, double, double,
  // 		int, double);  
  double default_v_ini(double, double);
  double set_v_ini(double, int, double);
  
  // double Iadap(double, double, double, double,
  // 		 double, double, double, double,
  // 	       double, int);
  // double Idep(double, double, double, double, int);
  double exp_cum(double, double, double);
  double monod(double, double, double, double, double);
  double mgblock(double);

  void evolve_ode();
  void update( Time const&, const long, const long );

  // make dynamics function quasi-member
  friend int migliore_ode_dynamics( double, const double*, double*, void* );

  
  // The next two classes need to be friends to access the State_ class/member
  friend class DynamicRecordablesMap< migliore_ode >;
  friend class DynamicUniversalDataLogger< migliore_ode >;
  friend class DataAccessFunctor< migliore_ode >;

  // ----------------------------------------------------------------

  /**
   * Independent parameters of the model.
   */
  struct Parameters_
  {
    double E_L_;
    double Vres_;
    // double vtm_;
    double Cm_;
    double I_th_;
    double tao_m_;
    double sc_;
    double bet_;
    double delta1_;
    double cost_idep_ini_;
    double Idep_ini_vr_;
    double psi1_;
    double a_;
    double b_;
    double c_;
    double alp_;
    double Iadap_start_;
    double istim_min_spikinig_exp_;
    double istim_max_spikinig_exp_;
    double corrcostatratti_;
    double corrcost_;
    double Delta_T; //!< Slope factor in ms
    bool firstSpikeFlag_;
    bool has_connections_;
    double tau_syn_NMDA_;
    double NMDA_ratio_;
    double mg_;
    double mg_ref_;
    double mgb_k_;
    double mgb_shift_;
    double mincurr_;

    /** External DC current */
    double I_e_;

    /** Threshold */
    double V_th_;

    /** Lower bound */
    double V_min_;

    /** Time constants of synaptic currents in ms. */
    std::vector< double > tau_syn_;

    /** Refractory period in ms. */
    double t_ref_;

    size_t n_receptors_() const; //!< Returns the size of tau_syn_

    Parameters_(); //!< Sets default parameter values

    void get( DictionaryDatum& ) const;             //!< Store current values in dictionary
    void set( const DictionaryDatum&, Node* node ); //!< Set values from dicitonary
  };

  // ----------------------------------------------------------------

  /**
   * State variables of the model.
   */
public:
  struct State_
  {
    double V_m_; // membrane potential
    double init_sign_; // Adaptation state variable
    double Iadap_ini_; // Adaptation state variable
    double Idep_ini_; // Adaptation state variable
    double sis_; // Adaptation state variable
    double I_inj_; // input current
    int r_ref_;       //!< Absolute refractory counter (no membrane potential propagation)
    int counter_;
    double current_;
    

    /** Accumulate spikes arriving during refractory period, discounted for
        decay until end of refractory period.
    */

    //! Symbolic indices to the elements of the state vector y_
    enum StateVecElems {
      V_M = 0,
      I_ve,    // 1
      I_SYN, // 2
      Idep,
      Iadap,
      STATE_VEC_SIZE			// This is the minimum state vector size
    };
    //! state vector, must be C-array for GSL solver

    static const size_t NUMBER_OF_FIXED_STATES_ELEMENTS = I_SYN; // V_M, I
    static const size_t NUM_STATE_ELEMENTS_PER_RECEPTOR = 1;     // I_SYN

    double y_[STATE_VEC_SIZE];//  - an array of all the state variables undergoing

    std::vector< double > i_syn_;
    std::vector< double > i_syn_fast_;
    std::vector< double > i_syn_slow_;
    /** Accumulate spikes arriving during refractory period, discounted for
        decay until end of refractory period.
    */

    State_( const Parameters_& ); //!< Default initialization
    // State_( const State_& );

    // State_& operator=( const State_& );

    void get( DictionaryDatum&, const Parameters_& ) const;
    void set( const DictionaryDatum&, const Parameters_&, Node* );
  };

  // ----------------------------------------------------------------

  /**
   * Buffers of the model.
   */
public:
  struct Buffers_
  {
    /**
     * Buffer for recording
     */
    Buffers_( migliore_ode& );
    Buffers_( const Buffers_&, migliore_ode& );
    // UniversalDataLogger< migliore_ode > logger_;

    /** buffers and sums up incoming spikes/currents */
    std::vector< RingBuffer > spikes_;
    RingBuffer currents_;

    /* GSL ODE stuff */
    gsl_odeiv_step* s_;    //!< stepping function
    gsl_odeiv_control* c_; //!< adaptive stepsize control function
    gsl_odeiv_evolve* e_;  //!< evolution function
    gsl_odeiv_system sys_; //!< struct describing system

    // Since IntergrationStep_ is initialized with step_, and the resolution
    // cannot change after nodes have been created, it is safe to place both
    // here.
    double step_;            //!< step size in ms
    double IntegrationStep_; //!< current integration time step, updated by GSL
    double model_alpha;

    inline nest::RingBuffer &get_currents() { return currents; }
    nest::RingBuffer currents;
    //!< Buffer incoming Buffers through delay, as sum
    ;
    double currents_last_value_;

    //! Logger for all analog data
    DynamicUniversalDataLogger< migliore_ode > logger_;

  };

  // ----------------------------------------------------------------

  /**
   * Internal variables of the model.
   */
private:
  struct Variables_
  {
    // time evolution operator
    std::vector< double > P11_syn_fast_;
    std::vector< double > P11_syn_slow_;

    int RefractoryCounts_;

    unsigned int receptor_types_size_;

    std::vector<double> out;
    double mincurr, V_star_min_, alpha_neg_;
    double t_spk;
    double time_scale_, d_dt, dt, t_step;
    double H, Vconvfact_, vrm_;

  };

  // Access functions for UniversalDataLogger -----------------------

  //! Read out the membrane potential
  double
  get_V_m_() const
  {
    return S_.y_[ State_::V_M ];
  }
  //! Read out the injected current
public:
  double
  get_I_inj_() const
  {
    return S_.I_inj_;
  }
  //! Read out the adaptation variable
private:
  double
  get_Iadap_ini_() const
  {
    return S_.y_[ State_::Iadap ];
  }
  //! Read out the adaptation variable
  double
  get_Idep_ini_() const
  {
    return S_.y_[ State_::Idep ];
  }
  //! Read out the init_sign_ variable
  double
  get_init_sign_() const
  {
    return S_.init_sign_;
  }
  //! Read out the sis_ variable
  double
  get_sis_() const
  {
    return S_.sis_;
  }

  // ----------------------------------------------------------------

public:
  Parameters_ P_;
private:
  State_ S_;
  Variables_ V_;
public:
  Buffers_ B_;

  //! Mapping of recordables names to access functions
  DynamicRecordablesMap< migliore_ode > recordablesMap_;

    // Data Access Functor getter
  DataAccessFunctor< migliore_ode > get_data_access_functor( size_t elem );
  inline double
  get_state_element( size_t elem )
  {
    if ( elem == State_::V_M)
    {
      return S_.V_m_;
    }
    else if ( elem == State_::I_ve )
    {
      return S_.current_;
    }
    else if ( elem == State_::Iadap )
    {
      return S_.Iadap_ini_;
    }
    else if ( elem == State_::Idep )
    {
      return S_.Idep_ini_;
    }
    else
    {
      return S_.i_syn_[ elem - S_.NUMBER_OF_FIXED_STATES_ELEMENTS ];
    }
  };

  // Utility function that inserts the synaptic conductances to the
  // recordables map
  Name get_i_syn_name( size_t elem );
  void insert_current_recordables( size_t first = 0 );
  // //! Mapping of recordables names to access functions
  // static RecordablesMap< migliore > recordablesMap_;
/** @} */
};

inline size_t
migliore_ode::Parameters_::n_receptors_() const
{
  return tau_syn_.size();
}

inline size_t
migliore_ode::send_test_event( Node& target, size_t receptor_type, synindex, bool )
{
  SpikeEvent e;
  e.set_sender( *this );

  return target.handles_test_event( e, receptor_type );
}

// inline port
// migliore::handles_test_event( SpikeEvent&, rport receptor_type )
// {
//   if ( receptor_type != 0 )
//   {
//     throw UnknownReceptorType( receptor_type, get_name() );
//   }
//   return 0;
// }

inline size_t
migliore_ode::handles_test_event( CurrentEvent&, size_t receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return 0;
}

inline size_t
migliore_ode::handles_test_event( DataLoggingRequest& dlr, size_t receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return B_.logger_.connect_logging_device( dlr, recordablesMap_ );
}

inline void
migliore_ode::get_status( DictionaryDatum& d ) const
{
  P_.get( d );
  S_.get( d, P_ );
  ArchivingNode::get_status( d );
  ( *d )[ names::recordables ] = recordablesMap_.get_list();
}

// inline void
// migliore_ode::set_status( const DictionaryDatum& d )
// {
//   Parameters_ ptmp = P_;     // temporary copy in case of errors
//   ptmp.set( d, this );       // throws if BadProperty
//   State_ stmp = S_;          // temporary copy in case of errors
//   stmp.set( d, ptmp, this ); // throws if BadProperty

//   // We now know that (ptmp, stmp) are consistent. We do not
//   // write them back to (P_, S_) before we are also sure that
//   // the properties to be set in the parent class are internally
//   // consistent.
//   ArchivingNode::set_status( d );

//   // if we get here, temporaries contain consistent set of properties
//   P_ = ptmp;
//   S_ = stmp;
// }

} // namesopace nest

#endif /* #ifndef MIGLIORE_ODE_H */
#endif /* HAVE GSL */
