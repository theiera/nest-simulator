/*
 *  migliore.h
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

#ifndef MIGLIORE_H
#define MIGLIORE_H

// Includes from nestkernel:
#include "archiving_node.h"
#include "connection.h"
#include "event.h"
#include "nest_types.h"
#include "recordables_map.h"
#include "ring_buffer.h"
#include "universal_data_logger.h"

namespace nest
{

/* BeginUserDocs: neuron, integrate-and-fire

Short description
+++++++++++++++++

Migliore neuron model

Description
+++++++++++

Implementation of the spiking neuron model introduced by Migliore
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

class migliore : public ArchivingNode
{

public:
  migliore();
  migliore( const migliore& );

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
  // friend class RecordablesMap< migliore >;
  // friend class UniversalDataLogger< migliore >;

  void init_buffers_() override;
  void pre_run_hook() override;

  double tagliorette(double);
  double migliV(double, double, double, 
		double, double, double, 
		double, double, double,
		int, double);  
  double default_v_ini(double, double, double);
  double set_v_ini(double, int, double);
  
  double Iadap(double, double, double, double,
		 double, double, double, double,
	       double, int);
  double Idep(double, double, double, double, int);
  double exp_cum(double, double, double);
  double monod(double, double, double, double, double);
  double mgblock(double);
  
  void update( const Time&, const long, const long ) override;

  // The next two classes need to be friends to access the State_ class/member
  friend class DynamicRecordablesMap< migliore >;
  friend class DynamicUniversalDataLogger< migliore >;
  friend class DataAccessFunctor< migliore >;

  // ----------------------------------------------------------------

  /**
   * Independent parameters of the model.
   */
  struct Parameters_
  {
    /** Lower bound */
    double V_min_;
    double E_L_;
    double Vres_;
//    double vtm_;
    /** Threshold */
    double V_th_;
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
    // double tau_syn_fast_rise_;
    // double tau_syn_fast_decay_;
    // double tau_syn_slow_rise_;
    // double tau_syn_slow_decay_;
    double NMDA_ratio_;
    double mg_;
    double mg_ref_;
    double mgb_k_;
    double mgb_shift_;
    double mincurr_;
    double coeffInf_;
    double constInf_;
    double coeffSup_;
    double constSup_;
    double aglif_p_;
    double vinc_inf_;
    double vinc_sup_;
    double zeta_;
    double eta_;
    double rho_;
    double csi_;
    bool plotit_;
    
    /** External DC current */
    double I_e_;

    /** Refractory period in ms. */
    double t_ref_;

    /** Time constants of synaptic currents in ms. */
    // std::vector< double > tau_syn_;
    std::vector< double >  tau_syn_fast_rise_;
    std::vector< double >  tau_syn_fast_decay_;
    std::vector< double >  tau_syn_slow_rise_;
    std::vector< double >  tau_syn_slow_decay_;

    // boolean flag which indicates whether the neuron has connections
    bool has_connections_;

    size_t n_receptors_() const; //!< Returns the size of tau_syn_

    Parameters_(); //!< Sets default parameter values

    void get( DictionaryDatum& ) const;             //!< Store current values in dictionary
    void set( const DictionaryDatum&, Node* node ); //!< Set values from dicitonary
  };

  // ----------------------------------------------------------------

  /**
   * State variables of the model.
   */
  struct State_
  {
    double V_m_; // membrane potential
    double init_sign_; // Adaptation state variable
    double Iadap_ini_; // Adaptation state variable
    double Idep_ini_; // Adaptation state variable
    double sis_; // Adaptation state variable
    double I_inj_; // input current
    int r_ref_;       //!< Absolute refractory counter (no membrane potential propagation)
    double current_; //!< This is the current in a time step.
    
    enum StateVecElems
    {
      V_M = 0,
      I_ve,    // 1
      I_SYN, // 2
      Iadap,
      Idep
    };

    static const size_t NUMBER_OF_FIXED_STATES_ELEMENTS = I_SYN; // V_M, I
    static const size_t NUM_STATE_ELEMENTS_PER_RECEPTOR = 1;     // I_SYN

    std::vector< double > i_syn_;
    std::vector< double > i_syn_fast_rise_A_;
    std::vector< double > i_syn_fast_decay_B_;
    std::vector< double > i_syn_slow_rise_A_;
    std::vector< double > i_syn_slow_decay_B_;
    /** Accumulate spikes arriving during refractory period, discounted for
        decay until end of refractory period.
    */
    State_(); //!< Default initialization

    void get( DictionaryDatum&, const Parameters_& ) const;
    void set( const DictionaryDatum&, const Parameters_&, Node* );
  };

  // ----------------------------------------------------------------

  /**
   * Buffers of the model.
   */
  struct Buffers_
  {
    /**
     * Buffer for recording
     */
    Buffers_( migliore& );
    Buffers_( const Buffers_&, migliore& );
    //UniversalDataLogger< migliore > logger_;

    /** buffers and sums up incoming spikes/currents */
    std::vector< RingBuffer > spikes_;
    RingBuffer currents_;

    //! Logger for all analog data
    DynamicUniversalDataLogger< migliore > logger_;

  };

  // ----------------------------------------------------------------

  /**
   * Internal variables of the model.
   */
  struct Variables_
  {

    // time evolution operator
    std::vector< double > P11_syn_fast_rise_;
    std::vector< double > P11_syn_fast_decay_;
    std::vector< double > P11_syn_slow_rise_;
    std::vector< double > P11_syn_slow_decay_;
//    std::vector< double > P21_syn_;
//    double P20_;
//    double P22_;

    std::vector< double > fast_tp_;
    std::vector< double > slow_tp_;
    std::vector< double > syn_fast_factor_;
    std::vector< double > syn_slow_factor_;

    bool blockActive;

    int RefractoryCounts_;

    unsigned int receptor_types_size_;

    double time_scale_, d_dt, dt, beta2, t_step;
    double V_star_min_, alpha_neg_;
    double H, Vconvfact, vrm, psi1;
    double t_spk;
    double I_th_, Cm_, sc_;
    
    // double GG;
    double C;
    // double JJ, JJ_1, JJ_2, JJ_3, JJ_4, JJ_5, JJ_6, JJ_7;
    double VV_1, VV_2, VV_3, VV_4, VV_5, VV_6, VV_7, VV_8, VV_9, VV_9b, VV_10, VV_11, VV_12, VV_13, VV_14; 
    double pdelta;
    double AA_1, AA_2, AA_3, AA_4, AA_5, AA_6, AA_7, AA_8, AA_9, AA_10, AA_11, AA_12;
    double DD_1;
    std::vector<double> out;
  };

  // Access functions for UniversalDataLogger -----------------------

  //! Read out the membrane potential
  double
  get_V_m_() const
  {
    return S_.V_m_;
  }
  //! Read out the injected current
  double
  get_I_inj_() const
  {
    return S_.I_inj_;
  }
  //! Read out the adaptation variable
  double
  get_Iadap_ini_() const
  {
    return S_.Iadap_ini_;
  }
  //! Read out the adaptation variable
  double
  get_Idep_ini_() const
  {
    return S_.Idep_ini_;
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
  //! Read out the I_syn variable
  double
  get_I_syn_() const
  {
    return S_.I_SYN;
  }

  // ----------------------------------------------------------------
  Parameters_ P_;
  State_ S_;
  Variables_ V_;
  Buffers_ B_;

  //! Mapping of recordables names to access functions
  DynamicRecordablesMap< migliore > recordablesMap_;

  // Data Access Functor getter
  DataAccessFunctor< migliore > get_data_access_functor( size_t elem );
  inline double
  get_state_element( size_t elem )
  {
    if ( elem == State_::V_M )
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
migliore::Parameters_::n_receptors_() const
{
  return tau_syn_fast_decay_.size();
}
  
inline size_t
migliore::send_test_event( Node& target, size_t receptor_type, synindex, bool )
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
migliore::handles_test_event( CurrentEvent&, size_t receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return 0;
}

inline size_t
migliore::handles_test_event( DataLoggingRequest& dlr, size_t receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return B_.logger_.connect_logging_device( dlr, recordablesMap_ );
}

inline void
migliore::get_status( DictionaryDatum& d ) const
{
  P_.get( d );
  S_.get( d, P_ );
  ArchivingNode::get_status( d );
  ( *d )[ names::recordables ] = recordablesMap_.get_list();
}

} // namesopace nest

#endif /* #ifndef MIGLIORE_H */
