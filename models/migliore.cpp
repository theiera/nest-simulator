/*
 *  migliore.cpp
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2.0004 The NEST Initiative
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
 *  Adapted from iaf_psc_exp_multisynapse model by Sergio Solinas in 2023
 */

#include "migliore.h"

// C++ includes:
#include <limits>

// Includes from libnestutil:
#include "dict_util.h"
#include "exceptions.h"
#include "iaf_propagator.h"
#include "kernel_manager.h"
#include "numerics.h"
#include "universal_data_logger_impl.h"

// Includes from sli:
#include "dict.h"
#include "dictutils.h"

/* ----------------------------------------------------------------
 * Recordables map
 * ---------------------------------------------------------------- */

//nest::RecordablesMap< nest::migliore > nest::migliore::recordablesMap_;

namespace nest
{
  // Override the create() method with one call to RecordablesMap::insert_()
  // for each quantity to be recorded.
  template <>
  void
  DynamicRecordablesMap< migliore >::create( migliore& host )
  {
    // use standard names wherever you can for consistency!
    insert( names::V_m, host.get_data_access_functor( migliore::State_::V_M ) );
    //    insert( names::I, host.get_data_access_functor( migliore::State_::I ) );
    insert( names::I_syn, host.get_data_access_functor( migliore::State_::I_ve ) );

    insert( names::Iadap_ini, host.get_data_access_functor( migliore::State_::Iadap ) );
    insert( names::Idep_ini, host.get_data_access_functor( migliore::State_::Idep ) );

    host.insert_current_recordables();
  }
  // void
  // RecordablesMap< migliore >::create()
  // {
  //   // use standard names whereever you can for consistency!
  //   insert_( names::V_m, &migliore::get_V_m_ );
  //   insert_( names::I, &migliore::get_I_inj_ );
  //   insert_( names::init_sign, &migliore::get_init_sign_ );
  //   insert_( names::Iadap_ini, &migliore::get_Iadap_ini_ );
  //   insert_( names::Idep_ini, &migliore::get_Idep_ini_ );
  //   insert_( names::sis, &migliore::get_sis_ );
  //   insert_( names::I_syn, &migliore::get_I_syn_ );
  // }


  Name
  migliore::get_i_syn_name( size_t elem )
  {
    std::stringstream i_syn_name;
    i_syn_name << "I_syn_" << elem + 1;
    return Name( i_syn_name.str() );
  }

  void
  migliore::insert_current_recordables( size_t first )
  {
    for ( size_t receptor = first; receptor < P_.tau_syn_fast_decay_.size(); ++receptor )
      {
	size_t elem = migliore::State_::I_SYN
	  + receptor * migliore::State_::NUM_STATE_ELEMENTS_PER_RECEPTOR;
	Name a = get_i_syn_name( receptor );
	recordablesMap_.insert( a , this->get_data_access_functor( elem ) );
      }
  }

  DataAccessFunctor< migliore >
  migliore::get_data_access_functor( size_t elem )
  {
    return DataAccessFunctor< migliore >( *this, elem );
  }

  /* ----------------------------------------------------------------
   * Default constructors defining default parameters and state
   * ---------------------------------------------------------------- */

  migliore::Parameters_::Parameters_()
    : V_min_( -90.0 ) //-std::numeric_limits< double >::max() ) // mV
    , E_L_( -72.5 )
    , Vres_( -65.0 )
    , V_th_( -52.0 )// mV vtm in pyhon
    , Cm_( 2047.4164432004916 )
    , I_th_( 300.9987901902274 )
    , tao_m_( 3310.462136574088 )
    , sc_( 4526.328798037026 )
    , bet_( 0.24522251335200956 )
    , delta1_( 0.009906254244852036 )
    , cost_idep_ini_( 0.017625482908326662 )
    , Idep_ini_vr_( 1.0122905259090516 )
    , psi1_( 0.1975362978159442 )
    , a_( 14.2787 )
    , b_( -2.10966 )
    , c_( 0.0608809 )
    , alp_( 225.491 )
    , Iadap_start_ ( 0.0 )
    , istim_min_spikinig_exp_( 400 )
    , istim_max_spikinig_exp_( 1000 )
    , I_e_( 0.0 ) // pA
    , t_ref_( 2.0 ) // in ms
    , tau_syn_fast_rise_ ( 0.1, 0.1 ) // in ms
    , tau_syn_fast_decay_ ( 3.0, 10.0 ) // in ms
    , tau_syn_slow_rise_ ( 0.1, 0.1 ) // in ms
    , tau_syn_slow_decay_ ( 3.0, 10.0 ) // in ms
    , has_connections_( false )
    , NMDA_ratio_( 1.22 ) // 2 mM in the Johnston et al. 2010, extracellula [MgCl2] = 1 mM in Edelman et al. 2015
    , mg_( 2.0) // 2 mM in the Johnston et al. 2010, extracellula [MgCl2] = 1 mM in Edelman et al. 2015
    , mgb_k_ ( 0.062 ) // (/mV) Johnston et al. 2010
    , mg_ref_ ( 3.57 ) // (mM) Johnston et al. 2010
    , mgb_shift_ ( -10.0 ) // (mM) Johnston et al. 2010
    , vinc_inf_ ( 700.0 ) // Default value for the NEURON model
    , vinc_sup_ ( std::numeric_limits<double>::infinity() ) // Default value for the NEURON model
    , coeffInf_ ( 0.68 ) // Default value for the NEURON model
    , constInf_ ( -190.0 ) // Default value for the NEURON model
    , coeffSup_ ( -0.028 )  // Default value for the NEURON model
    , constSup_ ( 76.86 ) // Default value for the NEURON model
    , mincurr_  ( -185.0 ) // Default value for the NEURON model
    , aglif_p_ ( 1.0 )
    , zeta_ ( 3.5e-3 )     // Neuron equilibrium params    
    , eta_ ( 2.5e-3 )
    , rho_ ( 1e-3 )
    , csi_ ( 3.5e-3 )
    , plotit_ ( false )
  {
  }

  migliore::State_::State_()
    : V_m_( -72.5 )       // membrane potential  (E_L_+(1-exp(-cor[i]/1000))*(V_th_-E_L_))/V_.Vconvfact  TOBEFIXED
    , init_sign_( 0.0 ) // membrane adaptation variable
    , Iadap_ini_( 0.0 ) // membrane adaptation variable
    , Idep_ini_( 0.0 ) // membrane dependent variable
    , sis_( 0.0 ) // membrane dependent variable
    , I_inj_( 0.0 )         // input current
    , r_ref_( 0.0 )
    , current_ ( 0.0 )
  {
    i_syn_.clear();
    i_syn_fast_rise_A_.clear();
    i_syn_fast_decay_B_.clear();
    i_syn_slow_rise_A_.clear();
    i_syn_slow_decay_B_.clear();
  }

  /* ----------------------------------------------------------------
   * Parameter and state extractions and manipulation functions
   * ---------------------------------------------------------------- */

  void
  migliore::Parameters_::get( DictionaryDatum& d ) const
  {
    def< double >( d, names::I_e, I_e_ ); // External current injection
    def< double >( d, names::V_th, V_th_ ); // threshold value
    def< double >( d, names::V_min, V_min_ );
    def< double >( d, names::E_L, E_L_ );
    def< double >( d, names::Vres, Vres_ );
    //    def< double >( d, names::vtm, vtm_);
    def< double >( d, names::Cm, Cm_);
    def< double >( d, names::I_th, I_th_);
    def< double >( d, names::tao_m, tao_m_);
    def< double >( d, names::sc, sc_);
    def< double >( d, names::bet, bet_);
    def< double >( d, names::delta1, delta1_);
    def< double >( d, names::cost_idep_ini, cost_idep_ini_);
    def< double >( d, names::Idep_ini_vr, Idep_ini_vr_);
    def< double >( d, names::psi1, psi1_);  // To be removed as it is calculated as Variable.
    def< double >( d, names::a, a_);
    def< double >( d, names::b, b_);
    def< double >( d, names::c, c_);
    def< double >( d, names::alp, alp_);
    def< double >( d, names::Iadap_start, Iadap_start_);
    def< double >( d, names::istim_min_spikinig_exp, istim_min_spikinig_exp_);
    def< double >( d, names::istim_max_spikinig_exp, istim_max_spikinig_exp_);
    def< double >( d, names::NMDA_ratio, NMDA_ratio_);
    def< double >( d, names::t_ref, t_ref_ );
    def< double >( d, names::mg, mg_ );
    def< double >( d, names::mg_ref, mg_ref_ );
    def< double >( d, names::mgb_k, mgb_k_ );
    def< double >( d, names::mgb_shift, mgb_shift_ );

    def< double >( d, names::coeffInf, coeffInf_ );
    def< double >( d, names::constInf, constInf_ );
    def< double >( d, names::coeffSup, coeffSup_ );
    def< double >( d, names::constSup, constSup_ );
    def< double >( d, names::aglif_p, aglif_p_ );
    def< double >( d, names::vinc_inf, vinc_inf_ );
    def< double >( d, names::vinc_sup, vinc_sup_ );
    def< double >( d, names::mincurr, mincurr_ );
    def< double >( d, names::zeta, zeta_ );     // Neuron equilibrium params    
    def< double >( d, names::eta,  eta_ );
    def< double >( d, names::rho, rho_ );
    def< double >( d, names::csi, csi_ );
    def< bool >( d, names::plotit, plotit_ );

    def< int >( d, names::n_synapses, n_receptors_() );
    def< bool >( d, names::has_connections, has_connections_ );

    ArrayDatum tau_syn_fr_ad( tau_syn_fast_rise_ );
    ArrayDatum tau_syn_fd_ad( tau_syn_fast_decay_ );
    ArrayDatum tau_syn_sr_ad( tau_syn_slow_rise_ );
    ArrayDatum tau_syn_sd_ad( tau_syn_slow_decay_ );
    def< ArrayDatum >( d, names::tau_syn_fast_rise, tau_syn_fr_ad );
    def< ArrayDatum >( d, names::tau_syn_fast_decay, tau_syn_fd_ad );
    def< ArrayDatum >( d, names::tau_syn_slow_rise, tau_syn_sr_ad );
    def< ArrayDatum >( d, names::tau_syn_slow_decay, tau_syn_sd_ad );
  }

  void
  migliore::Parameters_::set( const DictionaryDatum& d, Node* node )
  {
    updateValueParam< double >( d, names::V_th, V_th_, node );
    updateValueParam< double >( d, names::V_min, V_min_, node );
    updateValueParam< double >( d, names::I_e, I_e_, node ); // External current injection
    updateValueParam< double >( d, names::E_L, E_L_, node );
    updateValueParam< double >( d, names::Vres, Vres_, node );
    //    updateValueParam< double >( d, names::vtm, vtm_, node );
    updateValueParam< double >( d, names::Cm, Cm_, node );
    updateValueParam< double >( d, names::I_th, I_th_, node );
    updateValueParam< double >( d, names::tao_m, tao_m_, node );
    updateValueParam< double >( d, names::sc, sc_, node );
    updateValueParam< double >( d, names::bet, bet_, node );
    updateValueParam< double >( d, names::delta1, delta1_, node );
    updateValueParam< double >( d, names::cost_idep_ini, cost_idep_ini_, node );
    updateValueParam< double >( d, names::Idep_ini_vr, Idep_ini_vr_, node );
    updateValueParam< double >( d, names::psi1, psi1_, node );
    updateValueParam< double >( d, names::a, a_, node );
    updateValueParam< double >( d, names::b, b_, node );
    updateValueParam< double >( d, names::c, c_, node );
    updateValueParam< double >( d, names::alp, alp_, node );
    updateValueParam< double >( d, names::Iadap_start, Iadap_start_, node );
    updateValueParam< double >( d, names::istim_min_spikinig_exp, istim_min_spikinig_exp_, node );
    updateValueParam< double >( d, names::istim_max_spikinig_exp, istim_max_spikinig_exp_, node );
    updateValueParam< double >( d, names::NMDA_ratio, NMDA_ratio_, node );
    updateValueParam< double >( d, names::t_ref, t_ref_, node );
    updateValueParam< double >( d, names::mg, mg_, node );
    updateValueParam< double >( d, names::mg_ref, mg_ref_, node );
    updateValueParam< double >( d, names::mgb_k, mgb_k_, node );
    updateValueParam< double >( d, names::mgb_shift, mgb_shift_, node );

    updateValueParam< double >( d, names::coeffInf, coeffInf_, node );
    updateValueParam< double >( d, names::constInf, constInf_, node );
    updateValueParam< double >( d, names::coeffSup, coeffSup_, node );
    updateValueParam< double >( d, names::constSup, constSup_, node );
    updateValueParam< double >( d, names::aglif_p, aglif_p_, node );
    updateValueParam< double >( d, names::vinc_inf, vinc_inf_, node );
    updateValueParam< double >( d, names::vinc_sup, vinc_sup_, node );
    updateValueParam< double >( d, names::mincurr, mincurr_, node );
    updateValueParam< double >( d, names::zeta, zeta_, node );
    updateValueParam< double >( d, names::eta, eta_, node );
    updateValueParam< double >( d, names::rho, rho_, node );
    updateValueParam< double >( d, names::csi, csi_, node );
    updateValueParam< bool >( d, names::plotit, plotit_, node );

    if ( t_ref_ < 0 )
      {
    	throw BadProperty( "Refractory time must not be negative." );
      }

    const size_t old_n_receptors = this->n_receptors_();
    if ( updateValue< std::vector< double > >( d, "tau_syn_fast_rise", tau_syn_fast_rise_ ) )
      {
    	if ( this->n_receptors_() != old_n_receptors && has_connections_ == true )
	  {
	    throw BadProperty( "The neuron has connections, therefore the number of ports cannot be reduced." );
	  }
	for ( size_t i = 0; i < tau_syn_fast_rise_.size(); ++i )
	  {
	    if ( tau_syn_fast_rise_[ i ] <= 0 )
	      {
	    	throw BadProperty( "All synaptic time constants must be strictly positive." );
	      }
	    if ( tau_syn_fast_rise_[ i ] == tao_m_ )
	      {
	    	throw BadProperty( "Membrane and synapse time constant(s) must differ. See note in documentation." );
	      }
	  }
      }
    else
      {
    	if ( this->n_receptors_() != old_n_receptors )
	  {
	    throw BadProperty( "All synaptic time constants must defined. Set time constants for tau_syn_fast_rise." );
	  }
      }
    if ( updateValue< std::vector< double > >( d, "tau_syn_fast_decay", tau_syn_fast_decay_ ) )
      {
    	if ( this->n_receptors_() != old_n_receptors && has_connections_ == true )
	  {
	    throw BadProperty( "The neuron has connections, therefore the number of ports cannot be reduced." );
	  }
	for ( size_t i = 0; i < tau_syn_fast_decay_.size(); ++i )
	  {
	    if ( tau_syn_fast_decay_[ i ] <= 0 )
	      {
	    	throw BadProperty( "All synaptic time constants must be strictly positive." );
	      }
	    if ( tau_syn_fast_decay_[ i ] == tao_m_ )
	      {
	    	throw BadProperty( "Membrane and synapse time constant(s) must differ. See note in documentation." );
	      }
	  }
      }
    else
      {
    	if ( this->n_receptors_() != old_n_receptors )
	  {
	    throw BadProperty( "All synaptic time constants must defined. Set time constants for tau_syn_fast_decay." );
	  }
      }
    if ( updateValue< std::vector< double > >( d, "tau_syn_slow_rise", tau_syn_slow_rise_ ) )
      {
    	if ( this->n_receptors_() != old_n_receptors && has_connections_ == true )
	  {
	    throw BadProperty( "The neuron has connections, therefore the number of ports cannot be reduced." );
	  }
	for ( size_t i = 0; i < tau_syn_slow_rise_.size(); ++i )
	  {
	    if ( tau_syn_slow_rise_[ i ] <= 0 )
	      {
	    	throw BadProperty( "All synaptic time constants must be strictly positive." );
	      }
	    if ( tau_syn_slow_rise_[ i ] == tao_m_ )
	      {
	    	throw BadProperty( "Membrane and synapse time constant(s) must differ. See note in documentation." );
	      }
	  }
      }
    else
      {
    	if ( this->n_receptors_() != old_n_receptors )
	  {
	    throw BadProperty( "All synaptic time constants must defined. Set time constants for tau_syn_slow_rise." );
	  }
      }
    if ( updateValue< std::vector< double > >( d, "tau_syn_slow_decay", tau_syn_slow_decay_ ) )
      {
    	if ( this->n_receptors_() != old_n_receptors && has_connections_ == true )
	  {
	    throw BadProperty( "The neuron has connections, therefore the number of ports cannot be reduced." );
	  }
	for ( size_t i = 0; i < tau_syn_slow_decay_.size(); ++i )
	  {
	    if ( tau_syn_slow_decay_[ i ] <= 0 )
	      {
	    	throw BadProperty( "All synaptic time constants must be strictly positive." );
	      }
	    if ( tau_syn_slow_decay_[ i ] == tao_m_ )
	      {
	    	throw BadProperty( "Membrane and synapse time constant(s) must differ. See note in documentation." );
	      }
	  }
      }
    else
      {
    	if ( this->n_receptors_() != old_n_receptors )
	  {
	    throw BadProperty( "All synaptic time constants must defined. Set time constants for tau_syn_slow_decay." );
	  }
      }
  }

  void
  migliore::State_::get( DictionaryDatum& d, const Parameters_& ) const
  {
    def< double >( d, names::init_sign, init_sign_ ); // Membrane potential adaptationariable
    def< double >( d, names::Iadap_ini, Iadap_ini_ ); // Membrane potential adaptationariable
    def< double >( d, names::Idep_ini , Idep_ini_ ); // Membrane potential adaptationariable
    def< double >( d, names::sis , sis_ ); // Membrane potential adaptationariable
    def< double >( d, names::I, I_inj_ ); // Membrane potential
    def< double >( d, names::V_m, V_m_ ); // Membrane potential
  }

  void
  migliore::State_::set( const DictionaryDatum& d, const Parameters_&, Node* node )
  {
    updateValueParam< double >( d, names::init_sign, init_sign_, node );
    updateValueParam< double >( d, names::Iadap_ini, Iadap_ini_, node );
    updateValueParam< double >( d, names::Idep_ini, Idep_ini_, node );
    updateValueParam< double >( d, names::sis, sis_, node );
    updateValueParam< double >( d, names::V_m, V_m_, node );
    updateValueParam< double >( d, names::I, I_inj_, node );
  }

  migliore::Buffers_::Buffers_( migliore& n )
    : logger_( n )
  {
  }

  migliore::Buffers_::Buffers_( const Buffers_&, migliore& n )
    : logger_( n )
  {
  }

  /* ----------------------------------------------------------------
   * Default and copy constructor for node
   * ---------------------------------------------------------------- */

  migliore::migliore()
    : ArchivingNode()
    , P_()
    , S_()
    , B_( *this )
  {
    recordablesMap_.create( *this );
  }

  migliore::migliore( const migliore& n )
    : ArchivingNode( n )
    , P_( n.P_ )
    , S_( n.S_ )
    , B_( n.B_, *this )
  {
    recordablesMap_.create( *this );
  }

  /* ----------------------------------------------------------------
   * Node initialization functions
   * ---------------------------------------------------------------- */

  void
  migliore::init_buffers_()
  {
    B_.spikes_.clear();   // includes resize
    B_.currents_.clear(); // includes resize
    B_.logger_.reset();   // includes resize
    ArchivingNode::clear_history();
  }

  void
  nest::migliore::pre_run_hook()
  {
    B_.logger_.init();

    // t_ref_ specifies the length of the absolute refractory period as
    // a double in ms. The grid based iaf_psc_exp can only handle refractory
    // periods that are integer multiples of the computation step size (h).
    // To ensure consistency with the overall simulation scheme such conversion
    // should be carried out via objects of class nest::Time. The conversion
    // requires 2 steps:
    //     1. A time object r is constructed, defining representation of
    //        t_ref_ in tics. This representation is then converted to computation
    //        time steps again by a strategy defined by class nest::Time.
    //     2. The refractory time in units of steps is read out get_steps(), a
    //        member function of class nest::Time.
    //
    // Choosing a t_ref_ that is not an integer multiple of the computation time
    // step h will lead to accurate (up to the resolution h) and self-consistent
    // results. However, a neuron model capable of operating with real valued
    // spike time may exhibit a different effective refractory time.

    const double h = Time::get_resolution().get_ms();
    
    V_.blockActive = false;
    
    V_.P11_syn_fast_rise_.resize( P_.n_receptors_() );
    V_.P11_syn_fast_decay_.resize( P_.n_receptors_() );
    V_.fast_tp_.resize( P_.n_receptors_() );
    V_.syn_fast_factor_.resize( P_.n_receptors_() );

    V_.P11_syn_slow_rise_.resize( P_.n_receptors_() );
    V_.P11_syn_slow_decay_.resize( P_.n_receptors_() );
    V_.slow_tp_.resize( P_.n_receptors_() );
    V_.syn_slow_factor_.resize( P_.n_receptors_() );
    
    S_.i_syn_.resize( P_.n_receptors_() );
    S_.i_syn_fast_rise_A_.resize( P_.n_receptors_() );
    S_.i_syn_fast_decay_B_.resize( P_.n_receptors_() );
    S_.i_syn_slow_rise_A_.resize( P_.n_receptors_() );
    S_.i_syn_slow_decay_B_.resize( P_.n_receptors_() );

    S_.V_m_ = P_.E_L_;

    B_.spikes_.resize( P_.n_receptors_() );

    // Initialize states
    S_.Iadap_ini_ = P_.Iadap_start_;
    
    // V_.P22_ = std::exp( -h / P_.tao_m_ );
    // V_.P20_ = P_.tao_m_ / P_.Cm_ * ( 1.0 - V_.P22_ );

    for ( size_t i = 0; i < P_.n_receptors_(); i++ )
      {
	V_.P11_syn_fast_rise_[ i ] = std::exp( -h / P_.tau_syn_fast_rise_[ i ] );
	V_.P11_syn_fast_decay_[ i ] = std::exp( -h / P_.tau_syn_fast_decay_[ i ] );
	if (i == 0)
	  {
	    V_.P11_syn_slow_rise_[ i ] = std::exp( -h / P_.tau_syn_slow_rise_[ i ] );
	    V_.P11_syn_slow_decay_[ i ] = std::exp( -h / P_.tau_syn_slow_decay_[ i ] );
	  }
	else
	  {
	    V_.P11_syn_slow_rise_[ i ] = 0;
	    V_.P11_syn_slow_decay_[ i ] = 0;
	  }
	// these are determined according to a numeric stability criterion
	// V_.P21_syn_[ i ] = propagator_32( P_.tau_syn_fast_decay_[ i ], P_.tao_m_, P_.Cm_, h );
	B_.spikes_[ i ].resize();
	assert( P_.tau_syn_fast_rise_[ i ] < P_.tau_syn_fast_decay_[ i ] );
	assert( P_.tau_syn_slow_rise_[ i ] < P_.tau_syn_slow_decay_[ i ] );

	if (P_.tau_syn_fast_rise_[ i ] / P_.tau_syn_fast_decay_[ i ] > 0.9999) {
	  P_.tau_syn_fast_rise_[ i ] = 0.9999 * P_.tau_syn_fast_rise_[ i ];
	}
	if (P_.tau_syn_fast_rise_[ i ] / P_.tau_syn_fast_decay_[ i ] < 1e-9) {
	  P_.tau_syn_fast_rise_[ i ] = P_.tau_syn_fast_decay_[ i ] * 1e-9;
	}
    
	V_.fast_tp_[ i ] = (P_.tau_syn_fast_rise_[ i ] * P_.tau_syn_fast_decay_[ i ]) /
	  (P_.tau_syn_fast_decay_[ i ] - P_.tau_syn_fast_rise_[ i ]) *
	  log(P_.tau_syn_fast_decay_[ i ] / P_.tau_syn_fast_rise_[ i ]);
	V_.syn_fast_factor_[ i ] = -exp(-V_.fast_tp_[ i ] / P_.tau_syn_fast_rise_[ i ]) +
	  exp(-V_.fast_tp_[ i ] / P_.tau_syn_fast_decay_[ i ]);
	V_.syn_fast_factor_[ i ] = 1 / V_.syn_fast_factor_[ i ];
	
	V_.slow_tp_[ i ] = (P_.tau_syn_slow_rise_[ i ] * P_.tau_syn_slow_decay_[ i ]) /
	  (P_.tau_syn_slow_decay_[ i ] - P_.tau_syn_slow_rise_[ i ]) *
	  log(P_.tau_syn_slow_decay_[ i ] / P_.tau_syn_slow_rise_[ i ]);
	V_.syn_slow_factor_[ i ] = -exp(-V_.slow_tp_[ i ] / P_.tau_syn_slow_rise_[ i ]) +
	  exp(-V_.slow_tp_[ i ] / P_.tau_syn_slow_decay_[ i ]);
	V_.syn_slow_factor_[ i ] = 1 / V_.syn_slow_factor_[ i ];
	
      }
    V_.RefractoryCounts_ = Time( Time::ms( P_.t_ref_ ) ).get_steps();
    // since t_ref_ >= 0, this can only fail in error
    assert( V_.RefractoryCounts_ >= 0 );

    V_.time_scale_ = 1 / (-P_.sc_ / (P_.Cm_ * P_.E_L_));
    V_.d_dt = h;
    V_.dt = V_.d_dt/V_.time_scale_;
    V_.t_spk = -3 * V_.d_dt;

    V_.beta2 = pow(P_.bet_,2);
    V_.psi1 = pow(-4.0 * P_.bet_ + (pow(1+P_.delta1_,2.0)),0.5);
    V_.t_step = V_.dt;

    V_.I_th_ = P_.I_th_ * P_.aglif_p_;
    V_.Cm_ = P_.Cm_ * P_.aglif_p_;
    V_.sc_ = P_.sc_ * P_.aglif_p_;
      
    V_.V_star_min_ = -P_.V_min_ / P_.E_L_;
    V_.alpha_neg_ = P_.mincurr_ / V_.sc_;

    
    // V_.H = (90+P_.E_L_)*V_.sc_*(P_.bet_-P_.delta1_)/(P_.E_L_*(P_.mincurr_));
    V_.Vconvfact = -P_.E_L_;
    V_.vrm = P_.Vres_/V_.Vconvfact;
    // V_.psi1 = pow((-4)*P_.bet_+pow(1+P_.delta1_,2.0),0.5); // This should substitute P_.psi1_
    // // t_step = t - t0
    // V_.VV_1 = 0.5 / ((P_.bet_ -P_.delta1_) * (pow(P_.bet_,2.0) + (P_.bet_-1.0) * P_.delta1_) * (4.0 * P_.bet_ - pow((1.0 + P_.delta1_), 2.0))) * P_.psi1_;
    // V_.VV_2 = 2.0 * exp(-V_.t_step * P_.bet_) * (P_.bet_-1.0) * P_.bet_ * (P_.bet_ - P_.delta1_) * P_.psi1_;
    // V_.VV_3 = (pow(P_.bet_,2) + ((- 1.0) + P_.bet_) * P_.delta1_) * P_.psi1_;
    // V_.VV_4 = exp((1.0 / 2.0) * V_.t_step * (-1.0 + P_.delta1_ -P_.psi1_));
    // V_.VV_5 = P_.bet_ * (P_.bet_ -P_.delta1_) * (-1.0 -P_.delta1_ + P_.bet_ * (3.0 + P_.delta1_ -P_.psi1_) + P_.psi1_);
    // V_.VV_6 = (pow(P_.bet_,2) -P_.delta1_ + P_.bet_ * P_.delta1_);
    // V_.VV_7 = (1.0 + (-2.0) * P_.bet_ + P_.delta1_ -P_.psi1_);
    // V_.VV_8 = exp((1.0 / 2.0) * V_.t_step * (-1.0 + P_.delta1_ + P_.psi1_));
    // V_.VV_9 = P_.bet_ * (P_.bet_-P_.delta1_) * (-1.0 -P_.delta1_ -P_.psi1_ + P_.bet_ * (3.0 + P_.delta1_ + P_.psi1_));
    // V_.VV_10 = (pow(P_.bet_,2) - P_.delta1_ + P_.bet_ * P_.delta1_);
    // V_.VV_11 = (1.0 + (-2.0) * P_.bet_ + P_.delta1_ + P_.psi1_);
    // V_.C = P_.bet_ - P_.delta1_;

    // V_.pdelta = pow((1 + P_.delta1_), 2);

    // V_.t_step = t - t0;
    double beta = P_.bet_;
    double delta = P_.delta1_;
    double Psi = P_.psi1_;
    V_.VV_1 = 0.5 / (beta - delta) / (pow(beta,2.0) + (beta - 1.0) * delta) / (4.0 * beta - pow((1.0 + delta),2.0)) * Psi;
    V_.VV_2 = (beta - 1.0) * beta * (beta - delta) * Psi;
    V_.VV_3 = (pow(beta,2.0) + (beta - 1.0) * delta) * Psi;
    V_.VV_4 = 0.5 * (-1.0 + delta - Psi);
    V_.VV_5 = beta * (beta - delta) * (-1.0 - delta + beta * (3.0 + delta - Psi) + Psi);
    V_.VV_6 = (pow(beta,2.0) - delta + beta * delta);
    V_.VV_7 = (1.0 -2.0 * beta + delta - Psi);
    V_.VV_8 = beta * (beta - delta) * (-1.0 - delta - Psi + beta * (3.0 + delta + Psi));
    V_.VV_9 = (pow(beta,2.0) - delta + beta * delta);
    V_.VV_10 = (1.0 - 2.0 * beta + delta + Psi);


    V_.AA_1 = -4.0 * pow(P_.bet_,3.0) + pow(P_.bet_,2.0) * pow((-1.0+P_.delta1_),2.0) - P_.delta1_ * pow((1.0 + P_.delta1_),2.0) + P_.bet_ * P_.delta1_ * (5 + 2.0 * P_.delta1_ + pow(P_.delta1_,2.0));
    V_.AA_2 = 2.0 * exp(-V_.t_step * P_.bet_) * P_.bet_ * (4.0 * pow(P_.bet_,2.0) + P_.delta1_ * pow((1.0+P_.delta1_),2.0) - P_.bet_ * (1.0 + 6 * P_.delta1_ + pow(P_.delta1_,2.0)));
    V_.AA_3 = exp((1.0 / 2.0)*V_.t_step*(-1.0+P_.delta1_+P_.psi1_));
    V_.AA_4 = -1.0-2.0*P_.delta1_ - pow(P_.delta1_,2.0) - P_.psi1_ + P_.delta1_ * P_.psi1_+2.0*P_.bet_*(2.0+P_.psi1_);
    V_.AA_5 = (pow(P_.bet_,2.0)-P_.delta1_+P_.bet_*P_.delta1_);
    V_.AA_6 = (1.0-4.0*P_.bet_+2.0*P_.delta1_+pow(P_.delta1_,2.0)+P_.psi1_-P_.delta1_*P_.psi1_);
    V_.AA_7 = (1.0+P_.delta1_)*(-1.0-P_.delta1_+P_.psi1_);
    V_.AA_8 = exp((-1.0)*(1.0 / 2.0) * V_.t_step * (1.0-P_.delta1_+P_.psi1_));
    V_.AA_9 = P_.bet_ * (P_.bet_-P_.delta1_)*(1.0+2.0*P_.delta1_+ pow(P_.delta1_,2.0)-P_.psi1_+P_.delta1_*P_.psi1_+2.0*P_.bet_*(-2.0+P_.psi1_));
    V_.AA_10 = (pow(P_.bet_,2.0)-P_.delta1_+P_.bet_*P_.delta1_);
    V_.AA_11 = (1.0-4.0*P_.bet_+2.0*P_.delta1_ + pow(P_.delta1_,2.0) - P_.psi1_+P_.delta1_*P_.psi1_);
    V_.AA_12 = (2.0*(P_.bet_-P_.delta1_)*(pow(P_.bet_,2.0)+(-1.0+P_.bet_)*P_.delta1_)*(4.0*P_.bet_-pow((1.0+P_.delta1_),2.0)));

    V_.DD_1 = exp(-V_.t_step * P_.bet_);
  }

  /* ----------------------------------------------------------------
   * Update and spike handling functions
   */

  double
  migliore::tagliorette(double corr)
  {
    double dur_sign = std::numeric_limits<double>::infinity();
    if (corr < P_.vinc_inf_ && corr >= 0) // ValSupLineInf
      {
    	dur_sign = P_.coeffInf_ * corr + P_.constInf_;
      }
    if (corr > P_.vinc_sup_) // THIS CANNOT HAPPEN!! ValInfLineSup
      {
    	dur_sign = P_.coeffSup_ * corr + P_.constSup_;
      }
    return dur_sign;
  }

// # block lines
// def tagliorette(corr,retteParParsed):
//     vinc_sup = retteParParsed[0]
//     coeffSup = retteParParsed[1]
//     constSup = retteParParsed[2]
//     vinc_inf = retteParParsed[3]
//     coeffInf = retteParParsed[4]
//     constInf = retteParParsed[5]
    
//     dur_sign=np.inf

//     if corr<vinc_inf and corr>0:
//         dur_sign = coeffInf*corr + constInf
    
//     if corr>vinc_sup:
//         dur_sign = coeffSup*corr + constSup
//     return dur_sign

       // retteOrig = [float(x) for x in list(pyramidalNeuronsDatabase[nomeNeurone]['block_line_params'].values())]
       //                  rette = [0,0,0,0,0,0]
       //                  for it in range(6):
       //                      if retteOrig[it]>1000000:                                
       //                          rette[it] = np.inf
       //                      elif retteOrig[it]<-1000000:                                
       //                          rette[it] = -np.inf
  
  double
  migliore::migliV(double t, double delta, double Psi, 
		   double alpha, double beta, double IaA0, 
		   double IdA0, double t0, double V0,
		   int r_ref_, double vrm)
  {
    // double to_return = 0.5 / (beta -delta) / (pow(beta,2.0) + (-1.0 + beta) * delta) / (4 * beta + (- 1.0) * pow((1.0 + delta),2.0)) * Psi *      (2.0 * exp((-t + t0) * beta) * IdA0 * (-1.0 + beta) * beta * (beta -delta) * Psi + (-2.0) * (alpha + (-1.0) * beta + delta) * (pow(beta,2.0) + ((-1.0) + beta) * delta) * Psi + exp(0.5 * (t - t0) * (-1.0 + delta - Psi)) * (IdA0 * beta * (beta -delta) * (-1.0 -delta + beta * (3.0 + delta - Psi) + Psi) - (pow(beta,2.0) -delta + beta * delta) * (alpha * (1.0 + (-2.0) * beta + delta - Psi) + (beta -delta) * ((-1.0) + 2.0 * IaA0 * beta -delta + Psi + V0 * (-1.0 -delta + Psi)))) + exp(0.5 * (t + (-1.0) * t0) * (-1.0 + delta + Psi)) * (- IdA0 * beta * (beta - delta) * (-1.0 -delta - Psi + beta * (3.0 + delta + Psi)) + (pow(beta,2.0) -delta + beta * delta) * (alpha * (1.0 + (-2.0) * beta + delta + Psi) + (beta -delta) * (-1.0 + 2.0 * IaA0 * beta - delta - Psi - V0 * (1.0 + delta + Psi)))));

    double to_return = 0.5 / (beta - delta) / (pow(beta,2.0) + (-1.0 + beta) * delta) / (4.0 * beta - pow((1.0 + delta),2.0)) * Psi * (2.0 * exp((t0 - t ) * beta) * IdA0 * (-1.0 + beta) * beta * (beta - delta) * Psi -2.0 * (alpha - beta + delta) * (pow(beta,2.0) + (- 1.0 + beta) * delta) * Psi + exp(0.5 * (t - t0) * (-1.0 + delta - Psi)) * (IdA0 * beta * (beta - delta) * (-1.0 - delta + beta * (3.0 + delta - Psi) + Psi) - (pow(beta,2.0) - delta + beta * delta) * (alpha * (1.0 -2.0 * beta + delta - Psi) + (beta - delta) * (-1.0 + 2.0 * IaA0 * beta - delta + Psi + V0 * (-1.0 - delta + Psi)))) + exp(0.5 * (t - t0) * (-1.0 + delta + Psi)) * (-1.0 * IdA0 * beta * (beta - delta) * (-1.0 - delta - Psi + beta * (3.0 + delta + Psi)) + (pow(beta,2.0) - delta + beta * delta) * (alpha * (1.0 -2.0 * beta + delta + Psi) + (beta - delta) * (-1.0 + 2.0 * IaA0 * beta - delta - Psi - V0 * (1.0 + delta + Psi)))));
    
    //    double to_return = V_.VV_1 * (2.0 * exp((t0 - t) * beta) * IdA0 * V_.VV_2 - 2.0 * (alpha - beta + delta) * V_.VV_3 + exp(V_.VV_4 * (t - t0)) * (IdA0 * V_.VV_5 - V_.VV_6 * (alpha * V_.VV_7 + (beta - delta) * (-1.0 + 2.0 * IaA0 * beta - delta + Psi + V0 * (-1.0 - delta + Psi)))) + exp(0.5 * (t - t0) * (-1.0 + delta + Psi)) * (- IdA0 * V_.VV_8 + V_.VV_9 * (alpha * V_.VV_10 + (beta - delta) * (-1.0 + 2.0 * IaA0 * beta - delta - Psi - V0 * (1.0 + delta + Psi)))));

    // double to_return = 0.5 / (beta - delta) / (pow(beta,2.0) + (beta - 1.0) * delta) / (4 * beta - pow((1.0 + delta),2.0)) *
    // 		Psi * (2.0 * exp((t0 - t) * beta) * IdA0 * (beta - 1.0) * beta * (beta - delta) * Psi -2.0 * (alpha - beta + delta) *
    // 				(pow(beta,2.0) + (beta - 1.0) * delta) * Psi + exp(0.5 * (t - t0) * (-1.0 + delta - Psi)) *
    // 					(IdA0 * beta * (beta - delta) * (-1.0 -delta + beta * (3.0 + delta -Psi) + Psi) - (pow(beta,2.0)
    // 							-delta + beta * delta) * (alpha * (1.0 -2.0 * beta + delta -Psi) + (beta -delta) *
    // 									(-1.0 + 2.0 * IaA0 * beta -delta + Psi + V0 * (-1.0 -delta + Psi)))) +
    // 									exp(0.5 * (t -t0) * (-1.0 + delta + Psi)) *
    // 									(-IdA0 * beta * (beta-delta) * (-1.0 -delta -Psi + beta * (3.0 + delta + Psi)) +
    // 											(pow(beta,2.0) -delta +beta * delta) * (alpha * (1.0 -2.0 * beta + delta + Psi) +
    // 													(beta -delta) * (-1.0 + 2 * IaA0 * beta - delta - Psi - V0 * (1 + delta + Psi)))));

    // double to_return = V_.VV_1 * (2.0 * exp(-(t-t0) * beta) * IdA0 * V_.VV_3 -2.0 * (alpha - beta + delta) * V_.VV_4  + exp(0.5 * (t-t0) * (-1.0 + delta - Psi)) * (IdA0 * V_.VV_6 - V_.VV_7 * (alpha * V_.VV_8 + (beta - delta) * (-1.0 + 2.0 * IaA0 * beta + V_.VV_9 + V0 * V_.VV_9b))) + exp(0.5 * (t-t0) * (-1.0 + delta + Psi)) * (-IdA0 * V_.VV_11 + V_.VV_12 * (alpha * V_.VV_13 + V_.VV_14 * (-1.0 + 2 * IaA0 * beta - delta - Psi - V0 * (1 + delta + Psi)))));
    // double to_return = V_.VV_1 * (V_.VV_2 * IdA0 * V_.VV_3 -2.0 * (alpha - beta + delta) * V_.VV_4  + V_.VV_5 * (IdA0 * V_.VV_6 - V_.VV_7 * (alpha * V_.VV_8 + (beta - delta) * (-1.0 + 2.0 * IaA0 * beta + +  V_.VV_9 + V0 * V_.VV_9b))) + V_.VV_10 * (-IdA0 * V_.VV_11 + V_.VV_12 * (alpha * V_.VV_13 + V_.VV_14 * (-1.0 + 2 * IaA0 * beta - delta - Psi - V0 * (1 + delta + Psi)))));

    //double to_return = V_.VV_1 * (V_.VV_2 * IdA0 + -2.0 * (alpha -beta + delta) * V_.VV_3 + V_.VV_4 * (IdA0 * V_.VV_5 - V_.VV_6 * (alpha * V_.VV_7 + (beta -delta) * (-1.0 + 2.0 * IaA0 * beta -delta + Psi + V0 * (-1.0 -delta + Psi)))) + V_.VV_8 * (-IdA0 * V_.VV_9 + V_.VV_10 * (alpha * V_.VV_11 + (beta -delta) * (-1.0 + 2.0 * IaA0 * beta-delta -Psi -V0 * (1.0 + delta + Psi)))));
    if (to_return < 1e-11 && to_return > -1e-11) {to_return = 0;}
    // double G = (V_.beta2 - delta + beta * delta) * (alpha * (1 + (-2) * beta + delta + Psi) + V_.C * (-1 + 2 * IaA0 * beta - delta - Psi - V0 * (1 + delta + Psi)));
    // double F = -IdA0 * V_.GG + G;
    // double D = V_.JJ * IdA0 - 2 * (alpha - beta + delta) * V_.JJ_3 + V_.JJ_1 * (IdA0 * beta * V_.C * V_.JJ_2 - V_.JJ_4 * (alpha * V_.JJ_5 + V_.C * (-1 + 2 * IaA0 * beta - delta + Psi + V0 * (-1 - delta + Psi)))) + V_.JJ_6 * F;
    //  double D = (2 * exp(((-1) * t + t0) * beta) * IdA0 * (beta - 1) * beta * C * Psi -2 * (alpha - beta + delta) * (beta2 + (beta - 1) * delta) * Psi + exp(0.5 * (t - t0) * ((-1) + delta - Psi)) * (IdA0 * beta * C * ((-1) - delta + beta * (3 + delta - Psi) + Psi) - (beta2 - delta + beta * delta) * (alpha * (1 + (-2) * beta + delta - Psi) + C * ((-1) + 2 * IaA0 * beta - delta + Psi + V0 * ((-1) - delta + Psi)))) + exp(0.5 * (t - t0) * ((-1) + delta + Psi)) * ((-1) * IdA0 * beta * (beta+(-1) * delta) * ((-1) - delta - Psi + beta * (3 + delta + Psi)) + (beta2 - delta+beta * delta) * (alpha * (1 + (-2) * beta + delta + Psi) + C * ((-1) + 2 * IaA0 * beta+(-1) * delta - Psi - V0 * (1 + delta + Psi)))));
    if (r_ref_ == 0)
      {
    	//return 0.5 / C / (beta2 + (beta - 1) * delta) / (4 * beta + (- 1) * pow((1 + delta),2)) * Psi * D;
    	// return V_.JJ_7 * D;
    	return to_return;
      }
    else
      {
    	return vrm;
      }
  }


  double
  migliore::default_v_ini(double cor_i, double zeta_eta_csi, double rho)
  {
    //v_ini = ((EL + (1 - np.exp(-(zeta*1000*cor[i] - rho*1000*ith)/1000) )*(vtm - EL))/(-EL))
    double to_return = ((P_.E_L_ + (1 - exp(-(zeta_eta_csi * 1000 * cor_i - rho * 1000 * V_.I_th_) / 1000) )
			 *(P_.V_th_ - P_.E_L_))/(-P_.E_L_));
    return set_v_ini(to_return, S_.r_ref_, V_.vrm);
  }
  
  double
  migliore::set_v_ini(double to_v_ini,
		      int r_ref_, double vrm)
  {
    if (r_ref_ == 0)
      {
    	return to_v_ini;
      }
    else
      {
    	return vrm;
      }
  }



  double
  migliore::Iadap(double t, double delta, double Psi, double alpha,
		  double beta, double IaA0, double IdA0, double t0,
		  double V0, int r_ref_)
  {
    double to_return = ( - 2.0 * alpha * ( - 4 * pow(beta, 3.0) + pow(beta, 2.0) * pow(( - 1.0 + delta), 2.0) - delta * pow((1.0 + delta), 2.0) + beta * delta * (5 + 2.0 * delta + pow(delta, 2.0))) + 2.0 * exp((( - 1.0) * t + t0) * beta) * IdA0 * beta * (4 * pow(beta, 2.0) + delta * pow((1.0 + delta), 2.0) - beta * (1.0 + 6 * delta + pow(delta, 2.0))) + exp((1.0 / 2.0) * (t - t0) * ( - 1.0 + delta + Psi)) * ( - IdA0 * beta * (beta - delta) * ( - 1.0 + ( - 2.0) * delta - pow(delta, 2.0) - Psi + delta * Psi + 2.0 * beta * (2.0 + Psi)) + (pow(beta, 2.0) - delta + beta * delta) * (alpha * (1.0 + ( - 4) * beta + 2.0 * delta + pow(delta, 2.0) + Psi - delta * Psi) + (beta - delta) * (4 * IaA0 * beta - 2.0 * (1.0 + V0) * Psi + IaA0 * (1.0 + delta) * ( - 1.0 - delta + Psi)))) + exp(( - 1.0) * (1.0 / 2.0) * (t - t0) * (1.0 - delta + Psi)) * (IdA0 * beta * (beta - delta) * (1.0 + 2.0 * delta + pow(delta, 2.0) - Psi + delta * Psi + 2.0 * beta * ( - 2.0 + Psi)) + (pow(beta, 2.0) - delta + beta * delta) * (alpha * (1.0 - 4 * beta + 2.0 * delta + pow(delta, 2.0) - Psi + delta * Psi) - (beta - delta) * ( - 4 * IaA0 * beta - 2.0 * (1.0 + V0) * Psi + IaA0 * (1.0 + delta) * (1.0 + delta + Psi))))) / (2.0 * (beta - delta) * (pow(beta, 2.0) + ( - 1.0 + beta) * delta) * (4 * beta - pow((1.0 + delta), 2.0)));
    // double to_return = (-2.0 * alpha * (-4.0 * pow(beta,3.0) + pow(beta,2.0) * pow((-1.0 + delta),2.0) - delta * pow((1.0 + delta),2.0) + beta * delta * (5 + 2.0 * delta + pow(delta,2.0))) + 2.0 * exp(((-1.0) * t  +  t0)  *  beta) * IdA0 * beta * (4.0 * pow(beta,2.0) + delta * pow((1.0 + delta),2.0) - beta * (1.0 + 6 * delta + pow(delta,2.0))) + exp(0.5 * (t-t0) * (-1.0 + delta + Psi)) * (-IdA0 * beta * (beta-delta) * (-1.0 + (-2.0) * delta - pow(delta,2.0) - Psi + delta * Psi + 2.0 * beta * (2.0 + Psi)) + (pow(beta,2.0) - delta + beta * delta) * (alpha * (1.0 + (-4.0) * beta + 2.0 * delta + pow(delta,2.0) + Psi-delta * Psi) + (beta-delta) * (4.0 * IaA0 * beta-2.0 * (1.0 + V0) * Psi + IaA0 * (1.0 + delta) * (-1.0-delta + Psi)))) + exp((-1.0) * (1.0 / 2.0)  *  (t-t0)  *  (1.0-delta + Psi)) * (IdA0 * beta * (beta-delta) * (1.0 + 2.0 * delta + pow(delta,2.0)-Psi + delta * Psi + 2.0 * beta * (-2.0 + Psi)) + (pow(beta,2.0) - delta + beta * delta) * (alpha * (1.0 - 4.0 * beta + 2.0 * delta + pow(delta,2.0) - Psi + delta * Psi)-(beta-delta) * (-4.0 * IaA0 * beta-2.0 * (1.0 + V0) * Psi + IaA0 * (1.0 + delta) * (1.0 + delta + Psi)))))/(2.0 * (beta-delta) * (pow(beta,2.0) + (-1.0 + beta) * delta) * (4.0 * beta - pow((1.0 + delta),2.0)));
    // double to_return = (-2.0*alpha*(-4.0*pow(beta,3.0)+pow(beta,2.0)*pow((-1.0+delta),2.0)-delta*pow((1.0+delta),2.0)+beta*delta*(5+2.0*delta+pow(delta,2.0)))+2.0*exp(((-1.0)*t + t0) * beta)*IdA0*beta*(4.0*pow(beta,2.0)+delta*pow((1.0+delta),2.0)-beta*(1.0+6*delta+pow(delta,2.0)))+exp((1.0 / 2.0)*(t-t0)*(-1.0+delta+Psi))*(-IdA0*beta*(beta-delta)*(-1.0+(-2.0)*delta-pow(delta,2.0)-Psi+delta*Psi+2.0*beta*(2.0+Psi))+(pow(beta,2.0)-delta+beta*delta)*(alpha*(1.0+(-4.0)*beta+2.0*delta+pow(delta,2.0)+Psi-delta*Psi)+(beta-delta)*(4.0*IaA0*beta-2.0*(1.0+V0)*Psi+IaA0*(1.0+delta)*(-1.0-delta+Psi))))+exp((-1.0)*(1.0 / 2.0) * (t-t0) * (1.0-delta+Psi))*(IdA0*beta*(beta-delta)*(1.0+2.0*delta+pow(delta,2.0)-Psi+delta*Psi+2.0*beta*(-2.0+Psi))+(pow(beta,2.0)-delta+beta*delta)*(alpha*(1.0-4.0*beta+2.0*delta+pow(delta,2.0)-Psi+delta*Psi)-(beta-delta)*(-4.0*IaA0*beta-2.0*(1.0+V0)*Psi+IaA0*(1.0+delta)*(1.0+delta+Psi)))))/(2.0*(beta-delta)*(pow(beta,2.0)+(-1.0+beta)*delta)*(4.0*beta-pow((1.0+delta),2.0)));

    // double to_return = (-2.0 * alpha * V_.AA_1 + IdA0 * V_.AA_2 + V_.AA_3 * (-IdA0 * beta * (beta - delta) * V_.AA_4 + V_.AA_5 * (alpha * V_.AA_6 + (beta - delta) * (4.0 * IaA0 * beta -2.0 * (1.0 + V0) * Psi + IaA0 * V_.AA_7))) + V_.AA_8 * (IdA0 * V_.AA_9 + V_.AA_10 * (alpha * V_.AA_11 - (beta-delta) * (-4.0 * IaA0 * beta - 2.0*(1.0+V0)*Psi+IaA0*(1.0+delta)*(1.0+delta+Psi)))))/ V_.AA_12;
    // double AA = 2.0 * exp(t0 * (-1.0 + beta + delta) + 0.5 * t * (-1.0 + delta + Psi)) * beta * V_.C * (4.0 * beta - V_.pdelta);
    // double AA_2 = 2.0 * exp(t0 * (-1.0 + delta) + 0.5 * t * (-1.0 + 2.0 * beta + delta + Psi)) * (V_.beta2 + (-1.0 + beta) * delta) * (-4.0 * beta + V_.pdelta);
    // double AA_3 = exp(0.5 * t0 * (-1.0 + delta - Psi) + t * (-1.0 + beta + delta + Psi));
    // double AA_8 = exp(t * (-1.0 + beta + delta)+0.5 * t0 * (-1.0 + delta + Psi));
    // double AA_13 = 0.5 * exp(t0 - t0 * delta -0.5 * t * (-1.0 + 2.0 * beta + delta + Psi)) / (beta - delta) / (V_.beta2 + (-1.0 + beta) * delta) / (4.0 * beta - V_.pdelta);
    
    // double D = AA * IdA0 -alpha * AA_2 + AA_3 * (-IdA0 * V_.AA_4 + V_.AA_5 * (alpha * V_.AA_6 + V_.C * (4.0 * IaA0 * beta -2.0 * (1.0 + V0) * Psi + IaA0 * V_.AA_7))) +  AA_8 * (IdA0 * V_.AA_9 + V_.AA_10 * (alpha * V_.AA_11 + V_.C * (4.0 * IaA0 * beta + 2.0 * (1.0 + V0) * Psi - IaA0 * V_.AA_12)));
    // double to_return = AA_13 * D;
    if (r_ref_ == 0)
      {
    	return to_return;
      }
    else
      {
    	return IaA0;
      }
  }

  double
  migliore::Idep(double t, double beta, double IdA0, double t0, int r_ref_)
  {
    if (r_ref_ == 0)
      {
    	return exp(-(t - t0) * P_.bet_)  * IdA0;
      }
    else
      {
    	return IdA0;
      }
  }

  double
  migliore::exp_cum(double x, double a, double b)
  {
    return a * (1.0 - exp(-b * x));
  }

  double
  migliore::monod(double x, double a, double b, double c, double alp)
  {
    double to_return = c + (a * exp(b) * x) / (alp + x);
    return to_return;
  }

  double
  migliore::mgblock(double v)
  {
    double block;
    block = 1 / (1 + exp(P_.mgb_k_ * -(v - P_.mgb_shift_)) * (P_.mg_ / P_.mg_ref_));
    return block;
  }
  
  void
  migliore::update( Time const& origin, const long from, const long to )
  {
    //assert( to >= 0 && ( delay ) from < kernel().connection_manager.get_min_delay() );
    assert( from < to );

    double t0_val = origin.get_ms() / V_.time_scale_;
    double local_time_step = V_.dt;
    double t_final = t0_val + local_time_step;
    std::cout << "from " << from << " to " << to << " V_.d_dt " << V_.d_dt << "\n";
	  
    double v_ini = set_v_ini(S_.V_m_ / V_.Vconvfact,S_.r_ref_, V_.vrm);
    double vini_prec = v_ini;
    double teta;
    double c_aux;
    double currCoeff;
    // int counter = S_.sis_;
    // double timer = origin.get_ms() ;
    // double current;
    // double vmss, timess;
    double paramL_;
    double input_spike;
	  
    // evolve from timestep 'from' to timestep 'to' with steps of h each
    for ( long lag = from; lag < to; ++lag )
      {
	double corpre = S_.current_; // + S_.I;

	// set new input current
	S_.current_ = 0.0;
	S_.I_inj_ = B_.currents_.get_value( lag );
	S_.current_ = S_.I_inj_; // + S_.I;
		  	
        for ( size_t i = 0; i < P_.n_receptors_(); i++ )
          {
            S_.current_ =  S_.current_ + S_.i_syn_[ i ];
          }
	S_.current_ =  S_.current_ + P_.I_e_;
	//Update synaptic currents
	for ( size_t i = 0; i < P_.n_receptors_(); i++ )
	  {
	    // exponential decaying PSCs
	    S_.i_syn_fast_rise_A_[ i ] *= V_.P11_syn_fast_rise_[ i ];
	    S_.i_syn_fast_decay_B_[ i ] *= V_.P11_syn_fast_decay_[ i ];
	    S_.i_syn_slow_rise_A_[ i ] *= V_.P11_syn_slow_rise_[ i ];
	    S_.i_syn_slow_decay_B_[ i ] *= V_.P11_syn_slow_decay_[ i ];
	    // collect spikes
	    input_spike = B_.spikes_[ i ].get_value( lag );
	    S_.i_syn_fast_rise_A_[ i ] += input_spike * V_.syn_fast_factor_[ i ]; // not sure about this
	    S_.i_syn_fast_decay_B_[ i ] += input_spike * V_.syn_fast_factor_[ i ]; // not sure about this
	    if ( i == 0 )
	      {
		S_.i_syn_slow_rise_A_[ i ] += input_spike *
		  V_.syn_slow_factor_[ i ] *
		  P_.NMDA_ratio_ *
		  mgblock(S_.V_m_); // not sure about this
		S_.i_syn_slow_decay_B_[ i ] += input_spike *
		  V_.syn_slow_factor_[ i ] *
		  P_.NMDA_ratio_ *
		  mgblock(S_.V_m_); // not sure about this
	      }
	    S_.i_syn_[ i ] = S_.i_syn_fast_decay_B_[ i ] - S_.i_syn_fast_rise_A_[ i ] + S_.i_syn_slow_decay_B_[ i ] - S_.i_syn_slow_rise_A_[ i ];
	  }
	// current = S_.current_;
	// vmss = S_.V_m_;
	// timess = t_final * V_.time_scale_;
	// t0_val = t_final;
	// t_final = t0_val + local_time_step ;//Time::step( origin.get_steps() + lag + 1 );
	if ((t_final-S_.init_sign_)*V_.time_scale_ >= nest::migliore::tagliorette(S_.current_) and V_.blockActive)
	  {
	    if (S_.current_ > V_.I_th_)
	      {
		if (corpre < V_.I_th_ || S_.sis_ == 0){
		  S_.init_sign_ = t_final;
		  S_.Idep_ini_ = std::max(P_.Idep_ini_vr_, P_.cost_idep_ini_*(S_.current_ - V_.I_th_));
		  S_.Iadap_ini_ = 0;

		  V_.blockActive = false;
		  
		  // currCoeff = (S_.current_ - V_.I_th_)/S_.current_; removed from AGLIF_040 when introduced the copies
		  v_ini = default_v_ini(S_.current_, P_.zeta_, P_.rho_);
		  v_ini = migliV(t_final, P_.delta1_, V_.psi1,
				 S_.current_/V_.sc_, P_.bet_, S_.Iadap_ini_, S_.Idep_ini_, t0_val, v_ini, S_.r_ref_, V_.vrm);
		  S_.Iadap_ini_ = Iadap(t_final, P_.delta1_, V_.psi1, S_.current_ / V_.sc_, P_.bet_, S_.Iadap_ini_, S_.Idep_ini_, t0_val, v_ini, S_.r_ref_);
		  S_.Idep_ini_ = Idep(t_final, P_.bet_, S_.Idep_ini_, t0_val, S_.r_ref_);
		}
	      }
	    if (corpre == 0) {
	      v_ini = set_v_ini(vini_prec, S_.r_ref_, V_.vrm);
	    } else
	      {
		if (S_.current_ < V_.I_th_ && S_.current_ > 0) {
		  v_ini = default_v_ini(S_.current_, P_.eta_, 0 );
		} else if (S_.current_<=0) {
		  v_ini = default_v_ini(S_.current_, P_.csi_, 0 );
		} else {
		  v_ini = default_v_ini(S_.current_, P_.zeta_, P_.rho_ );
		}
	      }
	    vini_prec = v_ini;
	    V_.out.push_back(v_ini);
	    S_.V_m_ = v_ini * V_.Vconvfact;
	    while (V_.out.size() > 3)
	      {
		V_.out.erase(V_.out.begin());
	      }
	    // Count down for refractory period
	    if (S_.r_ref_ > 0) {--S_.r_ref_;}

	    // voltage logging
	    B_.logger_.record_data( origin.get_steps() + lag );

	  } else {
	  vini_prec = v_ini;
	  if ((S_.current_ < V_.I_th_ && S_.current_ >= 0) || S_.sis_ == 0)
	    {
	      v_ini = default_v_ini(S_.current_, P_.eta_, 0 );
	    } else
	    {
	      if ( V_.out.size() > 2 && S_.current_ < corpre && S_.current_ > 0 && ((V_.t_spk + 2 * V_.d_dt) < t_final * V_.time_scale_)) {
		teta = (V_.out[V_.out.size()-1] / (corpre / V_.sc_)) * (1/V_.dt-P_.delta1_)
		  -(V_.out[V_.out.size()-2] / ((corpre / V_.sc_)*V_.dt))
		  -P_.delta1_ / (corpre / V_.sc_) -1;
		if (teta < 0) {teta = 0;}
		S_.Idep_ini_ = S_.Iadap_ini_ + teta * (S_.current_/ V_.sc_) / P_.bet_;
		v_ini = migliV(t_final, P_.delta1_, V_.psi1,
			       S_.current_/V_.sc_, P_.bet_,
			       S_.Iadap_ini_, S_.Idep_ini_,
			       t0_val, v_ini, S_.r_ref_, V_.vrm);
		S_.Iadap_ini_ = Iadap(t_final, P_.delta1_, V_.psi1,
				      S_.current_ / V_.sc_, P_.bet_, S_.Iadap_ini_,
				      S_.Idep_ini_, t0_val, v_ini, S_.r_ref_);
		S_.Idep_ini_ = Idep(t_final, P_.bet_, S_.Idep_ini_, t0_val, S_.r_ref_);
	      } else {
		if (S_.current_ > 0) {
		  v_ini = migliV(t_final, P_.delta1_, V_.psi1,
				 S_.current_/V_.sc_, P_.bet_,
				 S_.Iadap_ini_, S_.Idep_ini_,
				 t0_val, v_ini, S_.r_ref_, V_.vrm);
		  S_.Iadap_ini_ = Iadap(t_final, P_.delta1_, V_.psi1,
					S_.current_ / V_.sc_, P_.bet_, S_.Iadap_ini_,
					S_.Idep_ini_, t0_val, v_ini, S_.r_ref_);
		  S_.Idep_ini_ = Idep(t_final, P_.bet_, S_.Idep_ini_, t0_val, S_.r_ref_);
	      }
	    }
	    if (corpre != S_.current_ && (S_.current_ < 0 && S_.current_ > P_.mincurr_))
	      {
		v_ini = set_v_ini(vini_prec, S_.r_ref_, V_.vrm);
	      }
	    if (S_.current_ < 0 && S_.current_ > P_.mincurr_)
	      {
		v_ini = default_v_ini(S_.current_, P_.csi_, 0 );
	      }
	    if (corpre != S_.current_  && S_.current_ <= P_.mincurr_){
	      S_.Iadap_ini_ = -P_.V_min_ / P_.E_L_ + 1;
	      S_.Idep_ini_ = 0;
	      v_ini = default_v_ini(S_.current_, P_.csi_, 0 );
	    }
	    if (S_.current_ <= P_.mincurr_) {
	      v_ini = V_.V_star_min_;
	    }
	  }
	  if (v_ini * V_.Vconvfact < P_.V_min_){
	    v_ini = P_.V_min_ / V_.Vconvfact;
	    S_.Iadap_ini_ = P_.Iadap_start_;
	  }


	// Count down for refractory period
	if (S_.r_ref_ > 0) {--S_.r_ref_;}

	if (S_.current_ > V_.I_th_) {
	  if (corpre < V_.I_th_) {
	    
	    V_.blockActive = false;

	    S_.init_sign_ = t_final;
	    S_.Idep_ini_ = std::max(P_.Idep_ini_vr_, P_.cost_idep_ini_*(S_.current_ - V_.I_th_));
	    S_.Iadap_ini_ = 0;                
	    v_ini = default_v_ini(S_.current_, P_.zeta_, P_.rho_ );
	    if (S_.current_<1e-11) {
	      v_ini = -1;
	    }
	  }
	}

	V_.out.push_back(v_ini);
	S_.V_m_ = v_ini * V_.Vconvfact;
	while (V_.out.size() > 3)
	  {
	    V_.out.erase(V_.out.begin());
	  }
	// lower bound of membrane potential REMOVED in 041 version
	// S_.V_m_ = ( S_.V_m_ < P_.V_min_ ? P_.V_min_ : S_.V_m_ );

	// Plot VM for debugging latency divergence
	if ( P_.plotit_ ) {
	  if ( t_final * V_.time_scale_ > 1322 ) {
	    std::cout << "t_final = " << t_final << ", t = " << t_final * V_.time_scale_ << ", v_ini = " << v_ini * V_.Vconvfact << ", V_m = " << S_.V_m_ << ", th = " << P_.V_th_ << "\n";
	  }
	}

	// threshold crossing
	if ( S_.V_m_ >= P_.V_th_ )
	  {
	    // S_.V_m_ = P_.V_th_; This is for output graphic so that alla spikes have the same height ASK MICHELE!!!
	    // voltage logging
	    B_.logger_.record_data( origin.get_steps() + lag );

	    S_.r_ref_ = V_.RefractoryCounts_; // Inizialize refractory

	    V_.t_spk = t_final * V_.time_scale_;
	    // f.write(str(round(t_spk, 3)) + ' \n')
	    S_.V_m_ = P_.Vres_; // -65
	    v_ini = V_.vrm;

	    V_.blockActive = true;
	    
	    c_aux = P_.c_;
	    if (S_.current_ < P_.istim_min_spikinig_exp_ || S_.current_ > P_.istim_max_spikinig_exp_)
	      {
		c_aux = P_.c_;
		S_.Iadap_ini_ = monod((t_final-S_.init_sign_) * V_.time_scale_,
				      P_.a_, P_.b_ * S_.current_/1000, P_.c_, P_.alp_);
	      } else {
	      // orig parameters
	      S_.Iadap_ini_ = monod((t_final - S_.init_sign_) * V_.time_scale_, P_.a_, P_.b_ * S_.current_ / 1000, P_.c_, P_.alp_);
	      // print('Iadap_ini: ',Iadap_ini)
	      if (S_.Iadap_ini_<0) {
		// print('monod negativa')
		//sinapt
		paramL_ = S_.Iadap_ini_;
		if (P_.a_ > 0) {
		  c_aux = P_.c_ - paramL_;
		} else {
		  c_aux = -P_.a_ * exp(P_.b_ * S_.current_ / 1000);
		}
		S_.Iadap_ini_ = monod((t_final-S_.init_sign_) * V_.time_scale_, P_.a_, P_.b_ * S_.current_/1000, c_aux, P_.alp_);
	      }
	    }
			  
	    if (S_.current_ < V_.I_th_) {
	      S_.Idep_ini_ = 0;
	      S_.Iadap_ini_ = P_.Iadap_start_;
	    } else {
	      S_.Idep_ini_ = P_.Idep_ini_vr_;
	    }

	    // compute spike time
	    set_spiketime( Time::step( origin.get_steps() + lag + 1 ) );

	    SpikeEvent se;
	    kernel().event_delivery_manager.send( *this, se, lag );
	  }
	else
	  {
	  // voltage logging
	  B_.logger_.record_data( origin.get_steps() + lag );
	  }
	}
	t0_val = t_final;
	t_final = t0_val + local_time_step;
      }
    S_.sis_++;
  }


  size_t
  nest::migliore::handles_test_event( SpikeEvent&, size_t receptor_type )
  {
    if ( receptor_type <= 0 or receptor_type > P_.n_receptors_() )
      {
	throw IncompatibleReceptorType( receptor_type, get_name(), "SpikeEvent" );
      }

    P_.has_connections_ = true;
    return receptor_type;
  }

  void
  migliore::handle( SpikeEvent& e )
  {
    assert( e.get_delay_steps() > 0 );
    B_.spikes_[ e.get_rport() - 1 ].add_value(
					      e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ), e.get_weight() * e.get_multiplicity() );
    //    B_.spikes_.add_value(
    //			 e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ), e.get_weight() * e.get_multiplicity() );
  }

  void
  migliore::handle( CurrentEvent& e )
  {
    assert( e.get_delay_steps() > 0 );

    const double c = e.get_current();
    const double w = e.get_weight();
    B_.currents_.add_value( e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ), w * c );
  }

  void
  migliore::handle( DataLoggingRequest& e )
  {
    B_.logger_.handle( e );
  }


  inline void
  migliore::set_status( const DictionaryDatum& d )
  {
    Parameters_ ptmp = P_;     // temporary copy in case of errors
    ptmp.set( d, this );       // throws if BadProperty
    State_ stmp = S_;          // temporary copy in case of errors
    stmp.set( d, ptmp, this ); // throws if BadProperty

    // We now know that (ptmp, stmp) are consistent. We do not
    // write them back to (P_, S_) before we are also sure that
    // the properties to be set in the parent class are internally
    // consistent.
    ArchivingNode::set_status( d );

    /*
     * Here is where we must update the recordablesMap_ if new receptors
     * are added!
     */
    if ( ptmp.tau_syn_fast_decay_.size() > P_.tau_syn_fast_decay_.size() ) // Number of receptors increased
      {
	for ( size_t i_syn = P_.tau_syn_fast_decay_.size(); i_syn < ptmp.tau_syn_fast_decay_.size(); ++i_syn )
	  {
	    size_t elem = migliore::State_::I_SYN
	      + i_syn * migliore::State_::NUM_STATE_ELEMENTS_PER_RECEPTOR;
	    recordablesMap_.insert( get_i_syn_name( i_syn ), get_data_access_functor( elem ) );
	  }
      }
    else if ( ptmp.tau_syn_fast_decay_.size() < P_.tau_syn_fast_decay_.size() )
      { // Number of receptors decreased
	for ( size_t i_syn = ptmp.tau_syn_fast_decay_.size(); i_syn < P_.tau_syn_fast_decay_.size(); ++i_syn )
	  {
	    recordablesMap_.erase( get_i_syn_name( i_syn ) );
	  }
      }

    // if we get here, temporaries contain consistent set of properties
    P_ = ptmp;
    S_ = stmp;
  }
  
} // namespace
