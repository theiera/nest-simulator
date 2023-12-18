/*
 *  migliore_ode.cpp
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

#include "migliore_ode.h"

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

//nest::RecordablesMap< nest::migliore_ode > nest::migliore_ode::recordablesMap_;

namespace nest
{
  // Override the create() method with one call to RecordablesMap::insert_()
  // for each quantity to be recorded.
  template <>
  void
  DynamicRecordablesMap< migliore_ode >::create( migliore_ode& host )
  {
    // use standard names wherever you can for consistency!
    insert( names::V_m, host.get_data_access_functor( nest::migliore_ode::State_::V_M ) );
    //    insert( names::I, host.get_data_access_functor( nest::migliore_ode::State_::I ) );
    insert( names::I_syn, host.get_data_access_functor( nest::migliore_ode::State_::I_ve ) );

    insert( names::Iadap_ini, host.get_data_access_functor( nest::migliore_ode::State_::Iadap ) );
    insert( names::Idep_ini, host.get_data_access_functor( nest::migliore_ode::State_::Idep ) );

    host.insert_current_recordables();
  }
  // void
  // RecordablesMap< migliore_ode >::create()
  // {
  //   // use standard names whereever you can for consistency!
  //   insert_( names::V_m, &nest::migliore_ode::get_V_m_ );
  //   insert_( names::I, &nest::migliore_ode::get_I_inj_ );
  //   insert_( names::init_sign, &nest::migliore_ode::get_init_sign_ );
  //   insert_( names::Iadap_ini, &nest::migliore_ode::get_I_adap );
  //   insert_( names::Idep_ini, &nest::migliore_ode::get_I_dep );
  //   insert_( names::sis, &nest::migliore_ode::get_sis_ );
  //   insert_( names::I_syn, &nest::migliore_ode::get_I_syn_ );
  // }


  Name
  nest::migliore_ode::get_i_syn_name( size_t elem )
  {
    std::stringstream i_syn_name;
    i_syn_name << "I_syn_" << elem + 1;
    return Name( i_syn_name.str() );
  }

  void
  nest::migliore_ode::insert_current_recordables( size_t first )
  {
    for ( size_t receptor = first; receptor < P_.tau_syn_.size(); ++receptor )
      {
	size_t elem = nest::migliore_ode::State_::I_SYN
	  + receptor * nest::migliore_ode::State_::NUM_STATE_ELEMENTS_PER_RECEPTOR;
	Name a = get_i_syn_name( receptor );
	recordablesMap_.insert( a , this->get_data_access_functor( elem ) );
      }
  }

  DataAccessFunctor< migliore_ode >
  nest::migliore_ode::get_data_access_functor( size_t elem )
  {
    return DataAccessFunctor< migliore_ode >( *this, elem );
  }

  /* ----------------------------------------------------------------
   * Default constructors defining default parameters and state
   * ---------------------------------------------------------------- */

  nest::migliore_ode::Parameters_::Parameters_()
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
    , tau_syn_( 3.0, 2.0 )  // in ms
    , has_connections_( false )
    , mg_( 2.0) // 2 mM in the Johnston et al. 2010, extracellula [MgCl2] = 1 mM in Edelman et al. 2015
    , mgb_k_ ( 0.062 ) // (/mV) Johnston et al. 2010
    , mg_ref_ ( 3.57 ) // (mM) Johnston et al. 2010
    , mgb_shift_ ( 0.0 ) // (mM) Johnston et al. 2010
    , mincurr_ (-185.0 )  // THIS SHOULD BE A PARAMETER!!!
  {
  }

  nest::migliore_ode::State_::State_( const Parameters_& p )
    : V_m_( p.E_L_ )       // membrane potential  (E_L_+(1-exp(-cor[i]/1000))*(V_th_-E_L_))/V_.Vconvfact_  TOBEFIXED
    , init_sign_( 0.0 ) // membrane adaptation variable
    , Iadap_ini_( 0.0 ) // membrane adaptation variable
    , Idep_ini_( 0.0 ) // membrane dependent variable
    , sis_( 0.0 ) // membrane dependent variable
    , I_inj_( 0.0 )         // input current
    , r_ref_( 0.0 )
    , current_ ( 0.0 )
  {
    i_syn_.clear();
    i_syn_fast_.clear();
    i_syn_slow_.clear();
    y_[ V_M ] = p.E_L_; // initialize to reversal potential
    for ( size_t i = 1; i < STATE_VEC_SIZE; ++i )
      {
	y_[ i ] = 0;
      }
}

  /* ----------------------------------------------------------------
   * Parameter and state extractions and manipulation functions
   * ---------------------------------------------------------------- */

  void
  nest::migliore_ode::Parameters_::get( DictionaryDatum& d ) const
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
    def< double >( d, names::tau_syn_NMDA, tau_syn_NMDA_);
    def< double >( d, names::NMDA_ratio, NMDA_ratio_);
    def< double >( d, names::t_ref, t_ref_ );
    def< double >( d, names::mg, mg_ );
    def< double >( d, names::mg_ref, mg_ref_ );
    def< double >( d, names::mgb_k, mgb_k_ );
    def< double >( d, names::mgb_shift, mgb_shift_ );
    def< double >( d, names::mincurr, mincurr_ );
    def< int >( d, names::n_synapses, n_receptors_() );
    def< bool >( d, names::has_connections, has_connections_ );

    ArrayDatum tau_syn_ad( tau_syn_ );
    def< ArrayDatum >( d, names::tau_syn, tau_syn_ad );
  }

  void
  nest::migliore_ode::Parameters_::set( const DictionaryDatum& d, Node* node )
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
    updateValueParam< double >( d, names::tau_syn_NMDA, tau_syn_NMDA_, node );
    updateValueParam< double >( d, names::NMDA_ratio, NMDA_ratio_, node );
    updateValueParam< double >( d, names::t_ref, t_ref_, node );
    updateValueParam< double >( d, names::mg, mg_, node );
    updateValueParam< double >( d, names::mg_ref, mg_ref_, node );
    updateValueParam< double >( d, names::mgb_k, mgb_k_, node );
    updateValueParam< double >( d, names::mgb_shift, mgb_shift_, node );
    updateValueParam< double >( d, names::mincurr, mincurr_, node );
    if ( t_ref_ < 0 )
      {
    	throw BadProperty( "Refractory time must not be negative." );
      }

    const size_t old_n_receptors = this->n_receptors_();
    if ( updateValue< std::vector< double > >( d, "tau_syn", tau_syn_ ) )
      {
    	if ( this->n_receptors_() != old_n_receptors && has_connections_ == true )
    		{
    			throw BadProperty( "The neuron has connections, therefore the number of ports cannot be reduced." );
    		}
	for ( size_t i = 0; i < tau_syn_.size(); ++i )
	  {
	    if ( tau_syn_[ i ] <= 0 )
	      {
	    	throw BadProperty( "All synaptic time constants must be strictly positive." );
	      }
	    if ( tau_syn_[ i ] == tao_m_ )
	      {
	    	throw BadProperty( "Membrane and synapse time constant(s) must differ. See note in documentation." );
	      }
	  }
      }
  }

  void
  nest::migliore_ode::State_::get( DictionaryDatum& d, const Parameters_& ) const
  {
    def< double >( d, names::init_sign, init_sign_ ); // Membrane potential adaptationariable
    def< double >( d, names::Iadap_ini, Iadap_ini_ ); // Membrane potential adaptationariable
    def< double >( d, names::Idep_ini , Idep_ini_ ); // Membrane potential adaptationariable
    def< double >( d, names::sis , sis_ ); // Membrane potential adaptationariable
    def< double >( d, names::I, I_inj_ ); // Membrane potential
    def< double >( d, names::V_m, V_m_ ); // Membrane potential
  }

  void
  nest::migliore_ode::State_::set( const DictionaryDatum& d, const Parameters_&, Node* node )
  {
    updateValueParam< double >( d, names::init_sign, init_sign_, node );
    updateValueParam< double >( d, names::Iadap_ini, Iadap_ini_, node );
    updateValueParam< double >( d, names::Idep_ini, Idep_ini_, node );
    updateValueParam< double >( d, names::sis, sis_, node );
    updateValueParam< double >( d, names::V_m, V_m_, node );
    updateValueParam< double >( d, names::I, I_inj_, node );
  }

  nest::migliore_ode::Buffers_::Buffers_( migliore_ode& n )
    : logger_( n )
    , s_( 0 )
    , c_( 0 )
    , e_( 0 )
  {
  }

  nest::migliore_ode::Buffers_::Buffers_( const Buffers_&, migliore_ode& n )
    : logger_( n )
    , s_( 0 )
    , c_( 0 )
    , e_( 0 )
  {
  }

  extern "C" inline int
  migliore_ode_dynamics(double, const double y[], double f[],
			void *pnode) {

    //const double tIe = nest::kernel().simulation_manager.get_time().get_ms();

    typedef nest::migliore_ode::State_ State_;
    // get access to node so we can almost work as in a member function
    assert(pnode);
    const nest::migliore_ode& node =
      *(reinterpret_cast<nest::migliore_ode* >(pnode));

    // y[] here is---and must be---the state vector supplied by the integrator,
    // not the state vector in the node, node.S_.y_[].


    // // Synaptic current: I_syn = - sum_k g_k (V - E_rev_k).
    // double I_syn = 0.0;
    // I_syn += y[ State_::G1] * ( node.get_E_rev1() - y[State_::V_M] );
    // I_syn += y[ State_::G2] * ( node.get_E_rev2() - y[State_::V_M] );
    // I_syn += y[ State_::G3] * ( node.get_E_rev3() - y[State_::V_M] );


    // Total current
    // double I_tot = y[State_::Idep] - y[State_::Iadap] + node.get_I_() + node.B_.currents_last_value_; // + I_syn;

    // std::cout << "PArams alpha " << node.P_.alp_ << " beta " << node.P_.bet_ << " delta " << node.P_.delta1_ << "\n";
    // Model currents
    f[State_::Idep] = (-node.P_.bet_) * y[State_::Idep];
    // f[State_::Iadap] = ((-node.get_k2())) * y[State_::Iadap] + node.get_kadap()*(y[State_::V_M] - node.get_E_L());
    f[State_::Iadap] = 1.0 - y[State_::Iadap] + y[State_::V_M];
    if (f[State_::Iadap] < 1e-11) { f[State_::Iadap] = 0.0;}
    // f[State_::Iadap] = f[State_::Iadap] < 1e-11 ? 0.0 : f[State_::Iadap];

 
    // f[State_::V_M] =
    //     ((1)) / node.get_tau_m() * (y[State_::V_M] - node.get_E_L()) +
    //     1 / node.get_C_m() * I_tot ;
    f[State_::V_M] = node.B_.model_alpha + node.P_.bet_ * ( y[State_::Idep] - y[State_::Iadap] ) + node.P_.delta1_ * (1 + y[State_::V_M]);
    // std::cout << node.P_.alp_ << " " << node.P_.bet_ * ( y[State_::Idep] - y[State_::Iadap] )  << " " <<  node.P_.delta1_ * (1 + y[State_::V_M]);
    // std::cout << " f[State_::V_M]  " << f[State_::V_M] << " f[State_::Iadap]  " << f[State_::Iadap] << " f[State_::Idep]  " << f[State_::Idep] << "\n";
  

    // // Conductance dynamics
    //  // Synaptic conductance derivative dG/dt
    //  f[State_::DG1] = -y[State_::DG1] / node.get_tau_syn1();
    //  f[State_::G1] = y[ State_::DG1] - y[ State_::G1] / node.get_tau_syn1();

    //  f[State_::DG2] = -y[State_::DG2] / node.get_tau_syn2();
    //  f[State_::G2] = y[ State_::DG2] - y[ State_::G2] / node.get_tau_syn2();

    //  f[State_::DG3] = -y[State_::DG3] / node.get_tau_syn3();
    //  f[State_::G3] = y[ State_::DG3] - y[ State_::G3] / node.get_tau_syn3();

    return GSL_SUCCESS;
  }

  /* ----------------------------------------------------------------
   * Default and copy constructor for node
   * ---------------------------------------------------------------- */

  nest::migliore_ode::migliore_ode()
    : ArchivingNode()
    , P_()
    , S_( P_ )
    , B_( *this )
  {
    recordablesMap_.create( *this );
  }

  nest::migliore_ode::migliore_ode( const migliore_ode& n )
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
  nest::migliore_ode::init_buffers_()
  {
    B_.spikes_.clear();   // includes resize
    B_.currents_.clear(); // includes resize
    B_.logger_.reset();   // includes resize
    ArchivingNode::clear_history();
    int state_size = State_::STATE_VEC_SIZE;

    B_.step_ = Time::get_resolution().get_ms();
    B_.IntegrationStep_ = B_.step_;
  
  if ( not B_.s_ )
  {
    B_.s_ = gsl_odeiv_step_alloc( gsl_odeiv_step_rkf45, State_::STATE_VEC_SIZE );
  }
  else
  {
    gsl_odeiv_step_reset( B_.s_ );
  }

  if ( not B_.c_ )
  {
    B_.c_ = gsl_odeiv_control_y_new( 1e-3, 0.0 );
  }
  else
  {
    gsl_odeiv_control_init( B_.c_, 1e-3, 0.0, 1.0, 0.0 );
  }

  if ( not B_.e_ )
  {
    B_.e_ = gsl_odeiv_evolve_alloc( State_::STATE_VEC_SIZE );
  }
  else
  {
    gsl_odeiv_evolve_reset( B_.e_ );
  }


    B_.sys_.function = migliore_ode_dynamics;
    B_.sys_.jacobian = nullptr;
    B_.sys_.dimension = state_size;
    B_.sys_.params = reinterpret_cast<void *>(this);
  }

  void
  nest::migliore_ode::pre_run_hook()
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

    V_.P11_syn_fast_.resize( P_.n_receptors_() );
    V_.P11_syn_slow_.resize( P_.n_receptors_() );
    // V_.P21_syn_.resize( P_.n_receptors_() );

    S_.i_syn_.resize( P_.n_receptors_() );
    S_.i_syn_fast_.resize( P_.n_receptors_() );
    S_.i_syn_slow_.resize( P_.n_receptors_() );

    B_.spikes_.resize( P_.n_receptors_() );

    // Initialize states
    S_.Iadap_ini_ = P_.Iadap_start_;
    
    // V_.P22_ = std::exp( -h / P_.tao_m_ );
    // V_.P20_ = P_.tao_m_ / P_.Cm_ * ( 1.0 - V_.P22_ );

    for ( size_t i = 0; i < P_.n_receptors_(); i++ )
      {
		V_.P11_syn_fast_[ i ] = std::exp( -h / P_.tau_syn_[ i ] );
		if (i == 0)
		  {
		    V_.P11_syn_slow_[ i ] = std::exp( -h / P_.tau_syn_NMDA_ );
		  }
		else
		  {
		    V_.P11_syn_slow_[ i ] = 0;
		  }
		// these are determined according to a numeric stability criterion
		// V_.P21_syn_[ i ] = propagator_32( P_.tau_syn_[ i ], P_.tao_m_, P_.Cm_, h );
		B_.spikes_[ i ].resize();
      }

    V_.RefractoryCounts_ = Time( Time::ms( P_.t_ref_ ) ).get_steps();
    // since t_ref_ >= 0, this can only fail in error
    assert( V_.RefractoryCounts_ >= 0 );

    V_.time_scale_ = 1 / (-P_.sc_ / (P_.Cm_ * P_.E_L_));
    V_.d_dt = h;
    V_.dt = V_.d_dt/V_.time_scale_;
    V_.t_spk = -3 * V_.d_dt;

    // V_.beta2 = pow(P_.bet_,2);
    // V_.psi1 = pow(-4.0 * P_.bet_ + (pow(1+P_.delta1_,2.0)),0.5);
    V_.t_step = V_.dt;

    V_.V_star_min_ = -P_.V_min_ / P_.E_L_;
    V_.alpha_neg_ = V_.mincurr / P_.sc_;

    
    // V_.H = (90+P_.E_L_)*P_.sc_*(P_.bet_-P_.delta1_)/(P_.E_L_*(V_.mincurr));
    V_.Vconvfact_ = -P_.E_L_;
    V_.vrm_ = P_.Vres_/V_.Vconvfact_;
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
    // V_.VV_1 = 0.5 / (beta - delta) / (pow(beta,2.0) + (beta - 1.0) * delta) / (4.0 * beta - pow((1.0 + delta),2.0)) * Psi;
    // V_.VV_2 = (beta - 1.0) * beta * (beta - delta) * Psi;
    // V_.VV_3 = (pow(beta,2.0) + (beta - 1.0) * delta) * Psi;
    // V_.VV_4 = 0.5 * (-1.0 + delta - Psi);
    // V_.VV_5 = beta * (beta - delta) * (-1.0 - delta + beta * (3.0 + delta - Psi) + Psi);
    // V_.VV_6 = (pow(beta,2.0) - delta + beta * delta);
    // V_.VV_7 = (1.0 -2.0 * beta + delta - Psi);
    // V_.VV_8 = beta * (beta - delta) * (-1.0 - delta - Psi + beta * (3.0 + delta + Psi));
    // V_.VV_9 = (pow(beta,2.0) - delta + beta * delta);
    // V_.VV_10 = (1.0 - 2.0 * beta + delta + Psi);


    // V_.AA_1 = -4.0 * pow(P_.bet_,3.0) + pow(P_.bet_,2.0) * pow((-1.0+P_.delta1_),2.0) - P_.delta1_ * pow((1.0 + P_.delta1_),2.0) + P_.bet_ * P_.delta1_ * (5 + 2.0 * P_.delta1_ + pow(P_.delta1_,2.0));
    // V_.AA_2 = 2.0 * exp(-V_.t_step * P_.bet_) * P_.bet_ * (4.0 * pow(P_.bet_,2.0) + P_.delta1_ * pow((1.0+P_.delta1_),2.0) - P_.bet_ * (1.0 + 6 * P_.delta1_ + pow(P_.delta1_,2.0)));
    // V_.AA_3 = exp((1.0 / 2.0)*V_.t_step*(-1.0+P_.delta1_+P_.psi1_));
    // V_.AA_4 = -1.0-2.0*P_.delta1_ - pow(P_.delta1_,2.0) - P_.psi1_ + P_.delta1_ * P_.psi1_+2.0*P_.bet_*(2.0+P_.psi1_);
    // V_.AA_5 = (pow(P_.bet_,2.0)-P_.delta1_+P_.bet_*P_.delta1_);
    // V_.AA_6 = (1.0-4.0*P_.bet_+2.0*P_.delta1_+pow(P_.delta1_,2.0)+P_.psi1_-P_.delta1_*P_.psi1_);
    // V_.AA_7 = (1.0+P_.delta1_)*(-1.0-P_.delta1_+P_.psi1_);
    // V_.AA_8 = exp((-1.0)*(1.0 / 2.0) * V_.t_step * (1.0-P_.delta1_+P_.psi1_));
    // V_.AA_9 = P_.bet_ * (P_.bet_-P_.delta1_)*(1.0+2.0*P_.delta1_+ pow(P_.delta1_,2.0)-P_.psi1_+P_.delta1_*P_.psi1_+2.0*P_.bet_*(-2.0+P_.psi1_));
    // V_.AA_10 = (pow(P_.bet_,2.0)-P_.delta1_+P_.bet_*P_.delta1_);
    // V_.AA_11 = (1.0-4.0*P_.bet_+2.0*P_.delta1_ + pow(P_.delta1_,2.0) - P_.psi1_+P_.delta1_*P_.psi1_);
    // V_.AA_12 = (2.0*(P_.bet_-P_.delta1_)*(pow(P_.bet_,2.0)+(-1.0+P_.bet_)*P_.delta1_)*(4.0*P_.bet_-pow((1.0+P_.delta1_),2.0)));

    // V_.DD_1 = exp(-V_.t_step * P_.bet_);
  }

  /* ----------------------------------------------------------------
   * Update and spike handling functions
   */

  double
  nest::migliore_ode::tagliorette(double corr)
  {
    double dur_sign = std::numeric_limits<double>::infinity();
    double vinc_inf = 700;
    if (corr < vinc_inf && corr >= 0)
      {
    	dur_sign = 0.68 * corr - 190.0;
      }
    double vinc_sup = std::numeric_limits<double>::infinity();
    if (corr > vinc_sup) // THIS CANNOT HAPPEN!!
      {
    	dur_sign = 76.86-0.028*corr;
      }
    return dur_sign;
  }

  // double
  // nest::migliore_ode::migliV(double t, double delta, double Psi, 
  // 			 double alpha, double beta, double IaA0, 
  // 			 double IdA0, double t0, double V0,
  // 			 int r_ref_, double vrm)
  // {
  //   double to_return = V_.VV_1 * (2.0 * exp((t0 - t) * beta) * IdA0 * V_.VV_2 - 2.0 * (alpha - beta + delta) * V_.VV_3 + exp(V_.VV_4 * (t - t0)) * (IdA0 * V_.VV_5 - V_.VV_6 * (alpha * V_.VV_7 + (beta - delta) * (-1.0 + 2.0 * IaA0 * beta - delta + Psi + V0 * (-1.0 - delta + Psi)))) + exp(0.5 * (t - t0) * (-1.0 + delta + Psi)) * (- IdA0 * V_.VV_8 + V_.VV_9 * (alpha * V_.VV_10 + (beta - delta) * (-1.0 + 2.0 * IaA0 * beta - delta - Psi - V0 * (1.0 + delta + Psi)))));

  //   if (to_return < 1e-11 && to_return > -1e-11) {to_return = 0;}
  //   if (r_ref_ == 0)
  //     {
  //   	//return 0.5 / C / (beta2 + (beta - 1) * delta) / (4 * beta + (- 1) * pow((1 + delta),2)) * Psi * D;
  //   	// return V_.JJ_7 * D;
  //   	return to_return;
  //     }
  //   else
  //     {
  //   	return vrm;
  //     }
  // }

  double
  nest::migliore_ode::default_v_ini(double currCoeff, double cor_i)
  {
    double to_return = (P_.E_L_ + (1 - exp(-(2.5 + currCoeff)*cor_i/1000) )*(P_.V_th_ - P_.E_L_))/(-P_.E_L_);
    return set_v_ini(to_return, S_.r_ref_, V_.vrm_);
  }
  
  double
  nest::migliore_ode::set_v_ini(double to_v_ini,
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



  // double
  // nest::migliore_ode::Iadap(double t, double delta, double Psi, double alpha,
  // 			double beta, double IaA0, double IdA0, double t0,
  // 			double V0, int r_ref_)
  // {
  //   double to_return = (-2.0 * alpha * (-4.0 * pow(beta,3.0) + pow(beta,2.0) * pow((-1.0 + delta),2.0) - delta * pow((1.0 + delta),2.0) + beta * delta * (5 + 2.0 * delta + pow(delta,2.0))) + 2.0 * exp(((-1.0) * t  +  t0)  *  beta) * IdA0 * beta * (4.0 * pow(beta,2.0) + delta * pow((1.0 + delta),2.0) - beta * (1.0 + 6 * delta + pow(delta,2.0))) + exp(0.5 * (t-t0) * (-1.0 + delta + Psi)) * (-IdA0 * beta * (beta-delta) * (-1.0 + (-2.0) * delta - pow(delta,2.0) - Psi + delta * Psi + 2.0 * beta * (2.0 + Psi)) + (pow(beta,2.0) - delta + beta * delta) * (alpha * (1.0 + (-4.0) * beta + 2.0 * delta + pow(delta,2.0) + Psi-delta * Psi) + (beta-delta) * (4.0 * IaA0 * beta-2.0 * (1.0 + V0) * Psi + IaA0 * (1.0 + delta) * (-1.0-delta + Psi)))) + exp((-1.0) * (1.0 / 2.0)  *  (t-t0)  *  (1.0-delta + Psi)) * (IdA0 * beta * (beta-delta) * (1.0 + 2.0 * delta + pow(delta,2.0)-Psi + delta * Psi + 2.0 * beta * (-2.0 + Psi)) + (pow(beta,2.0) - delta + beta * delta) * (alpha * (1.0 - 4.0 * beta + 2.0 * delta + pow(delta,2.0) - Psi + delta * Psi)-(beta-delta) * (-4.0 * IaA0 * beta-2.0 * (1.0 + V0) * Psi + IaA0 * (1.0 + delta) * (1.0 + delta + Psi)))))/(2.0 * (beta-delta) * (pow(beta,2.0) + (-1.0 + beta) * delta) * (4.0 * beta - pow((1.0 + delta),2.0)));
  //   if (r_ref_ == 0)
  //     {
  //   	return to_return;
  //     }
  //   else
  //     {
  //   	return IaA0;
  //     }
  // }

  // double
  // nest::migliore_ode::Idep(double t, double beta, double IdA0, double t0, int r_ref_)
  // {
  //   if (r_ref_ == 0)
  //     {
  //   	return exp(-(t - t0) * P_.bet_)  * IdA0;
  //     }
  //   else
  //     {
  //   	return IdA0;
  //     }
  // }

  double
  nest::migliore_ode::exp_cum(double x, double a, double b)
  {
    return a * (1.0 - exp(-b * x));
  }

  double
  nest::migliore_ode::monod(double x, double a, double b, double c, double alp)
  {
    double to_return = c + (a * exp(b) * x) / (alp + x);
    return to_return;
  }

  double
  nest::migliore_ode::mgblock(double v)
  {
    double block;
    block = 1 / (1 + exp(P_.mgb_k_ * -(v - P_.mgb_shift_)) * (P_.mg_ / P_.mg_ref_));
    return block;
  }
  
  void
  nest::migliore_ode::evolve_ode()
  {
    double t = 0;
    std::cout << "evolve 100\n";
    
    // std::cout << "Evolve ode " << S_.y_[State_::V_M] << "\n";
    while (t < B_.step_) {
      std::cout << "evolve 200\n";
      std::cout << "Evolving ode " << t << " " << B_.step_ << " " << " " << B_.IntegrationStep_ << " " << S_.y_[State_::V_M] << "\n";
      std::cout << "Evolving ode " << "Inj " << S_.current_ <<  " " << S_.r_ref_ <<  " I dep " << S_.y_[State_::Idep] << " Iadap " << S_.y_[State_::Iadap] << " I_ve " << S_.y_[State_::I_ve] << "\n";
      if (S_.y_[State_::Iadap] < 1e-11) { S_.y_[State_::Iadap] = 0.0;}
      
      const int status =
	gsl_odeiv_evolve_apply(B_.e_, B_.c_, B_.s_,
			       &B_.sys_,          // system of ODE
			       &t,                // from t
			       B_.step_,             // to t <= step
			       &B_.IntegrationStep_, // integration step size
			       S_.y_);
      if (status != GSL_SUCCESS) {
	std::cout << "evolve 300\n";
	throw nest::GSLSolverFailure(get_name(), status);
      }
    }
    // std::cout << "Evolved ode " << S_.y_[State_::V_M] << "\n";

  }

  void
  nest::migliore_ode::update( Time const& origin, const long from, const long to )
  {
    //assert( to >= 0 && ( delay ) from < kernel().connection_manager.get_min_delay() );
	  assert( from < to );

	  B_.step_ = V_.dt;
	  B_.IntegrationStep_ = V_.dt;

	  //std::cout << "100\n";
	  double t0_val = origin.get_ms() / V_.time_scale_;
	  double local_time_step = V_.dt;// / (to - 1);
	  double t_final = t0_val + local_time_step;

	  
	  // double v_ini = set_v_ini(S_.V_m_ / V_.Vconvfact_,S_.r_ref_, V_.vrm_);
	  S_.y_[State_::V_M] = S_.y_[ State_::V_M ] / V_.Vconvfact_;

	  double vini_prec = S_.y_[State_::V_M];
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
		  S_.I_inj_ = B_.currents_.get_value( lag );
		  S_.current_ = S_.I_inj_; // + S_.I;


		  std::cout <<  "Current " << S_.current_ << "\n";
		  // Update synaptic currents
		  // if ( S_.r_ref_ == 0 )
		  //   {
		      // neuron not refractory, so evolve add synaptic currents
		      for ( size_t i = 0; i < P_.n_receptors_(); i++ )
			{
			  // S_.V_m_ += V_.P21_syn_[ i ] * S_.i_syn_[ i ];
			  S_.current_ += S_.i_syn_[ i ];
			}
		    // }
		  for ( size_t i = 0; i < P_.n_receptors_(); i++ )
		    {
		      // exponential decaying PSCs
		      S_.i_syn_fast_[ i ] *= V_.P11_syn_fast_[ i ];
		      S_.i_syn_slow_[ i ] *= V_.P11_syn_slow_[ i ];
		      // collect spikes
		      input_spike = B_.spikes_[ i ].get_value( lag );
		      S_.i_syn_fast_[ i ] += input_spike; // not sure about this
		      if ( i == 0 )
			{
			  if (S_.i_syn_slow_[ i ] < 0) {
			    std::cout << S_.i_syn_slow_[ i ] << "\n";
			  }
			  S_.i_syn_slow_[ i ] += input_spike * P_.NMDA_ratio_ * mgblock(S_.V_m_); // not sure about this
			}
		      S_.i_syn_[ i ] = S_.i_syn_fast_[ i ] + S_.i_syn_slow_[ i ];
		    }
		  // current = S_.current_;
		  // vmss = S_.V_m_;
		  // timess = t_final * V_.time_scale_;
		  t0_val = t_final;
		  t_final = t0_val + local_time_step ;//Time::step( origin.get_steps() + lag + 1 );
		  if ((t_final-S_.init_sign_)*V_.time_scale_ >= nest::migliore_ode::tagliorette(S_.current_))
		  {
//		    std::cout << "300\n";
		    if (S_.current_ > P_.I_th_)
		      {
//			std::cout << "320\n";
		      if (corpre < P_.I_th_ || S_.sis_ == 0){
//			std::cout << "340\n";
			S_.init_sign_ = t_final;
			S_.y_[State_::Idep] = std::max(P_.Idep_ini_vr_, P_.cost_idep_ini_*(S_.current_ - P_.I_th_));
			S_.y_[State_::Iadap] = 0;

			currCoeff = (S_.current_ - P_.I_th_)/S_.current_;
			S_.y_[State_::V_M] = default_v_ini(currCoeff, S_.current_);
			B_.model_alpha = S_.current_ / P_.sc_;
			if (S_.r_ref_ == 0){
//			  std::cout << "360\n";
			  evolve_ode();
			}
			// v_ini = migliV(t_final, P_.delta1_, V_.psi1,
			// 	       S_.current_/P_.sc_, P_.bet_, S_.Iadap_ini_, S_.Idep_ini_, t0_val, v_ini, S_.r_ref_, V_.vrm_);
			// S_.Iadap_ini_ = Iadap(t_final, P_.delta1_, V_.psi1, S_.current_ / P_.sc_, P_.bet_, S_.Iadap_ini_, S_.Idep_ini_, t0_val, v_ini, S_.r_ref_);
			// S_.Idep_ini_ = Idep(t_final, P_.bet_, S_.Idep_ini_, t0_val, S_.r_ref_);
		      }
		    }
		    if (corpre == 0) {
		      S_.y_[State_::V_M] = set_v_ini(vini_prec, S_.r_ref_, V_.vrm_);
		    } else
		      {
			if (S_.current_ < P_.I_th_ && S_.current_ > 0) {currCoeff = 0;}
			else if (S_.current_<=0) {currCoeff = 1;}
			else
			  { currCoeff = (S_.current_ - P_.I_th_)/S_.current_;}
			S_.y_[State_::V_M] = default_v_ini(currCoeff, S_.current_);
		      }
		    vini_prec = S_.y_[State_::V_M];
		  } else {
//		    std::cout << "400\n";
		    vini_prec = S_.y_[State_::V_M];
		    if ((S_.current_ < P_.I_th_ && S_.current_ >= 0) || S_.sis_ == 0)
		      {
			currCoeff = 0;
			S_.y_[State_::V_M] = default_v_ini(currCoeff, S_.current_);
		    } else{
		      if ( V_.out.size() > 2 && S_.current_ < corpre && S_.current_ > 0 && ((V_.t_spk + 2 * V_.d_dt) < t_final * V_.time_scale_)) {
			teta = (V_.out[V_.out.size()-1] / (corpre / P_.sc_)) * (1/V_.dt-P_.delta1_)
			  -(V_.out[V_.out.size()-2] / ((corpre / P_.sc_)*V_.dt))
			  -P_.delta1_ / (corpre / P_.sc_) -1;
			if (teta < 0) {teta = 0;}
			// S_.Idep_ini_ = S_.Iadap_ini_ + teta * (S_.current_/ P_.sc_) / P_.bet_;
			S_.y_[State_::Idep] = S_.y_[State_::Iadap] + teta * (S_.current_/ P_.sc_) / P_.bet_;
			// v_ini = migliV(t_final, P_.delta1_, V_.psi1,
			// 	       S_.current_/P_.sc_, P_.bet_,
			// 	       S_.Iadap_ini_, S_.Idep_ini_,
			// 	       t0_val, v_ini, S_.r_ref_, V_.vrm_);
			// S_.Iadap_ini_ = Iadap(t_final, P_.delta1_, V_.psi1,
			// 		      S_.current_ / P_.sc_, P_.bet_, S_.Iadap_ini_,
			// 		      S_.Idep_ini_, t0_val, v_ini, S_.r_ref_);
			// S_.Idep_ini_ = Idep(t_final, P_.bet_, S_.Idep_ini_, t0_val, S_.r_ref_);
			B_.model_alpha = S_.current_ / P_.sc_;
			if (S_.r_ref_ == 0){
			  evolve_ode();
			}
		      } else {
			if (S_.current_ > 0) {
			  // v_ini = migliV(t_final, P_.delta1_, V_.psi1,
			  // 		 S_.current_/P_.sc_, P_.bet_,
			  // 		 S_.Iadap_ini_, S_.Idep_ini_,
			  // 		 t0_val, v_ini, S_.r_ref_, V_.vrm_);
			  // S_.Iadap_ini_ = Iadap(t_final, P_.delta1_, V_.psi1,
			  // 			S_.current_ / P_.sc_, P_.bet_, S_.Iadap_ini_,
			  // 			S_.Idep_ini_, t0_val, v_ini, S_.r_ref_);
			  // S_.Idep_ini_ = Idep(t_final, P_.bet_, S_.Idep_ini_, t0_val, S_.r_ref_);
			  B_.model_alpha = S_.current_ / P_.sc_;
			  if (S_.r_ref_ == 0){
			    evolve_ode();
			  }
			}
		      }
		      if (corpre != S_.current_ && (S_.current_ < 0 && S_.current_ > V_.mincurr))
			{
			  S_.y_[State_::V_M] = set_v_ini( vini_prec, S_.r_ref_, V_.vrm_);
			}
		      if (S_.current_ < 0 && S_.current_ > V_.mincurr)
			{
			  currCoeff = 1;
			  S_.y_[State_::V_M] = default_v_ini(currCoeff, S_.current_);
			}
		      if (corpre != S_.current_  && S_.current_ <= V_.mincurr){
			// S_.Iadap_ini_ = -P_.V_min_ / P_.E_L_ + 1;
			// S_.Idep_ini_ = 0;
		    	  std::cout << "Minucurr " << S_.y_[State_::V_M] << "\n";
			S_.y_[State_::Iadap] = -P_.V_min_ / P_.E_L_ + 1;
			S_.y_[State_::Idep] = 0;
			currCoeff = 1;
			S_.y_[State_::V_M] = default_v_ini(currCoeff, S_.current_);
		      }
		      if (S_.current_ <= V_.mincurr) {
			S_.y_[State_::V_M] = V_.V_star_min_;
		      }
		    }
		    if (S_.y_[State_::V_M] * V_.Vconvfact_ < P_.V_min_){
		    	std::cout << S_.y_[State_::V_M] << "\n";
		      S_.y_[State_::V_M] = P_.V_min_ / V_.Vconvfact_;
		      // S_.Iadap_ini_ = P_.Iadap_start_;
		      S_.y_[State_::Iadap] = P_.Iadap_start_;
		    }
		  }
		  if (S_.current_ > P_.I_th_) {
		    std::cout << "500\n";
		    if (corpre < P_.I_th_) {
		      S_.init_sign_ = t_final;
		      // S_.Idep_ini_ = std::max(P_.Idep_ini_vr_, P_.cost_idep_ini_*(S_.current_ - P_.I_th_));
		      // S_.Iadap_ini_ = 0;
		      S_.y_[State_::Idep] = std::max(P_.Idep_ini_vr_, P_.cost_idep_ini_*(S_.current_ - P_.I_th_));
		      S_.y_[State_::Iadap] = 0;
		      
		      currCoeff = (S_.current_ - P_.I_th_)/S_.current_;
		      S_.y_[State_::V_M] = default_v_ini(currCoeff, S_.current_);
		      if (S_.current_<1e-11) {
			S_.y_[State_::V_M] = -1;
		      }
		    }
		  }

		  V_.out.push_back(S_.y_[State_::V_M]);
		  std::cout <<  "Current " << S_.current_ << "\n";
		  S_.y_[State_::V_M] = S_.y_[State_::V_M] * V_.Vconvfact_;
		  S_.V_m_ = S_.y_[State_::V_M];
		  while (V_.out.size() > 3)
		  {
			  V_.out.erase(V_.out.begin());
		  }
		  // lower bound of membrane potential REMOVED in 041 version
		  // S_.V_m_ = ( S_.V_m_ < P_.V_min_ ? P_.V_min_ : S_.V_m_ );

		  // Count down for refractory period
		  if (S_.r_ref_ > 0) {--S_.r_ref_;}



		  
		  // threshold crossing
		  std::cout <<  "V_m mV " << S_.y_[State_::V_M] << "\n";
		  if ( S_.y_[State_::V_M] >= P_.V_th_ )
		  {
		    std::cout << "600\n";
			  S_.r_ref_ = V_.RefractoryCounts_; // Inizialize refractory

			  V_.t_spk = t_final * V_.time_scale_;
			  // f.write(str(round(t_spk, 3)) + ' \n')
			  S_.y_[State_::V_M] = P_.Vres_; // -65
			  // v_ini = V_.vrm_;
			  c_aux = P_.c_;
			  if (S_.current_ < P_.istim_min_spikinig_exp_ || S_.current_ > P_.istim_max_spikinig_exp_)
			    {
			      c_aux = P_.c_;
			      S_.y_[State_::Iadap] = monod((t_final-S_.init_sign_) * V_.time_scale_,
						    P_.a_, P_.b_ * S_.current_/1000, P_.c_, P_.alp_);
			    } else {
			    // orig parameters
			    S_.y_[State_::Iadap] = monod((t_final - S_.init_sign_) * V_.time_scale_, P_.a_, P_.b_ * S_.current_ / 1000, P_.c_, P_.alp_);
			    // print('Iadap_ini: ',Iadap_ini)
			    if (S_.y_[State_::Iadap] < 0) {
			      // print('monod negativa')
			      //sinapt
			      paramL_ = S_.y_[State_::Iadap];
			      if (P_.a_ > 0) {
				c_aux = P_.c_ - paramL_;
			      } else {
                                c_aux = -P_.a_ * exp(P_.b_ * S_.current_ / 1000);
			      }
			      S_.y_[State_::Iadap] = monod((t_final-S_.init_sign_) * V_.time_scale_, P_.a_, P_.b_ * S_.current_/1000, c_aux, P_.alp_);
			    }
			  }
			  
			  if (S_.current_ < P_.I_th_) {
			    S_.y_[State_::Idep] = 0;
			    S_.y_[State_::Iadap] = P_.Iadap_start_;
			  } else {
			    S_.y_[State_::Idep] = P_.Idep_ini_vr_;
			  }

			  // compute spike time
			  set_spiketime( Time::step( origin.get_steps() + lag + 1 ) );

			  SpikeEvent se;
			  kernel().event_delivery_manager.send( *this, se, lag );
		  }
		  // voltage logging
		  B_.logger_.record_data( origin.get_steps() + lag );
	  }
	  S_.sis_++;
  }


  size_t
  nest::migliore_ode::handles_test_event( SpikeEvent&, size_t receptor_type )
  {
    if ( receptor_type <= 0 or receptor_type > P_.n_receptors_() )
      {
	throw IncompatibleReceptorType( receptor_type, get_name(), "SpikeEvent" );
      }

    P_.has_connections_ = true;
    return receptor_type;
  }

  void
  nest::migliore_ode::handle( SpikeEvent& e )
  {
	  assert( e.get_delay_steps() > 0 );
	  B_.spikes_[ e.get_rport() - 1 ].add_value(
			  e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ), e.get_weight() * e.get_multiplicity() );
	  //    B_.spikes_.add_value(
	  //			 e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ), e.get_weight() * e.get_multiplicity() );
  }

  void
  nest::migliore_ode::handle( CurrentEvent& e )
  {
	  assert( e.get_delay_steps() > 0 );

	  const double c = e.get_current();
	  const double w = e.get_weight();
	  B_.currents_.add_value( e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ), w * c );
  }

  void
  nest::migliore_ode::handle( DataLoggingRequest& e )
  {
	  B_.logger_.handle( e );
  }


inline void
nest::migliore_ode::set_status( const DictionaryDatum& d )
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
  if ( ptmp.tau_syn_.size() > P_.tau_syn_.size() ) // Number of receptors increased
  {
    for ( size_t i_syn = P_.tau_syn_.size(); i_syn < ptmp.tau_syn_.size(); ++i_syn )
    {
      size_t elem = nest::migliore_ode::State_::I_SYN
        + i_syn * nest::migliore_ode::State_::NUM_STATE_ELEMENTS_PER_RECEPTOR;
      recordablesMap_.insert( get_i_syn_name( i_syn ), get_data_access_functor( elem ) );
    }
  }
  else if ( ptmp.tau_syn_.size() < P_.tau_syn_.size() )
  { // Number of receptors decreased
    for ( size_t i_syn = ptmp.tau_syn_.size(); i_syn < P_.tau_syn_.size(); ++i_syn )
    {
      recordablesMap_.erase( get_i_syn_name( i_syn ) );
    }
  }

  // if we get here, temporaries contain consistent set of properties
  P_ = ptmp;
  S_ = stmp;
}
  
} // namespace
