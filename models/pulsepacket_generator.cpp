/*
 *  pulsepacket_generator.cpp
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

#include "pulsepacket_generator.h"

// C++ includes:
#include <algorithm>

// Includes from libnestutil:
#include "dict_util.h"
#include "numerics.h"

// Includes from nestkernel:
#include "event_delivery_manager_impl.h"
#include "exceptions.h"
#include "kernel_manager.h"

// Includes from sli:
#include "dict.h"
#include "dictutils.h"
#include "doubledatum.h"


/* ----------------------------------------------------------------
 * Default constructors defining default parameters and variables
 * ---------------------------------------------------------------- */

nest::pulsepacket_generator::Parameters_::Parameters_()
  : pulse_times_()
  , a_( 0 )
  , sdev_( 0.0 )
  , sdev_tolerance_( 10.0 )
{
}

nest::pulsepacket_generator::Variables_::Variables_()
  : start_center_idx_( 0 )
  , stop_center_idx_( 0 )
  , tolerance( 0.0 )
{
}

/* ----------------------------------------------------------------
 * Parameter extraction and manipulation functions
 * ---------------------------------------------------------------- */

void
nest::pulsepacket_generator::Parameters_::get( DictionaryDatum& d ) const
{
  ( *d )[ names::pulse_times ] = DoubleVectorDatum( new std::vector< double >( pulse_times_ ) );
  ( *d )[ names::activity ] = a_;
  ( *d )[ names::sdev ] = sdev_;
}

void
nest::pulsepacket_generator::Parameters_::set( const DictionaryDatum& d, pulsepacket_generator& ppg, Node* node )
{
  // We cannot use a single line here since short-circuiting may stop evaluation
  // prematurely. Therefore, neednewpulse must be second arg on second line.
  bool neednewpulse = updateValueParam< long >( d, names::activity, a_, node );
  neednewpulse = updateValueParam< double >( d, names::sdev, sdev_, node ) or neednewpulse;
  if ( a_ < 0 )
  {
    throw BadProperty( "The activity cannot be negative." );
  }
  if ( sdev_ < 0 )
  {
    throw BadProperty( "The standard deviation cannot be negative." );
  }


  if ( updateValue< std::vector< double > >( d, "pulse_times", pulse_times_ ) or neednewpulse )
  {
    std::sort( pulse_times_.begin(), pulse_times_.end() );
    ppg.B_.spiketimes_.clear();
  }
}

/* ----------------------------------------------------------------
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

nest::pulsepacket_generator::pulsepacket_generator()
  : StimulationDevice()
  , P_()
{
}

nest::pulsepacket_generator::pulsepacket_generator( const pulsepacket_generator& ppg )
  : StimulationDevice( ppg )
  , P_( ppg.P_ )
{
}

/* ----------------------------------------------------------------
 * Node initialization functions
 * ---------------------------------------------------------------- */

void
nest::pulsepacket_generator::init_state_()
{
  StimulationDevice::init_state();
}

void
nest::pulsepacket_generator::init_buffers_()
{
  StimulationDevice::init_buffers();
}

void
nest::pulsepacket_generator::pre_run_hook()
{
  StimulationDevice::pre_run_hook();
  assert( V_.start_center_idx_ <= V_.stop_center_idx_ );

  if ( P_.sdev_ > 0.0 )
  {
    V_.tolerance = P_.sdev_ * P_.sdev_tolerance_;
  }
  else
  {
    V_.tolerance = 1.0;
  }

  const double now = ( kernel().simulation_manager.get_time() ).get_ms();

  V_.start_center_idx_ = 0;
  V_.stop_center_idx_ = 0;


  // determine pulse-center times that lie within
  // a window sdev*sdev_tolerance around the current time
  while (
    V_.stop_center_idx_ < P_.pulse_times_.size() and P_.pulse_times_.at( V_.stop_center_idx_ ) - now <= V_.tolerance )
  {
    if ( std::abs( P_.pulse_times_.at( V_.stop_center_idx_ ) - now ) > V_.tolerance )
    {
      V_.start_center_idx_++;
    }
    V_.stop_center_idx_++;
  }
}


void
nest::pulsepacket_generator::update( Time const& T, const long, const long to )
{
  if ( ( V_.start_center_idx_ == P_.pulse_times_.size() and B_.spiketimes_.empty() )
    or ( not StimulationDevice::is_active( T ) ) )
  {
    return; // nothing left to do
  }

  // determine next pulse-center times (around sdev*tolerance window)
  if ( V_.stop_center_idx_ < P_.pulse_times_.size() )
  {
    while ( V_.stop_center_idx_ < P_.pulse_times_.size()
      and ( Time( Time::ms( P_.pulse_times_.at( V_.stop_center_idx_ ) ) ) - T ).get_ms() <= V_.tolerance )
    {
      V_.stop_center_idx_++;
    }
  }

  if ( V_.start_center_idx_ < V_.stop_center_idx_ )
  {
    RngPtr rng = get_vp_specific_rng( get_thread() );

    bool needtosort = false;

    while ( V_.start_center_idx_ < V_.stop_center_idx_ )
    {
      for ( int i = 0; i < P_.a_; i++ )
      {
        double x = P_.sdev_ * V_.normal_dist_( rng ) + P_.pulse_times_.at( V_.start_center_idx_ );
        if ( Time( Time::ms( x ) ) >= T )
        {
          B_.spiketimes_.push_back( Time( Time::ms( x ) ).get_steps() );
        }
      }
      needtosort = true;
      V_.start_center_idx_++;
    }
    if ( needtosort )
    {
      std::sort( B_.spiketimes_.begin(), B_.spiketimes_.end() );
    }
  }

  int n_spikes = 0;

  // Since we have an ordered list of spiketimes,
  // we can compute the histogram on the fly.
  while ( not B_.spiketimes_.empty() and B_.spiketimes_.front() < ( T.get_steps() + to ) )
  {
    n_spikes++;
    long prev_spike = B_.spiketimes_.front();
    B_.spiketimes_.pop_front();

    if ( n_spikes > 0 and ( B_.spiketimes_.empty() or prev_spike != B_.spiketimes_.front() ) )
    {
      SpikeEvent se;
      se.set_multiplicity( n_spikes );
      kernel().event_delivery_manager.send( *this, se, prev_spike - T.get_steps() );
      n_spikes = 0;
    }
  }
}

/* ----------------------------------------------------------------
 * Other functions
 * ---------------------------------------------------------------- */

void
nest::pulsepacket_generator::set_data_from_stimulation_backend( std::vector< double >& input_param )
{
  Parameters_ ptmp = P_; // temporary copy in case of errors

  // For the input backend
  if ( not input_param.empty() )
  {
    if ( input_param.size() < 3 )
    {
      throw BadParameterValue(
        "The size of the data for the pulse_generator needs to be higher than 3 "
        "[activity, sdev, all the pulse times]." );
    }
    DictionaryDatum d = DictionaryDatum( new Dictionary );
    ( *d )[ names::activity ] = DoubleDatum( input_param[ 0 ] );
    ( *d )[ names::sdev ] = DoubleDatum( input_param[ 1 ] );
    input_param.erase( input_param.begin(), input_param.begin() + 2 );
    ( *d )[ names::pulse_times ] = DoubleVectorDatum( input_param );
    ptmp.set( d, *this, this );
  }

  // if we get here, temporary contains consistent set of properties
  P_ = ptmp;
}
