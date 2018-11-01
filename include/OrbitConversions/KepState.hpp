
// MIT License

// Copyright (c) 2018 Benjamin Bercovici

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
#ifndef KEPSTATE_HEADER 
#define KEPSTATE_HEADER

#include "State.hpp"
#include "CartState.hpp"

namespace OC{

	class CartState;
	class KepState : public State{

	public: 

		/**
		Constructor
		@param state 6x1 vector of orbital elements ordered like so :
		- sma : semi-major axis [L]
		- e : eccentricity [-]
		- i : inclination in [0,pi] [rad]
		- Omega : right-ascension of ascending node in [0,2 pi] [rad] 
		- omega : longitude of perigee [0,2 pi] [rad] 
		- M0 : mean anomaly at epoch [rad]
		@param mu standard gravitational parameter of central body [L^3/T^2]
		*/
		KepState(arma::vec state,double mu);
		KepState();


		/**
		Returns orbit energy
		@return energy (J)
		*/
		virtual double get_energy() const;

		/**
		Returns orbit sma
		@return sma (m)
		*/
		virtual double get_a() const;


		/**
		Returns orbit sma
		@return sma (m)
		*/
		virtual double get_eccentricity() const;

		/**
		Returns orbit momentum
		@return orbit momentum (m^2/s)
		*/
		virtual double get_momentum() const;

		double get_inclination() const;
		double get_Omega() const;
		double get_omega() const;
		double get_M0() const;

		/**
		Returns orbit speed 
		@param true anomaly
		@return orbit speed (m/s)
		*/
		double get_speed(double f) const;


		/**
		Returns orbit radius 
		@param true anomaly
		@return orbit radius (m)
		*/
		double get_radius(double f) const;

		/* 
		Returns the cartesian state corresponding to the 
		keplerian state
		@param delta_T time since epoch
		@return cartesian state vector (x, y, z, x_dot, y_dot, z_dot)
		*/

		CartState convert_to_cart(double delta_T) const;

		
	protected:

	};

}
#endif