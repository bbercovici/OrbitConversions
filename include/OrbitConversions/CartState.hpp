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

#ifndef CARTSTATE_HEADER 
#define CARTSTATE_HEADER

#include "OrbitConversions/State.hpp"
#include "OrbitConversions/KepState.hpp"

namespace OC{

	class KepState;

	class CartState : public State{

	public: 
		CartState( arma::vec state,double mu);
		CartState();

		virtual double get_momentum() const;
		virtual double get_energy() const;
		virtual double get_a() const;
		virtual double get_eccentricity() const;

		double get_speed() const;
		double get_radius() const;
		
		arma::vec::fixed<3> get_position_vector() const;
		arma::vec::fixed<3> get_velocity_vector() const;
		arma::vec::fixed<3> get_momentum_vector() const ;
		arma::vec::fixed<3> get_eccentricity_vector() const ;

		/* 
		Returns the keplerian orbital elements state corresponding to the 
		cartesian state
		@param delta_T time since epoch
		@return keplerian state vector of orbital_elements (a, e, i, Omega, omega, M0)
		*/

		KepState convert_to_kep(double delta_T) const;

	protected:

	};

}

#endif