
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

#include "KepState.hpp"
#include "CartState.hpp"



namespace OC{

	class CartState;
	class KepState : public State{

	public: 
		KepState(arma::vec state,double mu);
		KepState();

		virtual double get_energy() const;
		virtual double get_a() const;
		virtual double get_eccentricity() const;
		virtual double get_momentum() const;

		double get_inclination() const;
		double get_Omega() const;
		double get_omega() const;
		double get_M0() const;
		double get_speed(double f) const;
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