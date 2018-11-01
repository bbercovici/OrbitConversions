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

#include "OrbitConversions/KepState.hpp"
#include <RigidBodyKinematics.hpp>

namespace OC{

	KepState::KepState( arma::vec state,double mu) : State(state,mu){
	}

	KepState::KepState() : State(arma::zeros<arma::vec>(6),1){
	}


	double KepState::get_energy() const{
		return - this -> mu / (2 * this -> get_a());
	}

	double KepState::get_momentum() const{
		return std::sqrt(this -> mu * this -> get_parameter());
	}

	double KepState::get_speed(double f) const{
		return std::sqrt(this -> mu * (2./(this -> get_radius(f)) - 1./(this -> get_a())));
	}

	double KepState::get_radius(double f) const{
		return this -> get_parameter() / (1 + this -> get_eccentricity() * std::cos(f));
	}

	double KepState::get_a() const{
		return this -> state(0);
	}

	double KepState::get_eccentricity() const{
		return this -> state(1);
	}

	double KepState::get_inclination() const{
		return this -> state(2);
	}

	double KepState::get_Omega() const{
		return this -> state(3);
	}	

	double KepState::get_omega() const{
		return this -> state(4);
	}

	double KepState::get_M0() const{
		return this -> state(5);
	}


	CartState KepState::convert_to_cart(double dt) const{

		double M = this -> get_M0() + this -> get_n() * dt;
		double f = State::f_from_M(M,this -> get_eccentricity());
		

		double r = this -> get_radius(f);
		double r_dot = this -> get_momentum() / this -> get_parameter() * this -> get_eccentricity() * std::sin(f);

		arma::vec cartesian_state(6);

		arma::mat DCM_ON = RBK::M3(this -> get_omega()+ f) * RBK::M1(this -> get_inclination()) * RBK::M3(this -> get_Omega());
		arma::vec angular_velocity = {0,0,this -> get_momentum() / std::pow(r, 2)};
		arma::vec pos = {r,0,0};
		arma::vec vel = {r_dot,0,0};

		cartesian_state.subvec(0,2) = DCM_ON.t() * pos;
		cartesian_state.subvec(3,5) = DCM_ON.t() *(vel + arma::cross(angular_velocity,pos));

		return CartState(cartesian_state,this -> mu);


	}

}
