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

#include "CartState.hpp"
#include <RigidBodyKinematics.hpp>

namespace OC{

	CartState::CartState(arma::vec state,double mu) : State(state,mu){

	}

	CartState::CartState() : State(arma::zeros<arma::vec>(6),1){

	}


	arma::vec::fixed<3> CartState::get_position_vector() const{
		return this -> state.subvec(0,2);
	}


	arma::vec::fixed<3> CartState::get_velocity_vector() const{
		return this -> state.subvec(3,5);
	}

	double CartState::get_speed() const{
		return arma::norm(this -> get_velocity_vector());
	}

	double CartState::get_radius() const{
		return arma::norm(this -> get_position_vector());
	}

	double CartState::get_momentum() const{
		return arma::norm(this -> get_momentum_vector());
	}

	double CartState::get_energy() const{
		return std::pow(this -> get_speed(), 2) / 2 - this -> mu /(this -> get_radius());
	}

	double CartState::get_parameter() const{
		return this -> get_a() * (1 - std::pow(this -> get_eccentricity(),2));
	}

	double CartState::get_a() const{
		double energy = this -> get_energy();
		return  - this -> mu/(2 * energy);
	}

	arma::vec::fixed<3> CartState::get_momentum_vector() const {
		arma::vec::fixed<3> h_vector = arma::cross(this -> get_position_vector(), this -> get_velocity_vector());
		return h_vector;
	}

	arma::vec::fixed<3> CartState::get_eccentricity_vector() const {
		arma::vec::fixed<3> ecc_vector = (arma::cross(this -> get_velocity_vector(),this -> get_momentum_vector())/(this -> mu) 
			- this -> get_position_vector() / arma::norm(this -> get_position_vector()));
		return ecc_vector;
	}

	double CartState::get_eccentricity() const{
		return arma::norm(this -> get_eccentricity_vector());
	}

	KepState CartState::convert_to_kep(double delta_T) const{

    // semi major axis
		double a = this -> get_a();

    // spacecraft's angular momentum
		arma::vec::fixed<3> h_vector = this -> get_momentum_vector();
		double h = this -> get_momentum();

    // eccentricity
		arma::vec::fixed<3> ecc_vector = this -> get_eccentricity_vector();
		double e = this -> get_eccentricity();

    // conic parameter
		double p = this -> get_parameter();

    // line of nodes
		arma::vec::fixed<3> Z_axis = {0,0,1};
		arma::vec::fixed<3> nodal_line = arma::cross(Z_axis,h_vector /h);

    // orbit DCM
		arma::mat::fixed<3,3> ON = arma::zeros<arma::mat>(3,3);
		ON.row(0) = ecc_vector.t() / e;
		ON.row(1) = arma::cross(h_vector / h,ecc_vector / e).t();
		ON.row(2) = h_vector / h;

		double Omega = std::atan2(ON(2,0),-ON(2,1));
		double i = std::acos(ON(2,2));
		double omega = std::atan2(ON(0,2),ON(1,2));

		double f = std::acos(1./e * (p / this -> get_radius() - 1));

		if (arma::dot(this -> get_position_vector(),this -> get_velocity_vector()) < 0){
			f = 2 * arma::datum::pi - f;
		}

		double ecc,M,H;
		if (e < 1){
        // eccentric anomaly
			ecc = State::ecc_from_f(f,e);

        // mean anomaly
			M = ecc - e * std::sin(ecc);
		}

		else {

        // Hyperbolic anomaly
			H = State::H_from_f(f,e);


        // mean anomaly
			M = e * std::sinh(H) - H;
		}


    // mean motion
		double n = std::sqrt(this -> mu / std::pow(std::abs(a) , 3));

    // initial mean anomaly
		double M0 = M - n * delta_T;

		arma::vec kep_state = {a,e,i,Omega,omega,M0};

		return KepState(kep_state,this -> mu);

	}
}