#include "KepState.hpp"
#include <RigidBodyKinematics.hpp>


	KepState::KepState(arma::vec state,double mu) : State(state,mu){
	}

	KepState::KepState() : State(arma::zeros<arma::vec>(6),1){
	}


	double KepState::get_energy() const{
		return this -> mu / (2 * this -> get_a());
	}

	double KepState::get_momentum() const{
		return std::sqrt(this -> mu * this -> get_parameter());
	}

	double KepState::get_speed(double dt) const{
		return std::sqrt(this -> mu * (2./(this -> get_radius(dt)) - 1./(this -> get_a())));
	}

	double KepState::get_parameter() const{
		return this -> get_a() * (1 - std::pow(this -> get_eccentricity(),2));
	}

	double KepState::get_radius(double dt) const{
		double M = this -> get_M0() + this -> get_n() * dt;
		return this -> get_parameter() / (1 + this -> get_eccentricity() * std::cos(OC::f_from_M(M,this -> get_eccentricity())));
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


	