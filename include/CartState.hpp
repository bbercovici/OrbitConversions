#ifndef CARTSTATE_HEADER 
#define CARTSTATE_HEADER

#include "State.hpp"


	class CartState : public State{

	public: 
		CartState(arma::vec state,double mu);
		CartState();


		virtual double get_momentum() const;
		virtual double get_energy() const;
		virtual double get_parameter() const;
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

#endif