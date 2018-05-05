

#ifndef KEPSTATE_HEADER 
#define KEPSTATE_HEADER

#include "KepState.pp"


class KepState : public State{

public: 
	KepState(arma::vec state,double mu);
	KepState();

	virtual double get_energy() const;
	virtual double get_a() const;
	virtual double get_eccentricity() const;
	virtual double get_momentum() const;
	virtual double get_parameter() const;

	double get_inclination() const;
	double get_Omega() const;
	double get_omega() const;
	double get_M0() const;
	double get_speed(double dt) const;
	double get_radius(double dt) const;

		/* 
		Returns the cartesian state corresponding to the 
		keplerian state
		@param delta_T time since epoch
		@return cartesian state vector (x, y, z, x_dot, y_dot, z_dot)
		*/

	CartState convert_to_cart(double delta_T) const;

protected:

};

#endif