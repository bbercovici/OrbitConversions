#include "State.hpp"
#include <RigidBodyKinematics.hpp>

State::State(arma::vec state,double mu){
	this -> state = state;
	this -> mu = mu;

}

double State::f_from_ecc(const double & ecc,const double & e){

	double f = 2 * std::atan(std::sqrt((1 + e)/ (1 - e)) * std::tan(ecc/2));

	if (ecc < 0 && f > 0){
		f -= 2 * arma::datum::pi;
	}
	else if (ecc > 0 && f < 0){
		f += 2 * arma::datum::pi;
	}

	return f;
}


double State::get_n() const{
	return std::sqrt(this -> mu / std::pow(this -> get_a(),3));
}


double State::ecc_from_f(const double & f,const double & e){

	double ecc = 2 * std::atan(std::sqrt((1 - e)/ (1 + e)) * std::tan(f/2));

	if (ecc < 0 && f > 0){
		ecc += 2 * arma::datum::pi;
	}
	else if (ecc > 0 && f < 0){
		ecc -= 2 * arma::datum::pi;
	}

	return ecc;
}

double State::H_from_f(const double & f,const double & e){

	return 2 * std::atan(std::sqrt( (e - 1) / (1 + e) ) * std::tan( f / 2 ));
}

double State::f_from_H(const  double & H,const  double & e){

	double f = 2 * std::atan(std::sqrt( (1 + e) / (e - 1) ) * std::tanh( H / 2 ));

	if (H < 0 && f > 0){
		f -= 2 * arma::datum::pi;
	}
	else if (H > 0 && f < 0){
		f += 2 * arma::datum::pi;
	}

	return f;
}


double State::f_from_M(const double & M,const double & e){
	if (e < 1){
		return State::f_from_ecc(State::ecc_from_M(M,e),e);
	}
	else{
		return State::f_from_H(State::H_from_M(M,e),e) ;
	}
}

double State::ecc_from_M(const double & M,const double & e,const bool & pedantic){

    // The eccentric anomaly is found
	bool converge = false;
	double ecc = M;

	while (!converge){
		ecc = ecc - (State::M_from_ecc(ecc,e) - M)/(1 - e * std::cos(ecc));
		if (pedantic) {
			std::cout <<  "Eccentric anomaly: " <<  ecc << std::endl;
		}

		double error = std::abs(State::M_from_ecc(ecc,e) - M);

		if (error < 1e-9){
			converge = true;
		}
		else{
			if (pedantic){
				std::cout << "Residual: " <<  error << std::endl;
			}
		}
	}
	return ecc;

}


double State::H_from_M(const  double & M,const  double & e,const bool & pedantic){

    // The eccentric anomaly is found
	bool converge = false;

	double H = std::atan(M);

	while (!converge){

		double max_dH = 1;

		double dH = (State::M_from_H(H,e) - M)/(e * std::cosh(H) - 1);

		int sign;
		if (dH >= 0){
			sign = 1;
		}
		else{
			sign = -1;
		}

		dH = sign * std::min(max_dH,std::abs(dH));

		H = H - dH;

		if (pedantic){
			std::cout <<  "H: " <<  H << std::endl;
		}


		double error = std::abs(M - State::M_from_H(H,e));

		if (error < 1e-9){
			converge = true;
		}

		else{
			if (pedantic){
				std::cout <<  "Residual: " <<  error << std::endl;
			}
		}
	}
	return H;

}


double State::M_from_ecc(const double & ecc,const double & e){

	return ecc - e * std::sin(ecc);
}


double State::M_from_H(const double & H,const double & e){

	return e * std::sinh(H) - H;
}


double State::M_from_f(const double & f,const double & e){

	if (e < 1){
		return State::M_from_ecc(State::ecc_from_f(f,e),e);
	}
	else{
		return State::M_from_H(State::H_from_f(f,e),e);
	}
}
