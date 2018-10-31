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

#include "OrbitConversions/State.hpp"
#include <RigidBodyKinematics.hpp>

namespace OC{

	State::State(arma::vec state,double mu){
		this -> state = state;
		this -> mu = mu;

	}

	
	arma::vec State::get_state() const{
		return this -> state;
	}

	
	void State::set_state(arma::vec state) {
		this -> state = state;
	}


	double State::get_mu() const{
		return this -> mu;
	}

	void State::set_mu(double mu){
		this -> mu = mu;
	}

	double State::get_parameter() const{
		return this -> get_a() * (1 - std::pow(this -> get_eccentricity(),2));
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



	double State::get_n() const{
		return std::sqrt(this -> mu / std::pow(std::abs(this -> get_a()),3));
	}


	double State::H_from_f(const double & f,const double & e){

		return 2 * std::atanh(std::sqrt( (e - 1) / (1 + e) ) * std::tan( f / 2 ));
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

		double ecc = M;

		if (pedantic){
			std::cout << "Initial guess: " << ecc << " from M : " << M <<  " , e: " << e << std::endl;

		}

		for (unsigned int i = 0; i < 1000; ++i){

			double max_decc = 0.1;

			double decc = (State::M_from_ecc(ecc,e) - M)/(1 - e * std::cos(ecc));
			
			int sign;
			if (decc >= 0){
				sign = 1;
			}
			else{
				sign = -1;
			}

			decc =  sign * std::min(max_decc,std::abs(decc));

			ecc = ecc - decc;

			double error = std::abs(State::M_from_ecc(ecc,e) - M);
			if (pedantic) {
				std::cout << "decc : " << decc<<  " , ecc : " <<  ecc << " , Residual : " <<  error << std::endl;
			}

			if (error < 1e-13 || std::abs(decc) < 1e-12 || std::abs(ecc) < 1e-10){
				break;
			}

			if (i == 999){
				std::cout << "State::ecc_from_M did not converge\n";
			}
			
		}

		return ecc;

	}


	double State::H_from_M(const  double & M,const  double & e,const bool & pedantic){

		double H = std::atan(M);
		
		double damp = 1;
		if (pedantic){
			std::cout << "Initial guess: " << H << " from M : " << M <<  " , e: " << e << std::endl;
		}

		for (unsigned int i = 0; i < 1000; ++i){

			double max_dH = 0.1;

			double dH = (State::M_from_H(H,e) - M)/(e * std::cosh(H) - 1);

			int sign;
			if (dH >= 0){
				sign = 1;
			}
			else{
				sign = -1;
			}

			dH = damp * sign * std::min(max_dH,std::abs(dH));

			H = H - dH;

			double error = std::abs(M - State::M_from_H(H,e));

			if (pedantic){
				std::cout << "dH : " << dH <<  " , H : " <<  H <<" , Residual : " << error << std::endl;
			}


			if (error < 1e-13 || std::abs(dH) < 1e-12){
				break;
			}

			if (i == 999){
				std::cout << "State::H_from_M did not converge\n";
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


}
