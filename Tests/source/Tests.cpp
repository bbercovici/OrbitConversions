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

#include "Tests.hpp"
#include <RigidBodyKinematics.hpp>
#include <OrbitConversions.hpp>
#include <cassert>

namespace Tests{
	void run_tests(int N){
		

		Tests::test_H_from_M(N);
		Tests::test_f_from_H(N);

		Tests::test_ecc_from_M(N);
		Tests::test_f_from_ecc(N);

		Tests::test_cart_to_kep(N);
		Tests::test_kep_to_cart(N);
		Tests::test_cart_to_kep_to_cart(N);

	}


	void test_f_from_H(int N){

		std::cout <<  "\n- Running test_f_from_H... \n";

		arma::arma_rng::set_seed(N);

		for (int i = 0; i < N; ++i){

			arma::vec rands = arma::randu<arma::vec>(2);

			double e = 3 * rands(0) + 1;
			double H = 4 * arma::datum::pi * ( 0.5 - rands(1));

			double error = std::abs( OC::State::H_from_f(OC::State::f_from_H(H,e),e) - H);
			assert(error < 1e-8);
		}

		std::cout <<  "- test_f_from_H() passed" << std::endl;


	}

	void test_f_from_ecc(int N){

		std::cout <<  "\n- Running test_f_from_ecc... \n";

		arma::arma_rng::set_seed(N);

		for (int i = 0; i < N; ++i){

			arma::vec rands = arma::randu<arma::vec>(2);

			double e = rands(0);
			double f = 4 * arma::datum::pi * ( 0.5 - rands(1));

			double error = std::abs( OC::State::f_from_ecc(OC::State::ecc_from_f(f,e),e) - f);
			assert(error < 1e-8);
		}

		std::cout <<  "- test_f_from_ecc() passed" << std::endl;

		
	}



	void test_H_from_M(int N){

		std::cout <<  "\n- Running test_H_from_M... \n";

		arma::arma_rng::set_seed(N);

		for (int i = 0; i < N; ++i){

			arma::vec rands = arma::randu<arma::vec>(2);

			double e = 3 * rands(0) + 1;
			double M = 2 * arma::datum::pi * ( 0.5 - rands(1));

			double error = std::abs( OC::State::M_from_H(OC::State::H_from_M(M,e),e) - M);
			assert(error < 1e-8);
		}

		std::cout <<  "- test_H_from_M() passed" << std::endl;

	}



	void test_ecc_from_M(int N){

		std::cout <<  "\n- Running test_ecc_from_M... \n" ;
		arma::arma_rng::set_seed(N);
		for (int i = 0; i < N; ++i){

			arma::vec rands = arma::randu<arma::vec>(2);

			double e =rands(0);
			double M =  2 * arma::datum::pi * rands(1);

			double error = std::abs( OC::State::M_from_ecc(OC::State::ecc_from_M(M,e),e) - M);

			assert(error < 1e-8);
		}
		std::cout << "- test_ecc_from_M() passed\n";

	}

	void test_ecc_from_f(int N){

		std::cout <<  "\n- Running test_ecc_from_f... \n" ;
		arma::arma_rng::set_seed(N);
		for (int i = 0; i < N; ++i){

			arma::vec rands = arma::randu<arma::vec>(2);

			double e =rands(0);
			double M =  2 * arma::datum::pi * rands(1);
			double error = std::abs( OC::State::ecc_from_f(OC::State::f_from_ecc(M,e),e) - M);

			assert(error < 1e-8);
		}
		std::cout << "- test_ecc_from_M() passed\n";

	}

	void test_cart_to_kep_to_cart(int N){

		std::cout <<  "\n- Running test_cart_to_kep_to_cart... \n" ;

		arma::arma_rng::set_seed(N);

		for (int i = 0; i < N; ++i){
			arma::vec rands = arma::randu<arma::vec>(2);
			double dt = rands(0);
			double mu = rands(1) + 1;

			OC::CartState cart(arma::randn<arma::vec>(6),mu); 

			OC::KepState kep = cart.convert_to_kep(dt);

			OC::CartState cart_from_kep = kep.convert_to_cart(dt);

			double error = arma::norm(cart_from_kep.get_state() - cart.get_state())/arma::norm(cart.get_state());
			
			assert(error < 1e-7);
		}
		std::cout << "- test_cart_to_kep_to_cart() passed\n";

	}


	void  test_cart_to_kep(int N){
		std::cout <<  "\n- Running test_cart_to_kep\n";

		arma::arma_rng::set_seed(N);

		for (int i = 0; i < N; ++i){

			arma::vec rands = arma::randu<arma::vec>(2);
			double dt = rands(0);
			double mu = 1 + rands(1);

			OC::CartState cart(arma::randn<arma::vec>(6),mu);
			OC::KepState kep = cart.convert_to_kep(dt);

			double f = OC::State::f_from_M(kep.get_M0() + dt * kep.get_n(),kep.get_eccentricity());

			double error_radius = std::abs(kep.get_radius(f) - cart.get_radius())/cart.get_radius();
			assert(error_radius < 1e-7);

			double error_momentum_norm = std::abs(cart.get_momentum() - kep.get_momentum())/kep.get_momentum();
			assert(error_momentum_norm < 1e-7);

			double error_v = std::abs(kep.get_speed(f) - cart.get_speed()) / kep.get_speed(f);
			assert(error_v < 1e-7);

			double error_energy = std::abs(kep.get_energy() - cart.get_energy())/std::abs(cart.get_energy());
			assert(error_energy < 1e-7);



		}
		std::cout <<  "- test_cart_to_kep() passed" << std::endl;
	}

	void test_kep_to_cart(int N){
		
		std::cout << "\n- Running test_kep_to_cart \n" ;

		arma::arma_rng::set_seed(N);

		for (int i = 0; i < N; ++i){

			arma::vec kep_state_vec = arma::zeros<arma::vec>(6);
			arma::vec rands = arma::randu<arma::vec>(8);

			kep_state_vec(1) = 2 * rands(1);

			if (kep_state_vec(1) > 1){
				kep_state_vec(0) = - (rands(0) + 0.1);
			}
			else{
				kep_state_vec(0) = (rands(0) + 0.1);
			}

			kep_state_vec(2) = arma::datum::pi * rands(2);
			kep_state_vec(3) = 2 * arma::datum::pi * rands(3);

			kep_state_vec(4) = 2 * arma::datum::pi * rands(4);
			kep_state_vec(5) = 3 * (0.5 - rands(5));

			double dt = rands(6);
			double mu = 1 + rands(7);

			OC::KepState kep(kep_state_vec,mu);
			double M = kep.get_M0() + dt * kep.get_n();
			
			OC::CartState cart = kep.convert_to_cart(dt);

			double f = OC::State::f_from_M(M,kep.get_eccentricity());

			assert(std::abs(kep.get_radius(f) - cart.get_radius() )< 1e-8);

			double error_momentum_norm = std::abs(cart.get_momentum()- kep.get_momentum())/cart.get_momentum();
			assert(error_momentum_norm < 1e-8);

			double error_v = std::abs(cart.get_speed()- kep.get_speed(f))/kep.get_speed(f);
			assert(error_v < 1e-8);

			double error_energy = std::abs(cart.get_energy()- kep.get_energy())/kep.get_energy();
			assert(error_energy < 1e-8);
		}

		std::cout <<  "- test_kep_to_cart() passed\n";

	}




}

