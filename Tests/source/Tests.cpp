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
		Tests::test_ecc_from_M(N);
		Tests::test_cart_to_kep_to_cart(N);
		Tests::test_cart_to_kep(N);
		Tests::test_kep_to_cart(N);

	}


	void test_H_from_M(int N){

		std::cout <<  "\n- Running test_H_from_M... \n";

		arma::arma_rng::set_seed(N);

		for (int i = 0; i < N; ++i){

			arma::vec rands = arma::randu<arma::vec>(2);

			double e = rands(0) + 1;
			double M = arma::datum::pi * ( 0.5 - rands(1));

			double error = std::abs( OC::State::M_from_H(OC::State::H_from_M(M,e),e) - M);
			assert(error < 5e-8);
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

			assert(error < 5e-8);
		}
		std::cout << "- test_ecc_from_M() passed\n";

	}


	void test_cart_to_kep_to_cart(int N){

		std::cout <<  "\n- Running test_cart_to_kep_to_cart... \n" ;

		arma::arma_rng::set_seed(N);

		for (int i = 0; i < N; ++i){

			arma::vec rands = arma::randn<arma::vec>(1);
			double dt = rands(0);
			arma::vec rand_state = arma::randu<arma::vec>(6);
			OC::CartState cart_(rand_state,1); 
			std::cout << "before returning\n";
			std::cout << cart_.get_state() << std::endl;
			std::cout << "after returning\n";

			throw;

			OC::CartState cart(rand_state,1); 

			std::cout << cart.get_state().t() << std::endl;

			OC::KepState kep = cart.convert_to_kep(dt);
			std::cout << kep.get_state().t() << std::endl;

			OC::CartState cart_from_kep = kep.convert_to_cart(dt);
			std::cout << cart_from_kep.get_state().t() << std::endl;


			double error = arma::norm(cart_from_kep.get_state() - cart.get_state())/arma::norm(cart.get_state());

			assert(error < 1e-7);

		}
		std::cout << "- test_cart_to_kep_to_cart() passed\n";

	}


	void  test_cart_to_kep(int N){
		std::cout <<  "\n- Running test_cart_to_kep\n";

		arma::arma_rng::set_seed(N);

		for (int i = 0; i < N; ++i){

			
			arma::vec rand_v = arma::randu<arma::vec>(1);
			double dt = rand_v(0);

			OC::CartState cart(arma::randu<arma::vec>(6),1);
			OC::KepState kep = cart.convert_to_kep(0);

			assert(std::abs(kep.get_radius(dt) - cart.get_radius())/kep.get_radius(0) < 5e-8);

			double error_momentum_norm = abs(cart.get_momentum() - kep.get_momentum())/kep.get_momentum();
			assert(error_momentum_norm < 5e-8);

			double error_v = std::abs(kep.get_speed(dt) - cart.get_speed()) / kep.get_speed(dt);
			assert(error_v < 5e-8);

			double error_energy = (kep.get_energy() - cart.get_energy()) ;
			assert(error_energy < 5e-8);



		}
		std::cout <<  "- test_cart_to_kep() passed" << std::endl;
	}

	void test_kep_to_cart(int N){
		
		std::cout << "\n- Running test_kep_to_cart \n" ;

		arma::arma_rng::set_seed(N);

		for (int i = 0; i < N; ++i){

			arma::vec kep_state_vec = arma::zeros<arma::vec>(6);
			arma::vec rands = arma::randu<arma::vec>(7);

			kep_state_vec(1) = 2 * rands(1);

			if (kep_state_vec(1) > 1){
				kep_state_vec(0) = - rands(0);
			}
			else{
				kep_state_vec(0) = rands(0);
			}

			kep_state_vec(2) = arma::datum::pi * rands(2);
			kep_state_vec(3) = 2 * arma::datum::pi * rands(3);

			kep_state_vec(4) = 2 * arma::datum::pi * rands(4);
			kep_state_vec(5) = 10 * (0.5 - rands(5));

			double dt = rands(6);

			OC::KepState kep(kep_state_vec,1);

			OC::CartState cart = kep.convert_to_cart(dt);

			assert(std::abs(kep.get_radius(dt) - cart.get_radius() )< 5e-8);

			double error_momentum_norm = std::abs(cart.get_momentum()- kep.get_momentum())/cart.get_momentum();
			assert(error_momentum_norm < 5e-8);

			double error_v = std::abs(cart.get_speed()- kep.get_speed(dt))/kep.get_speed(dt);
			assert(error_v < 5e-8);

			double error_energy = std::abs(cart.get_energy()- kep.get_energy())/kep.get_energy();
			assert(error_energy < 5e-8);
		}

		std::cout <<  "- test_cart_to_kep() passed\n";

	}




	// void test_kep_to_cart_to_kep(N){
	// 	std::cout << "\n- Running test_kep_to_cart_to_kep..." ;

	// 	arma::arma_rng::set_seed(N);

	// 	for (int i = 0; i < N; ++i){

	// 		arma::vec kep_state_vec = arma::zeros<arma::vec>(6);
	// 		arma::vec rands = arma::randu<arma::vec>(6);

	// 		kep_state_vec(1) = 2 * rands(1);

	// 		if (kep_state_vec(1) > 1){
	// 			kep_state_vec(0) = - rands(0);
	// 		}
	// 		else{
	// 			kep_state_vec(0) = rands(0);
	// 		}

	// 		kep_state_vec(2) = arma::datum::pi * rands(2);
	// 		kep_state_vec(3) = 2 * arma::datum::pi * rands(3);

	// 		kep_state_vec(4) = 2 * arma::datum::pi * rands(4);
	// 		kep_state_vec(5) = 10 * (0.5 - rands(5));

	// 		KepState kep_state(kep_state,1);
	// 		CartState cart_from_kep = kep_state.convert_to_cart();
	// 		KepState kep_from_cart = cart_from_kep.convert_to_kep();

	// 		if kep[-1] < 0 and kep_from_cart[-1] > 0:
	// 			kep_from_cart[-1] -= 2 * arma::datum::pi

	// 		if kep[-2] < 0 and kep_from_cart[-2] > 0:
	// 			kep_from_cart[-2] -= 2 * arma::datum::pi

	// 		if kep[-3] < 0 and kep_from_cart[-3] > 0:
	// 			kep_from_cart[-3] -= 2 * arma::datum::pi

	// 		if kep[-1] > 0 and kep_from_cart[-1] < 0:
	// 			kep_from_cart[-1] += 2 * arma::datum::pi

	// 		if kep[-2] > 0 and kep_from_cart[-2] < 0:
	// 			kep_from_cart[-2] += 2 * arma::datum::pi

	// 		if kep[-3] > 0 and kep_from_cart[-3] < 0:
	// 			kep_from_cart[-3] += 2 * arma::datum::pi

	// 		error_vec = kep_from_cart - kep

	// 		error = np.linalg.norm(error_vec)/np.linalg.norm(kep)


	// 		assert(error < 5e-8)
	// 	}
	// 	print "- test_kep_to_cart_to_kep() passed"

	// }




















}

