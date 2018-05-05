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

			double error = std::abs( OC::M_from_H(OC::H_from_M(M,e),e) - M);
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

			double error = std::abs( OC::M_from_ecc(OC::ecc_from_M(M,e),e) - M);

			assert(error < 5e-8);
		}
		std::cout << "- test_ecc_from_M() passed\n";

	}


	void test_cart_to_kep_to_cart(int N){

		std::cout <<  "\n- Running test_cart_to_kep_to_cart... \n" ;

		arma::arma_rng::set_seed(N);

		for (int i = 0; i < N; ++i){

			OC::CartState cart(arma::randu<arma::vec>(6),1); 
			OC::KepState kep = cart.convert_to_kep(0);
			OC::CartState cart_from_kep = kep.convert_to_cart(0);

			double error = arma::norm(cart_from_kep.get_state() - cart.get_state())/np.linalg.norm(cart.get_state());

			assert(error < 1e-7);

		}
		std::cout << "- test_cart_to_kep_to_cart() passed\n";

	}


	void  test_cart_to_kep(int N){
		print "\n- Running test_cart_to_kep\n";

		arma::arma_rng::set_seed(N);

		for (int i = 0; i < N; ++i){

			OC::CartState cart(arma::randu<arma::vec>(6),1);
			OC::KepState kep = cart.convert_to_kep();

			double f = OC::f_from_M(kep.get_M0(),kep.get_eccentricity()); // Assume dt = 0
			double r_kep = kep.get_a() * ( 1- std::pow(kep.get_eccentricity(), 2)) / (1 + kep.get_eccentricity() * std::cos(f));

			assert(std::abs(r_kep - cart.get_radius())/r_kep < 5e-8);

			double h_vector =  cart.get_momentum_vector();

			double h =  cart.get_momentum();

			double error_momentum_norm = abs(h - std::sqrt(kep.get_a() * ( 1- std::pow(kep.get_eccentricity(), 2))))/h;
			assert(error_momentum_norm < 5e-8)

			double v_norm_kep = std::sqrt(2./kep.get_radius() - 1./kep.get_a());
			double error_v = std::abs(v_norm_kep - cart.get_speed()) / v_norm_kep;
			assert(error_v < 5e-8);

			double error_energy = (kep.get_energy() - cart.get_energy()) ;
			assert(error_energy < 5e-8);

			arma::vec h_dir_cart = h_vector / h;

			arma::mat DCM_ON = RBK::M3(f + kep.get_omega()) * RBK::M1(kep.get_inclination()) * RBK::M3(kep.get_Omega());
			arma::vec Z_vec = {0,0,1};
			arma::vec h_dir_kep = DCM_ON.t() * Z_vec;

			double error_h_dir = arma::norm(arma::cross(h_dir_kep,h_dir_cart));
			assert(error_h_dir < 5e-8);


		}
		std::cout <<  "- test_cart_to_kep() passed" << std::endl;
	}

	void test_kep_to_cart(int N){
		
		std::cout << "\n- Running test_kep_to_cart \n" ;

		arma::arma_rng::set_seed(N);

		for (int i = 0; i < N; ++i){

			arma::vec kep_state_vec = arma::zeros<arma::vec>(6);
			arma::vec rands = arma::randu<arma::vec>(6);

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

			KepState kep_state(kep_state,1);

			CartState cart = kep_state.convert_to_cart();

			double f = OC::f_from_M(kep_state.get_M0(),kep_state.get_eccentricity()); // Assume dt = 0

			double r_kep = kep.get_radius();

			assert(std::abs(r_kep - cart.get_radius() )< 5e-8);


			arma::vec h_vec_cart = cart.get_momentum_vector();
			double h_cart = cart.get_momentum();

			double error_momentum_norm = std::abs(cart.get_momentum()- kep.get_momentum());
			assert(error_momentum_norm < 5e-8)

			double error_v = std::abs(cart.get_speed()- kep.get_speed());
			assert(error_v < 5e-8);

			double error_energy = std::abs(cart.get_energy()- kep.get_energy());
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

