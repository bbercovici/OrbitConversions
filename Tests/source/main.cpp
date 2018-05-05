#include <Tests.hpp>
#include <armadillo>
#include<OrbitConversions.hpp>



int main(){

	Tests::run_tests(1000);



	// std::cout << OC::State::H_from_f(OC::State::f_from_H(-2.88515,3.30171),3.30171) << std::endl;

	return 0;



}