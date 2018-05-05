#ifndef STATE_HEADER 
#define STATE_HEADER

#include <armadillo>

	class State{

	public:

		State(arma::vec state,double mu);

		static double f_from_ecc(const double & ecc,const double & e);
		static double ecc_from_f(const double & f,const double & e);
		static double H_from_f(const double & f,const double & e);
		static double f_from_H(const  double & H,const  double & e);
		static double f_from_M(const double & M,const double & e);
		static double ecc_from_M(const double & M,const double & e,const bool & pedantic = false);
		static double H_from_M(const  double & M,const  double & e,const bool & pedantic = false);
		static double M_from_ecc(const double & ecc,const double & e);
		static double M_from_H(const double & H,const double & e);
		static double M_from_f(const double & f,const double & e);

		virtual double get_momentum() const = 0;
		virtual double get_radius() const = 0;
		virtual double get_energy() const = 0;
		virtual double get_parameter() const = 0;
		virtual double get_a() const = 0;
		virtual double get_eccentricity() const = 0;

		double get_n() const;


	protected:
		arma::vec::fixed<6> state;
		double mu;

	};

#endif