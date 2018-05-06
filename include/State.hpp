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

#ifndef STATE_HEADER 
#define STATE_HEADER

#include <armadillo>

namespace OC{

	class State{

	public:

		State(arma::vec state,double mu);

		static double f_from_ecc(const double & ecc,const double & e);
		static double ecc_from_f(const double & f,const double & e);
		static double H_from_f(const double & f,const double & e);
		static double f_from_H(const  double & H,const  double & e);
		static double f_from_M(const double & M,const double & e);
		static double ecc_from_M(const double & M,const double & e,const bool & pedantic = true);
		static double H_from_M(const  double & M,const  double & e,const bool & pedantic = true);
		static double M_from_ecc(const double & ecc,const double & e);
		static double M_from_H(const double & H,const double & e);
		static double M_from_f(const double & f,const double & e);

		virtual double get_momentum() const = 0;
		virtual double get_energy() const = 0;
		virtual double get_a() const = 0;
		virtual double get_eccentricity() const = 0;

		arma::vec get_state() const;
		double get_n() const;
		double get_parameter() const;
		double get_mu() const;
		void set_mu(double mu);


	protected:
		arma::vec state;
		double mu;

	};

}

#endif