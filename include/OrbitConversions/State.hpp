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


		/**
		Computes true anomaly from eccentric anomaly 
		@param ecc eccentric anomaly
		@param e eccentricity (0 =< e < 1)
		@return true anomaly
		*/
		static double f_from_ecc(const double & ecc,const double & e);

		/**
		Computes eccentric anomaly eccentric from true anomaly
		@param f true anomaly
		@param e eccentricity (0 =< e < 1)
		@return eccentric anomaly
		*/
		static double ecc_from_f(const double & f,const double & e);

		/**
		Computes hyperbolic anomaly from true anomaly 
		@param f true anomaly
		@param e eccentricity (1 < e)
		@return hyperbolic anomaly
		*/
		static double H_from_f(const double & f,const double & e);

		/**
		Computes true anomaly from hyperbolic anomaly
		@param H hyperbolic anomaly
		@param e eccentricity (1 < e)
		@return true anomaly
		*/
		static double f_from_H(const  double & H,const  double & e); 

		/**
		Computes true anomaly from mean anomaly
		@param M mean anomaly
		@param e eccentricity (0 =< e < 1)
		@return true anomaly
		*/
		static double f_from_M(const double & M,const double & e);

		/**
		Computes eccentric anomaly eccentric from mean anomaly
		@param M mean anomaly
		@param e eccentricity (0 =< e < 1)
		@param pedantic if true, will print out convergence details
		@return eccentric anomaly
		*/
		static double ecc_from_M(const double & M,const double & e,const bool & pedantic = false);
		

		/**
		Computes hyperbolic anomaly eccentric from mean anomaly
		@param M mean anomaly
		@param e eccentricity (1 < e)
		@param pedantic if true, will print out convergence details
		@return hyperbolic anomaly
		*/
		static double H_from_M(const  double & M,const  double & e,const bool & pedantic = false);
		

		/**
		Computes mean anomaly from eccentric anomaly 
		@param ecc eccentric anomaly
		@param e eccentricity (0 =< e < 1)
		@return mean anomaly
		*/
		static double M_from_ecc(const double & ecc,const double & e);

		/**
		Computes mean anomaly from hyperbolic anomaly 
		@param H hyperbolic anomaly
		@param e eccentricity (1 < e)
		@return mean anomaly
		*/
		static double M_from_H(const double & H,const double & e);

		/**
		Computes mean anomaly from true anomaly 
		@param f true anomaly
		@param e eccentricity (0 =< e < 1)
		@return mean anomaly
		*/
		static double M_from_f(const double & f,const double & e);

		virtual double get_momentum() const = 0;
		virtual double get_energy() const = 0;
		virtual double get_a() const = 0;
		virtual double get_eccentricity() const = 0;


		/**
		Get the state vector
		@return 6x1 state vector (either cartesian state or keplerian state) 
		*/
		arma::vec get_state() const;

		/**
		Get mean motion
		@return mean motion (rad/s)
		*/
		double get_n() const;


		/**
		Get ellipse parameter
		@return ellipse parameter (m)
		*/
		double get_parameter() const;


		/**
		Get standard gravitational parameter
		@return standard gravitational parameter (kg^3/s^2)
		*/
		double get_mu() const;


		/**
		Set standard gravitational parameter
		@param mu standard gravitational parameter (kg^3/s^2)
		*/
		void set_mu(double mu);

		/**
		Sets the state to the prescribed value
		@param state 6x1 state
		*/
		void set_state(arma::vec state) ;



	protected:
		arma::vec state;
		double mu;

	};

}

#endif