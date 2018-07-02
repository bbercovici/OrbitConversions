# OrbitConversions

A collection of orbit conversion routines. For now, the following parametrizations are supported:
 1. Cartesian coordinates (x,y,z,x_dot,y_dot,z_dot)
 2. Keplerian elements (a,e,i,Omega,omega,M0) where M0 is the mean anomaly at epoch

## Requires
1. Armadillo
2. CMake
3. [RigidBodyKinematics](https://github.com/bbercovici/RigidBodyKinematics)

## Installation: 

### Mac users

OrbitConversions can be retrieved from Homebrew:

    brew tap bbercovici/self
    brew update
    brew install orbit-conversions

### Unix users (Mac and Linux)

    git clone https://github.com/bbercovici/OrbitConversions.git
    cd OrbitConversions/build
    cmake ..
    make
    make install

## Getting updates

    git pull
    cd build
    cmake ..
    make
    make install
    
## Usage

Remember that the classical Keplerian elements are singular at zero eccentricity and inclinations!

### From Cartesian coordinates to Keplerian elements

    arma::vec cart_state_vec = {5.4970e+03  , 3.5750e+03   ,6.2943e+02,  -4.2249e+00,   5.6701e+00,   3.7519e+00};
    // Standard gravitational parameter of the earth (kg^3/s^2)
    double mu = 398600 ; 
    OC::CartState cart(cart_state_vec,mu);
    
    double dt = 300; // what is the keplerian state 300 seconds since epoch?
    OC::KepState kep = cart.convert_to_kep(dt);
    std::cout << kep.get_state().t() << std::endl;
    // cart.get_state() returns (7.0001e+03   6.0007e-02   5.0000e-01   4.0000e-01   4.0011e-01   5.7831e+00)
    

### From Keplerian coordinates to Cartesian elements

    // Keplerian elements sma, eccentricity, inclination, right-ascension of ascending node, longitude of perigee, true anomaly at epoch
    arma::vec kep_state_vec = {7000,0.06,0.5,0.4,0.4,-0.5};
    // Standard gravitational parameter of the earth (kg^3/s^2)
    double mu = 398600 ; 
    OC::KepState kep(kep_state_vec,mu);
    
    double dt = 300; // what is the cartesian state 300 seconds since epoch?
    OC::CartState cart = kep.convert_to_cart(dt);
    std::cout << cart.get_state().t() << std::endl;
    // cart.get_state() returns (5.4970e+03   3.5750e+03   6.2943e+02  -4.2249e+00   5.6701e+00   3.7519e+00)




## License

[This software is distributed under the MIT License](https://choosealicense.com/licenses/mit/)




