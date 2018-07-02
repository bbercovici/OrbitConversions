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
    

## License

[This software is distributed under the MIT License](https://choosealicense.com/licenses/mit/)




