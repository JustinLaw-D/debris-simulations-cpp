// implementation of the NASA standard breakup model

#include <random>
#include <cmath>
#include <array>
#pragma once

using namespace std;

// whether or not a collision is catestrophic
bool is_catastrophic(const double m_s, const double L, const double AM, const double v);

double calc_M(double m_t, double m_d, double v); // calculate M factor for collisions

double find_A(const double L); // finds the area given characteristic length

// calculate total number of debris generated by collisions/explosions
double calc_Ntot(double M, double Lmin, double Lmax, char typ, double C);

// cummulative distribution function in characteristic length for collisions/explosions
double L_cdf(double L, double L_min, double L_max, char typ);

// cummulative distribution function in chi for collisions/explosions
double X_cdf(double x, double x_min, double x_max, double L, char typ);

// cummulative distribution assuming that L < 8cm
double _X_cdf_8(double x, double x_min, double x_max, double L);

// cummulative distribution assuming that L > 11cm
double _X_cdf_11(double x, double x_min, double x_max, double L, char typ);

// cummulative distribution in log10(Delta v)
double v_cdf(double v, double x, char typ);

// cummulative distribution for post-collisions speed V
double vprime_cdf(double V, double v0, double theta, double phi, double x, char typ);

template <class URNG>
array<double, 3> rand_dir(URNG &generator) {
    /*
    generates a random 3-d unit vector, uniformly

    Input(s)
    generator : uniform random number generator

    Output(s):
    dir : random direction in 3d
    */
    uniform_real_distribution<double> udist = uniform_real_distribution<double>(0.0, 1.0);
    double x; double y; double z;
    do { // randomly generate numbers until the vector is inside the unit sphere
        x = udist(generator); y = udist(generator); z = udist(generator);
    } while (x*x + y*y + z*z > 1.0);
    double norm = sqrt(x*x + y*y + z*z); // get the magintude of the vector to normalize
    return array<double, 3>({x/norm, y/norm, z/norm});
}