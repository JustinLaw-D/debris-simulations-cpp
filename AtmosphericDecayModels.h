// contains models for calculating atmospheric decay lifetimes

#include <cstddef>

#pragma once

double density(const double alt, const double t, const double mo0); // function for calculating atmospheric density

double dadt(const double alt, const double t, const double m0, const double a_over_m, const double CD); // computes rate of decay

// calculates drag lifetime
double drag_lifetime(const double alt_i, const double alt_f, const double a_over_m, const double CD, double dt, 
                     const double m0, const double mindt, const double maxdt, const double dtfactor, const double tmax);
                    
// same thing but with default values
double drag_lifetime_default(const double alt_i, const double alt_f, const double a_over_m, const double m0);