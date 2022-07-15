// implementation of functions in AtmosphericDecayModels.h

#include "Constants.h"
#include <cmath>
#include <iostream>

using namespace std;

double density(const double alt, double t, const double mo0) {
    /*
    calculates atmospheric density using CIRA-2012 data

    Input(s):
    alt : altitude (km)
    t : time since arbitrary start point (yr)
    m0 : starting month in the solar cycle

    Output(s):
    rho : atmospheric density at the given altitude and time (kg/m^3)
    */

    size_t i = static_cast<size_t>((alt-100.0)/20.0); // calculate index of the altitude bin
    if (i > num_alts - 2) {i = num_alts - 2;}
    if (i < 0) {i = 0;}

    double logalt = log10(alt) + 3; // log(alt) in m
    double mo_frac = t*12 + mo0; // fractional month index
    double mo_int; // will contain the integer part of mo_frac
    double frac_part = modf(mo_frac, &mo_int); // split the two parts
    double mo = static_cast<double>(static_cast<size_t>(mo_int) % num_months) + frac_part; // calculate month (normalized by cycle)
    size_t moID = static_cast<size_t>(mo); // get integer part of month

    size_t moID1 = moID + 1; // get next month bin
    if (moID1 > num_months-1) {moID1 = 0;} 
    double F107 = f107_mo[moID] + (f107_mo[moID1]-f107_mo[moID])*(mo-moID); // get interpolated F107

    double rho = 0; // atmospheric density

    if (F107 <= 65) { // interpolate to get density value
        rho = pow(10, (logdenL[i]+(logdenL[i+1]-logdenL[i])/(logz[i+1]-logz[i])*(logalt-logz[i])));
    } else if (F107 <= 140) {
        double d0 = pow(10, (logdenL[i]+(logdenL[i+1]-logdenL[i])/(logz[i+1]-logz[i])*(logalt-logz[i])));
        double d1 = pow(10, (logdenM[i]+(logdenM[i+1]-logdenM[i])/(logz[i+1]-logz[i])*(logalt-logz[i])));
        rho = d0 + (d1-d0)*(F107-65.0)/75.0;
    } else if (F107 <= 250) {
        double d0 = pow(10, (logdenM[i]+(logdenM[i+1]-logdenM[i])/(logz[i+1]-logz[i])*(logalt-logz[i])));
        double d1 = pow(10, (logdenHL[i]+(logdenHL[i+1]-logdenHL[i])/(logz[i+1]-logz[i])*(logalt-logz[i])));
        rho = d0 + (d1-d0)*(F107-140.0)/110.0;
    } else {
        rho = pow(10, (logdenHL[i]+(logdenHL[i+1]-logdenHL[i])/(logz[i+1]-logz[i])*(logalt-logz[i])));
    }

    return rho;
}

double dadt(const double alt, const double t, const double m0, const double a_over_m, const double CD) {
    /*
    calculates the rate of change in the altitude of a circular orbit

    Input(s):
    alt : altitude of the orbit (km)
    t : time passed since the start of the solar cycle (yr)
    m0 : starting month in the solar cycle
    a_over_m : area-to-mass ratio of the object (m^2/kg)
    CD : drag coefficient of the object

    Output(s):
    dadt value (km/yr)
    */

    return -1.0*(CD*density(alt, t, m0)*a_over_m*sqrt(G*Me*(alt + Re)*1e3))*60.0*60.0*24.0*365.25*1e-3;
}

double drag_lifetime(const double alt_i, const double alt_f, const double a_over_m, const double CD, double dt, 
                     const double m0, const double mindt, const double maxdt, const double dtfactor, const double tmax) {
    /*
    estimates the drag lifetime of an object at altitude alt_i to degrade to altitude alt_f

    Input(s):
    alt_i : initial altitude of the object (km)
    alt_f : desired final altitude of the object (km)
    a_over_m : area-to-mass ratio of the object (m^2/kg)
    CD : drag coefficient of the object
    dt : initial time step of the integration (yr)
    m0 : starting month in the solar cycle
    mindt : minimum time step for integration (yr)
    maxdt : maximum time step of the integration (yr)
    dtfactor : fraction of altitude/rate of change to take as dt
    tmax : maximum time to search to (yr)

    Output(s):
    tau : drag lifetime, possibly infinite (yr)

    Note(s): outputs tmax if the computed value is larger than tmax
    */

    // initialize variables
    double time = 0.0;
    double alt = alt_i;
    double dadt0 = 0.0;
    double alt1 = 0.0;
    double dadt1 = 0.0;
    double ave_dadt = 0.0;

    // integrate using predictor-corrector method
    while (alt > alt_f) {
        dadt0 = dadt(alt, time, m0, a_over_m, CD);
        alt1 = alt + dadt0*dt;
        dadt1 = dadt(alt1, time + dt, m0, a_over_m, CD);
        ave_dadt = (dadt0 + dadt1)/2;
        alt += ave_dadt*dt;
        time += dt;
        dt = -(alt/ave_dadt)*dtfactor;
        if (dt < mindt) {
            cout << "WARNING: Problem is possibly too stiff for integrator." << endl;
            dt = mindt;
        }
        else if (dt > maxdt) {dt = maxdt;}
        if (time > tmax) {return tmax;} // give up
    }

    return time;
}

double drag_lifetime_default(const double alt_i, const double alt_f, const double a_over_m, const double m0) {
    /*
    estimates the drag lifetime of an object at altitude alt_i to degrade to altitude alt_f, assuming that
    CD = 2.2, dt=100 secdons, m0=0, mindt=0, maxdt=100, dtfactor=1/100, tmax=1000

    Input(s):
    alt_i : initial altitude of the object (km)
    alt_f : desired final altitude of the object (km)
    a_over_m : area-to-mass ratio of the object (m^2/kg)
    m0 : starting month in the solar cycle

    Output(s):
    tau : drag lifetime (yr)

    Note(s): outputs tmax if the computed value is larger than tmax
    */
    return drag_lifetime(alt_i, alt_f, a_over_m, 2.2, 100.0/(60.0*60.0*24.0*365.25), m0, 0.0, 100.0, 1.0/100.0, 1000.0);
}