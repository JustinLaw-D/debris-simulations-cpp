// file for Object/Event structs and related functions
// TODO add default event functions

#include "Arrays.h"
#include <vector>
#pragma once

using namespace std;

struct Satellite
{
    vector<double> S; // number of live satellites at each time
    vector<double> S_d; // number of de-orbiting satellites at each time
    vector<double> D; // number of derelict satellites at each time
    double m; // mass of the satellite (kg)
    double sigma; // cross-section of the satellite (m^2)
    double lam; // satellite launch rate (1/yr)
    double del_t; // live satellite lifetime (yr)
    double tau_do; // satellite de-orbiting time (yr)
    double target_alt; // satellite target altitude (km)
    double up_time; // time it takes satellite to ascend through the cell (yr)
    double alphaS; double alphaD; double alphaN; double alphaR; // failed avoidance chance for each object type
    double P; // chance of successful de-orbit
    double AM; // area-to-mass ratio of satellite (m^2/kg)
    double tau; // atmospheric drag lifetime (yr)
    double C; // fit constant for explosions
    double expl_rate_L; // number of explosions per 100 live satellites per year
    double expl_rate_D; // number of explosions per 100 derelict satellites per year
};

inline Satellite make_sat(vector<double> S, vector<double> S_d, vector<double> D, double m, double sigma, double lam,
                          double del_t, double tau_do, double target_alt, double up_time, double alphaS, double alphaD,
                          double alphaN, double alphaR, double P, double AM, double tau, double C, double expl_rate_L,
                          double expl_rate_D) {
    /*
    basic constructor for Satellite

    Input(s):
    S : number of live satellites at each time
    S_d : number of de-orbiting satellites at each time
    D : number of derelict satellites at each time
    m : mass of the satellite (kg)
    sigma : cross-section of the satellite (m^2)
    lam : satellite launch rate (1/yr)
    del_t : live satellite lifetime (yr)
    tau_do : satellite de-orbiting time (yr)
    target_alt : satellite target altitude (km)
    up_time : time it takes satellite to ascend through the cell (yr)
    alphaS : failed avoidance chance for live satellites 
    alphaD : failed avoidance chance for derelict satellites 
    alphaN : failed avoidance chance for trackable debris 
    alphaR : failed avoidance chance for rocket bodies
    P : chance of successful de-orbit
    AM : area-to-mass ratio of satellite (m^2/kg)
    tau : atmospheric drag lifetime (yr)
    C : fit constant for explosions
    expl_rate_L : number of explosions per 100 live satellites per year
    expl_rate_D : number of explosions per 100 derelict satellites per year

    Output(s):
    sat : Satellite object
    */
    Satellite sat = {
        S, S_d, D, m, sigma, lam, del_t, tau_do, target_alt, up_time, alphaS, alphaD, alphaN, alphaR,
        P, AM, tau, C, expl_rate_L, expl_rate_D
    };
    return sat;
}

struct RocketBody
{
    vector<double> R; // number of rocket bodies at each given time
    double m; // mass of the rocket body (kg)
    double sigma; // cross-section of the rocket body (m^2)
    double lam; // rocket body launch rate (1/yr)
    double AM; // area-to-mass ratio of the rocket body (m^2/kg)
    double tau; // atmospheric drag lifetime (yr)
    double C; // fit constant for explosions
    double expl_rate; // number of explosions per 100 rocket bodies per year
};

inline RocketBody make_rb(vector<double> R, double m, double sigma, double lam, double AM, double tau, 
                          double C, double expl_rate) {
    /*
    default constructor for rocket bodies

    Input(s):
    R : number of rocket bodies at each given time
    m : mass of the rocket body (kg)
    sigma : cross-section of the rocket body (m^2)
    lam : rocket body launch rate (1/yr)
    AM : area-to-mass ratio of the rocket body (m^2/kg)
    tau : atmospheric drag lifetime (yr)
    C : fit constant for explosions
    expl_rate : number of explosions per 100 rocket bodies per year

    Output:
    rb : RocketBody instance
    */
    RocketBody rb = {R, m, sigma, lam, AM, tau, C, expl_rate};
    return rb;
}

struct Coll
{
    double m1; // mass of first object in the collision (kg)
    double m2; // mass of the second object (kg)
    char typ; // type of collision ('s' for satellite-satellite, 'm' for satellite-rocket, 'r' for rocket-rocket)
    size_t num; // number of collisions of this type
};

inline Coll make_coll(double m1, double m2, char typ, size_t num) {
    /*
    constructor for collision type

    Input(s):
    m1 : mass of first object in the collision (kg)
    m2 : mass of the second object (kg)
    typ : type of collision ('s' for satellite-satellite, 'm' for satellite-rocket, 'r' for rocket-rocket)
    num : number of collisions of this type

    Output(s):
    coll : Coll instance
    */
    Coll coll = {m1, m2, typ, num};
    return coll;
}

struct Expl
{
    double C; // fit constant of the exploding body
    char typ; // type of the exploding body ('s' for satellite, 'r' for rocket)
    size_t num; // number of explosions of this type
};

inline Expl make_expl(double C, char typ, size_t num) {
    /*
    constructor for explosion type

    Input(s):
    C : fit constant of the exploding body
    typ : type of the exploding body ('s' for satellite, 'r' for rocket)
    num : number of explosions of this type

    Output(s):
    expl : Expl instance
    */
    Expl expl = {C, typ, num};
    return expl;
}

/*
type for function used to run events
signature is the following

Input(s):
S : number of live satellites of each type
S_d : number of de-orbiting satellites of each type
D : number of derelict satellites of each type
num_sat_types : the number of satellite types
R : number of rocket bodies of each type
num_rb_types : the number of rocket body types
N : binned number of debris
logL_edges : bin edges in logL
num_L : number of bins in L
chi_edges : bin edges in chi
num_chi : number of bins in chi
dS : change in S due to the event
dS_d : change in S_d due to the event
dD : change in D due to the event
dR : change in R due to the event
dN : change in N due to the event
coll_list : list of collisions occuring in the event
expl_list : list of explosions occuring in the event

Output(s): None
*/
typedef void (*event_func)(const Array2D<double>&, const Array2D<double>&, const Array2D<double>&, size_t, 
                           const Array2D<double>&, size_t, const Array2D<double>&, const double*, size_t, 
                           const double*, size_t, Array2D<double>&, Array2D<double>&, Array2D<double>&, 
                           Array2D<double>&, Array2D<double>&, vector<Coll>&, vector<Expl>&);

struct Event 
{
    vector<double> times; // pointer to the list of times to trigger the event
    double freq; // how frequently the event should trigger, 0 is never
    double alt; // altitude of the event
    event_func run_event; // pointer to function that simulates the event
};

inline Event make_event(vector<double> times, double freq, double alt, event_func run_event) {
    /*
    Event constructor

    Input(s):
    times : list of times to trigger the event (yr)
    freq : frequency of the event (1/yr, 0 is never)
    alt : altitude of the event (km)
    run_event : pointer to function that simulates the event

    Output(s):
    event : Event instance
    */
    Event event = {times, freq, alt, run_event};
    return event;
}

void default_event(const Array2D<double>& S, const Array2D<double>& S_d, const Array2D<double>& D, size_t num_sat_types, 
                   const Array2D<double>& R, size_t num_rb_types, const Array2D<double>& N, const double* logL_bins, 
                   size_t num_L, const double* chi_bins, size_t num_chi, Array2D<double>& dS, Array2D<double>& dS_d, 
                   Array2D<double>& dD, Array2D<double>& dR, Array2D<double>& dN, vector<Coll>& coll_list, 
                   vector<Expl>& expl_list); // dummy default event