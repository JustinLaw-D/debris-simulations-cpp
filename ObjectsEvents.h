// file for Object/Event structs and related functions
// TODO add default event functions

#include "Arrays.h"
#include <vector>
#pragma once

using namespace std;
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
logL_bins : bin edges in logL
num_L : number of bins in L
chi_bins : bin edges in chi
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
                          double expl_rate_D); // constructor for Satellite

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
                          double C, double expl_rate); // constructor for RocketBody

struct Coll
{
    double m1; // mass of first object in the collision (kg)
    double m2; // mass of the second object (kg)
    char typ; // type of collision ('s' for satellite-satellite, 'm' for satellite-rocket, 'r' for rocket-rocket)
    size_t num; // number of collisions of this type
};

inline Coll make_coll(double m1, double m2, char typ, size_t num); // constructor for Coll

struct Expl
{
    double C; // fit constant of the exploding body
    char typ; // type of the exploding body ('s' for satellite, 'r' for rocket)
    size_t num; // number of explosions of this type
};

inline Expl make_expl(double C, char typ, size_t num); // constructor for Expl

struct Event 
{
    vector<double> times; // pointer to the list of times to trigger the event
    double period; // how frequently the event should trigger
    double alt; // altitude of the event
    event_func run_event; // pointer to function that simulates the event
};

inline Event make_event(vector<double> times, double period, double alt, event_func run_event); // Event constructor
