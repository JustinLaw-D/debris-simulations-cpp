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
    double fail_t; // failure lifetime of ascending satellites (yr)
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

struct Coll
{
    double m1; // mass of first object in the collision (kg)
    double m2; // mass of the second object (kg)
    char typ; // type of collision ('s' for satellite-satellite, 'm' for satellite-rocket, 'r' for rocket-rocket)
    size_t num; // number of collisions of this type

    Coll();
    Coll(const Coll &coll);
};

struct Expl
{
    double C; // fit constant of the exploding body
    char typ; // type of the exploding body ('s' for satellite, 'r' for rocket)
    size_t num; // number of explosions of this type

    Expl();
    Expl(const Expl &expl); // copy constructor
};

/*
type for function used to run events
signature is the following

Input(s):
S : number of live satellites of each type
S_d : number of de-orbiting satellites of each type
D : number of derelict satellites of each type
R : number of rocket bodies of each type
N : binned number of debris
logL_edges : bin edges in logL
chi_edges : bin edges in chi
dS : change in S due to the event
dS_d : change in S_d due to the event
dD : change in D due to the event
dR : change in R due to the event
dN : change in N due to the event
coll_list : list of collisions occuring in the event
expl_list : list of explosions occuring in the event

Output(s): None
*/
typedef void (*event_func)(const Array1D<double>&, const Array1D<double>&, const Array1D<double>&, 
                           const Array1D<double>&, const ArrayND<double,2>&, const Array1D<double>&, 
                           const Array1D<double>&, Array1D<double>&, Array1D<double>&, Array1D<double>&, 
                           Array1D<double>&, ArrayND<double,2>&, vector<Coll>&, vector<Expl>&);

struct Event 
{
    vector<double> * times; // pointer to the list of times to trigger the event, must be in reverse chronological order
    double freq; // how frequently the event should trigger, 0 is never
    double last_event; // last time the event triggered, frequency-based
    double alt; // altitude of the event
    event_func run_event; // pointer to function that simulates the event

    Event(const Event &event); // copy constructor
    ~Event(); // destructor
};

void default_event(const Array1D<double>& S, const Array1D<double>& S_d, const Array1D<double>& D, const Array1D<double>& R, 
                   const ArrayND<double,2>& N, const Array1D<double>& logL_bins, const Array1D<double>& chi_bins,
                   Array1D<double>& dS, Array1D<double>& dS_d, Array1D<double>& dD, Array1D<double>& dR, ArrayND<double,2>& dN,
                   vector<Coll>& coll_list, vector<Expl>& expl_list); // dummy default event function