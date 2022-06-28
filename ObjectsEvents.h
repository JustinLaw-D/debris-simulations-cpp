// file for Object/Event structs and related functions

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
    double P; // change of successful de-orbit
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

struct Event 
{
    double *times; // pointer to the list of times to trigger the event
    unsigned int num_times; // length of the list
    double period; // how frequently the event should trigger
    double alt; // altitude of the event
};
