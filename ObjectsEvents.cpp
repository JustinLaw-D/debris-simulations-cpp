// implementation of ObjectsEvents.h

#include "ObjectsEvents.h"
#include "Arrays.h"
#include <vector>
#pragma once

using namespace std;

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

inline Event make_event(vector<double> times, double period, double alt, event_func run_event) {
    /*
    Event constructor

    Input(s):
    times : list of times to trigger the event
    period : how frequently the event should trigger
    alt : altitude of the event
    run_event : pointer to function that simulates the event

    Output(s):
    event : Event instance
    */
    Event event = {times, period, alt, run_event};
    return event;
}