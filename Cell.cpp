// implementation of Cell.h

#include "Cell.h"
#include "Arrays.h"
#include "ObjectsEvents.h"
#include <vector>
#include <cmath>
#pragma once

const double G = 6.6743e-11; // Gravitation constant, standard units
const double Re = 6371; // radius of the Earth (km)
const double Me = 5.97219e24; // mass of Earth (kg)

using namespace std;

Cell::Cell(Satellite * satellites, RocketBody * rockets, Array2D<double> N_i, size_t num_sat_types,
           size_t num_rb_types, double * logL_edges, size_t num_L, double * chi_edges, size_t num_chi,
           Event * event_list, size_t num_events, double alt, double dh, Array2D<double> tau_N, double v) {
    /*
    detailed constructor for Cell class
    
    Input(s):
    satellites : list of satellite types with initial values
    rockets : list of rocket body types with initial values
    N_i : initial array of number of debris by L and A/M
    num_sat_types : number of satellite types
    num_rb_types : number of rocket body types
    logL_edges : bin edges in log10 of characteristic length (log10(m))
    num_L : number of L bins
    chi_edges : bin edges in log10(A/M) (log10(m^2/kg))
    num_chi : number of chi bins
    event_list : list of discrete events that occur in the cell
    num_events : number of discrete event types
    alt : altitude of the shell centre (km)
    dh : width of the shell (km)
    tau_N : array of atmospheric drag lifetimes for debris (yr)
    v : relative collision speed (km/s)

    Output(s):
    Cell instance
    */
    this->num_sat_types = num_sat_types;
    this->num_rb_types = num_rb_types;

    // setup satellites
    this->S = vector<Array2D<double>>();
    Array2D<double> S_i = Array2D<double>::zeroes_double(num_sat_types, 1);
    this->S_d = vector<Array2D<double>>();
    Array2D<double> S_di = Array2D<double>::zeroes_double(num_sat_types, 1);
    this->D = vector<Array2D<double>>();
    Array2D<double> D_i = Array2D<double>::zeroes_double(num_sat_types, 1);
    this->m_s = Array2D<double>::zeroes_double(num_sat_types, 1);
    this->sigma_s = Array2D<double>::zeroes_double(num_sat_types, 1);
    this->lam_s = Array2D<double>::zeroes_double(num_sat_types, 1);
    this->del_t = Array2D<double>::zeroes_double(num_sat_types, 1);
    this->tau_do = Array2D<double>::zeroes_double(num_sat_types, 1);
    this->target_alt = Array2D<double>::zeroes_double(num_sat_types, 1);
    this->up_time = Array2D<double>::zeroes_double(num_sat_types, 1);
    this->alphaS = Array2D<double>::zeroes_double(num_sat_types, 1); 
    this->alphaD = Array2D<double>::zeroes_double(num_sat_types, 1); 
    this->alphaN = Array2D<double>::zeroes_double(num_sat_types, 1); 
    this->alphaR = Array2D<double>::zeroes_double(num_sat_types, 1);
    this->P = Array2D<double>::zeroes_double(num_sat_types, 1);
    this->AM_s = Array2D<double>::zeroes_double(num_sat_types, 1);
    this->tau_s = Array2D<double>::zeroes_double(num_sat_types, 1);
    this->C_s = Array2D<double>::zeroes_double(num_sat_types, 1);
    this->expl_rate_L = Array2D<double>::zeroes_double(num_sat_types, 1);
    this->expl_rate_D = Array2D<double>::zeroes_double(num_sat_types, 1);
    for (size_t i = 0; i < num_sat_types; i++){
        Satellite sat = satellites[i];
        S_i.set(sat.S[0], i, 0); S_di.set(sat.S_d[0], i, 0); D_i.set(sat.D[0], i, 0);
        (this->m_s).set(sat.m, i, 0); (this->sigma_s).set(sat.sigma, i, 0); (this->lam_s).set(sat.lam, i, 0);
        (this->del_t).set(sat.del_t, i, 0); (this->tau_do).set(sat.tau_do, i, 0); 
        (this->target_alt).set(sat.target_alt, i, 0); (this->up_time).set(sat.up_time, i, 0);
        (this->alphaS).set(sat.alphaS, i, 0); (this->alphaD).set(sat.alphaD, i, 0); (this->alphaN).set(sat.alphaN, i, 0);
        (this->alphaR).set(sat.alphaR, i, 0); (this->P).set(sat.P, i, 0); (this->AM_s).set(sat.AM, i, 0);
        (this->tau_s).set(sat.tau, i, 0); (this->C_s).set(sat.C, i, 0); (this->expl_rate_L).set(sat.expl_rate_L, i, 0);
        (this->expl_rate_D).set(sat.expl_rate_D, i, 0);
    }
    (this->S).push_back(S_i); (this->S_d).push_back(S_di); (this->D).push_back(D_i);

    // setup rocket bodies
    this->R = vector<Array2D<double>>();
    Array2D<double> R_i = Array2D<double>::zeroes_double(num_sat_types, 1);
    this->m_rb = Array2D<double>::zeroes_double(num_sat_types, 1);
    this->sigma_rb = Array2D<double>::zeroes_double(num_sat_types, 1);
    this->lam_rb = Array2D<double>::zeroes_double(num_sat_types, 1);
    this->AM_rb = Array2D<double>::zeroes_double(num_sat_types, 1);
    this->tau_rb = Array2D<double>::zeroes_double(num_sat_types, 1);
    this->C_rb = Array2D<double>::zeroes_double(num_sat_types, 1);
    this->expl_rate_R = Array2D<double>::zeroes_double(num_sat_types, 1);
    for (size_t i = 0; i < num_rb_types; i++){
        RocketBody rb = rockets[i];
        R_i.set(rb.R[0], i, 0);
        (this->m_rb).set(rb.m, i, 0); (this->sigma_rb).set(rb.sigma, i, 0); (this->lam_rb).set(rb.lam, i, 0);
        (this->AM_rb).set(rb.AM, i, 0); (this->tau_rb).set(rb.tau, i, 0); (this->C_rb).set(rb.C, i, 0); 
        (this->expl_rate_R).set(rb.expl_rate, i, 0);
    }
    (this->R).push_back(R_i);

    // setup basic parameters
    this->N_bins = vector<Array2D<double>>(); (this->N_bins).push_back(N_i);
    this->logL_edges = logL_edges; this->num_L = num_L; this->chi_edges = chi_edges; this->num_chi = num_chi;
    this->event_list = event_list; this->num_events = num_events;
    this->alt = alt; this->dh = dh; this->v = v;
    this->v_orbit = sqrt(G*Me/((Re + alt)*1000));
    this->tau_N = tau_N;

    // setup other matrices
    this->lethal_sat_N = new Array2D<bool>[this->num_sat_types];
    this->lethal_rb_N = new Array2D<bool>[this->num_rb_types];
    this->ascending = new bool[this->num_sat_types];
    this->trackable = Array2D<bool>::fill(false, this->num_L, this->num_chi);
    for (size_t i = 0; i < this->num_sat_types; i++) {
        (this->lethal_sat_N)[i] = Array2D<bool>::fill(true, this->num_L, this->num_chi);
        if (this->target_alt.get(i, 0) > (this->alt + (this->dh)/2)) {(this->trackable).set(true, i, 0);}
    }
    for (size_t i = 0; i < this->num_rb_types; i++) {
        (this->lethal_rb_N)[i] = Array2D<bool>::fill(true, this->num_L, this->num_chi);
    }
    this->update_lethal_N();
}

void Cell::update_lethal_N() {return;} // TODO

Cell::~Cell() {
    // basic class destructor
    // Note : the Cell is not assumed to own the logL_edges and chi_edges arrays, and hence
    //        will not free them
    delete [] this->ascending; delete [] this->lethal_rb_N; 
    delete [] this->lethal_sat_N; delete [] this->event_list;
}