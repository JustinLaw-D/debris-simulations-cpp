// implementation of Cell.h

#include "Cell.h"
#include "Arrays.h"
#include "BreakupModel.h"
#include "ObjectsEvents.h"
#include "Constants.h"
#include <vector>
#include <cmath>
#include <math.h>

using namespace std;

Cell::Cell(Satellite * satellites, RocketBody * rockets, Array2D<double> & N_i, size_t num_sat_types,
           size_t num_rb_types, double * logL_edges, size_t num_L, double * chi_edges, size_t num_chi,
           Event * event_list, size_t num_events, double alt, double dh, Array2D<double> & tau_N, double v) {
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

    Note(s) : event_list must be a dynamically allocated array
    */
    
    this->num_sat_types = num_sat_types;
    this->num_rb_types = num_rb_types;

    // setup satellites
    
    this->S = new vector<Array2D<double> *>();
    Array2D<double> * S_i = Array2D<double>::zeroes_double_dyn(num_sat_types, 1);
    this->S_d = new vector<Array2D<double> *>();
    Array2D<double> * S_di = Array2D<double>::zeroes_double_dyn(num_sat_types, 1);
    this->D = new vector<Array2D<double> *>();
    Array2D<double> * D_i = Array2D<double>::zeroes_double_dyn(num_sat_types, 1);
    
    this->m_s = Array2D<double>::zeroes_double_dyn(num_sat_types, 1);
    this->sigma_s = Array2D<double>::zeroes_double_dyn(num_sat_types, 1);
    this->sigma_s_km = Array2D<double>::zeroes_double_dyn(num_sat_types, 1);
    this->lam_s = Array2D<double>::zeroes_double_dyn(num_sat_types, 1);
    this->del_t = Array2D<double>::zeroes_double_dyn(num_sat_types, 1);
    this->tau_do = Array2D<double>::zeroes_double_dyn(num_sat_types, 1);
    this->target_alt = Array2D<double>::zeroes_double_dyn(num_sat_types, 1);
    this->up_time = Array2D<double>::zeroes_double_dyn(num_sat_types, 1);
    this->alphaS = Array2D<double>::zeroes_double_dyn(num_sat_types, 1); 
    this->alphaD = Array2D<double>::zeroes_double_dyn(num_sat_types, 1); 
    this->alphaN = Array2D<double>::zeroes_double_dyn(num_sat_types, 1); 
    this->alphaR = Array2D<double>::zeroes_double_dyn(num_sat_types, 1);
    this->P = Array2D<double>::zeroes_double_dyn(num_sat_types, 1);
    this->AM_s = Array2D<double>::zeroes_double_dyn(num_sat_types, 1);
    this->tau_s = Array2D<double>::zeroes_double_dyn(num_sat_types, 1);
    this->C_s = Array2D<double>::zeroes_double_dyn(num_sat_types, 1);
    this->expl_rate_L = Array2D<double>::zeroes_double_dyn(num_sat_types, 1);
    this->expl_rate_D = Array2D<double>::zeroes_double_dyn(num_sat_types, 1);
    
    for (size_t i = 0; i < num_sat_types; i++){
        Satellite sat = satellites[i];
        S_i->set(sat.S[0], i, 0); S_di->set(sat.S_d[0], i, 0); D_i->set(sat.D[0], i, 0);
        (this->m_s)->set(sat.m, i, 0); (this->sigma_s)->set(sat.sigma, i, 0); (this->lam_s)->set(sat.lam, i, 0);
        (this->del_t)->set(sat.del_t, i, 0); (this->tau_do)->set(sat.tau_do, i, 0); 
        (this->target_alt)->set(sat.target_alt, i, 0); (this->up_time)->set(sat.up_time, i, 0);
        (this->alphaS)->set(sat.alphaS, i, 0); (this->alphaD)->set(sat.alphaD, i, 0); (this->alphaN)->set(sat.alphaN, i, 0);
        (this->alphaR)->set(sat.alphaR, i, 0); (this->P)->set(sat.P, i, 0); (this->AM_s)->set(sat.AM, i, 0);
        (this->tau_s)->set(sat.tau, i, 0); (this->C_s)->set(sat.C, i, 0); (this->expl_rate_L)->set(sat.expl_rate_L, i, 0);
        (this->expl_rate_D)->set(sat.expl_rate_D, i, 0); (this->sigma_s_km)->set(sat.sigma/1e6, i, 0);
    }
    (this->S)->push_back(S_i); (this->S_d)->push_back(S_di); (this->D)->push_back(D_i);
    
    // setup rocket bodies
    this->R = new vector<Array2D<double> *>();
    Array2D<double> * R_i = Array2D<double>::zeroes_double_dyn(num_rb_types, 1);
    
    this->m_rb = Array2D<double>::zeroes_double_dyn(num_rb_types, 1);
    this->sigma_rb = Array2D<double>::zeroes_double_dyn(num_rb_types, 1);
    this->sigma_rb_km = Array2D<double>::zeroes_double_dyn(num_rb_types, 1);
    this->lam_rb = Array2D<double>::zeroes_double_dyn(num_rb_types, 1);
    this->AM_rb = Array2D<double>::zeroes_double_dyn(num_rb_types, 1);
    this->tau_rb = Array2D<double>::zeroes_double_dyn(num_rb_types, 1);
    this->C_rb = Array2D<double>::zeroes_double_dyn(num_rb_types, 1);
    this->expl_rate_R = Array2D<double>::zeroes_double_dyn(num_rb_types, 1);
    
    for (size_t i = 0; i < num_rb_types; i++){
        RocketBody rb = rockets[i];
        R_i->set(rb.R[0], i, 0);
        (this->m_rb)->set(rb.m, i, 0); (this->sigma_rb)->set(rb.sigma, i, 0); (this->lam_rb)->set(rb.lam, i, 0);
        (this->AM_rb)->set(rb.AM, i, 0); (this->tau_rb)->set(rb.tau, i, 0); (this->C_rb)->set(rb.C, i, 0); 
        (this->expl_rate_R)->set(rb.expl_rate, i, 0); (this->sigma_rb_km)->set(rb.sigma/1e6, i, 0);
    }
    (this->R)->push_back(R_i);
    

    // setup basic parameters
    this->N_bins = new vector<Array2D<double> *>(); (this->N_bins)->push_back(new Array2D<double>(N_i));
    
    this->logL_edges = logL_edges; this->num_L = num_L; this->chi_edges = chi_edges; this->num_chi = num_chi;
    this->logL_ave = new double[num_L]; this->chi_ave = new double[num_chi];
    for (size_t i = 0; i < num_L; i++) {
        (this->logL_ave)[i] = ((this->logL_edges)[i] + (this->logL_edges)[i+1])/2;
    } for (size_t i = 0; i < num_chi; i++) {
        (this->chi_ave)[i] = ((this->chi_edges)[i] + (this->chi_edges)[i+1])/2;
    }
    this->event_list = event_list; this->num_events = num_events;
    this->alt = alt; this->dh = dh; this->v = v; this->vyr = v*60.0*60.0*24.0*365.25;
    this->v_orbit = sqrt(G*Me/((Re + alt)*1000));
    this->V = 4*M_PI*(Re + this->alt)*(Re + this->alt)*(this->dh);
    this->tau_N = new Array2D<double>(tau_N);

    // setup other matrices
    
    this->cat_sat_N = new vector<Array2D<bool> *>();
    this->cat_rb_N = new vector<Array2D<bool> *>();
    this->ascending = new bool[this->num_sat_types];
    this->trackable = Array2D<bool>::fill_dyn(false, this->num_L, this->num_chi);
    for (size_t i = 0; i < this->num_sat_types; i++) {
        (this->cat_sat_N)->push_back(Array2D<bool>::fill_dyn(true, this->num_L, this->num_chi));
        if (this->target_alt->get(i, 0) > (this->alt + (this->dh)/2)) {(this->trackable)->set(true, i, 0);}
    }
    for (size_t i = 0; i < this->num_rb_types; i++) {
        (this->cat_rb_N)->push_back(Array2D<bool>::fill_dyn(true, this->num_L, this->num_chi));
    }
    this->update_cat_N();
}

void Cell::update_cat_N() {
    // updates catestrophic table values
    bool is_cat;
    double ave_L;
    double ave_AM;
    for (size_t i = 0; i < this->num_L; i++) {
        ave_L = pow(10, (this->logL_ave)[i]); // get average length in m
        for (size_t j = 0; j < this->num_chi; j++) {
            ave_AM = pow(10, (this->chi_ave)[i]); // get average AM in m^2/kg
            for (size_t k = 0; k < this->num_sat_types; k++) {
                is_cat = is_catastrophic((this->m_s)->get(k, 0), ave_L, ave_AM, this->v); // determine if catestrophic
                (this->cat_sat_N->at(k))->set(is_cat, k, 0); // set the result
            } for (size_t k = 0; k < this->num_rb_types; k++) { // do the same thing for rockets
                is_cat = is_catastrophic((this->m_rb)->get(k, 0), ave_L, ave_AM, this->v); // determine if catestrophic
                (this->cat_rb_N->at(k))->set(is_cat, k, 0); // set the result
            }
        }
    }
}

void Cell::dxdt_cell(size_t time, Array2D<double> &dSdt, Array2D<double> &dS_ddt, Array2D<double> &dDdt, Array2D<double> &dRdt, 
                     Array2D<double> &S_out, Array2D<double> &S_dout, Array2D<double> &D_out, Array2D<double> &R_out,
                     Array2D<double> &N_out, Array2D<double> &D_dt, Array2D<double> &DR_dt, Array2D<double> &R_dt, 
                     Array2D<double> **CS_dt, Array2D<double> **RS_dt, Array2D<double> &expl_S, Array2D<double> &expl_R) {
    /*
    calculates the rate of collisions and decays from each debris bin, the rate
    of decaying/de-orbiting satellites, the rate of launches/deorbit starts of satallites, 
    and the rate of creation of derelicts at the given time, due only to events in the cell

    Inputs(s):
    time : index of the values to use

    Output(s) (by reference):
    dSdt : array of rate of change of the number of live satellites in the cell of each type due to only processes
            withing the cell (yr^(-1))
    dS_ddt : array of rate of change of the number of de-orbiting satellites in the cell of each type (yr^(-1))
    dDdt : array of rate of change of the number of derelict satellites in the cell of each type (yr^(-1))
    dRdt : array of rate of change of number of rocket bodies in the cell of each type (yr^(-1))
    S_out : array of rate of satellites ascending from the cell of each type (yr^(-1))
    S_dout : array of rate of satellites de-orbiting from the cell of each type (yr^(-1))
    D_out : array of rate of satellites decaying from the cell of each type (yr^(-1))
    R_out : array of rate of rocket bodies decaying from the cell of each type (yr^(-1))
    N_out : matrix with the rate of exiting debris from each bin (yr^(-1))
    D_dt : matrix with total rate of collisions between satellites (yr^(-1))
    DR_dt : matrix with total rate of collisions between satellites and rocket bodies (yr^(-1))
    R_dt : matrix with total rate of collisions between rocket bodies (yr^(-1))
    CS_dt : array of pointers to matrices with the rate of collisions from each bin with each satellite type (yr^(-1))
    CR_dt : array of pointers to matrices with the rate of collisions from each bin with each rocket body type (yr^(-1))
    expl_S : array of rate of explosions for satellites of each type (yr^(-1))
    expl_R : array of rate of explosions for rocket bodies of each type (yr^(-1))

    Note: Assumes that collisions with debris of L_cm < 10cm cannot be avoided, and that the given time input is valid
    */

    // get current N value
    Array2D<double> N = *(this->N_bins->back());

    // setup some temp parameters
    double sigma_loc_km0; double R1; double sigma_loc_km1; double sigma_comb; double n;
    
    // start with satellite collisions
    for (size_t i = 0; i < this->num_sat_types; i++) {

        // get current parameters
        double S0 = this->S->at(time)->get(i, 0);
        double SD0 = this->S_d->at(time)->get(i,0);
        double D0 = this->D->at(time)->get(i,0);
        sigma_loc_km0 = this->sigma_s_km->get(i,0);
        double alphaS0 = this->alphaS->get(i, 0);
        double alphaD0 = this->alphaD->get(i, 0);
        double alphaN0 = this->alphaN->get(i, 0);
        double alphaR0 = this->alphaR->get(i, 0);

        // setup temp parameters
        double dSdt_loc; double dS_ddt_loc; double dDdt_loc;
        double dSSdt_loc; double dSS_ddt_loc; double dSDdt_loc;
        double dSRdt_loc; double dS_dS_ddt_loc; double dS_dDdt_loc;
        double dS_dRdt_loc; double dDDdt_loc; double dDRdt_loc;
        double S1; double SD1; double D1; double alphaS1;
        Array2D<double> * CS_dt_loc = CS_dt[i]; // pointer to array for this satellite type

        // handle satellite-debris collisions
        for (size_t j = 0; j < this->num_L; j++) {
            for (size_t k = 0; k < this->num_chi; k++) {
                n = N.get(j,k)/(this->V); // calculate debris density
                if (this->trackable->get(j,k)) { // if avoidance is possible
                    dSdt_loc = alphaN0*n*sigma_loc_km0*(this->vyr)*S0;
                    dS_ddt_loc = alphaN0*n*sigma_loc_km0*(this->vyr)*SD0;
                } else { // if avoidance is impossible
                    dSdt_loc = n*sigma_loc_km0*(this->vyr)*S0;
                    dS_ddt_loc = n*sigma_loc_km0*(this->vyr)*SD0;
                }
                dDdt_loc = n*sigma_loc_km0*(this->vyr)*D0; // never avoids
                // update return values
                dSdt.set(dSdt.get(i,0) - dSdt_loc, i, 0);
                dS_ddt.set(dS_ddt.get(i,0) - dS_ddt_loc, i, 0);
                dDdt.set(dDdt.get(i,0) - dDdt_loc, i, 0);
                CS_dt_loc->set(CS_dt_loc->get(j,k) - dSdt_loc - dS_ddt_loc - dDdt_loc, j, k);
            }
        }

        // handle satellite-satellite collisions
        for (size_t j = 0; j < this->num_sat_types; j++) {
            // get local parameters for this satellite type
            S1 = this->S->at(time)->get(j, 0);
            SD1 = this->S_d->at(time)->get(j,0);
            D1 = this->D->at(time)->get(j,0);
            sigma_loc_km1 = this->sigma_s_km->get(j,0);
            alphaS1 = this->alphaS->get(j, 0);
            
            // calculate combined cross-section
            sigma_comb = sigma_loc_km0 + sigma_loc_km1 + 2*sqrt(sigma_loc_km0*sigma_loc_km1);

            // calculate collisions
            dSSdt_loc = alphaS0*alphaS1*sigma_comb*(this->vyr)*S0*S1/(this->V);
            dSS_ddt_loc = alphaS0*alphaS1*sigma_comb*(this->vyr)*S0*SD1/(this->V);
            dSDdt_loc = alphaD0*sigma_comb*(this->vyr)*S0*D1/(this->V);
            dS_dS_ddt_loc = alphaS0*alphaS1*sigma_comb*(this->vyr)*SD0*SD1/(this->V);
            dS_dDdt_loc = alphaS0*sigma_comb*(this->vyr)*SD0*D1/(this->V);
            dDDdt_loc = sigma_comb*(this->vyr)*D0*D1/(this->V);

            // update return values
            if (i <= j) { // avoid double counting
                D_dt.set(D_dt.get(i,j) + dSSdt_loc + dS_dS_ddt_loc + dDDdt_loc, i, j);
            }
            D_dt.set(D_dt.get(i,j) + dSS_ddt_loc + dSDdt_loc + dS_dDdt_loc, i, j); // these aren't double counted
            if (i == j) { // destroys two of the same type in one go
                dSdt.set(dSdt.get(i,0) - 2*dSSdt_loc - dSS_ddt_loc - dSDdt_loc, i, 0);
                dS_ddt.set(dS_ddt.get(i,0) - dSS_ddt_loc - 2*dS_dS_ddt_loc - dS_dDdt_loc, i, 0);
                dDdt.set(dDdt.get(i,0) - dSDdt_loc - dS_dDdt_loc - 2*dDDdt_loc, i, 0);
            } else { // doesnt, and have to treat first and second satellite types seperately sometimes
                dSdt.set(dSdt.get(i,0) - dSSdt_loc - dSS_ddt_loc - dSDdt_loc, i, 0);
                dS_ddt.set(dS_ddt.get(i,0) - dS_dS_ddt_loc - dS_dDdt_loc, i, 0);
                dDdt.set(dDdt.get(i,0) - dDDdt_loc, i, 0);
                dS_ddt.set(dS_ddt.get(j, 0) - dSS_ddt_loc, j, 0);
                dDdt.set(dDdt.get(j,0) - dSDdt_loc - dS_dDdt_loc, j, 0);
            }
        }

        // handle satellite-rocket collisions
        for (size_t j = 0; j < this->num_rb_types; j++) {
            // get local parameters for this rocket type
            R1 = this->R->at(time)->get(j, 0);
            sigma_loc_km1 = this->sigma_rb_km->get(j,0);
            
            // calculate combined cross-section
            sigma_comb = sigma_loc_km0 + sigma_loc_km1 + 2*sqrt(sigma_loc_km0*sigma_loc_km1);

            // calculate collisions
            dSRdt_loc = alphaR0*sigma_comb*(this->vyr)*S0*R1/(this->V);
            dS_dRdt_loc = alphaR0*sigma_comb*(this->vyr)*SD0*R1/(this->V);
            dDRdt_loc = sigma_comb*(this->vyr)*D0*R1/(this->V);

            // update return values
            DR_dt.set(DR_dt.get(i,j) + dSRdt_loc + dS_dRdt_loc + dDRdt_loc, i, j); // these aren't double counted
            dSdt.set(dSdt.get(i,0) - dSRdt_loc, i, 0);
            dS_ddt.set(dS_ddt.get(i,0) - dS_dRdt_loc, i, 0);
            dDdt.set(dDdt.get(i,0) - dDRdt_loc, i, 0);
            dRdt.set(dRdt.get(j,0) - dSRdt_loc - dS_dRdt_loc - dDRdt_loc, j, 0);
        }
    }
    // TODO ROCKET-ROCKET COLLISIONS, DECAYS AND UP-TIMES
}

Cell::~Cell() {
    
    // basic class destructor
    // Note : the Cell is not assumed to own the logL_edges and chi_edges arrays, and hence
    //        will not free them
    delete this->trackable; delete [] this->ascending; delete this->cat_rb_N;
    delete this->cat_sat_N; delete this->tau_N; delete [] this->event_list;
    delete [] this->logL_ave; delete [] this->chi_ave; delete this->sigma_s_km;
    delete this->N_bins; delete this->R; delete this->expl_rate_R;
    delete this->C_rb; delete this->tau_rb; delete this->AM_rb; delete this->lam_rb;
    delete this->sigma_rb; delete this->m_rb; delete this->D; delete this->S_d;
    delete this->S; delete this->expl_rate_D; delete this->expl_rate_L;
    delete this->C_s; delete this->tau_s; delete this->AM_s; delete this->P;
    delete this->alphaR; delete this->alphaN; delete this->alphaD; delete this->alphaS;
    delete this->up_time; delete this->target_alt; delete this->tau_do; delete this->del_t;
    delete this->lam_s; delete this->sigma_s; delete this->m_s; delete this->sigma_rb_km;
}