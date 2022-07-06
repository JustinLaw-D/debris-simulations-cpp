// Cell class

#include "Arrays.h"
#include "ObjectsEvents.h"
#include <vector>
#pragma once

using namespace std;

class Cell
{
    private:
        vector<Array1D<double> *> * S; // number of live satellites of each type at each time
        vector<Array1D<double> *> * S_d; // number of de-orbiting satellites at each time
        vector<Array1D<double> *> * D; // number of derelict satellites at each time
        Array1D<double> * m_s; // mass of the satellite types (kg)
        Array1D<double> * sigma_s; // cross-section of the satellite types (m^2)
        Array1D<double> * sigma_s_km; // cross-section of the satellite types (km^2)
        Array1D<double> * lam_s; // satellite launch rates (1/yr)
        Array1D<double> * del_t; // live satellite lifetimes (yr)
        Array1D<double> * fail_t; // ascending satellite failure lifetime (yr)
        Array1D<double> * tau_do; // satellite de-orbiting times (yr)
        Array1D<double> * target_alt; // satellite target altitudes (km)
        Array1D<double> * up_time; // time it takes satellites to ascend through the cell (yr)
        // failed avoidance chance for each object type
        Array1D<double> * alphaS; Array1D<double> * alphaD; Array1D<double> * alphaN; Array1D<double> * alphaR;
        Array1D<double> * P; // chance of successful de-orbit
        Array1D<double> * AM_s; // area-to-mass ratio of satellite (m^2/kg)
        Array1D<double> * tau_s; // atmospheric drag lifetime (yr)
        Array1D<double> * C_s; // fit constant for explosions
        Array1D<double> * expl_rate_L; // number of explosions per 100 live satellites per year
        Array1D<double> * expl_rate_D; // number of explosions per 100 derelict satellites per year
        vector<Array1D<double> *> * R; // number of rocket bodies of each type at each given time
        Array1D<double> * m_rb; // mass of the rocket bodies (kg)
        Array1D<double> * sigma_rb; // cross-section of the rocket bodies (m^2)
        Array1D<double> * sigma_rb_km; // cross-section of the rocket bodies (km^2)
        Array1D<double> * lam_rb; // rocket body launch rates (1/yr)
        Array1D<double> * AM_rb; // area-to-mass ratio of the rocket bodies (m^2/kg)
        Array1D<double> * tau_rb; // atmospheric drag lifetimes (yr)
        Array1D<double> * C_rb; // fit constants for explosions
        Array1D<double> * expl_rate_R; // number of explosions per 100 rocket bodies per year
        vector<ArrayND<double, 2> *> * N_bins; // amount of debris of each type
        size_t num_sat_types; // number of satellite types
        size_t num_rb_types; // number of rocket body types
        Array1D<double> * logL_edges; // edges of the length bins, logarithmic
        Array1D<double> * logL_ave; // middle of each bin, logarithmic
        size_t num_L; // number of length bins
        Array1D<double> * chi_edges; // edges of the chi bins
        Array1D<double> * chi_ave; // middle of each bin
        size_t num_chi; // number of chi bins
        vector<Event *> * event_list; // list of events that can occur in the cell
        size_t num_events; // number of events in the event_list
        double alt; // altitude of the shell centre (km)
        double dh; // width of the shell (km)
        ArrayND<double, 2> * tau_N; // decay lifetimes for each debris bin (yr)
        double v; // relative velocity of collisions (km/s)
        double vyr; // relative velocity of collisions (m/s)
        double v_orbit; // orbital velocity of the shell (km/s)
        double V; // volume of the shell (km^3)
        ArrayND<bool, 3> * cat_sat_N; // array of which bins have catestrophic collisions for satellite types
        ArrayND<bool, 3> * cat_rb_N; // array of which bins have catestrophic for rocket types
        Array1D<bool> * ascending; // whether or not each satellite type is ascending
        ArrayND<bool, 2> * trackable; // which bins are trackable
        vector<double> * C_l; // total number of catestrophic collisions in each time step
        vector<double> * C_nl; // total number of non-catestrophic collisions in each time step

    public:
        // full constructor
        Cell(Satellite * satellites, RocketBody * rockets, ArrayND<double,2> * N_i, size_t num_sat_types,
             size_t num_rb_types, Array1D<double> * logL_edges, size_t num_L, Array1D<double> * chi_edges, size_t num_chi,
             vector<Event *> * event_list, size_t num_events, double alt, double dh, ArrayND<double,2> * tau_N, double v,
             vector<double> * C_l, vector<double> * C_nl);
        // calculating local rates of change
        void dxdt_cell(size_t time, Array1D<double> &dSdt, Array1D<double> &dS_ddt, Array1D<double> &dDdt, Array1D<double> &dRdt,
                       double &dC_ldt, double &dC_nldt, Array1D<double> &S_out, Array1D<double> &S_dout, Array1D<double> &D_out, 
                       Array1D<double> &R_out, ArrayND<double,2> &N_out, ArrayND<double,2> &D_dt, ArrayND<double,2> &DR_dt, 
                       ArrayND<double,2> &R_dt, ArrayND<double,3> &CS_dt, ArrayND<double,3> &CR_dt, Array1D<double> &expl_S, 
                       Array1D<double> &expl_R);
        void update_cat_N(); // updates catestrophic debris tables
        ~Cell(); // destructor
};