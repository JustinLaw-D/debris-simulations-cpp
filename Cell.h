// Cell class

#include "Arrays.h"
#include "ObjectsEvents.h"
#include <vector>
#pragma once


using namespace std;

class Cell
{
    private:
        vector<Array2D<double> *> * S; // number of live satellites of each type at each time
        vector<Array2D<double> *> * S_d; // number of de-orbiting satellites at each time
        vector<Array2D<double> *> * D; // number of derelict satellites at each time
        Array2D<double> * m_s; // mass of the satellite types (kg)
        Array2D<double> * sigma_s; // cross-section of the satellite types (m^2)
        Array2D<double> * lam_s; // satellite launch rates (1/yr)
        Array2D<double> * del_t; // live satellite lifetimes (yr)
        Array2D<double> * tau_do; // satellite de-orbiting times (yr)
        Array2D<double> * target_alt; // satellite target altitudes (km)
        Array2D<double> * up_time; // time it takes satellites to ascend through the cell (yr)
        // failed avoidance chance for each object type
        Array2D<double> * alphaS; Array2D<double> * alphaD; Array2D<double> * alphaN; Array2D<double> * alphaR;
        Array2D<double> * P; // chance of successful de-orbit
        Array2D<double> * AM_s; // area-to-mass ratio of satellite (m^2/kg)
        Array2D<double> * tau_s; // atmospheric drag lifetime (yr)
        Array2D<double> * C_s; // fit constant for explosions
        Array2D<double> * expl_rate_L; // number of explosions per 100 live satellites per year
        Array2D<double> * expl_rate_D; // number of explosions per 100 derelict satellites per year
        vector<Array2D<double> *> * R; // number of rocket bodies of each type at each given time
        Array2D<double> * m_rb; // mass of the rocket bodies (kg)
        Array2D<double> * sigma_rb; // cross-section of the rocket bodies (m^2)
        Array2D<double> * lam_rb; // rocket body launch rates (1/yr)
        Array2D<double> * AM_rb; // area-to-mass ratio of the rocket bodies (m^2/kg)
        Array2D<double> * tau_rb; // atmospheric drag lifetimes (yr)
        Array2D<double> * C_rb; // fit constants for explosions
        Array2D<double> * expl_rate_R; // number of explosions per 100 rocket bodies per year
        vector<Array2D<double> *> * N_bins; // amount of debris of each type
        size_t num_sat_types; // number of satellite types
        size_t num_rb_types; // number of rocket body types
        double * logL_edges; // edges of the length bins, logarithmic
        size_t num_L; // number of length bins
        double * chi_edges; // edges of the chi bins
        size_t num_chi; // number of chi bins
        Event * event_list; // list of events that can occur in the cell
        size_t num_events; // number of events in the event_list
        double alt; // altitude of the shell centre (km)
        double dh; // width of the shell (km)
        Array2D<double> * tau_N; // decay lifetimes for each debris bin (yr)
        double v; // relative velocity of collisions (km/s)
        double v_orbit; // orbital velocity of the shell (km/s)
        vector<Array2D<bool> *> * lethal_sat_N; // array of which bins have lethal collisions for satellite types
        vector<Array2D<bool> *> * lethal_rb_N; // array of which bins have lethal collisions for rocket types
        bool * ascending; // whether or not each satellite type is ascending
        Array2D<bool> * trackable; // which bins are trackable

    public:
        // full constructor
        Cell(Satellite * satellites, RocketBody * rockets, Array2D<double> & N_i, size_t num_sat_types,
             size_t num_rb_types, double * logL_edges, size_t num_L, double * chi_edges, size_t num_chi,
             Event * event_list, size_t num_events, double alt, double dh, Array2D<double> & tau_N, double v);
        // calculating local rates of change
        void dxdt(size_t time, Array2D<double> &dSdt, Array2D<double> &dS_ddt, Array2D<double> &dDdt, Array2D<double> &dRdt, 
                  Array2D<double> &S_out, Array2D<double> &S_dout, Array2D<double> &D_out, Array2D<double> &R_out,
                  Array2D<double> &N_out, Array2D<double> &D_dt, Array2D<double> &RD_dt, Array2D<double> &R_dt, 
                  Array2D<double> &CS_dt, Array2D<double> &RS_dt, Array2D<double> &expl_S, Array2D<double> &expl_R);
        void update_lethal_N(); // updates lethal debris tables
        ~Cell(); // destructor
};