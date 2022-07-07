// declaration of NCell class

#include <cstddef>
#include <vector>
#include <iostream>
#include "ObjectsEvents.h"
#include "Cell.h"
#include "Arrays.h"
#include "AtmosphericDecayModels.h"
#include "BreakupModel.h"
#pragma once

using namespace std;

class NCell {

    private:
        Array1D<double> * alts;
        Array1D<double> * dhs;
        size_t num_L;
        size_t num_chi;
        double update_period;
        size_t time;
        size_t lupdate_time;
        vector<double> * t;
        Cell * cells;
        size_t num_cells;
        Array1D<double> * logL_edges;
        Array1D<double> * logL_ave;
        Array1D<double> * chi_edges;
        Array1D<double> * chi_ave;
        ArrayND<double, 4> * sat_coll_prob_tables;
        ArrayND<double, 4> * rb_coll_prob_tables;
        ArrayND<double, 4> * sat_expl_prob_tables;
        ArrayND<double, 4> * rb_expl_prob_tables;

    public:
        NCell(string &filepath); // constructor based on loading data from file
        void save(string &filepath, string &name, double gap); // function for saving data
        void add_event(Event * event, double alt); // adds event to the system
        void dxdt(size_t time, bool upper, Array1D<double> &dSdt, Array1D<double> &dS_ddt, Array1D<double> &dDdt,
                  Array1D<double> &dRdt, ArrayND<double,2> &dNdt, double &dC_ldt, double &dC_nldt); // calculate rates of change
        void run_sim_euler(double T, double dt, bool upper); // run simulation with euler method
        // run simulation using predictor-corrector method
        void run_sim_precor(double T, double dt_i, double dt_min, double dt_max, double tolerance, bool upper);
        void sim_colls(ArrayND<double,3> &dNdt, double rate, double m1, double m2, size_t index, char typ); // simulate collisions
        // same function for collisions between a satellite and rocket body
        void sim_colls_satrb(ArrayND<double,3> &dNdt, double rate, double m, size_t index, char typ);
        void sim_expl(ArrayND<double,3> &dNdt, double rate, double C, size_t index, char typ); // simulate explosions
        void sim_events(); // runs event handling for the system
        void parse_coll(vector<Coll> &coll_list, ArrayND<double,3> &dN, size_t index); // handles list of collision objects
        void parse_expl(vector<Expl> &expl_list, ArrayND<double,3> &dN, size_t index); // handles list of explosion objects
        void update_lifetimes(double t); // updates all the atmospheric lifetimes in the system
        void alt_to_index(double alt); // converts a given altitude to the index of the corresponding cell4
        ~NCell(); // destructor
};