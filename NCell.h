// declaration of NCell class

#include <cstddef>
#include <vector>
#include <iostream>
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
        vector<double> t;
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
                  Array1D<double> &dRdt, ArrayND<double,2> &dNdt, double &dC_ldt, double & dC_nldt); // calculate rates of change
        void run_sim_euler(double T, double dt, bool upper); // run simulation with euler method
        // run simulation using predictor-corrector method
        void run_sim_precor(double T, double dt_i, double dt_min, double dt_max, double tolerance, bool upper);
};