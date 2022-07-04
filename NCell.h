// declaration of NCell class

#include <cstddef>
#include <vector>
#include "Cell.h"
#pragma once

using namespace std;

class NCell {

    private:
        double * alts;
        double * dhs;
        double num_L;
        double num_chi;
        double update_period;
        size_t time;
        size_t lupdate_time;
        vector<double> t;
        Cell * cells;
        double num_cells;
        double * logL_edges;
        double * logL_ave;
        double * chi_edges;
        double * chi_ave;
};