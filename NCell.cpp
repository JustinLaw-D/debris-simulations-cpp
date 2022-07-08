// implementation of NCell.h

#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "Arrays.h"
#include "NCell.h"
#include "BreakupModel.h"
#include "AtmosphericDecayModels.h"
#include "Cell.h"

using namespace std;

NCell::NCell(string &filepath) {
    /*
    loads NCell object from saved data

    Input(s):
    filepath : string containing relative or absolute path to saved data

    Output(s):
    NCell instance

    Note(s): the intention is to create objects in Python, then pass them over to C++
             via this saving/reading method
    */
    ifstream param_file; // file containing basic NCell parameters
    param_file.open(filepath + string("params.csv"), ios::in); // open relevant file
    if (param_file.is_open()) {
        vector<string> row; string line; string word; // extract the data
        while(getline(param_file, line)) {
            stringstream str(line);
            while(getline(str, word, ',')) {
                word.erase(remove(word.begin(), word.end(), '"'), word.end());
                row.push_back(word);
            }
        }
        istringstream num_L_temp(row[0]); // need to make these into string streams to convert
        istringstream num_chi_temp(row[1]);
        istringstream num_cell_temp(row[2]);
        num_L_temp >> this->num_L; num_chi_temp >> this->num_chi; num_cell_temp >> this->num_cells;
        this->update_period = stod(row[3]);
    } else {
        throw invalid_argument("NCell params.csv file could not be opened");
    }
    param_file.close();

    // pull non-cell based arrays and data
    this->t = load_vec<double>(filepath + string("t.npy"));
    this->alts = new Array1D<double>(filepath + string("alts.npy"));
    this->dhs = new Array1D<double>(filepath + string("dhs.npy"));
    this->logL_edges = new Array1D<double>(filepath + string("logL.npy"));
    this->chi_edges = new Array1D<double>(filepath + string("chi.npy"));
    this->time = t->size() - 1; this->lupdate_time = this->time;
    this->logL_ave = new Array1D<double>(this->num_L);
    this->chi_ave = new Array1D<double>(this->num_chi);
    for (size_t i = 0; i < this->num_L; i++) { // setup average arrays
        this->logL_ave->at(i) = (this->logL_edges->at(i) + this->logL_edges->at(i+1))/2;
    } 
    for (size_t i = 0; i < this->num_chi; i++) {
        this->chi_ave->at(i) = (this->chi_edges->at(i) + this->chi_edges->at(i+1))/2;
    }
    for (size_t i = 0; i < this->num_cells; i++) {
        this->cells->push_back(new Cell(filepath + "cell" + to_string(i) + "/"));
    }
    // TODO: ADD PROBABILITY TABLES, UPDATE DRAG LIFETIMES
}

NCell::~NCell() {}