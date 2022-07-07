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

    // parse time array
    this->t = new vector<double>();
    ifstream t_file; t_file.open(filepath + string("t.npy"), ios::in | ios::binary);
    if (t_file.is_open()) {
        char waste[8]; // opening bytes of the file that are essentially useless
        unsigned short header_len; // for holding information on the length of the header
        t_file.read(waste, 8); t_file.read((char *) &header_len, 2);
        // header info is useless
        char * header_waste = new char[header_len]; t_file.read(header_waste, header_len); delete [] header_waste;
        double * data = new double[8]; // pull data in 64 byte chuncks
        t_file.read((char *)data, 64);
        while (!t_file.eof()) { // keep going until eof is set
            for (size_t i = 0; i < 8; i++) {t->push_back(data[i]);}
            t_file.read((char *)data, 64);
        }
        streamsize to_write = t_file.gcount()/8; // figure out how many numbers were read on the last call
        for (streamsize i = 0; i < to_write; i++) {t->push_back(data[i]);}
        for (size_t i = 0; i < t->size(); i++) { // temporary for tests
            cout << t->at(i) << endl;
        }
    } else {
        cout << "Ah this one" << endl;
        throw invalid_argument("NCell t.npy file could not be opened");
    }
    
}

NCell::~NCell() {}