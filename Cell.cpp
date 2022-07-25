// implementation of Cell.h

#include "Cells.h"
#include "Arrays.h"
#include "BreakupModel.h"
#include "ObjectsEvents.h"
#include "Constants.h"
#include <vector>
#include <cmath>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <experimental/filesystem>

using namespace std;
namespace fs = std::experimental::filesystem;

Cell::Cell(Satellite * satellites, RocketBody * rockets, ArrayND<double,2> * N_i, size_t num_sat_types,
           size_t num_rb_types, Array1D<double> * logL_edges, size_t num_L, Array1D<double> * chi_edges, size_t num_chi,
           vector<Event *> * event_list, size_t num_events, double alt, double dh, Array1D<double> * tau_N, double v,
           vector<double> * C_l, vector<double> * C_nl) {
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
    C_l : total number of catestrophic collisions in each time step
    C_nl : total number of non-catestrophic collisions in each time step

    Output(s):
    Cell instance

    Note(s): all pointers must be to dynamically allocated objects.
    */
    
    this->num_sat_types = num_sat_types;
    this->num_rb_types = num_rb_types;

    // setup satellites
    
    this->S = new vector<Array1D<double> *>();
    Array1D<double> * S_i = new Array1D<double>(0.0, num_sat_types);
    this->S_d = new vector<Array1D<double> *>();
    Array1D<double> * S_di = new Array1D<double>(0.0, num_sat_types);
    this->D = new vector<Array1D<double> *>();
    Array1D<double> * D_i = new Array1D<double>(0.0, num_sat_types);
    
    this->m_s = new Array1D<double>(0.0, num_sat_types);
    this->sigma_s = new Array1D<double>(0.0, num_sat_types);
    this->sigma_s_km = new Array1D<double>(0.0, num_sat_types);
    this->lam_s = new Array1D<double>(0.0, num_sat_types);
    this->del_t = new Array1D<double>(0.0, num_sat_types);
    this->fail_t = new Array1D<double>(0.0, num_sat_types);
    this->tau_do = new Array1D<double>(0.0, num_sat_types);
    this->target_alt = new Array1D<double>(0.0, num_sat_types);
    this->up_time = new Array1D<double>(0.0, num_sat_types);
    this->alphaS = new Array1D<double>(0.0, num_sat_types);
    this->alphaD = new Array1D<double>(0.0, num_sat_types); 
    this->alphaN = new Array1D<double>(0.0, num_sat_types); 
    this->alphaR = new Array1D<double>(0.0, num_sat_types);
    this->P = new Array1D<double>(0.0, num_sat_types);
    this->AM_s = new Array1D<double>(0.0, num_sat_types);
    this->tau_s = new Array1D<double>(0.0, num_sat_types);
    this->C_s = new Array1D<double>(0.0, num_sat_types);
    this->expl_rate_L = new Array1D<double>(0.0, num_sat_types);
    this->expl_rate_D = new Array1D<double>(0.0, num_sat_types);
    
    for (size_t i = 0; i < num_sat_types; i++){
        Satellite sat = satellites[i];
        S_i->at(i) = sat.S[0]; S_di->at(i) = sat.S_d[0]; D_i->at(i) = sat.D[0]; (this->m_s)->at(i) = sat.m;
        (this->sigma_s)->at(i) = sat.sigma; (this->lam_s)->at(i) = sat.lam; (this->del_t)->at(i) = sat.del_t;
        (this->fail_t)->at(i) = sat.fail_t; (this->tau_do)->at(i) = sat.tau_do; (this->target_alt)->at(i) = sat.target_alt;
        (this->up_time)->at(i) = sat.up_time; (this->alphaS)->at(i) = sat.alphaS; (this->alphaD)->at(i) = sat.alphaD;
        (this->alphaN)->at(i) = sat.alphaN; (this->alphaR)->at(i) = sat.alphaR; (this->P)->at(i) = sat.P;
        (this->AM_s)->at(i) = sat.AM; (this->tau_s)->at(i) = sat.tau; (this->C_s)->at(i) = sat.C;
        (this->expl_rate_L)->at(i) = sat.expl_rate_L; (this->expl_rate_D)->at(i) = sat.expl_rate_D;
        (this->sigma_s_km)->at(i) = sat.sigma/1e6;
    }
    (this->S)->push_back(S_i); (this->S_d)->push_back(S_di); (this->D)->push_back(D_i);
    
    // setup rocket bodies
    this->R = new vector<Array1D<double> *>();
    Array1D<double> * R_i = new Array1D<double>(0.0, num_rb_types);
    
    this->m_rb = new Array1D<double>(0.0, num_rb_types);
    this->sigma_rb = new Array1D<double>(0.0, num_rb_types);
    this->sigma_rb_km = new Array1D<double>(0.0, num_rb_types);
    this->lam_rb = new Array1D<double>(0.0, num_rb_types);
    this->AM_rb = new Array1D<double>(0.0, num_rb_types);
    this->tau_rb = new Array1D<double>(0.0, num_rb_types);
    this->C_rb = new Array1D<double>(0.0, num_rb_types);
    this->expl_rate_R = new Array1D<double>(0.0, num_rb_types);
    
    for (size_t i = 0; i < num_rb_types; i++){
        RocketBody rb = rockets[i];
        R_i->at(i) = rb.R[0]; (this->m_rb)->at(i) = rb.m; (this->sigma_rb)->at(i) = rb.sigma;
        (this->lam_rb)->at(i) = rb.lam; (this->AM_rb)->at(i) = rb.AM; (this->tau_rb)->at(i) = rb.tau;
        (this->C_rb)->at(i) = rb.C; (this->expl_rate_R)->at(i) = rb.expl_rate;
        (this->sigma_rb_km)->at(i) = rb.sigma/1e6;
    }
    (this->R)->push_back(R_i);

    // setup basic parameters
    this->N_bins = new vector<ArrayND<double, 2> *>(); this->N_bins->push_back(N_i);
    
    this->logL_edges = logL_edges; this->num_L = num_L; this->chi_edges = chi_edges; this->num_chi = num_chi;
    this->logL_ave = new Array1D<double>(0.0, num_L); this->chi_ave = new Array1D<double>(0.0, num_chi);
    for (size_t i = 0; i < num_L; i++) {
        (this->logL_ave)->at(i) = ((this->logL_edges)->at(i) + (this->logL_edges)->at(i+1))/2;
    } for (size_t i = 0; i < num_chi; i++) {
        (this->chi_ave)->at(i) = ((this->chi_edges)->at(i) + (this->chi_edges)->at(i+1))/2;
    }
    this->event_list = event_list; this->num_events = num_events;
    this->alt = alt; this->dh = dh; this->v = v; this->vyr = v*60.0*60.0*24.0*365.25;
    this->v_orbit = sqrt(G*Me/((Re + alt)*1000));
    this->V = 4*M_PI*(Re + this->alt)*(Re + this->alt)*(this->dh);
    this->tau_N = tau_N; this->C_l = C_l; this->C_nl = C_nl;

    // setup other matrices
    
    this->cat_sat_N = new ArrayND<bool,3>(true, array<size_t,3>({this->num_sat_types, this->num_L, this->num_chi}));
    this->cat_rb_N = new ArrayND<bool,3>(true, array<size_t,3>({this->num_rb_types, this->num_L, this->num_chi}));
    this->ascending = new Array1D<bool>(false, this->num_sat_types);
    this->trackable = new ArrayND<bool,2>(false, array<size_t,2>({this->num_L, this->num_chi}));
    for (size_t i = 0; i < this->num_sat_types; i++) {
        if (this->target_alt->at(i) > (this->alt + (this->dh)/2)) {(this->ascending)->at(i) = true;}
    }
    for (size_t i = 0; i < this->num_L; i++) {
        if (pow(10, logL_ave->at(i)) >= 1.0/10.0) {
            for (size_t j = 0; j < this->num_chi; j++) {
                (this->trackable)->at(array<size_t,2>({i,j})) = true;
            }
        }
    }
    this->update_cat_N();
}

Cell::Cell(const string &filepath) {
    /*
    loads Cell object from saved data

    Input(s):
    filepath : string containing relative or absolute path to saved data

    Output(s):
    NCell instance

    Note(s): the intention is to create objects in Python, then pass them over to C++
             via this saving/reading method
    */
    ifstream param_file; // file containing basic Cell parameters
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
        istringstream num_sat_types_temp(row[0]); // need to make these into string streams to convert
        istringstream num_rb_types_temp(row[1]);
        istringstream num_L_temp(row[5]);
        istringstream num_chi_temp(row[6]);
        num_L_temp >> this->num_L; num_chi_temp >> this->num_chi;
        num_rb_types_temp >> this->num_rb_types; num_sat_types_temp >> this->num_sat_types;
        this->alt = stod(row[2]); this->dh = stod(row[3]); this->v = stod(row[4]);
        // calculate related parameters
        this->vyr = (this->v)*60.0*60.0*24.0*365.25; this->v_orbit = sqrt(G*Me/((Re + this->alt)*1000.0))/1000.0;
        this->V = 4*M_PI*(Re + this->alt)*(Re + this->alt)*(this->dh);
    } else {
        throw invalid_argument("Cell params.csv file could not be opened");
    }
    param_file.close();

    // pull basic arrays and vectors associated with Cell
    this->C_l = load_vec<double>(filepath + string("dataCl.npy"));
    this->C_nl = load_vec<double>(filepath + string("dataCnl.npy"));
    this->logL_edges = new Array1D<double>(filepath + string("logL.npy"));
    this->chi_edges = new Array1D<double>(filepath + string("chi.npy"));
    // calculate related values
    this->logL_ave = new Array1D<double>(num_L); this->chi_ave = new Array1D<double>(num_chi);
    for (size_t i = 0; i < this->num_L; i++) { // setup average arrays
        this->logL_ave->at(i) = (this->logL_edges->at(i) + this->logL_edges->at(i+1))/2.0;
    } 
    for (size_t i = 0; i < this->num_chi; i++) {
        this->chi_ave->at(i) = (this->chi_edges->at(i) + this->chi_edges->at(i+1))/2.0;
    }
    this->trackable = new ArrayND<bool,2>(false, array<size_t,2>({this->num_L, this->num_chi}));
    for (size_t i = 0; i < this->num_L; i++) {
        if (pow(10.0, logL_ave->at(i)) >= 1.0/10.0) {
            for (size_t j = 0; j < this->num_chi; j++) {
                (this->trackable)->at(array<size_t,2>({i,j})) = true;
            }
        }
    }

    // load N_bin values
    this->N_bins = new vector<ArrayND<double, 2> *>;
    string N_bins_path = filepath + string("N_bins/");
    ArrayND<double, 2> * N_bin_loc;
    size_t i = 0; // count through the files
    bool keep_going = true; // whether to keep trying new files
    while (keep_going) {
        try {
            N_bin_loc = new ArrayND<double,2>(N_bins_path + to_string(i) + string(".npy")); // get file
            this->N_bins->push_back(N_bin_loc);
            i++;
        } catch (invalid_argument &e) { // throws this when you reach the end of the files
            keep_going = false;
        }
    }

    size_t num_data = this->N_bins->size(); // number of time data points

    // setup satellite parameters
    this->S = new vector<Array1D<double> *>();
    this->S_d = new vector<Array1D<double> *>();
    this->D = new vector<Array1D<double> *>();
    for (size_t i = 0; i < num_data; i++) { // fill with empty arrays
        this->S->push_back(new Array1D<double>(this->num_sat_types));
        this->S_d->push_back(new Array1D<double>(this->num_sat_types));
        this->D->push_back(new Array1D<double>(this->num_sat_types));
    }
    
    this->m_s = new Array1D<double>(this->num_sat_types);
    this->sigma_s = new Array1D<double>(this->num_sat_types);
    this->sigma_s_km = new Array1D<double>(this->num_sat_types);
    this->lam_s = new Array1D<double>(this->num_sat_types);
    this->del_t = new Array1D<double>(this->num_sat_types);
    this->fail_t = new Array1D<double>(this->num_sat_types);
    this->tau_do = new Array1D<double>(this->num_sat_types);
    this->target_alt = new Array1D<double>(this->num_sat_types);
    this->up_time = new Array1D<double>(this->num_sat_types);
    this->alphaS = new Array1D<double>(this->num_sat_types);
    this->alphaD = new Array1D<double>(this->num_sat_types); 
    this->alphaN = new Array1D<double>(this->num_sat_types); 
    this->alphaR = new Array1D<double>(this->num_sat_types);
    this->P = new Array1D<double>(this->num_sat_types);
    this->AM_s = new Array1D<double>(this->num_sat_types);
    this->tau_s = new Array1D<double>(this->num_sat_types);
    this->C_s = new Array1D<double>(this->num_sat_types);
    this->expl_rate_L = new Array1D<double>(this->num_sat_types);
    this->expl_rate_D = new Array1D<double>(this->num_sat_types);

    // load satellites
    for (size_t i = 0; i < this->num_sat_types; i++) {
        this->load_sat(filepath + string("Satellite") + to_string(i) + string("/"), i);
    }

    // calculate related variables
    this->ascending = new Array1D<bool>(this->num_sat_types);
    for (size_t i = 0; i < this->num_sat_types; i++) {
        if (this->target_alt->at(i) > (this->alt + (this->dh)/2.0)) {(this->ascending)->at(i) = true;}
        else {(this->ascending)->at(i) = false;}
    }

    // setup rocket body parameters
    this->R = new vector<Array1D<double> *>();
    for (size_t i = 0; i < num_data; i++) { // fill with empty arrays
        this->R->push_back(new Array1D<double>(this->num_rb_types));
    }

    this->m_rb = new Array1D<double>(this->num_rb_types);
    this->sigma_rb = new Array1D<double>(this->num_rb_types);
    this->sigma_rb_km = new Array1D<double>(this->num_rb_types);
    this->lam_rb = new Array1D<double>(this->num_rb_types);
    this->AM_rb = new Array1D<double>(this->num_rb_types);
    this->tau_rb = new Array1D<double>(this->num_rb_types);
    this->C_rb = new Array1D<double>(this->num_rb_types);
    this->expl_rate_R = new Array1D<double>(this->num_rb_types);

    // load rockets
    for (size_t i = 0; i < this->num_rb_types; i++) {
        this->load_rb(filepath + string("RocketBody") + to_string(i) + string("/"), i);
    }

    // set up blank arrays for NCell and Cell to handle
    this->tau_N = new Array1D<double>(this->num_chi);
    this->cat_sat_N = new ArrayND<bool,3>(true, array<size_t,3>({this->num_sat_types, this->num_L, this->num_chi}));
    this->cat_rb_N = new ArrayND<bool,3>(true, array<size_t,3>({this->num_rb_types, this->num_L, this->num_chi}));

    // update catestrophic arrays
    this->update_cat_N();
}

void Cell::save(string &filepath, Array1D<bool> &filter, size_t filter_len) {
    /*
    saves current Cell object

    Input(s):
    filepath : path to location to save the object
    filter : filter of which datapoints to keep
    filter_len : number of values kept by the filter

    Output(s): None
    */
    
    // create directory
    bool exists; string file_path;
    try {
        exists = fs::create_directory(fs::path(filepath)); // make the folder representing the object
    } catch (const std::exception& e) {
        cout << "From Cell save " << e.what();
        return;
    }
    if (!exists) {
        cout << "Cell save failed : directory already exists" << endl;
        return;
    }

    // write parameters
    ofstream param_file; // file containing basic Cell parameters
    param_file.open(filepath + string("params.csv"), ios::out); // open relevant file
    if (param_file.is_open()) {
        param_file << "\"" << to_string(this->num_sat_types) << "\"" << ","; // write values
        param_file << "\"" << to_string(this->num_rb_types) << "\"" << ",";
        param_file << "\"" << to_string(this->alt) << "\"" << ",";
        param_file << "\"" << to_string(this->dh) << "\"" << ",";
        param_file << "\"" << to_string(this->v) << "\"" << ",";
        param_file << "\"" << to_string(this->num_L) << "\"" << ",";
        param_file << "\"" << to_string(this->num_chi) << "\"";
    } else {
        throw invalid_argument("Cell params.csv file could not be opened");
    }
    param_file.close();

    // write easy arrays
    file_path = filepath + "dataCl.npy"; save_vec(file_path, this->C_l, filter, filter_len);
    file_path = filepath + "dataCnl.npy"; save_vec(file_path, this->C_nl, filter, filter_len);
    file_path = filepath + "logL.npy"; this->logL_edges->save(file_path);
    file_path = filepath + "chi.npy"; this->chi_edges->save(file_path);

    // write N_bins values
    string bin_path = filepath + "N_bins/";
    try {
        exists = fs::create_directory(fs::path(bin_path)); // make the folder representing the list
    } catch (const std::exception& e) {
        cout << "From N_bins save " << e.what();
        return;
    }
    if (!exists) {
        cout << "N_bins save failed : directory already exists" << endl;
        return;
    }
    size_t index = 0;
    for (size_t i = 0; i < this->N_bins->size(); i++) {
        if (filter.at(i)) {
            file_path = bin_path + to_string(index); this->N_bins->at(i)->save(file_path);
            index++;
        }
    }

    // write satellites and rockets
    for (size_t i = 0; i < this->num_sat_types; i++) {
        this->save_sat(filepath + string("Satellite") + to_string(i) + string("/"), i, filter, filter_len);
    } for (size_t i = 0; i < this->num_rb_types; i++) {
        this->save_rb(filepath + string("RocketBody") + to_string(i) + string("/"), i, filter, filter_len);
    }
}

void Cell::add_event(Event * event) {
    // adds the given event to the cell
    this->event_list->push_back(event); this->num_events++;
}

void Cell::load_sat(const string &filepath, size_t i) {
    /*
    loads a saved satellite object into the cell

    Input(s):
    filepath : relative or absolute filepath to the saved satellite object
    i : indicates that this is the "ith" satellite type in the cell

    Output(s): None

    Note(s): function assumes that all relavant variables have already been
             allocated to the correct size.
    */
    ifstream param_file; // file containing basic Satellite parameters
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
        this->m_s->at(i) = stod(row[0]); this->sigma_s->at(i) = stod(row[1]);
        this->lam_s->at(i) = stod(row[2]); this->del_t->at(i) = stod(row[3]);
        this->fail_t->at(i) = stod(row[4]); this->tau_do->at(i) = stod(row[5]);
        this->target_alt->at(i) = stod(row[6]); this->up_time->at(i) = stod(row[7]);
        this->alphaS->at(i) = stod(row[8]); this->alphaD->at(i) = stod(row[9]);
        this->alphaN->at(i) = stod(row[10]); this->alphaR->at(i) = stod(row[11]);
        this->P->at(i) = stod(row[12]); this->AM_s->at(i) = stod(row[13]);
        this->C_s->at(i) = stod(row[14]); this->expl_rate_L->at(i) = stod(row[15]);
        this->expl_rate_D->at(i) = stod(row[16]);
        // calculate related parameters
        this->sigma_s_km->at(i) = (this->sigma_s->at(i))/1e6;
    } else {
        throw invalid_argument("Satellite params.csv file could not be opened");
    }
    param_file.close();
    
    // get tracked values
    Array1D<double> S_temp = Array1D<double>(filepath + string("dataS.npy"));
    Array1D<double> SD_temp = Array1D<double>(filepath + string("dataSd.npy"));
    Array1D<double> D_temp = Array1D<double>(filepath + string("dataD.npy"));
    size_t tot_size = S_temp.get_tot_size();
    for (size_t j = 0; j < tot_size; j++) {
        this->S->at(j)->at(i) = S_temp.at(j);
        this->S_d->at(j)->at(i) = SD_temp.at(j);
        this->D->at(j)->at(i) = D_temp.at(j);
    }
}

void Cell::save_sat(const string &filepath, size_t i, Array1D<bool> &filter, size_t filter_len) {
    /*
    saves a single satellite object

    Input(s):
    filepath : path to location to save the object
    i : indicates that this is the "ith" satellite type in the cell
    filter : filter of which datapoints to keep
    filter_len : number of values kept by the filter

    Output(s): None
    */
    // create directory
    bool exists; string file_path;
    try {
        exists = fs::create_directory(fs::path(filepath)); // make the folder representing the object
    } catch (const std::exception& e) {
        cout << "From sat save " << e.what();
        return;
    }
    if (!exists) {
        cout << "sat save failed : directory already exists" << endl;
        return;
    }

    // write parameters
    ofstream param_file; // file containing basic Cell parameters
    param_file.open(filepath + string("params.csv"), ios::out); // open relevant file
    if (param_file.is_open()) {
        param_file << "\"" << to_string(this->m_s->at(i)) << "\"" << ","; // write values
        param_file << "\"" << to_string(this->sigma_s->at(i)) << "\"" << ",";
        param_file << "\"" << to_string(this->lam_s->at(i)) << "\"" << ",";
        param_file << "\"" << to_string(this->del_t->at(i)) << "\"" << ",";
        param_file << "\"" << to_string(this->fail_t->at(i)) << "\"" << ",";
        param_file << "\"" << to_string(this->tau_do->at(i)) << "\"" << ",";
        param_file << "\"" << to_string(this->target_alt->at(i)) << "\"" << ",";
        param_file << "\"" << to_string(this->up_time->at(i)) << "\"" << ",";
        param_file << "\"" << to_string(this->alphaS->at(i)) << "\"" << ",";
        param_file << "\"" << to_string(this->alphaD->at(i)) << "\"" << ",";
        param_file << "\"" << to_string(this->alphaN->at(i)) << "\"" << ",";
        param_file << "\"" << to_string(this->alphaR->at(i)) << "\"" << ",";
        param_file << "\"" << to_string(this->P->at(i)) << "\"" << ",";
        param_file << "\"" << to_string(this->AM_s->at(i)) << "\"" << ",";
        param_file << "\"" << to_string(this->C_s->at(i)) << "\"" << ",";
        param_file << "\"" << to_string(this->expl_rate_L->at(i)) << "\"" << ",";
        param_file << "\"" << to_string(this->expl_rate_D->at(i)) << "\"";
    } else {
        throw invalid_argument("Cell params.csv file could not be opened");
    }
    param_file.close();

    // write data arrays
    vector<double> S_loc = vector<double>(this->S->size()); // copy arrays for this satellite type
    vector<double> SD_loc = vector<double>(this->D->size());
    vector<double> D_loc = vector<double>(this->S->size());
    for (size_t j = 0; j < this->S->size(); j++) {
        S_loc[j] = this->S->at(j)->at(i);
        SD_loc[j] = this->S_d->at(j)->at(i);
        D_loc[j] = this->D->at(j)->at(i);
    }
    file_path = filepath + "dataS.npy"; save_vec(file_path, &S_loc, filter, filter_len);
    file_path = filepath + "dataSd.npy"; save_vec(file_path, &SD_loc, filter, filter_len);
    file_path = filepath + "dataD.npy"; save_vec(file_path, &D_loc, filter, filter_len);
}

void Cell::load_rb(const string &filepath, size_t i) {
    /*
    loads a saved rocket body object into the cell

    Input(s):
    filepath : relative or absolute filepath to the saved rocket body
    i : indicates that this is the "ith" rocket type in the cell

    Output(s): None

    Note(s): function assumes that all relavant variables have already been
             allocated to the correct size.
    */
    ifstream param_file; // file containing basic RocketBody parameters
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
        this->m_rb->at(i) = stod(row[0]); this->sigma_rb->at(i) = stod(row[1]);
        this->lam_rb->at(i) = stod(row[2]); this->AM_rb->at(i) = stod(row[3]);
        this->C_rb->at(i) = stod(row[4]); this->expl_rate_R->at(i) = stod(row[5]);
        // calculate related parameters
        this->sigma_rb_km->at(i) = (this->sigma_rb->at(i))/1e6;
    } else {
        throw invalid_argument("RocketBody params.csv file could not be opened");
    }
    param_file.close();
    
    // get tracked values
    Array1D<double> R_temp = Array1D<double>(filepath + string("data.npy"));
    size_t tot_size = R_temp.get_tot_size();
    for (size_t j = 0; j < tot_size; j++) {
        this->R->at(j)->at(i) = R_temp.at(j);
    }
}

void Cell::save_rb(const string &filepath, size_t i, Array1D<bool> &filter, size_t filter_len) {
    /*
    saves a single rocket body object

    Input(s):
    filepath : path to location to save the object
    i : indicates that this is the "ith" rocket type in the cell
    filter : filter of which datapoints to keep
    filter_len : number of values kept by the filter

    Output(s): None
    */
    // create directory
    bool exists; string file_path;
    try {
        exists = fs::create_directory(fs::path(filepath)); // make the folder representing the object
    } catch (const std::exception& e) {
        cout << "From rb save " << e.what();
        return;
    }
    if (!exists) {
        cout << "rb save failed : directory already exists" << endl;
        return;
    }

    // write parameters
    ofstream param_file; // file containing basic Cell parameters
    param_file.open(filepath + string("params.csv"), ios::out); // open relevant file
    if (param_file.is_open()) {
        param_file << "\"" << to_string(this->m_rb->at(i)) << "\"" << ","; // write values
        param_file << "\"" << to_string(this->sigma_rb->at(i)) << "\"" << ",";
        param_file << "\"" << to_string(this->lam_rb->at(i)) << "\"" << ",";
        param_file << "\"" << to_string(this->AM_rb->at(i)) << "\"" << ",";
        param_file << "\"" << to_string(this->C_rb->at(i)) << "\"" << ",";
        param_file << "\"" << to_string(this->expl_rate_R->at(i)) << "\"";
    } else {
        throw invalid_argument("Cell params.csv file could not be opened");
    }
    param_file.close();

    // write data arrays
    vector<double> R_loc = vector<double>(this->R->size()); // copy arrays for this rocket type
    for (size_t j = 0; j < this->R->size(); j++) {
        R_loc[j] = this->R->at(j)->at(i);
    }
    file_path = filepath + "data.npy"; save_vec(file_path, &R_loc, filter, filter_len);
}

void Cell::update_cat_N() {
    // updates catestrophic table values
    bool is_cat;
    double ave_L;
    double ave_AM;
    for (size_t i = 0; i < this->num_L; i++) {
        ave_L = pow(10.0, (this->logL_ave)->at(i)); // get average length in m
        for (size_t j = 0; j < this->num_chi; j++) {
            ave_AM = pow(10.0, (this->chi_ave)->at(j)); // get average AM in m^2/kg
            for (size_t k = 0; k < this->num_sat_types; k++) {
                is_cat = is_catastrophic((this->m_s)->at(k), ave_L, ave_AM, this->v); // determine if catestrophic
                (this->cat_sat_N)->at(array<size_t,3>({k,i,j})) = is_cat; // set the result
            } for (size_t k = 0; k < this->num_rb_types; k++) { // do the same thing for rockets
                is_cat = is_catastrophic((this->m_rb)->at(k), ave_L, ave_AM, this->v); // determine if catestrophic
                (this->cat_rb_N)->at(array<size_t,3>({k,i,j})) = is_cat; // set the result
            }
        }
    }
}

void Cell::dxdt_cell(size_t time, Array1D<double> &dSdt, Array1D<double> &dS_ddt, Array1D<double> &dDdt, Array1D<double> &dRdt, 
                     double &dC_ldt, double &dC_nldt, Array1D<double> &S_out, Array1D<double> &S_dout, Array1D<double> &D_out, 
                     Array1D<double> &R_out, ArrayND<double,2> &N_out, ArrayND<double,2> &D_dt, ArrayND<double,2> &DR_dt, 
                     ArrayND<double,2> &R_dt, ArrayND<double,3> &CS_dt, ArrayND<double,3> &CR_dt, Array1D<double> &expl_S, Array1D<double> &expl_R) {
    /*
    calculates the rate of collisions and decays from each debris bin, the rate
    of decaying/de-orbiting satellites, the rate of launches/deorbit starts of satallites, 
    and the rate of creation of derelicts at the given time, due only to events in the cell

    Inputs(s):
    time : index of the values to use

    Output(s) (by reference):
    dSdt : array of rate of change of the number of live satellites in the cell of each type due to only processes
            withing the cell (not including ascending satellites) (yr^(-1))
    dS_ddt : array of rate of change of the number of de-orbiting satellites in the cell of each type
             (not including de-orbiting satellites) (yr^(-1))
    dDdt : array of rate of change of the number of derelict satellites in the cell of each type
           (not including decaying derelicts) (yr^(-1))
    dRdt : array of rate of change of number of rocket bodies in the cell of each type (yr^(-1))
           (not including decaying rocket bodies)
    dC_ldt : rate of change of total number of catestrophic collisions (1/yr)
    dC_nldt : rate of change of total number of non-catestrophic collisions (1/yr)
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
    array<size_t,2> index_2d = {0,0}; array<size_t,3> index_3d = {0,0,0}; 
    // get current N value
    ArrayND<double,2> N = *(this->N_bins->back());

    // setup some temp parameters
    double sigma_loc_km0; double R1; double sigma_loc_km1; double sigma_comb; double n;
    
    // start with satellite collisions/events
    for (size_t i = 0; i < this->num_sat_types; i++) {

        // get current parameters
        double S0 = this->S->at(time)->at(i);
        double SD0 = this->S_d->at(time)->at(i);
        double D0 = this->D->at(time)->at(i);
        sigma_loc_km0 = this->sigma_s_km->at(i);
        double alphaS0 = this->alphaS->at(i);
        double alphaD0 = this->alphaD->at(i);
        double alphaN0 = this->alphaN->at(i);
        double alphaR0 = this->alphaR->at(i);

        // setup temp parameters
        double dSdt_loc; double dS_ddt_loc; double dDdt_loc;
        double dSSdt_loc; double dSS_ddt_loc; double dSDdt_loc;
        double dSRdt_loc; double dS_dS_ddt_loc; double dS_dDdt_loc;
        double dS_dRdt_loc; double dDDdt_loc; double dDRdt_loc;
        double S1; double SD1; double D1; double alphaS1;

        index_3d[0] = i;
        // handle satellite-debris collisions
        for (size_t j = 0; j < this->num_L; j++) {
            index_2d[0] = j; index_3d[1] = j;
            for (size_t k = 0; k < this->num_chi; k++) {
                index_2d[1] = k; index_3d[2] = k;
                n = N.at(index_2d)/(this->V); // calculate debris density
                if (this->trackable->at(index_2d)) { // if avoidance is possible
                    dSdt_loc = alphaN0*n*sigma_loc_km0*(this->vyr)*S0;
                    dS_ddt_loc = alphaN0*n*sigma_loc_km0*(this->vyr)*SD0;
                } else { // if avoidance is impossible
                    dSdt_loc = n*sigma_loc_km0*(this->vyr)*S0;
                    dS_ddt_loc = n*sigma_loc_km0*(this->vyr)*SD0;
                }
                dDdt_loc = n*sigma_loc_km0*(this->vyr)*D0; // never avoids
                // update return values
                dSdt.at(i) -= dSdt_loc; dS_ddt.at(i) -= dS_ddt_loc;
                if ((this->cat_sat_N)->at(index_3d) == false) { // collisions that create derelicts
                    dDdt.at(i) += dSdt_loc + dS_ddt_loc;
                } else { // collisions that destroy derelicts
                    dDdt.at(i) -= dDdt_loc;
                }
                CS_dt.at(array<size_t,3>({i,j,k})) += dSdt_loc + dS_ddt_loc + dDdt_loc;
                if (this->cat_sat_N->at(index_3d) == true) { // update collision count
                    dC_ldt += dSdt_loc + dS_ddt_loc + dDdt_loc;
                } else {
                    dC_nldt += dSdt_loc + dS_ddt_loc + dDdt_loc;
                }
            }
        }

        index_2d[0] = i;
        // handle satellite-satellite collisions
        for (size_t j = 0; j < this->num_sat_types; j++) {
            index_2d[1] = j;
            // get local parameters for this satellite type
            S1 = this->S->at(time)->at(j);
            SD1 = this->S_d->at(time)->at(j);
            D1 = this->D->at(time)->at(j);
            sigma_loc_km1 = this->sigma_s_km->at(j);
            alphaS1 = this->alphaS->at(j);
            
            // calculate combined cross-section
            sigma_comb = sigma_loc_km0 + sigma_loc_km1 + 2.0*sqrt(sigma_loc_km0*sigma_loc_km1);

            // calculate collisions
            dSSdt_loc = alphaS0*alphaS1*sigma_comb*(this->vyr)*S0*S1/(this->V);
            dSS_ddt_loc = alphaS0*alphaS1*sigma_comb*(this->vyr)*S0*SD1/(this->V);
            dSDdt_loc = alphaD0*sigma_comb*(this->vyr)*S0*D1/(this->V);
            dS_dS_ddt_loc = alphaS0*alphaS1*sigma_comb*(this->vyr)*SD0*SD1/(this->V);
            dS_dDdt_loc = alphaS0*sigma_comb*(this->vyr)*SD0*D1/(this->V);
            dDDdt_loc = sigma_comb*(this->vyr)*D0*D1/(this->V);

            // update return values
            if (i <= j) { // avoid double counting
                D_dt.at(index_2d) += dSSdt_loc + dS_dS_ddt_loc + dDDdt_loc;
                dC_ldt += dSSdt_loc + dS_dS_ddt_loc + dDDdt_loc;
            }
            D_dt.at(index_2d) += dSS_ddt_loc + dSDdt_loc + dS_dDdt_loc; // these aren't double counted
            dC_ldt += dSS_ddt_loc + dSDdt_loc + dS_dDdt_loc;
            if (i == j) { // destroys two of the same type in one go
                dSdt.at(i) -= 2.0*dSSdt_loc + dSS_ddt_loc + dSDdt_loc;
                dS_ddt.at(i) -= dSS_ddt_loc + 2.0*dS_dS_ddt_loc + dS_dDdt_loc;
                dDdt.at(i) -= dSDdt_loc + dS_dDdt_loc + 2.0*dDDdt_loc;
            } else { // doesnt, and have to treat first and second satellite types seperately sometimes
                dSdt.at(i) -= dSSdt_loc + dSS_ddt_loc + dSDdt_loc;
                dS_ddt.at(i) -= dS_dS_ddt_loc + dS_dDdt_loc;
                dDdt.at(i) -= dDDdt_loc;
                dS_ddt.at(j) -= dSS_ddt_loc;
                dDdt.at(j) -= dSDdt_loc + dS_dDdt_loc;
            }
        }

        // handle satellite-rocket collisions
        for (size_t j = 0; j < this->num_rb_types; j++) {
            index_2d[1] = j;
            // get local parameters for this rocket type
            R1 = this->R->at(time)->at(j);
            sigma_loc_km1 = this->sigma_rb_km->at(j);
            
            // calculate combined cross-section
            sigma_comb = sigma_loc_km0 + sigma_loc_km1 + 2.0*sqrt(sigma_loc_km0*sigma_loc_km1);

            // calculate collisions
            dSRdt_loc = alphaR0*sigma_comb*(this->vyr)*S0*R1/(this->V);
            dS_dRdt_loc = alphaR0*sigma_comb*(this->vyr)*SD0*R1/(this->V);
            dDRdt_loc = sigma_comb*(this->vyr)*D0*R1/(this->V);

            // update return values
            DR_dt.at(index_2d) += dSRdt_loc + dS_dRdt_loc + dDRdt_loc; // these aren't double counted
            dC_ldt += dSRdt_loc + dS_dRdt_loc + dDRdt_loc;
            dSdt.at(i) -= dSRdt_loc;
            dS_ddt.at(i) -= dS_dRdt_loc;
            dDdt.at(i) -= dDRdt_loc;
            dRdt.at(j) -= dSRdt_loc + dS_dRdt_loc + dDRdt_loc;
        }

        // handle decays/explosions
        double expl_S = (this->expl_rate_L->at(i))*S0/100.0; // calculate explosion rates
        double expl_S_d = (this->expl_rate_L->at(i))*SD0/100.0;
        double expl_D = (this->expl_rate_D->at(i))*D0/100.0;
        double decay_S = D0/(this->tau_s->at(i)); // satellites ascending/descending, switching to de-orbit
        double kill_S = 0.0;
        double ascend_S = 0.0;
        if (this->ascending->at(i)) { // failure on ascending
            kill_S = S0/(this->fail_t->at(i));
            ascend_S = S0/(this->up_time->at(i));
        } else { // end of lifetime
            kill_S = S0/(this->del_t->at(i));
        }
        double deorbit_S = SD0/(this->tau_do->at(i));

        // update return values
        double P_loc = this->P->at(i);
        dSdt.at(i) -= expl_S + kill_S;
        dS_ddt.at(i) += P_loc*kill_S - expl_S_d;
        dDdt.at(i) += kill_S*(1-P_loc) - expl_D;
        S_out.at(i) += ascend_S;
        S_dout.at(i) += deorbit_S;
        D_out.at(i) += decay_S;
    }

    // handle rocket-body only events
    for (size_t i = 0; i < this->num_rb_types; i++) {
        index_3d[0] = i;
        // setup temp local parameter
        double dRRdt_loc; double dRdt_loc;

        // get relevant local values
        double R0 = this->R->at(i)->at(time);
        sigma_loc_km0 = this->sigma_rb_km->at(i);

        // handle rocket-debris collisions
        for (size_t j = 0; j < this->num_L; j++) {
            index_3d[1] = j; index_2d[0] = j;
            for (size_t k = 0; k < this->num_chi; k++) {
                index_3d[2] = k; index_2d[1] = k;
                n = N.at(index_2d)/(this->V); // calculate debris density
                dRdt_loc = n*sigma_loc_km0*(this->vyr)*R0; // never avoids
                // update return values
                if ((this->cat_rb_N)->at(index_3d) == true) { // collisions that destroy rockets
                    dRdt.at(i) -= dRdt_loc;
                }
                CR_dt.at(index_3d) += dRdt_loc;
                if (this->cat_rb_N->at(index_3d) == true) { // update collision count
                    dC_ldt += dRdt_loc;
                } else {
                    dC_nldt += dRdt_loc;
                }
            }
        }

        index_2d[0] = i;
        // handle rocket-rocket collisions
        for (size_t j = 0; j < this->num_rb_types; j++) {
            index_2d[1] = j;
            // get local parameters for this rocket type
            R1 = this->R->at(i)->at(time);
            sigma_loc_km1 = this->sigma_rb_km->at(i);
            
            // calculate combined cross-section
            sigma_comb = sigma_loc_km0 + sigma_loc_km1 + 2.0*sqrt(sigma_loc_km0*sigma_loc_km1);

            // calculate collision rate
            dRRdt_loc = R0*R1*sigma_comb*(this->vyr)/(this->V);

            // update return values
            if (i <= j) { // avoid double counting
                R_dt.at(index_2d) += dRRdt_loc;
                dC_ldt += dRRdt_loc;
            }
            if (i == j) { // destroys two of the same type in one go
                dRdt.at(i) -= 2.0*dRRdt_loc;
            } else {
                dRdt.at(i) -= dRRdt_loc;
            }
        }

        // handle decays/explosions
        double expl_R = (this->expl_rate_R->at(i))*R0/100.0; // calculate explosion rates
        double decay_R = R0/(this->tau_rb->at(i)); // calculate decay rates

        // update return values
        dRdt.at(i) -= expl_R + decay_R;
        R_out.at(i) += decay_R;
    }

    // handle debris decays
    for (size_t i = 0; i < this->num_L; i++) {
        index_2d[0] = i;
        for (size_t j = 0; j < this->num_chi; j++) {
            index_2d[1] = j;
            N_out.at(index_2d) = N.at(index_2d)/(this->tau_N->at(j));
        }
    }
}

Cell::~Cell() {
    
    // basic class destructor
    // Note : the Cell is not assumed to own the logL_edges and chi_edges arrays, and hence
    //        will not free them
    delete this->trackable; delete this->ascending; delete this->cat_rb_N;
    delete this->cat_sat_N; delete this->tau_N; delete this->logL_ave; delete this->chi_ave; 
    delete this->sigma_s_km; delete this->expl_rate_R; delete this->C_rb; delete this->tau_rb; 
    delete this->AM_rb; delete this->lam_rb; delete this->sigma_rb; delete this->m_rb; 
    delete this->expl_rate_D; delete this->expl_rate_L; delete this->C_s; delete this->tau_s; 
    delete this->AM_s; delete this->P; delete this->alphaR; delete this->alphaN; delete this->alphaD; 
    delete this->alphaS; delete this->up_time; delete this->target_alt; delete this->tau_do; delete this->del_t;
    delete this->fail_t; delete this->lam_s; delete this->sigma_s; delete this->m_s; delete this->sigma_rb_km;
    delete this->C_l; delete this->C_nl;

    for (size_t i = 0; i < this->S->size(); i++) {
        delete (*this->S)[i]; delete (*this->S_d)[i]; delete (*this->D)[i]; delete (*this->N_bins)[i];
    }
    for (size_t i = 0; i < this->R->size(); i++) {delete (*this->R)[i];}
    for (size_t i = 0; i < this->num_events; i++) {delete (*this->event_list)[i];}
    delete this->S; delete this->S_d; delete this->D; delete this->R; delete this->N_bins; delete this->event_list;
}