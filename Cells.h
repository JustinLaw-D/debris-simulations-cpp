// declaration of NCell and Cell classes

#include <cstddef>
#include <vector>
#include <iostream>
#include <random>
#include <algorithm>
#include <experimental/filesystem>
#include "ObjectsEvents.h"
#include "Arrays.h"
#include "AtmosphericDecayModels.h"
#include "BreakupModel.h"
#include "Constants.h"
#pragma once

using namespace std;
namespace fs = std::experimental::filesystem;

class Cell
{
    template <class URNG>
    friend class NCell;

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
        Array1D<double> * tau_N; // decay lifetimes for each debris bin (yr)
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
             vector<Event *> * event_list, size_t num_events, double alt, double dh, Array1D<double> * tau_N, double v,
             vector<double> * C_l, vector<double> * C_nl);
        Cell(const string &filepath); // load cell from file
        void save(string &filepath, Array1D<bool> &filter, size_t filter_len); // function for saving data
        void add_event(Event * event); // adds event to the cell
        void load_sat(const string &filepath, size_t i); // loads a single satellite type into the cell
        void save_sat(const string &filepath, size_t i, Array1D<bool> &filter, size_t filter_len); // saves satellite type
        void load_rb(const string &filepath, size_t i); // loads a single rocket body into the cell
        void save_rb(const string &filepath, size_t i, Array1D<bool> &filter, size_t filter_len); // saves rocket body type
        // calculating local rates of change
        void dxdt_cell(size_t time, Array1D<double> &dSdt, Array1D<double> &dS_ddt, Array1D<double> &dDdt, Array1D<double> &dRdt,
                       double &dC_ldt, double &dC_nldt, Array1D<double> &S_out, Array1D<double> &S_dout, Array1D<double> &D_out, 
                       Array1D<double> &R_out, ArrayND<double,2> &N_out, ArrayND<double,2> &D_dt, ArrayND<double,2> &DR_dt, 
                       ArrayND<double,2> &R_dt, ArrayND<double,3> &CS_dt, ArrayND<double,3> &CR_dt, Array1D<double> &expl_S, 
                       Array1D<double> &expl_R);
        void update_cat_N(); // updates catestrophic debris tables
        ~Cell(); // destructor
};

template <class URNG>
class NCell {

    private:
        Array1D<double> * alts;
        Array1D<double> * dhs;
        size_t num_L;
        size_t num_chi;
        double update_period;
        double min_lifetime;
        size_t time;
        size_t lupdate_time;
        vector<double> * t;
        Cell ** cells;
        size_t num_cells;
        Array1D<double> * logL_edges;
        Array1D<double> * logL_ave;
        Array1D<double> * chi_edges;
        Array1D<double> * chi_ave;
        ArrayND<double, 4> * sat_coll_prob_tables;
        ArrayND<double, 4> * rb_coll_prob_tables;
        ArrayND<double, 4> * sat_expl_prob_tables;
        ArrayND<double, 4> * rb_expl_prob_tables;
        ArrayND<double, 2> * N_upper_init;

    public:
        NCell(string &filepath, size_t num_dir, URNG &generator); // constructor based on loading data from file
        void save(string &filepath, string &name, double gap); // function for saving data
        // calculates the probability tables
        void calc_prob_table(ArrayND<double, 4> * table, size_t indx, char etyp, char ttyp, size_t num_dir, URNG &generator);
        bool add_event(Event * event, double alt); // adds event to the system
        void dxdt(size_t time, bool upper, Array1D<double> *dSdt, Array1D<double> *dS_ddt, Array1D<double> *dDdt,
                  Array1D<double> *dRdt, ArrayND<double,3> &dNdt, double *dC_ldt, double *dC_nldt); // calculate rates of change
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
        size_t alt_to_index(double alt); // converts a given altitude to the index of the corresponding cell
        ~NCell(); // destructor
};

template <class URNG>
NCell<URNG>::NCell(string &filepath, size_t num_dir, URNG &generator) {
    /*
    loads NCell object from saved data

    Input(s):
    filepath : string containing relative or absolute path to saved data
    num_dir : number of random directions to sample
    generator : uniform random number generator used for direction randomization

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
        this->min_lifetime = stod(row[4]);
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

    this->cells = new Cell *[this->num_cells];
    for (size_t i = 0; i < this->num_cells; i++) {
        this->cells[i] = new Cell(filepath + "cell" + to_string(i) + "/");
    }

    // set initial N table for debris decaying into top cell
    this->N_upper_init = new ArrayND<double, 2>(*this->cells[this->num_cells-1]->N_bins->at(0));

    // initialize/calculate probability tables
    this->sat_coll_prob_tables = new ArrayND<double, 4>(array<size_t,4>({this->num_cells, this->num_cells, this->num_L, this->num_chi}));
    this->rb_coll_prob_tables = new ArrayND<double, 4>(array<size_t,4>({this->num_cells, this->num_cells, this->num_L, this->num_chi}));
    this->sat_expl_prob_tables = new ArrayND<double, 4>(array<size_t,4>({this->num_cells, this->num_cells, this->num_L, this->num_chi}));
    this->rb_expl_prob_tables = new ArrayND<double, 4>(array<size_t,4>({this->num_cells, this->num_cells, this->num_L, this->num_chi}));
    for (size_t i = 0; i < this->num_cells; i++) {
        this->calc_prob_table(this->sat_coll_prob_tables, i, 'c', 's', num_dir, generator);
        this->calc_prob_table(this->rb_coll_prob_tables, i, 'c', 'r', num_dir, generator);
        this->calc_prob_table(this->sat_expl_prob_tables, i, 'e', 's', num_dir, generator);
        this->calc_prob_table(this->rb_expl_prob_tables, i, 'e', 'r', num_dir, generator);
    }

    // calculate atmospheric decay lifetimes
    this->update_lifetimes(t->back());
}

template <class URNG>
void NCell<URNG>::save(string &filepath, string &name, double gap) {
    /*
    saves current NCell object

    Input(s):
    filepath : path to location to save the object
    name : name of the object
    gap : smallest gap in time between data points saved to accept (yr)

    Output(s): None
    */
    
    // create directory
    string true_path = filepath + name + string("/"); string file_path;
    bool exists;
    try {
        exists = fs::create_directory(fs::path(true_path)); // make the folder representing the object
    } catch (const std::exception& e) {
        cout << "From NCell save " << e.what();
        return;
    }
    if (!exists) {
        cout << "NCell save failed : directory already exists" << endl;
        return;
    }

    // write parameters
    ofstream param_file; // file containing basic NCell parameters
    param_file.open(true_path + string("params.csv"), ios::out); // open relevant file
    if (param_file.is_open()) {
        param_file << "\"" << to_string(this->num_L) << "\"" << ","; // write values
        param_file << "\"" << to_string(this->num_chi) << "\"" << ",";
        param_file << "\"" << to_string(this->num_cells) << "\"" << ",";
        param_file << "\"" << to_string(this->update_period) << "\"" << ",";
        param_file << "\"" << to_string(this->min_lifetime) << "\"";
    } else {
        throw invalid_argument("NCell params.csv file could not be opened");
    }
    param_file.close();

    // build the filter
    Array1D<bool> filter = Array1D<bool>(this->t->size());
    size_t filter_len = 0; // number of values that survive the filter
    if (filter.get_tot_size() > 0) {
        double prev_t = this->t->at(0);
        filter.at(0) = true;
        filter_len++;
        for (size_t i = 1; i < filter.get_tot_size(); i++) {
            if (this->t->at(i) - prev_t >= gap) {
                prev_t = this->t->at(i); filter.at(i) = true; filter_len++;
            } else {filter.at(i) = false;}
        }
    }

    // write easy arrays
    file_path = true_path + string("t.npy"); save_vec(file_path, this->t, filter, filter_len);
    file_path = true_path + string("alts.npy"); this->alts->save(file_path);
    file_path = true_path + string("dhs.npy"); this->dhs->save(file_path);
    file_path = true_path + string("logL.npy"); this->logL_edges->save(file_path);
    file_path = true_path + string("chi.npy"); this->chi_edges->save(file_path);

    for (size_t i = 0; i < this->num_cells; i++) {
        string cell_path = true_path + string("cell") + to_string(i) + "/";
        this->cells[i]->save(cell_path, filter, filter_len);
    }

}

template <class URNG>
void NCell<URNG>::calc_prob_table(ArrayND<double, 4> * table, size_t indx, char etyp, char ttyp, size_t num_dir, URNG &generator) {
    /*
    calculates probability table for explosion/collision debris generation in the system

    Input(s):
    table : probability table to fill
    indx : index of the cell to calculate probability table for
    etyp : event type, either collision ('c') or explosion ('e')
    ttyp : target object type, either satellite ('s') or rocket body ('r')
    num_dir : number of random directions to sample
    generator : uniform random number generator used for direction randomization

    Output(s): None

    Note(s): assumes that the probability table objects have been properly instantiated
    */
    double v0 = (this->cells)[indx]->v_orbit*1000.0; // orbital velocity in m/s
    double r = (this->cells)[indx]->alt; // in km
    double L_min = pow(10, this->logL_edges->at(0));
    double L_max = pow(10, this->logL_edges->at(this->num_L));
    double chi_min = this->chi_edges->at(0);
    double chi_max = this->chi_edges->at(this->num_chi);
    double * phi = new double[num_dir]; // random directions
    double * theta = new double[num_dir];
    uniform_real_distribution<double> dist = uniform_real_distribution<double>(0.0, 1.0); // standard uniform distribution
    double P_temp_theta; // used to hold cummulative probability generated for theta
    for (size_t i = 0; i < num_dir; i++) {
        phi[i] = dist(generator)*2*M_PI;
        P_temp_theta = dist(generator);
        theta[i] = acos(1.0-2.0*P_temp_theta);
    }
    for (size_t i = 0; i < this->num_cells; i++) { // iterate through cells
        Cell * curr_cell = (this->cells)[i];
        double alt_min = curr_cell->alt - curr_cell->dh/2.0; // in km
        double alt_max = curr_cell->alt + curr_cell->dh/2.0;
        double v_min2 = G*Me*(2.0/((Re + r)*1000.0) - 1.0/((Re + alt_min)*1000.0)); // minimum velocity squared (m/s)
        double v_max2 = G*Me*(2.0/((Re + r)*1000.0) - 1.0/((Re + alt_max)*1000.0)); // maximum velocity squared (m/s)
        for (size_t j = 0; j < this->num_L; j++) { // iterate through bins
            double bin_bot_L = this->logL_edges->at(j);
            double bin_top_L = this->logL_edges->at(j+1);
            // probability of L being in this bin
            double L_prob = L_cdf(pow(10, bin_top_L), L_min, L_max, etyp) - L_cdf(pow(10, bin_bot_L), L_min, L_max, etyp);
            double L_ave = pow(10, this->logL_ave->at(j));
            for (size_t k = 0; k < this->num_chi; k++) {
                double bin_bot_chi = this->chi_edges->at(k);
                double bin_top_chi = this->chi_edges->at(k+1);
                // total probability of being in this bin
                double curr_prob = L_prob*(X_cdf(bin_top_chi, chi_min, chi_max, L_ave, ttyp) - X_cdf(bin_bot_chi, chi_min, chi_max, L_ave, ttyp));
                double sum = 0.0; // total of all random directions sampled
                for (size_t l = 0; l < num_dir; l++) { // sample random directions
                    if ((v_min2 < 0) && (v_max2 < 0)) {sum += 0.0;}
                    else if (v_min2 < 0) {sum += curr_prob*(vprime_cdf(sqrt(v_max2), v0, theta[l], phi[l], this->chi_ave->at(k), etyp));}
                    else {sum += curr_prob*(vprime_cdf(sqrt(v_max2), theta[l], phi[l], v0, this->chi_ave->at(k), etyp) - vprime_cdf(sqrt(v_min2), v0, theta[l], phi[l], this->chi_ave->at(k), etyp));}
                }
                table->at(array<size_t, 4>({indx, i, j, k})) = sum/num_dir;
            }
        }
    }
}

template <class URNG>
bool NCell<URNG>::add_event(Event * event, double alt) {
    /* 
    adds event to the system at the given altitude

    Input(s):
    event : pointer to event to add
    alt : altitude to add the event at (km)

    Output(s):
    in : true if the given altitude is in the system, false otherwise

    Note: if false is returned, event is not added to the system
    */
    size_t indx = this->alt_to_index(alt); // get index
    if (indx == this->num_cells) {
        return false;
    } else {
        this->cells[indx]->add_event(event);
        return true;
    }
}

template <class URNG>
void NCell<URNG>::dxdt(size_t time, bool upper, Array1D<double> *dSdt, Array1D<double> *dS_ddt, Array1D<double> *dDdt,
                       Array1D<double> *dRdt, ArrayND<double,3> &dNdt, double *dC_ldt, double *dC_nldt) {
    /*
    calculates the total rate of change of all parameters in the system in at the given time

    Input(s):
    time : index of the time in the system rates of change are being caculated at
    upper : whether or not to including debris coming down into the top cell
    dSdt : list of rates of change in S for each cell of each type (1/yr)
    dS_ddt : list of rates of change in S_d for each cell of each type (1/yr)
    dDdt : list of rates of change in D for each cell of each type (1/yr)
    dRdt : list of rates of change in R for each cell of each type (1/yr)
    dNdt : rates of change in N for each cell, binned (1/yr)
    dC_ldt : list of rates of change in C_l for each cell (1/yr)
    dC_nldt : list of rates of change in C_nl for each cell (1/yr)

    Output(s): None

    Note(s): assumes that the decay lifetimes are correctly set
    */
    Cell * top_cell = this->cells[this->num_cells-1]; // get top cell for easy access
    size_t num_sat_types = top_cell->num_sat_types; // get number of types
    size_t num_rb_types = top_cell->num_rb_types;

    // initialize collision arrays
    ArrayND<double, 2> S_coll = ArrayND<double, 2>(array<size_t,2>({num_sat_types, num_sat_types})); // sat-sat
    ArrayND<double, 2> RS_coll = ArrayND<double, 2>(array<size_t,2>({num_sat_types, num_rb_types})); // sat-rb
    ArrayND<double, 2> R_coll = ArrayND<double, 2>(array<size_t,2>({num_rb_types, num_rb_types})); // rb-rb
    ArrayND<double, 3> NS_coll = ArrayND<double, 3>(array<size_t,3>({num_sat_types, this->num_L, this->num_chi})); // sat-debris
    ArrayND<double, 3> NR_coll = ArrayND<double, 3>(array<size_t,3>({num_rb_types, this->num_L, this->num_chi})); // rb-debris
    // intialize explosion arrays
    Array1D<double> NS_expl = Array1D<double>(num_sat_types); // sat
    Array1D<double> NR_expl = Array1D<double>(num_rb_types); // rb
    // initialize S_in, S_din, D_in, R_in, N_in values
    Array1D<double> ** S_in = new Array1D<double>*[this->num_cells+1];
    Array1D<double> ** S_din = new Array1D<double>*[this->num_cells+1];
    Array1D<double> ** D_in = new Array1D<double>*[this->num_cells+1];
    Array1D<double> ** R_in = new Array1D<double>*[this->num_cells+1];
    ArrayND<double, 2> ** N_in = new ArrayND<double, 2>*[this->num_cells+1];
    for (size_t i = 0; i < this->num_cells+1; i++) {
        S_in[i] = new Array1D<double>(0.0, num_sat_types);
        S_din[i] = new Array1D<double>(0.0, num_sat_types);
        D_in[i] = new Array1D<double>(0.0, num_sat_types);
        R_in[i] = new Array1D<double>(0.0, num_rb_types);
        N_in[i] = new ArrayND<double,2>(0.0, array<size_t,2>({this->num_L, this->num_chi}));
    } for (size_t i = 0; i < num_sat_types; i++) {
        S_in[0]->at(i) = top_cell->lam_s->at(i);
    } if (upper) {
        for (size_t i = 0; i < this->num_L; i++) {
            for (size_t j = 0; j < this->num_chi; j++) {
                array<size_t, 2> index = {i,j};
                N_in[this->num_cells]->at(index) = this->N_upper_init->at(index)/this->cells[this->num_cells-1]->tau_N->at(j);
            }
        }
    }

    // iterate through cells, from bottom to top
    for (size_t i = 0; i < this->num_cells; i++) {
        Cell * curr_cell = this->cells[i]; // get current cell
        // zero the re-used arrays
        S_coll.zero(); RS_coll.zero(); R_coll.zero(); NS_coll.zero(); NR_coll.zero(); NS_expl.zero(); NR_expl.zero();
        // get rates of change for the current cell 
        curr_cell->dxdt_cell(time, dSdt[i], dS_ddt[i], dDdt[i], dRdt[i], dC_ldt[i], dC_nldt[i], *S_in[i+1], *S_din[i], 
                             *D_in[i], *R_in[i], *N_in[i], S_coll, RS_coll, R_coll, NS_coll, NR_coll, NS_expl, NR_expl);
        // update catastrophic and non-catastrophic collisions
        dC_ldt[i] += S_coll.sum_arr() + R_coll.sum_arr() + RS_coll.sum_arr();
        for (size_t j = 0; j < this->num_L; j++) {
            for (size_t k = 0; k < this->num_chi; k++) {
                for (size_t l = 0; l < num_sat_types; l++) {
                    if (curr_cell->cat_sat_N->at(array<size_t,3>({l,j,k}))) {
                        dC_ldt[i] += NS_coll.at(array<size_t,3>({l,j,k}));
                    } else {
                        dC_nldt[i] += NS_coll.at(array<size_t,3>({l,j,k}));
                    }
                } for (size_t l = 0; l < num_rb_types; l++) {
                    if (curr_cell->cat_rb_N->at(array<size_t,3>({l,j,k}))) {
                        dC_ldt[i] += NR_coll.at(array<size_t,3>({l,j,k}));
                    } else {
                        dC_nldt[i] += NR_coll.at(array<size_t,3>({l,j,k}));
                    }
                }
            }
        }
        // simulate collisions and explosions
        for (size_t j = 0; j < num_sat_types; j++) { // iterate through satellite types
            double m_s1 = curr_cell->m_s->at(j); // mass of the first satellite
            double C = curr_cell->C_s->at(j); // explosion constant of first satellite

            for (size_t k = j+1; k < num_sat_types; k++) { // sat-sat collisions
                double m_s2 = curr_cell->m_s->at(k); // mass of the second satellite
                this->sim_colls(dNdt, S_coll.at(array<size_t,2>({j,k})) + S_coll.at(array<size_t,2>({k,j})),
                                m_s1, m_s2, i, 's');
            } this->sim_colls(dNdt, S_coll.at(array<size_t,2>({j,j})), m_s1, m_s1, i, 's'); // collision of same type

            for (size_t k = 0; k < num_rb_types; k++) { // sat-rb collision
                double m_rb2 = curr_cell->m_rb->at(k); // get rocket mass
                this->sim_colls_satrb(dNdt, RS_coll.at(array<size_t,2>({j,k})), m_s1, i, 's');
                this->sim_colls_satrb(dNdt, RS_coll.at(array<size_t,2>({j,k})), m_rb2, i, 'r');
            }

            for (size_t k = 0; k < this->num_L; k++) { // sat-debris collisions
                double ave_L = pow(10, this->logL_ave->at(k)); // get average L value for this bin
                double A = find_A(ave_L); // average surface area of the debris
                for (size_t l = 0; l < this->num_chi; l++) {
                    double ave_AM = pow(10, this->chi_ave->at(l)); // get average debris AM for this bin
                    double m_d = A/ave_AM;
                    this->sim_colls(dNdt, NS_coll.at(array<size_t,3>({j,k,l})), m_s1, m_d, i, 's');
                }
            }

            this->sim_expl(dNdt, NS_expl.at(j), C, i, 's'); // simulate explosions for sat
        } for (size_t j = 0; j < num_rb_types; j++) { // iterate through rb types
            double m_rb1 = curr_cell->m_rb->at(j); // mass of the first rocket
            double C = curr_cell->C_rb->at(j); // current rocket explosion constant

            for (size_t k = j+1; k < num_rb_types; k++) { // rb-rb collisions
                double m_rb2 = curr_cell->m_rb->at(k); // mass of the second rocket
                this->sim_colls(dNdt, R_coll.at(array<size_t,2>({j,k})) + R_coll.at(array<size_t,2>({k,j})),
                                m_rb1, m_rb2, i, 'r');
            } this->sim_colls(dNdt, R_coll.at(array<size_t,2>({j,j})), m_rb1, m_rb1, i, 'r'); // collision of same type

            for (size_t k = 0; k < this->num_L; k++) { // rb-debris collisions
                double ave_L = pow(10, this->logL_ave->at(k)); // get average L value for this bin
                double A = find_A(ave_L); // average surface area of the debris
                for (size_t l = 0; l < this->num_chi; l++) {
                    double ave_AM = pow(10, this->chi_ave->at(l)); // get average debris AM for this bin
                    double m_d = A/ave_AM;
                    this->sim_colls(dNdt, NR_coll.at(array<size_t,3>({j,k,l})), m_rb1, m_d, i, 'r');
                }
            }
            
            this->sim_expl(dNdt, NR_expl.at(j), C, i, 'r'); // simulate explosions for rockets
        }

        // add on debris lost to collisions
        for (size_t j = 0; j < this->num_L; j++) {
            for (size_t k = 0; k < this->num_chi; k++) {
                for (size_t l = 0; l < num_sat_types; l++) {
                    dNdt.at(array<size_t,3>({i,j,k})) -= NS_coll.at(array<size_t,3>({l,j,k}));
                } for (size_t l = 0; l < num_rb_types; l++) {
                    dNdt.at(array<size_t,3>({i,j,k})) -= NR_coll.at(array<size_t,3>({l,j,k}));
                }
            }
        }
    }

    // go through cells from bottom to top to correct values
    for (size_t i = 0; i < this->num_cells; i++) {
        for (size_t j = 0; j < num_sat_types; j++) {
            dSdt[i].at(j) += S_in[i]->at(j) - S_in[i+1]->at(j);
            dS_ddt[i].at(j) += S_din[i+1]->at(j) - S_din[i]->at(j);
            dDdt[i].at(j) += D_in[i+1]->at(j) - D_in[i]->at(j);
        } for (size_t j = 0; j < num_rb_types; j++) {
            dRdt[i].at(j) += R_in[i+1]->at(j) - R_in[i]->at(j);
        } for (size_t j = 0; j < this->num_L; j++) {
            for (size_t k = 0; k < this->num_chi; k++) {
                dNdt.at(array<size_t,3>({i,j,k})) += N_in[i+1]->at(array<size_t,2>({j,k})) - N_in[i]->at(array<size_t,2>({j,k}));
            }
        }
    }
}

template <class URNG>
void NCell<URNG>::run_sim_euler(double T, double dt, bool upper) {
    /*
    runs a simulation using the Euler method

    Input(s):
    T : time to run to (yr)
    dt : time step (yr)
    upper : whether or not to have debris flow into top of system

    Output(s): None
    */
    this->sim_events(); // events that happen at the start
    // get number of object types
    size_t num_sat_types = this->cells[0]->num_sat_types; size_t num_rb_types = this->cells[0]->num_rb_types;
    // setup rate of change arrays
    Array1D<double> * dSdt = new Array1D<double>[this->num_cells];
    Array1D<double> * dS_ddt = new Array1D<double>[this->num_cells];
    Array1D<double> * dDdt = new Array1D<double>[this->num_cells];
    Array1D<double> * dRdt = new Array1D<double>[this->num_cells];
    ArrayND<double,3> dNdt = ArrayND<double,3>(array<size_t,3>({this->num_cells,this->num_L,this->num_chi}));
    double * dC_ldt = new double[this->num_cells];
    double * dC_nldt = new double[this->num_cells];
    for (size_t i = 0; i < this->num_cells; i++) {
        dSdt[i] = Array1D<double>(num_sat_types);
        dS_ddt[i] = Array1D<double>(num_sat_types);
        dDdt[i] = Array1D<double>(num_sat_types);
        dRdt[i] = Array1D<double>(num_rb_types);
    }

    while (this->t->at(this->time) < T) {
        // check if we need to update lifetimes
        if (this->t->at(this->time) - this->t->at(this->lupdate_time) >= this->update_period) {
            this->update_lifetimes(this->t->at(this->time));
            this->lupdate_time = this->time;
        }

        // clear re-used arrays
        for (size_t i = 0; i < this->num_cells; i++) {
            dSdt[i].zero(); dS_ddt[i].zero(); dDdt[i].zero(); dRdt[i].zero();
            dC_ldt[i] = 0.0; dC_nldt[i] = 0.0;
        } dNdt.zero();

        // get current rates of change
        this->dxdt(this->time, upper, dSdt, dS_ddt, dDdt, dRdt, dNdt, dC_ldt, dC_nldt);
        
        for (size_t i = 0; i < this->num_cells; i++) { // iterate through cells and update values
            // get current cell and setup new values
            Cell * curr_cell = this->cells[i];
            Array1D<double> * S_new = new Array1D<double>(0.0, num_sat_types);
            Array1D<double> * S_dnew = new Array1D<double>(0.0, num_sat_types);
            Array1D<double> * D_new = new Array1D<double>(0.0, num_sat_types);
            Array1D<double> * R_new = new Array1D<double>(0.0, num_rb_types);
            ArrayND<double,2> * N_new = new ArrayND<double,2>(0.0, array<size_t,2>({this->num_L,this->num_chi}));
            // calculate new values
            for (size_t j = 0; j < num_sat_types; j++) {
                S_new->at(j) += curr_cell->S->at(this->time)->at(j) + dSdt[i].at(j)*dt;
                S_dnew->at(j) += curr_cell->S_d->at(this->time)->at(j) + dS_ddt[i].at(j)*dt;
                D_new->at(j) += curr_cell->D->at(this->time)->at(j) + dDdt[i].at(j)*dt;
            } for (size_t j = 0; j < num_rb_types; j++) {
                R_new->at(j) += curr_cell->R->at(this->time)->at(j) + dRdt[i].at(j)*dt;
            } for (size_t j = 0; j < this->num_L; j++) {
                for (size_t k = 0; k < this->num_chi; k++) {
                    array<size_t, 2> index = {j,k};
                    N_new->at(index) += curr_cell->N_bins->at(this->time)->at(index) + dNdt.at(array<size_t,3>({i,j,k}))*dt;
                }
            }
            // update values
            curr_cell->S->push_back(S_new);
            curr_cell->S_d->push_back(S_dnew);
            curr_cell->D->push_back(D_new);
            curr_cell->R->push_back(R_new);
            curr_cell->N_bins->push_back(N_new);
            curr_cell->C_l->push_back(curr_cell->C_l->at(this->time) + dC_ldt[i]*dt);
            curr_cell->C_nl->push_back(curr_cell->C_nl->at(this->time) + dC_nldt[i]*dt);
        }

        // update times and run events
        this->t->push_back(this->t->at(this->time) + dt);
        this->time += 1; this->sim_events();
    }
}

template <class URNG>
void NCell<URNG>::sim_colls(ArrayND<double,3> &dNdt, double rate, double m1, double m2, size_t index, char typ) {
    /*
    updates dNdt by distributing a rate of collisions between two objects of mass m_1, m_2 in
    the index'th cell
    
    Input(s):
    dNdt : current dNdt values (1/yr)
    rate : rate of collisions to simulate (1/yr)
    m1 : mass of the first object (kg)
    m2 : mass of the second object (kg)
    index : index of the cell the collision occurs in
    typ : object type of the main (first) object, either 's' (satellite) or 'r' (rocket body)

    Output(s): None
    */

    if (rate == 0.0) {return;} // just skip everything if you can
    double v_rel = this->cells[index]->v; // collision velocity (km/s)
    double M = calc_M(m1, m2, v_rel); // M factor
    // min and max characteristic lengths
    double Lmin = pow(10, this->logL_edges->at(0)); double Lmax = pow(10, this->logL_edges->at(this->num_L));
    double N_debris = calc_Ntot(M, Lmin, Lmax, 'c', 0.0)*rate; // total rate of debris creation
    ArrayND<double, 4> * prob_table; // pointer to probability table to use
    if (typ == 's') {prob_table = this->sat_coll_prob_tables;} // get right probability table
    else {prob_table = this->rb_coll_probability_tables;}
    for (size_t i = 0; i < this->num_cells; i++) { // iterate through cells to send debris to
        for (size_t j = 0; j < this->num_L; j++) { // iterate through bins
            for (size_t k = 0; k < this->num_chi; k++) {
                dNdt.at(array<size_t,3>({i,j,k})) += N_debris*(prob_table->at(array<size_t,4>({index,i,j,k})));
            }
        }
    }
}

template <class URNG>
void NCell<URNG>::sim_colls_satrb(ArrayND<double,3> &dNdt, double rate, double m, size_t index, char typ) {
    /*
    version of sim_coll used for the satellite-rocket body collisions workaround, where
    each object is simulated as having its own catastrophic collision
    
    Input(s):
    dNdt : current dNdt values (1/yr)
    rate : rate of collisions to simulate (1/yr)
    m : mass of the object (kg)
    index : index of the cell the collision occurs in
    typ : object type in the collision, either 's' (satellite) or 'r' (rocket body)

    Output(s): None
    */

    if (rate == 0.0) {return;} // just skip everything if you can
    // min and max characteristic lengths
    double Lmin = pow(10, this->logL_edges->at(0)); double Lmax = pow(10, this->logL_edges->at(this->num_L));
    double N_debris = calc_Ntot(m, Lmin, Lmax, 'c', 0.0)*rate; // total rate of debris creation
    ArrayND<double, 4> * prob_table; // pointer to probability table to use
    if (typ == 's') {prob_table = this->sat_coll_prob_tables;} // get right probability table
    else {prob_table = this->rb_coll_probability_tables;}
    for (size_t i = 0; i < this->num_cells; i++) { // iterate through cells to send debris to
        for (size_t j = 0; j < this->num_L; j++) { // iterate through bins
            for (size_t k = 0; k < this->num_chi; k++) {
                dNdt.at(array<size_t,3>({i,j,k})) += N_debris*(prob_table->at(array<size_t,4>({index,i,j,k})));
            }
        }
    }
}

template <class URNG>
void NCell<URNG>::sim_expl(ArrayND<double,3> &dNdt, double rate, double C, size_t index, char typ) {
    /*
    updates dNdt by distributing a rate of explosions for an object with constant C in
    the index'th cell
    
    Parameter(s):
    dNdt : current dNdt values (1/yr)
    rate : rate of explosions to simulate (1/yr)
    C : fit constant for the explosion
    index : index of the cell the collision occurs in
    typ : object type of the main object, either 's' (satellite) or 'r' (rocket body)

    Output(s): None
    */

    if (rate == 0.0) {return;} // just skip everything if you can
    // min and max characteristic lengths
    double Lmin = pow(10, this->logL_edges->at(0)); double Lmax = pow(10, this->logL_edges->at(this->num_L));
    double N_debris = calc_Ntot(0.0, Lmin, Lmax, 'e', C)*rate; // total rate of debris creation
    ArrayND<double, 4> * prob_table; // pointer to probability table to use
    if (typ == 's') {prob_table = this->sat_expl_prob_tables;} // get right probability table
    else {prob_table = this->rb_expl_probability_tables;}
    for (size_t i = 0; i < this->num_cells; i++) { // iterate through cells to send debris to
        for (size_t j = 0; j < this->num_L; j++) { // iterate through bins
            for (size_t k = 0; k < this->num_chi; k++) {
                dNdt.at(array<size_t,3>({i,j,k})) += N_debris*(prob_table->at(array<size_t,4>({index,i,j,k})));
            }
        }
    }
}

template <class URNG>
void NCell<URNG>::sim_events() {
    // sims all events that need to be run in the system at the current time

    // setup arrays to hold changes in debris
    ArrayND<double, 3> dN = ArrayND<double,3>(0.0, array<size_t,3>({this->num_cells,this->num_L,this->num_chi}));
    ArrayND<double, 2> dN_loc = ArrayND<double,2>(array<size_t,3>({this->num_L,this->num_chi})); // for non-collision/explosion sources
    // setup arrays to hold changes in everything else
    size_t num_sat_types = this->cells[0]->num_sat_types; size_t num_rb_types = this->cells[0]->num_rb_types;
    Array1D<double> dS = Array1D<double>(num_sat_types);
    Array1D<double> dS_d = Array1D<double>(num_sat_types);
    Array1D<double> dD = Array1D<double>(num_sat_types);
    Array1D<double> dR = Array1D<double>(num_rb_types);

    for (size_t i = 0; i < this->num_cells; i++) { // iterate through cells
        Cell * curr_cell = this->cells[i]; // get current cell
        // zero the re-used arrays
        dS.zero(); dS_d.zero(); dD.zero(); dR.zero(); dN_loc.zero();
        vector<Coll> coll_list = vector<Coll>();
        vector<Expl> expl_list = vector<Expl>();
        // get current values
        Array1D<double> * S = curr_cell->S->at(i);
        Array1D<double> * S_d = curr_cell->S_d->at(i);
        Array1D<double> * D = curr_cell->D->at(i);
        Array1D<double> * R = curr_cell->R->at(i);
        ArrayND<double,2> * N = curr_cell->N_bins->at(i);
        
        for (size_t j = 0; j < curr_cell->num_events; j++) { // iterate through events
            Event * event = curr_cell->event_list->at(j); // get current event
            
            // handle specific time events
            while ((event->times->size() != 0) && (event->times->back() <= this->t->at(this->time))) {
                event->run_event(*S, *S_d, *D, *R, *N, *this->logL_edges, *this->chi_edges, dS, dS_d, dD, dR, dN_loc, coll_list, expl_list);
                event->times->pop_back();
            }
            // handle frequency-based events
            if (event->freq != 0) {
                if (this->t->at(this->time) - event->last_event >= 1.0/(event->freq)) {
                    event->run_event(*S, *S_d, *D, *R, *N, *this->logL_edges, *this->chi_edges, dS, dS_d, dD, dR, dN_loc, coll_list, expl_list);
                    event->last_event = this->t->at(this->time);
                }
            }
        }

        // update values
        for (size_t j = 0; j < num_sat_types; j++) {
            curr_cell->S->at(this->time)->at(j) += dS.at(j);
            curr_cell->S_d->at(this->time)->at(j) += dS_d.at(j);
            curr_cell->D->at(this->time)->at(j) += dD.at(j);
        } for (size_t j = 0; j < num_rb_types; j++) {
            curr_cell->R->at(this->time)->at(j) += dR.at(j);
        } for (size_t j = 0; j < this->num_L; j++) {
            for (size_t k = 0; k < this->num_chi; k++) {
                array<size_t,2> index = {i,j};
                curr_cell->N_bins->at(this->time)->at(index) += dN_loc.at(index);
            }
        }

        // handle collisions and explosions
        this->parse_coll(coll_list, dN, i); this->parse_expl(expl_list, dN, i);
    }

    // update with debris from collisions/explosions
    for (size_t i = 0; i < this->num_cells; i++) {
        Cell * curr_cell = this->cells[i];
        for (size_t j = 0; j < this->num_L; j++) {
            for (size_t k = 0; k < this->num_chi; k++) {
                array<size_t,3> index = {i,j,k};
                curr_cell->N_bins->at(this->time)->at(array<size_t,2>({j,k})) += dN.at(index);
            }
        }
    }
}

template <class URNG>
void NCell<URNG>::parse_coll(vector<Coll> &coll_list, ArrayND<double,3> &dN, size_t index) {
    /*
    parses and runs discrete collision events, storing the debris generated in dN

    Input(s):
    coll_list : list of collisions occuring in the current cell
    dN : 3d matrix of changes in debris for each bin and cell
    index : index of the current cell

    Output(s): none
    */
    Coll coll; // collision object
    for (size_t i = 0; i < coll_list.size(); i++) { // iterate through collisions
        coll = coll_list[i];
        if (coll.typ == 's' || coll.typ == 'r') { // standard collisions
            this->sim_colls(dN, coll.num, coll.m1, coll.m2, index, coll.typ);
        } else { // satellite-rocket body collision
            this->sim_colls_satrb(dN, coll.num, coll.m1, index, 's');
            this->sim_colls_satrb(dN, coll.num, coll.m2, index, 'r');
        }
    }
}

template <class URNG>
void NCell<URNG>::parse_expl(vector<Expl> &expl_list, ArrayND<double,3> &dN, size_t index) {
    /*
    parses and runs discrete explosion events, storing the debris generated in dN

    Input(s):
    expl_list : list of explosions occuring in the current cell
    dN : 3d matrix of changes in debris for each bin and cell
    i : index of the current cell

    Output(s): None
    */
    Expl expl; // explosion object
    for (size_t i = 0; i < expl_list.size(); i++) { // iterate through explosions
        expl = expl_list[i];
        this->sim_expl(dN, expl.num, expl.C, index, expl.typ);
    }
}

template <class URNG>
void NCell<URNG>::update_lifetimes(double t) {
    /*
    updates decay lifetimes in all cells in the system

    Input(s):
    t : time since the start of the simulation (yr)

    Output(s): None
    */
    Cell * curr_cell; double alt; double dh; double AM;
    double m0 = t*12.0; // compute starting month
    for (size_t i = 0; i < this->num_cells; i++) { // iterate through cells
        curr_cell = (this->cells)[i]; alt = curr_cell->alt; dh = curr_cell->dh;
        for (size_t j = 0; j < curr_cell->num_sat_types; j++) { // handle satellites
            AM = curr_cell->AM_s->at(j);
            curr_cell->tau_s->at(j) = max(this->min_lifetime,drag_lifetime_default(alt + dh/2, alt - dh/2, AM, m0));
        }
        for (size_t j = 0; j < curr_cell->num_rb_types; j++) { // handle rockets
            AM = curr_cell->AM_rb->at(j);
            curr_cell->tau_rb->at(j) = max(this->min_lifetime, drag_lifetime_default(alt + dh/2, alt - dh/2, AM, m0));
        }
        for (size_t j = 0; j < this->num_chi; j++) {// handle debris
            AM = pow(10, this->chi_ave->at(j));
            curr_cell->tau_N->at(j) = max(this->min_lifetime, drag_lifetime_default(alt + dh/2, alt - dh/2, AM, m0));
        }
    }
}

template <class URNG>
size_t NCell<URNG>::alt_to_index(double h) {
    /*
    converts altitude to index of corresponding cell

    Input(s):
    h : altitude (km)

    Output(s):
    indx : cell index

    Note(s): returns this->num_cells if the alt isn't in any cell
    */
    double alt; double dh;
    for (size_t i = 0; i < this->num_cells; i++) {
        alt = this->alts->at(i); dh = this->dh->at(i);
        if ((alt - dh/2 <= h) && (alt + dh/2 >= h)) {return i;}
    }
    return this->num_cells;
}

template <class URNG>
NCell<URNG>::~NCell() {
    // destructor
    delete alts; delete dhs; delete t;
    for (size_t i = 0; i < this->num_cells; i++) {delete this->cells[i];}
    delete logL_edges; delete chi_edges; delete logL_ave; delete chi_ave;
    delete sat_coll_prob_tables; delete sat_expl_prob_tables;
    delete rb_coll_prob_tables; delete rb_expl_prob_tables;
}