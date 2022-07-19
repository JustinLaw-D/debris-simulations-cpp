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

    public:
        NCell(string &filepath, size_t num_dir, URNG &generator); // constructor based on loading data from file
        void save(string &filepath, string &name, double gap); // function for saving data
        // calculates the probability tables
        void calc_prob_table(ArrayND<double, 4> * table, size_t indx, char etyp, char ttyp, size_t num_dir, URNG &generator);
        bool add_event(Event * event, double alt); // adds event to the system
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
        param_file << "\"" << to_string(this->update_period) << "\"";
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
void NCell<URNG>::sim_colls(ArrayND<double,3> &dNdt, double rate, double m1, double m2, size_t index, char typ) {
    /*
    updates dNdt by distributing a rate of collisions between two objects of mass m_1, m_2 in
    the index'th cell
    
    Input(s):
    dNdt : current dNdt values (1/yr)
    rate : rate of collisions to simulate (1/yr)
    m_1 : mass of the first object (kg)
    m_2 : mass of the second object (kg)
    index : index of the cell the collision occurs in
    typ : object type of the main (first) object, either 's' (satellite) or 'r' (rocket body)

    Output(s): None
    */

    if (rate == 0.0) {return;} // just skip everything if you can
    double v_rel = this->cells[index]->v; // collision velocity (km/s)
    double M = calc_M(m_1, m_2, v_rel); // M factor
    // min and max characteristic lengths
    double L_min = pow(10, this->logL_edges->at(0)); double L_max = pow(10, this->logL_edges->at(this->num_L));
    double N_debris = calc_Ntot(M, Lmin, Lmax, 'c', 0.0)*rate; // total rate of debris creation
    ArrayND<double, 4> * prob_table; // pointer to probability table to use
    if (typ == 's') {prob_table = this->sat_coll_prob_tables;} // get right probability table
    else {prob_table = self.rb_coll_probability_tables;}
    for (size_t i = 0; i < this->num_cells; i++) { // iterate through cells to send debris to
        for (size_t j = 0; j < this->num_L; j++) { // iterate through bins
            for (size_t k = 0; k < this->num_chi; k++) {
                dNdt->at(array<size_t,3>({i,j,k})) += N_debris*(prob_table->at(array<size_t,4>({index,i,j,k})));
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
    double L_min = pow(10, this->logL_edges->at(0)); double L_max = pow(10, this->logL_edges->at(this->num_L));
    double N_debris = calc_Ntot(m, Lmin, Lmax, 'c', 0.0)*rate // total rate of debris creation
    ArrayND<double, 4> * prob_table; // pointer to probability table to use
    if (typ == 's') {prob_table = this->sat_coll_prob_tables;} // get right probability table
    else {prob_table = self.rb_coll_probability_tables;}
    for (size_t i = 0; i < this->num_cells; i++) { // iterate through cells to send debris to
        for (size_t j = 0; j < this->num_L; j++) { // iterate through bins
            for (size_t k = 0; k < this->num_chi; k++) {
                dNdt->at(array<size_t,3>({i,j,k})) += N_debris*(prob_table->at(array<size_t,4>({index,i,j,k})));
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
    double L_min = pow(10, this->logL_edges->at(0)); double L_max = pow(10, this->logL_edges->at(this->num_L));
    double N_debris = calc_Ntot(0.0, Lmin, Lmax, 'e', C)*rate // total rate of debris creation
    ArrayND<double, 4> * prob_table; // pointer to probability table to use
    if (typ == 's') {prob_table = this->sat_expl_prob_tables;} // get right probability table
    else {prob_table = self.rb_expl_probability_tables;}
    for (size_t i = 0; i < this->num_cells; i++) { // iterate through cells to send debris to
        for (size_t j = 0; j < this->num_L; j++) { // iterate through bins
            for (size_t k = 0; k < this->num_chi; k++) {
                dNdt->at(array<size_t,3>({i,j,k})) += N_debris*(prob_table->at(array<size_t,4>({index,i,j,k})));
            }
        }
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
            curr_cell->tau_s->at(j) = drag_lifetime_default(alt + dh/2, alt - dh/2, AM, m0);
        }
        for (size_t j = 0; j < curr_cell->num_rb_types; j++) { // handle rockets
            AM = curr_cell->AM_rb->at(j);
            curr_cell->tau_rb->at(j) = drag_lifetime_default(alt + dh/2, alt - dh/2, AM, m0);
        }
        for (size_t j = 0; j < this->num_chi; j++) {// handle debris
            AM = pow(10, this->chi_ave->at(j));
            curr_cell->tau_N->at(j) = drag_lifetime_default(alt + dh/2, alt - dh/2, AM, m0);
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