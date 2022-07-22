// this is just here for basic testing

#include "Arrays.h"
#include "AtmosphericDecayModels.h"
#include "ObjectsEvents.h"
#include "Cells.h"
#include "BreakupModel.h"
#include <array>
#include <vector>
#include <iostream>
#include <time.h>
using namespace std;

void array_test() {
    const array<size_t,3> dim = {3,3,3};
    ArrayND<double, 3> arr = ArrayND<double, 3>(0.0, dim);
    arr.print();
    cout << endl;
    ArrayND<double, 2> * arr2 = new ArrayND<double, 2>(0.0, array<size_t,2>({3,2}));
    arr2->print();
    cout << endl;
    arr2->at(array<size_t,2>({1,0})) = 10.0;
    arr2->copy_sum(571.0);
    arr2->print();
    cout << endl;
    string filepath = "./test.npy";
    arr2->save(filepath);
    ArrayND<bool, 2> arr3 = ArrayND<bool, 2>(true, array<size_t, 2>({6,2}));
    arr3.print();
    cout << endl;
    ArrayND<bool, 2> * arr4 = new ArrayND<bool, 2>(arr3);
    arr4->print();
    cout << endl;
    ArrayND<float, 2> * arr5 = new ArrayND<float, 2>(array<size_t, 2>({2,2}));
    for (size_t i = 0; i < 2; i++) {
        for (size_t j = 0; j < 2; j++) {
            arr5->at(array<size_t, 2>({i,j})) = 0.0;
        }
    }
    arr5->print();
    cout << endl;
    string test_arr_loc = string("./test_save/t.npy");
    ArrayND<double, 1> arr6 = ArrayND<double, 1>(test_arr_loc);
    arr6.print();
    cout << endl;
}

void array1D_test() {
    Array1D<double> arr = Array1D<double>(0.0, 5);
    arr.print();
    cout << endl;
    Array1D<int> * arr2 = new Array1D<int>(0.0, 3);
    arr2->print();
    cout << endl;
    arr2->at(1) = 10;
    arr2->copy_sum(571);
    arr2->print();
    cout << endl;
    Array1D<bool> arr3 = Array1D<bool>(true, 6);
    arr3.print();
    cout << endl;
    Array1D<bool> * arr4 = new Array1D<bool>(arr3);
    arr4->print();
    cout << endl;
    Array1D<float> * arr5 = new Array1D<float>(10);
    for (size_t i = 0; i < 10; i++) {
        arr5->at(i) = 0.0;
    }
    arr5->print();
    cout << endl;
    string test_arr_loc = string("./test_save/t.npy");
    Array1D<double> arr6 = Array1D<double>(test_arr_loc);
    arr6.print();
    cout << endl;
}


void cell_test() {
    vector<double> S = {5}; vector<double> S_d = {2}; vector<double> D = {0};
    Satellite sat = {S, S_d, D, 250.0, 10.0, 200.0, 5.0, 2.0, 550.0, 1/12, 0.0, 0.2, 0.2, 0.2,
                     0.95, 1/(20*2.2), 5.0, 5.0, 1.0, 0.0, 0.0};
    vector<double> R = {10};
    RocketBody rb = {R, 1000.0, 20.0, 0.0, 1/(20*2.2), 5.0, 1.0, 1.0};
    size_t num_sat_types = 1; size_t num_rb_types = 1;
    size_t num_L = 10; size_t num_chi = 10;
    Array1D<double> * chi_edges = new Array1D<double>(0.0, 11);
    Array1D<double> * logL_edges = new Array1D<double>(0.0, 11);
    for (size_t i = 0; i < num_L + 1; i++) {
        chi_edges->at(i) = -2.0 + 0.35*(static_cast<double>(i));
        logL_edges->at(i) = -3 + 0.3*(static_cast<double>(i));
    }
    ArrayND<double,2> * N_i = new ArrayND<double,2>(0.0, array<size_t,2>({num_L, num_chi}));
    Array1D<double> * tau_N = new Array1D<double>(5.0, num_chi);
    vector<Event *> * event_list = new vector<Event *>();
    size_t num_events = 0;
    double alt = 550.0; double dh = 50.0; double v = 10.0;
    vector<double> * C_l = new vector<double>({0});
    vector<double> * C_nl = new vector<double>({0});
    cout << "Here" << endl;
    Cell test_cell = Cell(&sat, &rb, N_i, num_sat_types, num_rb_types, logL_edges, num_L, chi_edges,
                          num_chi, event_list, num_events, alt, dh, tau_N, v, C_l, C_nl);
    cout << "Cell made!" << endl;
}

void load_test() {
    string filepath = string("./test_save/");
    cout << "Hi" << endl;
    NCell test = NCell(filepath, 1000);
    string temp[2] = {string("./"), string("test_resave")};
    cout << "Here" << endl;
    test.save(temp[0], temp[1], 0.0);
    cout << "Made it!" << endl;
}

void atmospheric_test() {
    cout << drag_lifetime_default(600.0+25.0/2.0, 600.0-25.0/2.0, 1.0/40.0, 1.0) << endl;
}

void run_test() {
    time_t org_time;
    time_t dummy;
    string filepath = string("./test_save/");
    cout << "Hi" << endl;
    NCell test = NCell(filepath, 1000);
    cout << "Running simulation" << endl;
    cout << time(&org_time) << endl;
    test.run_sim_euler(30.0, 0.001, true);
    cout << time(&dummy)-org_time << endl;
    cout << "Finished simulation" << endl;
    string temp[2] = {string("./"), string("test_resave")};
    test.save(temp[0], temp[1], 0.0);
    cout << "Made it!" << endl;
}

int main() {
    //array_test();
    //array1D_test();
    //cell_test();
    //load_test();
    //atmospheric_test();
    run_test();
    return 0;
}