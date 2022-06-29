// this is just here for basic testing

#include "Arrays.h"
#include "ObjectsEvents.h"
#include "Cell.h"
#include <iostream>
using namespace std;

void array_test() {
    Array2D<double> arr = Array2D<double>::zeroes_double(3, 3);
    arr.print();
    cout << endl;
    arr.set(10.0,1,2);
    arr.print();
    cout << endl;
    Array2D<double> arr2 = Array2D<double>(arr);
    arr2.set(50.0,0,0);
    arr2.print();
    cout << endl;
    arr.print();
    cout << endl;
    arr2.get_range(0, 2, 0, 2).print();
    cout << endl;
    Array2D<double> sum = arr + arr2;
    sum.print();
    cout << endl;
    arr.copy_sum(arr2);
    arr.print();
    cout << endl;
    Array2D<int> vec = Array2D<int>::zeroes_int(1,5);
    vec.print();
    cout << endl;
    vec.swap_dim();
    vec.print();
    cout << endl;
}


void cell_test() {
    vector<double> S = {5}; vector<double> S_d = {2}; vector<double> D = {0};
    Satellite sat = {S, S_d, D, 250.0, 10.0, 200.0, 5.0, 2.0, 550.0, 1/12, 0.0, 0.2, 0.2, 0.2,
                     0.95, 1/(20*2.2), 5.0, 1.0, 0.0, 0.0};
    vector<double> R = {10};
    RocketBody rb = {R, 1000.0, 20.0, 0.0, 1/(20*2.2), 5.0, 1.0, 1.0};
    size_t num_sat_types = 1; size_t num_rb_types = 1;
    size_t num_L = 10; size_t num_chi = 10;
    double chi_edges[11]; double logL_edges[11];
    for (size_t i = 0; i < num_L + 1; i++) {
        chi_edges[i] = -2.0 + 0.35*(static_cast<double>(i));
        logL_edges[i] = -3 + 0.3*(static_cast<double>(i));
    }
    Array2D<double> N_i = Array2D<double>::zeroes_double(num_L, num_chi);
    Array2D<double> tau_N = Array2D<double>::fill(5.0, num_L, num_chi);
    Event * event_list = new Event[1]();
    size_t num_events = 0;
    double alt = 550.0; double dh = 50.0; double v = 10.0;
    cout << "Here" << endl;
    Cell test_cell = Cell(&sat, &rb, N_i, num_sat_types, num_rb_types, logL_edges, num_L, chi_edges,
                          num_chi, event_list, num_events, alt, dh, tau_N, v);
    cout << "Cell made!" << endl;
}

int main() {
    //array_test();
    cell_test();
    return 0;
}