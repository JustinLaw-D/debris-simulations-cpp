#include <iostream>
#include <ctime>
#include "Cells.h"

using namespace std;

int main() {
    string filename = "./data/NCellDeorbitFlowTest/";
    NCell atmosphere = NCell(filename, 1000);
    time_t start_time = time(NULL);
    atmosphere.run_sim_precor(50.0, 1.0, 0.001, 1.0, 1.0, true);
    cout << time(NULL) - start_time << endl;
    string folder = "./data/"; string name = "NCellDeorbitFlowTest";
    atmosphere.save(folder, name, 0.1);
}
