#include "AtmosphericDecayModels.h"
#include <iostream>

using namespace std;

void atmospheric_test() {
    cout << drag_lifetime_default(600.0+25.0/2.0, 600.0-25.0/2.0, 1.0/40.0, 1.0) << endl;
}

int main() {
    atmospheric_test();
    return 0;
}