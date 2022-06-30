// implementation of BreakupModel.h

#include <cmath>
#include "BreakupModel.h"

bool is_catastrophic(const double m, const double L, const double AM, const double v) {
    /*
    determines if a collision between debris and object is catestrophic

    Input(s):
    m : mass of the object (kg)
    L : characteristic length of the debris (m)
    AM : area to mass ratio of the debris (m^2/kg)
    v : relative velocity of the objects (km/s)

    Output(s):
    cat : true if the collision is catastrophic, false otherwise
    */
    
    double v_loc = v*1000; // convert to m/s
    double m_d = find_A(L)/AM; // mass of the debris (on average)
    double k_d = 0.5*m_d*(v_loc*v_loc); // relative kinetic energy of the debris
    double dec_fact = (k_d/m)/1000; // factor for making the decision (J/g)
    return dec_fact >= 40.0;
}

double find_A(const double L) {
    /*
    calculates the average cross-sectional area of an object with a given characteristic length

    Input(s):
    L : characteristic length (m)

    Output(s):
    A : average cross-sectional area (m^2)
    */

    if (L < 0.00167) {return 0.540424*(L*L);}
    else {return 0.556945*(pow(L,2.0047077));}
}