// implementation of ObjectsEvents.h

#include "Arrays.h"
#include "ObjectsEvents.h"
#include <vector>

void default_event(const Array2D<double>& S, const Array2D<double>& S_d, const Array2D<double>& D, size_t num_sat_types, 
                   const Array2D<double>& R, size_t num_rb_types, const Array2D<double>& N, const double* logL_bins, 
                   size_t num_L, const double* chi_bins, size_t num_chi, Array2D<double>& dS, Array2D<double>& dS_d, 
                   Array2D<double>& dD, Array2D<double>& dR, Array2D<double>& dN, vector<Coll>& coll_list, 
                   vector<Expl>& expl_list) {
    /*
    dummy default event function, does nothing

    Input(s):
    S : number of live satellites of each type
    S_d : number of de-orbiting satellites of each type
    D : number of derelict satellites of each type
    num_sat_types : the number of satellite types
    R : number of rocket bodies of each type
    num_rb_types : the number of rocket body types
    N : binned number of debris
    logL_edges : bin edges in logL
    num_L : number of bins in L
    chi_edges : bin edges in chi
    num_chi : number of bins in chi
    dS : change in S due to the event
    dS_d : change in S_d due to the event
    dD : change in D due to the event
    dR : change in R due to the event
    dN : change in N due to the event
    coll_list : list of collisions occuring in the event
    expl_list : list of explosions occuring in the event

    Output(s): None
    */
    return;
}