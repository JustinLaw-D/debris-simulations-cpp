// implementation of ObjectsEvents.h

#include "Arrays.h"
#include "ObjectsEvents.h"
#include <vector>

Coll::Coll() {
    // basic constructor
    this->m1 = 0.0; this->m2 = 0.0; this->num = 0; this->typ = 's';
}

Coll::Coll(const Coll &coll) {
    // copy constructor
    this->m1 = coll.m1; this->m2 = coll.m2;
    this->num = coll.num; this->typ = coll.typ;
}

Expl::Expl() {
    // basic constructor
    this->C = 0.0; this->num = 0; this->typ = 's';
}

Expl::Expl(const Expl &expl) {
    // copy constructor
    this->C = expl.C; this->num = expl.num; this->typ = expl.typ;
}

Event::Event(const Event &event) {
    // copy constructor
    this->alt = event.alt; this->freq = event.freq;
    this->times = new vector<double>(*event.times);
    this->run_event = run_event;
}

Event::~Event() {
    // destructor
    delete this->times;
}

void default_event(const Array1D<double>& S, const Array1D<double>& S_d, const Array1D<double>& D, const Array1D<double>& R, 
                   const ArrayND<double,2>& N, const Array1D<double>& logL_bins, const Array1D<double>& chi_bins, 
                   Array1D<double>& dS, Array1D<double>& dS_d, Array1D<double>& dD, Array1D<double>& dR, ArrayND<double,2>& dN, 
                   vector<Coll>& coll_list, vector<Expl>& expl_list) {
    /*
    dummy default event function, does nothing

    Input(s):
    S : number of live satellites of each type
    S_d : number of de-orbiting satellites of each type
    D : number of derelict satellites of each type
    R : number of rocket bodies of each type
    N : binned number of debris
    logL_edges : bin edges in logL
    chi_edges : bin edges in chi
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