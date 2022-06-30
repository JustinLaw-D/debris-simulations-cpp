// implementation of the NASA standard breakup model

#pragma once

// whether or not a collision is catestrophic
bool is_catastrophic(const double m_s, const double L, const double AM, const double v);

double find_A(const double L); // finds the area given characteristic length