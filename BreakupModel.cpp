// implementation of BreakupModel.h

#include <cmath>
#include "BreakupModel.h"
#include "Constants.h"

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
    
    double v_loc = v*1000.0; // convert to m/s
    double m_d = find_A(L)/AM; // mass of the debris (on average)
    double k_d = 0.5*m_d*(v_loc*v_loc); // relative kinetic energy of the debris
    double dec_fact = (k_d/m)/1000.0; // factor for making the decision (J/g)
    return (dec_fact >= 40.0);
}

double calc_M(double m_t, double m_d, double v) {
    /*
    calculates the M factor used for L distribution calculation
    
    Input(s):
    m_t : target mass (kg)
    m_d : mass of the debris/projectile (kg)
    v : collision velocity (km/s)

    Output(s):
    M : value of M parameter (variable units)
    */

    double E_p = (0.5*m_d*((v*1000.0)*(v*1000.0))/m_t)/1000.0; // E_p in J/g

    if (E_p >= 40.0) {return m_t + m_d;} // catestrophic collision
    else {return m_d*v;} // non-catestrophic collision
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

double calc_Ntot(double M, double Lmin, double Lmax, char typ, double C) {
    /*
    calculates the total number of debris produced with characteristic length
    between Lmin and Lmax

    Input(s):
    M : fit parameter given by calc_M, ignored for explosions (variable units)
    Lmin : minimum characteristic length (m)
    Lmax : maximum characteristic length (m)
    typ : one of 'c' (collision) or 'e' (explosion)
    C : fit parameter for explosions, ignored for collisions

    Output(s):
    N : total number of fragments of Lmax > size > Lmin
    */

    if (typ == 'c') {
        return 0.1*pow(M, 0.75)*(pow(Lmin, -1.71) - pow(Lmax, -1.71));
    } else {
        return 6.0*C*(pow(Lmin,-1.6) - pow(Lmax,-1.6));
    }
}

double L_cdf(double L, double L_min, double L_max, char typ) {
    /*
    calculates the cumulative distribution function for characteristic lengths
    from a collision/explosion at length L, assuming the distribution is truncated 
    at L_min and L_max

    Input(s):
    L : characteristic length (m)
    L_min : minimum characteristic length to consider (m)
    L_max : maximum characteristic length to consider (m)
    typ : one of 'c' (collision) or 'e' (explosion)

    Output(s):
    P : value of CDF at L
    */
    double beta;
    if (typ == 'c') {beta = -1.71;}
    else {beta = -1.6;}
    return (pow(L_min,beta) - pow(L,beta))/(pow(L_min,beta) - pow(L_max,beta));
}

double X_cdf(double x, double x_min, double x_max, double L, char typ) {
    /*
    calculates the cumulative distribution function for log10(A/M) from a
    collision/explosion at value x and length L, assuming the distribution
    is truncated at x_min and x_max

    Input(s):
    x : log10(A/M) (log10(m^2/kg))
    x_min : minimum log10(A/M) value to consider (log10(m^2/kg))
    x_max : maximum log10(A/M) value to consider (log10(m^2/kg))
    L : characteristic length of the debris (m)
    typ : one of 's' (satellite) or 'r' (rocket body)

    Output(s):
    P : value of CDF at x, L
    */
    if (L >= 11.0/100.0) {return _X_cdf_11(x, x_min, x_max, L, typ);}
    else if (L <= 8.0/100.0) {return _X_cdf_8(x, x_min, x_max, L);}
    else {
        double lam_min = log10(8.0/100.0);
        double lam_max = log10(11.0/100.0);
        double P = (log10(L)-lam_min)/(lam_max-lam_min);
        return P*_X_cdf_11(x, x_min, x_max, L, typ) + (1-P)*_X_cdf_8(x, x_min, x_max, L);
    }
}

double _X_cdf_8(double x, double x_min, double x_max, double L) {
    /*
    calculates the cumulative distribution function for log10(A/M) from both collisions
    and explosions at value x and length L<=8cm, assuming the distribution is truncated 
    at x_min and x_max

    Input(s):
    x : log10(A/M) (log10(m^2/kg))
    x_min : minimum log10(A/M) value to consider (log10(m^2/kg))
    x_max : maximum log10(A/M) value to consider (log10(m^2/kg))
    L : characteristic length of the debris (m)

    Output(s):
    P : value of CDF at x, L
    */
    double lam = log10(L); double mu; double sigma;

    // determine normal distribution parameters
    if (lam <= -1.75) {mu = -0.3;}
    else if (lam < -1.25) {mu = -0.3 - 1.4*(lam + 1.75);}
    else {mu = -1.0;}

    if (lam <= -3.5) {sigma = 0.2;}
    else {sigma = 0.2 + 0.1333*(lam + 3.5);}
    
    // compute normalization factor
    double C = 1.0/(erf((x_max-mu)/(SQRT2*sigma)) - erf((x_min-mu)/(SQRT2*sigma)));
    // compute total distribution
    return C*(erf((x-mu)/(SQRT2*sigma)) - erf((x_min-mu)/(SQRT2*sigma)));
}

double _X_cdf_11(double x, double x_min, double x_max, double L, char typ) {
    /*
    calculates the cumulative distribution function for log10(A/M) from collisions/explosions 
    at value x and length L>=11cm, assuming the distribution is truncated at x_min 
    and x_max

    Input(s):
    x : log10(A/M) (log10(m^2/kg))
    x_min : minimum log10(A/M) value to consider (log10(m^2/kg))
    x_max : maximum log10(A/M) value to consider (log10(m^2/kg))
    L : characteristic length of the debris (m)
    typ : one of 's' (satellite) or 'r' (rocket body)

    Output(s):
    P : value of CDF at x, L
    */
    double lam = log10(L); double mu1; double mu2;
    double sigma1; double sigma2; double alpha;

    // determine normal distribution parameters
    if (typ == 's') {
        if (lam <= -1.95) {alpha = 0.0;}
        else if (lam < 0.55) {alpha = 0.3 + 0.4*(lam + 1.2);}
        else {alpha = 1.0;}

        if (lam <= -1.1) {mu1 = -0.6;}
        else if (lam < 0.0) {mu1 = -0.6 - 0.318*(lam + 1.1);}
        else {mu1 = -0.95;}

        if (lam <= -1.3) {sigma1 = 0.1;}
        else if (lam < -0.3) {sigma1 = 0.1 + 0.2*(lam + 1.3);}
        else {sigma1 = 0.3;}

        if (lam <= -0.7) {mu2 = -1.2;}
        else if (lam < -0.1) {mu2 = -1.2 - 1.333*(lam + 0.7);}
        else {mu2 = -2.0;}

        if (lam <= -0.5) {sigma2 = 0.5;}
        else if (lam < -0.3) {sigma2 = 0.5 - (lam + 0.5);}
        else {sigma2 = 0.3;}
    } else {
        if (lam <= -1.4) {alpha = 1.0;}
        else if (lam < 0.0) {alpha = 1.0 - 0.3571*(lam + 1.4);}
        else {alpha = 0.5;}

        if (lam <= -0.5) {mu1 = -0.45;}
        else if (lam < 0.0) {mu1 = -0.45 - 0.9*(lam + 0.5);}
        else {mu1 = -0.9;}

        sigma1 = 0.55; mu2 = -0.9;

        if (lam <= -1.0) {sigma2 = 0.28;}
        else if (lam < 0.1) {sigma2 = 0.28 - 0.1636*(lam + 1.0);}
        else {sigma2 = 0.1;}
    }   
    
    // compute normalization factor
    double top = alpha*erf((x_max-mu1)/(SQRT2*sigma1)) + (1.0-alpha)*erf((x_max-mu2)/(SQRT2*sigma2));
    double bot = alpha*erf((x_min-mu1)/(SQRT2*sigma1)) + (1.0-alpha)*erf((x_min-mu2)/(SQRT2*sigma2));
    double C = 1.0/(top - bot);
    // compute total distribution
    double fac_one = erf((x-mu1)/(SQRT2*sigma1)) - erf((x_min-mu1)/(SQRT2*sigma1));
    double fac_two = erf((x-mu2)/(SQRT2*sigma2)) - erf((x_min-mu2)/(SQRT2*sigma2));
    return C*(alpha*fac_one + (1.0-alpha)*fac_two);
}

double v_cdf(double v, double x, char typ) {
    /*
    evaluates cdf for log10(Delta v) values at given v and x

    Input(s):
    v : log10(delta v) value to evaluate at (log10(m/s))
    x : log10(A/M) value of the debris (log10(m^2/kg))
    typ : one of 'c' (collision) or 'e' (explosion)

    Output(s):
    P : value of the CDF at v, x
    */
    double mu;
    if (typ == 'c') {mu = 0.9*x + 2.9;} // calculate normal distribution parameters
    else {mu = 0.2*x + 1.85;}
    double sigma_fac = 0.4*SQRT2;
    // calculate CDF value
    return 0.5*(erf((v-mu)/sigma_fac) + 1.0);
}

double vprime_cdf(double V, double v0, double theta, double phi, double x, char typ) {
    /*
    evaluates cdf for the post-collision speed V, given a pre-collision
    orbital speed v0, post-collision direction of (theta, phi), and x
    evaluates cdf for log10(Delta v) values at given v and x

    Input(s):
    V : post-collison speed to evaluate at (m/s)
    v0 : pre-collision orbital speed (m/s)
    theta : inclination angle (rad)
    phi : azimuthal angle (rad)
    x : log10(A/M) value of the debris (log10(m^2/kg))
    typ : one of 'c' (collision) or 'e' (explosion)

    Output(s):
    P : value of the CDF at V
    */
    double descriminate = pow((v0*sin(theta)*cos(phi)), 2.0) - (v0*v0-V*V);
    if (descriminate < 0.0) {return 0.0;} // cannot get a post-collision velocity this low
    double del_v_max = -v0*sin(theta)*cos(phi) + sqrt(descriminate);
    double del_v_min = -v0*sin(theta)*cos(phi) - sqrt(descriminate);
    if (del_v_max < 0.0) {return 0.0;} // cannot get a post-collision velocity this low
    else if (del_v_min <= 0.0) {return v_cdf(log10(del_v_max), x, typ);}
    else {return v_cdf(log10(del_v_max), x, typ) - v_cdf(log10(del_v_min), x, typ);}
}