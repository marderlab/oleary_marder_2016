// OU input
#ifndef OUINPUT
#define OUINPUT
#include "conductance.h"
#include <cstdlib>
#include <ctime>
#include <cmath>

//inherit conductance class spec
class OUInput: public conductance {
protected:
    double tau_inv, sd, g_inf;
    int rmax;
    double r1, r2, out, dg;
    
    double dw(double);
public:
    
    OUInput(double g_in, double e_in, double tau_in, double sd_in)
    {
        g_inf = g_in;
        g = g_in;
        gbar = g_in;
        e_rev = e_in;
        tau_inv = 1.0/tau_in;
        sd = sd_in;
        
        // other initialization
        state_dim = 0;
        srand (time(NULL));
        rmax = RAND_MAX;
    }
    
    double get_g(void);
    double i(void);
    int get_state_dim(void);
    void integrate(double dt); //integrate = do nothing
    void connect(compstate *s_in);
    void state2double(double*);
    void double2state(double*);    
};

//over-write to prevent negative conductance going to exp-integration 'sum_g'
double OUInput::get_g(void)
{
    return g;
}

int OUInput::get_state_dim(void)
{
    return 0;
}

double OUInput::i(void)
{
    return g*(e_rev - compartment_state->get_v());
}

// can maybe speed things up by overwriting get_g to return gbar & making integrate empty.
void OUInput::integrate(double dt)
{
    dg = tau_inv*dt*(g_inf-g) + sd*dw(dt);
    
    if ((g+dg) < 0.0)
    {
        g -= dg;
    }
    else
    {
        g += dg;
    }
    
    gbar = g;
}

void OUInput::connect(compstate *s_in)
{
    compartment_state = s_in;
}

void OUInput::state2double(double* st) { }

void OUInput::double2state(double*) { }

// wiener increment
double OUInput::dw(double dt)
{
    r1 = ((double) rand())/((double) rmax);
    r2 = ((double) rand())/((double) rmax);
    
    out = sqrt(-2.0*log(r1)*dt)*cos(6.283185307179586*r2);
    
    if (r1 > 1.0/((double) rmax))
    {
        out = sqrt(-2.0*log(r1)*dt)*cos(6.283185307179586*r2);
    }
    else
    {
        //mexPrintf("r1 = %f\nr2 = %f\n", r1, r2);
        out = 0.0;
    }
    return out;
}

#endif
