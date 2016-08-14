// LEAK CONDUCTANCE
#ifndef LEAK
#define LEAK
#include "conductance.h"

//inherit conductance class spec
class leak: public conductance {
public:
    //constructors - must set gbar!    
    // gbar only
    leak(double g_in)//: conductance(g_in)
    {
        gbar = g_in;
        state_dim = 0;
        // default if unspecified
        e_rev = -50.0;
    }
    
    //specify both gbar and erev
    leak(double g_in, double e_in)//: conductance(g_in)
    {
        gbar = g_in;
        state_dim = 0;
        e_rev = e_in;
    }
    
    leak(double g_in, double e_in, double gq10, double reftemp)//: conductance(g_in)
    {
        gbar = g_in;
        state_dim = 2;
        e_rev = e_in;
        g_q10 = gq10;
        ref_temp = reftemp;
    }
    
    double i(void);
    int get_state_dim(void);
    void integrate(double dt); //integrate = do nothing
    void connect(compstate *s_in);
    void state2double(double*);
    void double2state(double*);    
};

int leak::get_state_dim(void)
{
    return 0;
}

double leak::i(void)
{
    return g*(e_rev - compartment_state->get_v());
}

// can maybe speed things up by overwriting get_g to return gbar & making integrate empty.
void leak::integrate(double dt)
{
    g = gbar*pow(g_q10,((get_T()- ref_temp)/10.0));
}

void leak::connect(compstate *s_in)
{
    compartment_state = s_in;
}

void leak::state2double(double* st)
{
}

void leak::double2state(double*) { }

#endif
