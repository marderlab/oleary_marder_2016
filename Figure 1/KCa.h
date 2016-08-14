// LEAK CONDUCTANCE
#ifndef KCA
#define KCA
#include "conductance.h"

//inherit conductance class spec
class KCa: public conductance {
protected:
    double m;
public:
    //constructors - must set gbar!    
    // gbar only
    KCa(double g_in)//: conductance(g_in)
    {
        gbar = g_in;
        state_dim = 1;
        // default if unspecified
        e_rev = -80.0;
    }
    
    //specify both gbar and erev
    KCa(double g_in, double e_in)//: conductance(g_in)
    {
        gbar = g_in;
        state_dim = 1;
        e_rev = e_in;
    }
    
    KCa(double g_in, double e_in, double gq10, double actq10, double inactq10, double reftemp)//: conductance(g_in)
    {
        gbar = g_in;
        state_dim = 2;
        e_rev = e_in;
        g_q10 = gq10;
        act_q10 = actq10;
        inact_q10 = inactq10;
        ref_temp = reftemp;
    }
    
    //return current
    double i(void);
    int get_state_dim(void);
    void integrate(double dt);
    //let conductance access cell/compartment state
    void connect(compstate *s_in);
    double m_inf(void);
    double tau_m(void);
    void state2double(double*);
    void double2state(double*);
};

int KCa::get_state_dim(void)
{
    return 1;
}

//membrane current
double KCa::i(void)
{
    return g*(e_rev - v_m());
}

// can maybe speed things up by overwriting get_g to return gbar & making integrate empty.
void KCa::integrate(double dt)
{
    m = m_inf() + (m - m_inf())*exp(-dt/(tau_m()*pow(act_q10,((ref_temp-get_T())/10.0))));
    g = gbar*m*m*m*m*pow(g_q10,((get_T()- ref_temp)/10.0));
}

void KCa::connect(compstate *s_in)
{
    compartment_state = s_in;
    m = m_inf();
    //mexPrintf("state ptr val: %f\n",compartment_state->get_v());
    //mexPrintf("state ptr: %d\n",compartment_state);
}

double KCa::m_inf(void) 
{
    double c = compartment_state->get_Ca();
    return (c/(c+3.0))/(1.0+exp((v_m()+28.3)/-12.6));
}

double KCa::tau_m(void) {return 90.3 - 75.1/(1.0+exp((v_m()+46.0)/-22.7));}

void KCa::state2double(double* st)
{
    st[0] = m;
}

void KCa::double2state(double* st)
{
    m = st[0];
}

#endif
