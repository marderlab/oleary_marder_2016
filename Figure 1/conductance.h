//Abstract class for defining conductances
#ifndef CONDUCTANCE
#define CONDUCTANCE
#include "compstate.h"
#include <cmath>

class conductance {
protected:
    double gbar;
    double g;
    double e_rev;
    int state_dim;
    double g_q10;
    double act_q10;
    double inact_q10;
    double ref_temp;
    compstate *compartment_state;
public:
    conductance()
    {
        //null pointer for safety
        compartment_state = 0;
        state_dim = 0;
        ref_temp = 10.0;
        g_q10 = 1.0;
        act_q10 = 1.0;
        inact_q10 = 1.0;
    }
    
    ~conductance() {}
    
    double v_m(void); //return compartment voltage
    double get_T(void); //return compartment temperature (deg C)
    double get_gbar(void);
    double* get_pgbar(void); //overwritten for control of e_rev
    double get_g(void);
    double get_ge(void); // for Dayan-Abbot integration
    
    // all of these must be defined and will be unique
    // for each conductance in general
    virtual int get_state_dim(void) =0;
    virtual void integrate(double) =0;
    virtual double i(void) =0; //consider removing this
    virtual void connect(compstate*) =0;
    virtual void state2double(double*) =0;
    virtual void double2state(double*) =0;
};

double conductance::v_m(void)
{
    return compartment_state->get_v();
}

double conductance::get_T(void)
{
    return compartment_state->get_T();
}

double conductance::get_gbar(void)
{
    return gbar;
}

double* conductance::get_pgbar(void)
{
    return &gbar;
}

double conductance::get_g(void)
{
    return g;
}

double conductance::get_ge(void)
{
    return g*e_rev;
}

#endif
