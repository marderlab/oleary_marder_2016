//just run to get membrane trace etc, no regulation
#include <cmath>
#include <vector>
#include "mex.h"
#include "compartment.h"
#include "NaV.h"
#include "CaT.h"
#include "CaS.h"
#include "Ka.h"
#include "KCa.h"
#include "Kdr.h"
#include "Ih.h"
#include "leak.h"

using namespace std;

// simple max function
int max(int x, int y)
{
  return (x > y) ? x: y;
}

// usage [output] = f([dt tstop res],[gbars]
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *output;
    double *simparams, *gmax;
    double *inj_times, *inj_values;
    double *T_times, *T_values; //temperature ramp values
    int nits = 10, res = 1;
    double dt, tstop, v0 = -60.0;
    double gbar_leak, gbar_na, gbar_kd, gbar_ka, gbar_kca, gbar_cas, gbar_cat, gbar_h;
    double e_na, e_k, e_h, e_leak;
    int state_dim = 0, full_state_dim = 0;
    double *state, *full_state;
    int trows, tcols, vrows, vcols, nsteps;
    int inj_step = 0;
    int T_trows, T_tcols, T_rows, T_cols, T_nsteps; //temperature ramps
    int T_step = 0;
    double delta_T = 0.0; //Temp increment
    double T_old;
    
    //input + params
    if(!(nrhs==6)) {
        mexErrMsgTxt("Wrong number of inputs");
    }
    // 1st input - simulation parameters
    simparams = mxGetPr(prhs[0]);
    dt = simparams[0];
    tstop = simparams[1];
    res = simparams[2];
    v0 = simparams[3];

    // 2nd input - conductance parameters
    gmax = mxGetPr(prhs[1]);
    gbar_na = gmax[0];
    gbar_cat = gmax[1];
    gbar_cas = gmax[2];
    gbar_ka = gmax[3];
    gbar_kca = gmax[4];
    gbar_kd = gmax[5];
    gbar_h = gmax[6];
    gbar_leak = gmax[7];
    e_leak = gmax[8];
    e_na = gmax[9];
    e_k = gmax[10];
    e_h = gmax[11];
    
    // 3rd & 4th inputs - current injection
    // array of times and values
    inj_times = mxGetPr(prhs[2]);
    inj_values = mxGetPr(prhs[3]);
    
    trows = mxGetM(prhs[2]);
    tcols = mxGetN(prhs[2]);
    vrows = mxGetM(prhs[3]);
    vcols = mxGetN(prhs[3]);
    
    if( !((trows==1 || tcols==1) && (vrows==1 || vcols==1)) ||
        !(max(trows,tcols)==max(vrows,vcols)) ) {
        mexErrMsgTxt("Injection arrays must be 1xN or Nx1 and equal in length");
    }
    
    nsteps = max(trows,tcols);
    //mexPrintf("dim = %d\n", nsteps);
    
    // 3rd & 4th inputs - temperature ramp
    // array of times and values - linearly interp between them
    T_times = mxGetPr(prhs[4]);
    T_values = mxGetPr(prhs[5]);
    
    T_trows = mxGetM(prhs[4]);
    T_tcols = mxGetN(prhs[4]);
    T_rows = mxGetM(prhs[5]);
    T_cols = mxGetN(prhs[5]);
    
    if( !((T_trows==1 || T_tcols==1) && (T_rows==1 || T_cols==1)) ||
        !(max(T_rows,T_cols)==max(T_rows,T_cols)) ) {
        mexErrMsgTxt("Temperature arrays must be 1xN or Nx1 and equal in length");
    }
    
    T_nsteps = max(T_trows,T_tcols);
    //mexPrintf("dim = %d\n", nsteps);
    
    //make conductances
    NaV gna(gbar_na,e_na        ,gmax[12],gmax[13],gmax[14],10.0);
    CaT gcat(gbar_cat,120.0     ,gmax[15],gmax[16],gmax[17],10.0);
    CaS gcas(gbar_cas,120.0     ,gmax[18],gmax[19],gmax[20],10.0);
    Ka gka(gbar_ka,e_k          ,gmax[21],gmax[22],gmax[23],10.0);
    KCa gkca(gbar_kca,e_k       ,gmax[24],gmax[25],gmax[26],10.0);
    Kdr gkdr(gbar_kd,e_k        ,gmax[27],gmax[28],gmax[29],10.0);
    Ih gh(gbar_h,e_h            ,gmax[30],gmax[31],gmax[32],10.0);
    leak gleak(gbar_leak,e_leak ,gmax[33],10.0);
    
    //make compartment
    compartment cell(v0, 0.05, gmax[34], gmax[35], 10.0);
    cell.add_conductance(&gna);
    cell.add_conductance(&gcat);
    cell.add_conductance(&gcas);
    cell.add_conductance(&gka);
    cell.add_conductance(&gkca);
    cell.add_conductance(&gkdr);
    cell.add_conductance(&gh);
    cell.add_conductance(&gleak);
    nits = (int) floor(tstop/dt);
    
    state_dim = 3;
    state = new double[state_dim];
    //full_state_dim = cell.get_state_dim();
    //full_state = new double[full_state_dim];
    
    plhs[0] = mxCreateDoubleMatrix(state_dim, ((int)nits)/((int)res), mxREAL);
    //plhs[0] = mxCreateDoubleMatrix(state_dim, nits, mxREAL);
    output = mxGetPr(plhs[0]);
    
    if (inj_times[0] < dt)
    {
        cell.input(inj_values[0]);
        inj_step = 1;
    }
    
    cell.set_T(T_values[0]);

    if (T_times[0] < dt)
    {
        T_step = 1;
    }
    
    //integration loop
    for(int I = 0, i=0; i<nits; i++)
    //for(int i=0; i<nits; i++)
    {
        (void) cell.integrate(dt);
        if (i%res == 0)
        {
            //mexPrintf("temp = %f\n", cell.get_T());

            cell.get_vct(state);
            //cell.get_state(full_state);
            for(int j=0; j<state_dim; j++)
            {
                output[state_dim*I+j] = state[j];
            }
            //output[(state_dim+1)*I+state_dim] = full_state[full_state_dim - 1];
            
            if ((inj_step < nsteps) && (((double) i)*dt > inj_times[inj_step]))
            {
                cell.input(inj_values[inj_step]);
                //mexPrintf("dim = %f\n", inj_values[inj_step]);
                inj_step++;
            }
            
            if ((T_step < T_nsteps) && (((double) i)*dt > T_times[T_step]))
            {
                cell.set_T(T_values[T_step]);
                //mexPrintf("dim = %f\n", inj_values[inj_step]);
                T_step++;
            }
            else if (T_step < T_nsteps)
            {
                if (T_step > 0)
                {
                    T_old = cell.get_T();
                    delta_T = (T_values[T_step] - T_values[T_step-1])/(T_times[T_step] - T_times[T_step-1]);
                    cell.set_T(T_old + delta_T*dt*((double) res));
                }
                else
                {
                    cell.set_T(T_values[0]);
                }
            }
               
            I++;
        }
    }

}
