mex run_comp_inj_q10.cpp

% three examples with random Q10s
%load example_q10s_1.mat;
%load example_q10s_2.mat;
load example_q10s_3.mat;

% simulation params
dt = 0.1; % timestep
tstop = 2e4; % duration
dt_per_samp = 1; % sample resolution (must be integer)

% reversal potentials
e_na = 30;
e_leak = -50;
e_k = -80;
e_h = -20;

% initial memb potential
v0 = -60;

%max conductances
gs = [95.0174    0.7370    0.9627   10.1772   11.2365   18.0006    0.0110    0.0187];

t_inj = [0 100]; % time of current injection (ms)
i_inj = [0 0]; %current (nA)
t_Temp = [tstop/3 2*tstop/3]; % time of current injection (ms)
Temps = [5 25]; %temperature (deg C)

% run the simulation
v = run_comp_inj_q10([dt tstop dt_per_samp v0],g_params,t_inj,i_inj,t_Temp,Temps);
t = linspace(0,tstop,length(v));

% plot
figure(1);
subplot(2,1,1);
plot(t,v(1,:));
set(gca,'xlim',[tstop/4 3*tstop/4],'xcolor',[1 1 1]);
ylabel('Vm');
box off;
subplot(2,1,2);
plot(t,v(3,:),'r');
set(gca,'xlim',[tstop/4 3*tstop/4],'xcolor',[1 1 1]);
ylabel('temperature (^oC)');
box off;

%{
figure(2)

t_inj = [0 100]; % time of current injection (ms)
i_inj = [0 0]; %current (nA)
t_Temp = [tstop/3 2*tstop/3]; % time of current injection (ms)
Temps = [10 10]; %temperature (deg C)

% run the simulation
v = run_comp_inj_q10([dt tstop dt_per_samp v0],g_params,t_inj,i_inj,t_Temp,Temps);
t = linspace(0,tstop,length(v));

subplot(2,1,1);
plot(t,v(1,:));
set(gca,'xlim',[tstop/4 3*tstop/4],'xcolor',[1 1 1]);
ylabel('Vm');
box off;
title 'reference temperature'
%}
