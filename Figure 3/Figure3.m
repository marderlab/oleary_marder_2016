mex int_control_compartment_converge.cpp
mex run_comp_inj_q10.cpp

%note: different random examples are generated on each run

close all;
load('params.mat');

% simulation params
rang = 1.0e-1;
dt = 0.1;
taug = 1.e1;
taum = 10*taug;
tstop = 10000*taug;
tstop_trace = 1e4;
ms_per_samp = tstop/1000;

t_inj = [0 100]; % time of current injection (ms)
i_inj = [0 0]; %current (nA)
t_Temp = [tstop_trace/3 2*tstop_trace/3]; % time of temp ramp (ms)
Temps = [5 25]; %temperature (deg C)

e_leak = gparams_hi(9);
e_na = gparams_hi(10);
e_k = gparams_hi(11);
e_h = gparams_hi(12);

DC = zeros(5,5);
Ca = zeros(5,5);

allgs = [];

pal = 'Mrainbow'; %'Paired12';%'Accent8';%
colrs = othercolor(pal,8);
set(0,'DefaultAxesColorOrder',colrs);
tracecol = othercolor('Mdarkrainbow',5);

colrs = get(0,'defaultaxescolororder');

for j=1:5
    
    gbar_leak = 0.01 + rand*0.09;
    gs_init = [rand(1,7)*rang gbar_leak];
    v0 = e_leak + randn*2;
    [X, stopit] = int_control_compartment_converge([dt tstop floor(ms_per_samp/dt) v0],[gs_init e_leak e_na e_k e_h],simparams);
    X = X(:,1:stopit);
    tX = linspace(0,tstop,length(X)+1);
    gs = X(3:10,end)';
    
    allgs = [allgs' gs']';
    g_params = [gs(1:8) gparams_hi(9:36)];
    
    for p=1:5
        u = run_comp_inj_q10([dt tstop_trace 1 v0],g_params,t_inj,i_inj,t_Temp,[5*p 5*p]);
        b = burstParams2(u(1,end/2:end),dt,0);
        DC(j,p) = b(3);
        Ca(j,p) = mean(u(2,end/2:end));
    end
    
    % run the simulation
    v = run_comp_inj_q10([dt tstop_trace 1 v0],g_params,t_inj,i_inj,t_Temp,Temps);
    t = linspace(0,tstop_trace,length(v));
    
    % plot
    figure(1);
    subplot(6,1,j);
    plot(t,v(1,:),'linewidth',1.5,'color',tracecol(j,:));
    set(gca,'xlim',[tstop_trace/5 4*tstop_trace/5],'ylim',[-80 40],'ytick',[-60 0]);
    ylabel('Vm');
    box off;
end


figure(1);
subplot(6,1,6);
plot(t,v(3,:),'r','linewidth',2);
set(gca,'xlim',[tstop_trace/5 4*tstop_trace/5]);
ylabel('temperature (^oC)');
xlabel('time (ms)');
box off;

figure(2);
[H,AX,BigAx,P,PAx] = plotmatrix_spaced(allgs(:,1:7));
for i=1:7
    cla(PAx(i));
    set(PAx(i),'box','off','xcolor',[1 1 1],'ycolor',[1 1 1]);
    for j=1:7
        set(AX(i,j),'ticklength',[0.025 0.025]);
        if i<=j
            cla(AX(i,j));
            set(AX(i,j),'xcolor',[1 1 1],'ycolor',[1 1 1]);
            set(AX(i,j),'box','off');
        end
    end
end

figure(3);
for i=1:5
    hold on;
    plot(DC(:,i),'o-','color',tracecol(i,:));
end
set(gca,'xticklabel',(5:5:25),'ylim',[0 0.5]);
xlabel('Temperature');
ylabel('Duty cycle');
box off;

figure(4);
subplot(2,1,1);
vb = v(1,1:t_Temp(1)/dt);
[t0, ~] = burstTimes(vb,0);
t = linspace(1,t_Temp(1),length(vb));
vb = vb(t0(5):end);
t = t(1:end-t0(5)+1);
plot(t,vb,'linewidth',2);
set(gca,'xlim',[0 (t0(8)-t0(5))*dt],'ylim',[-80 40],'ytick',[-60 0]);
ylabel('Vm');
box off;

subplot(2,1,2);
va = v(1,t_Temp(2)/dt:end);
[t0, ~] = burstTimes(va,0);
t = linspace(1,t_Temp(1),length(va));
va = va(t0(5):end);
t = t(1:end-t0(5)+1);
plot(t,va,'linewidth',2);
set(gca,'xlim',[0 (t0(8)-t0(5))*dt],'ylim',[-80 40],'ytick',[-60 0]);
ylabel('Vm');
xlabel('time(ms)');
box off;
