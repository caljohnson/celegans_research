%Simple Springs model of c. elegans

%------------------- ODEs and relevant parameters ------------------

%spring parameters
xi = 10; %viscosity
k = 1; %spring constant
L0 = 1; %length of incompressible middle rod
L_0D = 1; %rest length of dorsal spring
L_0V = 1; %rest length of ventral spring

%resting torque equation - for spring ODE
pw_lin = @(A) 0*(A<=0) + A*(0<A & A<=1) + 1*(A>=1);
m0 = @(AD, AV) pw_lin(AD) - pw_lin(AV);

%ODE 1: -xi dm/dt = k(m-m0)
m_rhs = @(m, AD, AV) (-1/xi)*k*(m-m0(AD,AV));

%neural parameters
C = 1; %membrane capacitance, pF
G_mem = 500; %membrane conductance, pS
G_act = 20; %maximum conductance, pS
G_AVB = 150; %AVB gap-junction conductance, pS
G_GABA = 10; %max synaptic conductance for neural inhibiton, pS
G_SR = 20; %stretch receptor conductance, decreases linearly
k_GABA = 100; %inhibition function rate, 1/mV
k_act = 500; %activation function rate, 1/mV
V_rest = -72; %membrane reversal potential, mV
V_act = -60; %membrane activation potential, mV
V_AVB = -20.05;%-87.5; %AVB membrane potential, mV
V_0GABA = V_rest; %inhibition activation potential
I_bias = 8; %AVB ventral bias, pA

%neural functions
I_leak = @(V) G_mem*(V - V_rest);
I_act = @(V) G_act./(1 + exp(-k_act*(V - V_act)));
I_AVB = @(V) G_AVB*(V_AVB - V);
I_inh = @(V) G_GABA/(1 + exp(-k_GABA*(V - V_0GABA)));
I_sr = @(m) G_SR*(2*L0 + m - L_0D)/(2*L_0D);
% I_sr_V = @(m,A) G_SR*(2*L0 - m - L_0V)/(2*L_0V);

%ODE 2: dVD/dt = ...
VD_rhs = @(m, VD, VV, AD) (1/C)*(-I_leak(VD)+I_act(VD)+I_AVB(VD) + I_sr(m));

%ODE 3: dVV/dt = ...
VV_rhs = @(m, VD, VV, AV) (1/C)*(-I_leak(VV) + I_act(VV) + I_AVB(VV) +I_sr(m));%...
                            %+ I_bias + I_inh(VD) +  I_sr(-m));
                        
%muscle input parameters
w_ACH = 3; %synaptic weights for excitatory ACH synapses, unitless
w_GABA = -0.5; %synaptic weights for inhibitory GABA synapses, unitless
k_NMJ = 50; %muscle input activation rate, 1/mV
V_0NMJ = -22; %neuromuscular resting potential, mV
%muscle input function
I_NMJ = @(V1,V2) w_ACH/(1+exp(-k_NMJ*(V1 - V_0NMJ))) ...
               + w_GABA/(1+exp(-k_NMJ*(V2 - V_0NMJ)));

%muscle activity parameter
tau_m = 100; %time constant of muscle input integration, mS

%ODE 4: dAD/dt = ...
AD_rhs = @(AD, VD, VV) (1/tau_m)*(I_NMJ(VD,VV) - AD);

%ODE 5L dAV/dt = ...
AV_rhs = @(AV, VD, VV) (1/tau_m)*(I_NMJ(VV,VD) - AV);

system_rhs = @(t,x) [m_rhs(x(1), x(4), x(5)); VD_rhs(x(1), x(2), x(3), x(4)); ...
              VV_rhs(x(1), x(2), x(3), x(5)); AD_rhs(x(4), x(2), x(3)); AV_rhs(x(5), x(2), x(3));];

%------------------------- Simulation ---------------------------
m_init = 0;
VD_init = V_rest;
VV_init = V_act;
AD_init = 0;
AV_init = 1;
y_init = [m_init; VD_init; VV_init; AD_init; AV_init;];

tspan=[0 10];
[t,y] = ode45(system_rhs, tspan, y_init);

figure(1);
for i =[1, 4,5]
plot(y(:,i)); hold on
end; hold off
legend('m', 'AD', 'AV');
figure(2);
for i = [2,3]
    plot(y(:,i)); hold on
end; hold off
legend('VD', 'VV');
