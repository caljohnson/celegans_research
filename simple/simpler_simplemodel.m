%simpler_simplemodel
%2-spring unit oscillator

%clear everything
clear; clc;
tic

%timescales
t_f = 10; %mechanical force timescale
t_m = 100; %muscular activity timescale

%assume neural dynamics occur on a fast timescale compared to
%muscle+mechanics s.t. always at steady state - 0 or 1
SD_init = 0;
SV_init = 0;
I_AVB = 0.48; %driving AVB current
k_SR = 1; %stretch receptor weight
eps_h = 0.2; %hysteresis window

I = @(m) I_AVB - k_SR*m;
S = @(I,s) 1*(I>=0.5+eps_h*(0.5-s)) + 0*(I<=0.5+eps_h*(0.5-s));
S_D = @(m) S(I(m), SD_init);
S_V = @(m) S(I(-m), SV_init);

%driving torque- depends pw-linearly on muscle activities AD, AV
pw_lin = @(A) 0*(A<=0) + A*(0<A & A<=1) + 1*(A>=1);
m0 = @(AD, AV) pw_lin(AD) - pw_lin(AV);


%muscle and mechanical ODEs
m_rhs = @(m, AD, AV) -(1/t_f)*(m-m0(AD,AV));
AD_rhs = @(m, AD, AV) (1/t_m)*(S_D(m) - S_V(m) - AD);
AV_rhs = @(m, AD, AV) (1/t_m)*(S_V(m) - S_D(m) - AV);

system = @(t,x) [m_rhs(x(1), x(2), x(3)); AD_rhs(x(1), x(2), x(3)); AV_rhs(x(1), x(2), x(3));];

%integrate
dt = (1/t_m)*2;
tspan = 1:dt:1000;
m_init = 0.5;
AD_init= 0.5;
AV_init = 0.25;

y_temp = [m_init; AD_init; AV_init;];

%holders for ODE solve data
y=[y_temp];
sd = [SD_init];
sv = [SV_init];


for i=1:size(tspan,2)-1

    %solve ODE w/ fixed neural state - use forward Euler
    y_temp = y_temp + system(i*dt, y_temp)*dt;
    
    %save ODE data
    y = [y y_temp;];
    
    %update neural states
    SD_init = S_D(y_temp(1,end)); 
    SV_init = S_V(y_temp(1,end));
    sd = [sd SD_init;];
    sv = [sv SV_init;];
    
    %remake ODE equations
    S_D = @(m) S(I(m), SD_init);
    S_V = @(m) S(I(-m), SV_init);
    AD_rhs = @(m, AD, AV) (1/t_m)*(S_D(m) - S_V(m) - AD);
    AV_rhs = @(m, AD, AV) (1/t_m)*(S_V(m) - S_D(m) - AV);
    system = @(t,x) [m_rhs(x(1), x(2), x(3)); AD_rhs(x(1), x(2), x(3)); AV_rhs(x(1), x(2), x(3));];
end

figure(1);
for i=1:3
    plot(tspan,y(i,:)); hold on
end
xlabel('t'); hold off;
legend('m', 'AD', 'AV');

figure(2);
plot(tspan, sd, 'o'); hold on
plot(tspan,sv, 'o'); hold off;
xlabel('t')
legend('SD', 'SV');

toc



