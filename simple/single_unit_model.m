%single_unit_model
%2-spring unit oscillator


I_AVB = -100;

%timescales
t_f = 1; %mechanical force timescale
t_m = 10; %muscular activity timescale
t_n = 100; %neural activity timescale

a = 0.33; %steepness of fitzhugh nagumo curve
k_SR = 10; %stretch receptor weight
C = 10; %muscle activity torque/strength

%all ODEs
v1_rhs = @(v1, K) (1/t_n)*(v1 - a*v1^3+I_AVB - k_SR*K);
v2_rhs = @(v2, K) (1/t_n)*(v2 - a*v2^3+I_AVB + k_SR*K);
M_rhs = @(v1,v2,M) (1/t_m)*(2*v1 - 2*v2 - M);
K_rhs = @(M, K) -(1/t_f)*(K - C*M);

system_rhs = @(t,x) [v1_rhs(x(1), x(4)); v2_rhs(x(2), x(4)); M_rhs(x(1), x(2), x(3)); K_rhs(x(3), x(4));];

%integrate
v1_init = 0;
v2_init = 0;
M_init = 0;
K_init = 60;

y_init = [v1_init; v2_init; M_init; K_init;];

tspan=[0 200];
[t,y] = ode45(system_rhs, tspan, y_init);

figure(1);
plot(t, y(:,1));
ylabel('V1'); xlabel('t');

figure(2);
plot(t, y(:,2));
ylabel('V2'); xlabel('t');

figure(3);
plot(t, y(:,3));
ylabel('M'); xlabel('t');

figure(4);
plot(t, y(:,4));
ylabel('K'); xlabel('t');


