%Mech_coupling_Discrete_neurons_Sim
%simulates mechanical coupling on discrete neuron-system pair

%mechanical params
gamma = 1e-2; %parameter sweep this, calculate stable phase differences
mu = 5;  %]
k = 1;   %] - set so that tau = mu/k = 5
a = 0.1;     %]
delX = 0.5;  %]- set so that c= 1/(2*a*delX) = 10

% % Mechanically coupled system - no external fluid viscosity:
% kappa_dot = @(t,kappa, A) [(-k/mu)*(kappa(1) ...
%     - (A(1) -A(2))/(2*a*delX) ...
%     - delX/(8*a^2)*(kappa(1)+kappa(2))*(A(1) + A(2)));
%     (-k/mu)*(kappa(2) ...
%     - (A(3) -A(4))/(2*a*delX) ...
%     - delX/(8*a^2)*(kappa(1)+kappa(2))*(A(3) + A(4)));];
%    
% % ignoring "high-order" terms, no external fluid viscosity
% kappa_dot1 = @(t,kappa, A) [(-k/mu)*(kappa(1) ...
%     - (A(1) -A(2))/(2*a*delX));
%     (-k/mu)*(kappa(2) ...
%     - (A(3) -A(4))/(2*a*delX));]; 

% % with fluid viscosity, no higher order terms
% d = gamma*delX/(4*a^2*mu);
% RHS_matrix = [d 1; 1 d;];
% kappa_dot = @(t,kappa, A) RHS_matrix\[(-k/mu)*(kappa(2) ...
%     - (A(3) -A(4))/(2*a*delX));
%     (-k/mu)*(kappa(1) ...
%     - (A(1) -A(2))/(2*a*delX));];

%with fluid viscosity, with higher order terms
d = gamma*delX/(4*a^2*mu);
RHS_matrix = [d 1; 1 d;];
kappa_dot = @(t,kappa, A) RHS_matrix\[(-k/mu)*(kappa(2) ...
    - (A(3) -A(4))/(2*a*delX) ...
    - delX/(8*a^2)*(kappa(1)+kappa(2))*(A(3) + A(4)));
    (-k/mu)*(kappa(1) ...
    - (A(1) -A(2))/(2*a*delX) ...
    - delX/(8*a^2)*(kappa(1)+kappa(2))*(A(1) + A(2)));];
   
%neural params
eps = 2;   %%-] These two determine thresholds together
I = 0.01;  %%-]

%results in the thresholds
K_V_ON = eps/2-I;   %K_D_ON is negative of this
K_V_OFF = -eps/2-I; %K_D_OFF is negative of this

%simulation stuff
N = 1e6;
dt = 5e-4;
K = zeros(N,2);
deltaV = zeros(N,2);
A1 = zeros(2,1);
A2 = zeros(2,1);

% %IC - start antiphase
% SD(1) = 1;
% SV(1) = 0;
% SD(2) = 0;
% SV(2) = 1;
% K(1,1) = K_V_OFF;
% K(1,2) = -K_V_OFF;

% %IC - start in-phase
% SD(1) = 1;
% SV(1) = 0;
% SD(2) = 1;
% SV(2) = 0;
% K(1,1) = K_V_OFF;
% K(1,2) = K_V_OFF;

%IC - start w/ phase difference
SD(1) = 1;
SV(1) = 0;
SD(2) = 1;
SV(2) = 0;
K(1,1) = K_V_OFF;
K(1,2) = -K_V_ON;

%neural functions
state_v_1 = discrete_neural_state_init(SV(1), K(1,1), K_V_OFF, K_V_ON, 1);
state_d_1 = discrete_neural_state_init(SD(1), K(1,1), K_V_OFF, K_V_ON, 0);
state_v_2 = discrete_neural_state_init(SV(2), K(1,2), K_V_OFF, K_V_ON, 1);
state_d_2 = discrete_neural_state_init(SD(2), K(1,2), K_V_OFF, K_V_ON, 0);

%muscle eqns:
muscle_activity = @(t, K,A) [-A(1) + (state_d_1(K(1)) - state_v_1(K(1))); ...
                -A(2) + (state_v_1(K(1)) - state_d_1(K(1)));...
                -A(3) + (state_d_2(K(2)) - state_v_2(K(2))); ...
                -A(4) + (state_v_2(K(2)) - state_d_2(K(2)));];

ode_rhss = @(t,X) [kappa_dot(t,X(1:2),X(3:6)); muscle_activity(t,X(1:2),X(3:6));];
init_cond = [K(1,1); K(1,2); 0;0;0;0;];

[t,y] = ode113(ode_rhss,[0,1e2], init_cond);
            
% for i = 2:N
%     %update ?V
%     deltaV(i,1) = SD(1) - SV(1);
%     deltaV(i,2) = SD(2) - SV(2);
%     
%     %forward-Euler step
% %     A1 = A1 + dt*muscle_activity(SD(1), SV(1), A1);
% %     A2 = A2 + dt*musc`le_activity(SD(2), SV(2), A2);
% %     K(i,:) = K(i-1,:) + dt*kappa_dot(K(i-1,:), A1,A2)';
%   
%     %update states for unit 1
%     SV(1) = state_v_1(K(i,1));
%     SD(1) = state_d_1(K(i,1));
%     %update states for unit 2
%     SV(2) = state_v_2(K(i,2));
%     SD(2) = state_d_2(K(i,2));
% end

figure(1); clf;
% t = 0:dt:dt*(N-1);
subplot(2,2,1); plot(t,y(:,1), '-'); ylabel('\kappa_1'); xlabel('t');
subplot(2,2,2); plot(t,y(:,2), '-'); ylabel('\kappa_2'); xlabel('t');
subplot(2,2,3); plot(t,y(:,3)-y(:,4), '-'); ylabel('M_1'); xlabel('t');
subplot(2,2,4); plot(t,y(:,5)-y(:,6), '-'); ylabel('M_2'); xlabel('t');
% figure(2); clf;
% plot(K(:,1), K(:,2), '-');
% plot(K(end-2000:end,1), K(end-2000:end,2), '-');
% figure(3); clf;
% plot(K, Kp, '-'); hold on; ylabel('d \kappa / dt'); xlabel('\kappa');
% line([-c-1,c+1],[0,0], 'Color',[0 0 0])
% line([K_V_ON, K_V_ON], [-c-1,c+1], 'Color',[1 0 0])
% line([-K_V_OFF, -K_V_OFF], [-c-1,c+1], 'Color',[0 1 0])
% line([K_V_OFF, K_V_OFF], [-c-1,c+1], 'Color',[1 0 0])
% line([-K_V_ON, -K_V_ON], [-c-1,c+1], 'Color',[0 1 0])
% line([-c, -c], [0,0], 'Color',[0 0 0], 'Marker', 'o')
% line([c, c], [0,0], 'Color',[0 0 0], 'Marker', 'o')
% line([0,0],[-c-1,c+1], 'Color',[0 0 0])
% xlim([min(K), max(K)])
% ylim([min(Kp), max(Kp)])
