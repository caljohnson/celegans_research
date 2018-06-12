%non_dimlzed_Mech_coupling_Discrete_neurons_Sim
%simulates mechanical coupling on discrete neuron-system pair
%with "nondimensionalized" parameters

%simulation runtime
TF = 1e2;

%mechanical params
Gamma = 0; %parameter sweep this - Coupling strength
tau_f = 5; %mu/k
c_MA = 10; %muscle activity mechanical feedback strength
%timescale of muscle activation
tau_m = 1;

% % Mechanically coupled system 
LHS_matrix = tau_f.*[(Gamma/tau_f)+1, (Gamma/tau_f)+1; 1, -1;];
RHS_matrix = [-1, -1; -1, +1;];
%R*[theta1'; theta2;] = L[theta1; theta2;] - c_MA*[M1+M2; M1-M2;]
kappa_dot = @(t,kappa, A) LHS_matrix\(RHS_matrix*[kappa(1); kappa(2);]...
    -c_MA.*[(A(2)-A(1))+(A(4)-A(3));(A(2)-A(1))-(A(4)-A(3));]); 
%GAMMA = 0:
% kappa_dot = @(t,kappa,A) (1/tau_f).*[-kappa(1) - c_MA*(A(2)-A(1));...
%     -kappa(2) - c_MA*(A(4)-A(3));];
   
%neural params
eps = 2;   %%-] These two determine thresholds together
I = 0.01;  %%-]

%results in the thresholds
K_V_ON = eps/2-I;   %K_D_ON is negative of this
K_V_OFF = -eps/2-I; %K_D_OFF is negative of this

% %IC - start antiphase
% SD(1) = 1;
% SV(1) = 0;
% SD(2) = 0;
% SV(2) = 1;
% K(1) = K_V_OFF;
% K(2) = -K_V_OFF;
% A(1) = 0;
% A(2) = 0;
% A(3) = 0;
% A(4) = 0;
% % % 
% %IC - start in-phase
% SD(1) = 1;
% SV(1) = 0;
% SD(2) = 1;
% SV(2) = 0;
% K(1) = K_V_OFF;
% K(2) = K_V_OFF;
% A(1) = 1; % A_1^D
% A(2) = 0; % A_1^V
% A(3) = 1;% A_2^D
% A(4) = 0;% A_2^V

% %IC - start w/ phase difference
% SD(1) = 1;
% SV(1) = 0;
% SD(2) = 1;
% SV(2) = 0;
% K(1) = K_V_OFF;
% K(2) = 0;
% A(1) = 1;
% A(2) = 0;
% A(3) = 0.5;
% A(4) = 0;

%IC - only one is oscillating
SD(1) = 1;
SV(1) = 0;
SD(2) = 0;
SV(2) = 0;
K(1) = K_V_OFF;
K(2) = 0;
A(1) = 1; % A_1^D
A(2) = 0; % A_1^V
A(3) = 0;% A_2^D
A(4) = 0;% A_2^V


%neural functions
state_d_1 = discrete_neural_state_init(SD(1), K(1), K_V_OFF, K_V_ON, 0);
state_v_1 = discrete_neural_state_init(SV(1), K(1), K_V_OFF, K_V_ON, 1);
state_d_2 = discrete_neural_state_init(SD(2), K(2), K_V_OFF, K_V_ON, 0);
state_v_2 = discrete_neural_state_init(SV(2), K(2), K_V_OFF, K_V_ON, 1);

%muscle eqns:
muscle_activity = @(t, K,A) (1/tau_m).*[-A(1) + (state_d_1(K(1)) - state_v_1(K(1))); ...
                -A(2) + (state_v_1(K(1)) - state_d_1(K(1)));...
                -A(3) + (state_d_2(K(2)) - state_v_2(K(2))); ...
                -A(4) + (state_v_2(K(2)) - state_d_2(K(2)));];

ode_rhss = @(t,X) [kappa_dot(t,X(1:2),X(3:6)); muscle_activity(t,X(1:2),X(3:6));];
init_cond = [K(1); K(2); A(1);A(2);A(3);A(4);];

options = odeset('RelTol',1e-8,'AbsTol',1e-10,  'MaxStep', 1e-2);
[t,y] = ode23(ode_rhss,[0,TF], init_cond, options);
            
figure(1); clf;
subplot(3,2,1); plot(t,y(:,1), '-'); ylabel('\kappa_1'); xlabel('t');%ylim([-1 1]);
subplot(3,2,2); plot(t,y(:,2), '-'); ylabel('\kappa_2'); xlabel('t');%ylim([-1 1]);
subplot(3,2,3); plot(t,y(:,3), '-'); ylabel('A_1^D'); xlabel('t');
subplot(3,2,5); plot(t,y(:,4), '-'); ylabel('A_1^V'); xlabel('t');
subplot(3,2,4); plot(t,y(:,5), '-'); ylabel('A_2^D'); xlabel('t');
subplot(3,2,6); plot(t,y(:,6), '-'); ylabel('A_2^V'); xlabel('t');
% subplot(3,2,3); plot(t,y(:,4)-y(:,3), '-'); ylabel('A_1^V - A_1^D'); xlabel('t');
% subplot(3,2,4); plot(t,y(:,6)-y(:,5), '-'); ylabel('A_2^V - A_2^D'); xlabel('t');
figure(2); clf;
plot(y(:,1), y(:,2), '-');

