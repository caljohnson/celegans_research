%New_Mech_coupling_Discrete_neurons_Sim
%simulates mechanical coupling on discrete neuron-system pair

%simulation runtime
TF = 1e3;

%mechanical params
gamma = 1e0; %parameter sweep this, calculate stable phase differences
mu = 5;  %]
k = 1;   %] - set so that tau = mu/k = 5
a = 0.1;     %]
delX = 0.5;  %]- set so that c= 1/(2*a*delX) = 10
%timescale of muscle activation
tau_m = 1e1;

% % Mechanically coupled system 
thing1 = gamma*delX+4*a^2*mu;
thing2 = 4*a^2*mu;
RHS_matrix = (2/delX).*[thing1 thing1; thing2 -thing2;];
LHS_matrix = (-(2/delX)*4*a^2*k).*[1 1; 1 -1;];
%R*[theta1'; theta2;] = L[theta1; theta2;] - ak[M1+M2; M1-M2;]
kappa_dot = @(t,kappa, A) RHS_matrix\(LHS_matrix*[kappa(1); kappa(2);]...
    -a*k.*[(A(2)-A(1))+(A(4)-A(3));(A(2)-A(1))-(A(4)-A(3));]); 
   
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
% A(1,1) = 0;
% A(1,2) = 0;
% A(2,1) = 0;
% A(2,2) = 0;

%IC - start in-phase
SD(1) = 1;
SV(1) = 0;
SD(2) = 1;
SV(2) = 0;
K(1) = K_V_OFF;
K(2) = K_V_OFF;
A(1,1) = 1;
A(1,2) = 0;
A(2,1) = 1;
A(2,2) = 0;

% %IC - start w/ phase difference
% SD(1) = 1;
% SV(1) = 0;
% SD(2) = 1;
% SV(2) = 0;
% K(1) = K_V_OFF;
% K(2) = 0;
% A(1,1) = 1;
% A(1,2) = 0;
% A(2,1) = 0.5;
% A(2,2) = 0;

%neural functions
state_v_1 = discrete_neural_state_init(SV(1), K(1), K_V_OFF, K_V_ON, 1);
state_d_1 = discrete_neural_state_init(SD(1), K(1), K_V_OFF, K_V_ON, 0);
state_v_2 = discrete_neural_state_init(SV(2), K(2), K_V_OFF, K_V_ON, 1);
state_d_2 = discrete_neural_state_init(SD(2), K(2), K_V_OFF, K_V_ON, 0);

%muscle eqns:
muscle_activity = @(t, K,A) (1/tau_m).*[-A(1) + (state_d_1(K(1)) - state_v_1(K(1))); ...
                -A(2) + (state_v_1(K(1)) - state_d_1(K(1)));...
                -A(3) + (state_d_2(K(2)) - state_v_2(K(2))); ...
                -A(4) + (state_v_2(K(2)) - state_d_2(K(2)));];

ode_rhss = @(t,X) [kappa_dot(t,X(1:2),X(3:6)); muscle_activity(t,X(1:2),X(3:6));];
init_cond = [K(1); K(2); A(1,1);A(1,2);A(2,1);A(2,2);];

[t,y] = ode113(ode_rhss,[0,TF], init_cond);
            
figure(1); clf;
subplot(3,2,1); plot(t,y(:,1), '-'); ylabel('\kappa_1'); xlabel('t');ylim([-1 1]);
subplot(3,2,2); plot(t,y(:,2), '-'); ylabel('\kappa_2'); xlabel('t');ylim([-1 1]);
subplot(3,2,3); plot(t,y(:,3), '-'); ylabel('A_1^D'); xlabel('t');
subplot(3,2,5); plot(t,y(:,4), '-'); ylabel('A_1^V'); xlabel('t');
subplot(3,2,4); plot(t,y(:,5), '-'); ylabel('A_2^D'); xlabel('t');
subplot(3,2,6); plot(t,y(:,6), '-'); ylabel('A_2^V'); xlabel('t');
% subplot(3,2,3); plot(t,y(:,4)-y(:,3), '-'); ylabel('A_1^V - A_1^D'); xlabel('t');
% subplot(3,2,4); plot(t,y(:,6)-y(:,5), '-'); ylabel('A_2^V - A_2^D'); xlabel('t');
% figure(2); clf;
% plot(y(:,1), y(:,2), '-');
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
