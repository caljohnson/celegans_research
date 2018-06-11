%New_Mech_coupling_no_neurons_Sim
%simulates mechanical coupling with no neural-muscular inputs

%simulation runtime
TF = 1e4;

%mechanical params
gamma = 1e2; %parameter sweep this, calculate stable phase differences
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

% fixed muscle tension - only one side contracts
A_dot = @(t,A) (1/tau_m).*[-A(1)+1; -A(2); -A(3); -A(4)+1;];

%Initial conditions
K(1) = 0;
K(2) = 0;
A(1) = 0;
A(2) = 0;
A(3) = 0;
A(4) = 0;

ode_rhss = @(t,X) [kappa_dot(t,X(1:2),X(3:6)); A_dot(t,X(3:6));];
init_cond = [K(1); K(2);A(1); A(2); A(3); A(4);];

[t,y] = ode113(ode_rhss,[0,TF], init_cond);
            
figure(1); clf;
subplot(3,2,1); plot(t,y(:,1), '-'); ylabel('\kappa_1'); xlabel('t');%ylim([-1 1]);
subplot(3,2,2); plot(t,y(:,2), '-'); ylabel('\kappa_2'); xlabel('t');%ylim([-1 1]);
subplot(3,2,3); plot(t,y(:,3), '-'); ylabel('A_1^D'); xlabel('t');
subplot(3,2,5); plot(t,y(:,4), '-'); ylabel('A_1^V'); xlabel('t');
subplot(3,2,4); plot(t,y(:,5), '-'); ylabel('A_2^D'); xlabel('t');
subplot(3,2,6); plot(t,y(:,6), '-'); ylabel('A_2^V'); xlabel('t');
% subplot(3,2,3); plot(t,y(:,4)-y(:,3), '-'); ylabel('A_1^V - A_1^D'); xlabel('t');
% subplot(3,2,4); plot(t,y(:,6)-y(:,5), '-'); ylabel('A_2^V - A_2^D'); xlabel('t');
