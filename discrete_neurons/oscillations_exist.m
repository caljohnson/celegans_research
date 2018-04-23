%Oscillations_Exist
%Checks if parameters yield oscillations

%2nd order ODE for K
% tau K'' + (1+tau) K' + K = c ?V
%homogeneous solution form:
% K(t) = c1 exp(-t) + c2 exp(-t/tau)
% c2 = (K(0) + K'(0)) / (1-1/tau)
% c1 = K(0) - c2

%what should a0 be from find_alphas?
a0 = -1;

%params
eps = 2; %%-] These two determine thresholds together
I = -0.1;  %%-]

c = 1.3;  %%- Must be bigger than K_V_ON to get past first threshold
tau = 1e10;  %%- defines timescale of system

%results in
K_V_ON = eps/2-I;
K_V_OFF = -eps/2-I;

%worked-out oscillation checks
if (-K_V_OFF + K_V_ON*(tau-1))/tau >= c
    display('shouldnt be oscillations');
elseif (-K_V_OFF/c)^tau >= -K_V_OFF/K_V_ON
    display('shouldnt be oscillations');
else
    display('SHOULD get oscillations');
end

%visualize thresholds
n= 10;
ss = linspace(0,1,10);
l1 = (eps/2 + I)*ones(n,1);
l2 = (eps/2 - I)*ones(n,1);
l23 = linspace(eps/2 - I, -eps/2  - I, n);
l3 = (-eps/2  - I)*ones(n,1);
l4 = (-eps/2 + I)*ones(n,1);
l41 = linspace(-eps/2  + I, eps/2  + I, n);
figure(2); clf;
plot(ss, l1, '<-b'); hold on; text(-0.05, eps/2+I, '\kappa^D_{off}', 'Color', 'blue')
plot(ss, l2, '->r'); text(-0.08, eps/2-I, '\kappa^V_{on}', 'Color', 'red')
plot(ss, l3, '<-r'); text(-0.08, -eps/2-I, '\kappa^V_{off}', 'Color', 'red')
% plot(ss, l23, 'o-r');
plot(ss, l4, '->b'); text(-0.05, -eps/2+I, '\kappa^D_{on}', 'Color', 'blue')
% plot(ss, l41, '*-b'); hold off;
xlabel('S'); ylabel('K');

%IC
k0 = K_V_OFF;
kp0 = a0;
SD = 1;
SV = 0;

%simulation stuff
N = 1e5;
dt = 5e-3;
K = zeros(N,1);
Kp = zeros(N,1);
deltaV = zeros(N,1);
K(1) = k0;
Kp(1) = kp0;

for i = 2:N
   
    %update ?V
    deltaV(i) = SD - SV;
    
    %update k0, kp0
    k0 = K(i-1);
    kp0 = Kp(i-1);
    
    %update homog. soln. coeffs
    c2 = (k0 - c*deltaV(i) + kp0)/ (1-1/tau);
    c1 = k0 - c*deltaV(i) - c2;
    
    %step forward in ODE soln
    K(i) = c1*exp(-dt) + c2*exp(-dt/tau) + c*deltaV(i);
    Kp(i) = -c1*exp(-dt) - (c2/tau)*exp(-dt/tau);
    
   %update states
    if K(i) >= -K_V_OFF
       SD = 0;
    end
    if K(i) >= K_V_ON
       SV = 1;
    end
    if K(i) <= K_V_OFF
       SV = 0;
    end
    if K(i) <= -K_V_ON
       SD = 1;
   end
end

figure(1); clf;
t = 0:dt:dt*(N-1);
subplot(2,1,1); plot(t,K, '-'); ylabel('\kappa'); xlabel('t');
subplot(2,1,2); plot(t,Kp, '-'); ylabel('d\kappa /dt'); xlabel('t');
figure(3); clf;
plot(K, Kp, '-'); hold on; ylabel('d \kappa / dt'); xlabel('\kappa');
line([-c-1,c+1],[0,0], 'Color',[0 0 0])
line([K_V_ON, K_V_ON], [-c-1,c+1], 'Color',[1 0 0])
line([-K_V_OFF, -K_V_OFF], [-c-1,c+1], 'Color',[0 1 0])
line([K_V_OFF, K_V_OFF], [-c-1,c+1], 'Color',[1 0 0])
line([-K_V_ON, -K_V_ON], [-c-1,c+1], 'Color',[0 1 0])
line([-c, -c], [0,0], 'Color',[0 0 0], 'Marker', 'o')
line([c, c], [0,0], 'Color',[0 0 0], 'Marker', 'o')
line([0,0],[-c-1,c+1], 'Color',[0 0 0])
xlim([min(K), max(K)])
ylim([min(Kp), max(Kp)])
% xlim([-c-1, c+1]);
% ylim([-c-1, c+1]);

% hold on;
% plot(deltaV)