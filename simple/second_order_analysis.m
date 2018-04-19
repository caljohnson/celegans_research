%Second_Order_Analysis

%2nd order ODE for K
% tau K'' + (1+tau) K' + K = c ?V
%homogeneous solution form:
% K(t) = c1 exp(-t) + c2 exp(-t/tau)
% c1 = ((-1/tau)K(0) - K'(0)) / (1-1/tau)
% c2 = K(0) - c1

%params
eps = 1;
c = 1.5;
I = 0.1;
tau = 0.1;

n= 10;
ss = linspace(0,1,10);
% l1 = (eps/2 - 1/2 + I)*ones(n,1);
% l2 = (eps/2 + 1/2 - I)*ones(n,1);
% l23 = linspace(eps/2 + 1/2 - I, -eps/2 + 1/2 - I, n);
% l3 = (-eps/2 + 1/2 - I)*ones(n,1);
% l4 = (-eps/2 - 1/2 + I)*ones(n,1);
% l41 = linspace(-eps/2 - 1/2 + I, eps/2 - 1/2 + I, n);
l1 = (eps/2 + I)*ones(n,1);
l2 = (eps/2 - I)*ones(n,1);
l23 = linspace(eps/2 - I, -eps/2  - I, n);
l3 = (-eps/2  - I)*ones(n,1);
l4 = (-eps/2 + I)*ones(n,1);
l41 = linspace(-eps/2  + I, eps/2  + I, n);
figure(2); clf;
plot(ss, l1, '*-g'); hold on;
plot(ss, l2, 'o-r');
plot(ss, l3, 'o-r');
plot(ss, l23, 'o-r');
plot(ss, l4, '*-g'); 
plot(ss, l41, '*-g'); hold off;
xlabel('S'); ylabel('K');



%IC
k0 = -2;
kp0 = 10;
SD = 1;
SV = 0;

%simulation stuff
N = 1000;
dt = 0.005;
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
    c2 = (k0 + kp0)/ (1-1/tau);
    c1 = k0 - c1;
%     c1 = ((-1/tau)*k0 - kp0) / (1-1/tau);
%     c2 = k0 - c1;
    
    %step forward in ODE soln
    K(i) = c1*exp(-dt) + c2*exp(-dt/tau) + c*deltaV(i);
    Kp(i) = -c1*exp(-dt) - (c2/tau)*exp(-dt/tau);
    
   %update states
%    if K(i) >= eps/2 - 1/2 + I
    if K(i) >= eps/2  + I
       SD = 0;
   end
%    if K(i) >= eps/2 + 1/2 - I
    if K(i) >= eps/2 - I
       SV = 1;
   end
%    if K(i) <= -eps/2 + 1/2 - I
    if K(i) <= -eps/2 - I
       SV = 0;
   end
%    if K(i) <= -eps/2 - 1/2 + I
    if K(i) <= -eps/2 + I
       SD = 1;
   end
end

figure(1); clf;
plot(K); 
% hold on;
% plot(deltaV)