%Discrete Beam Eqn model

%parameters
% a = 0.1; %thickness of beam
% Ks = 10; %spring constant of discrete spring compartments in beam
% L_0 = 1/n; %length of neutral axis in each segment
n = 100; %number of discrete segments in beam
Kb = 10; %bending stiffness
xi = 10; %viscosity
dt = 0.5; %time step
delX = 0.5; %segment spacing

%initial segment displacements
x = linspace(0,n,n)';
y = sin(pi/100*x);
y = zeros(n,1);

%active moments
tau_D = zeros(n-2,1); %0 for now
% tau_V = zeros(n-2,1); %0 for now

%CPG cycle - add to tau_D in head
t = 0:dt:100;
cpg = 1*sin(0.05*t');
% cpg = 0*t';

%torque 
% w = a*(tau_D - tau_V) + 2*Ks*a^2*L_0*K;

%discrete beam PDE matrix
e = ones(n,1);
A = spdiags([e -2*e e], [0 1 2], n-2, n);

%solve via backward Euler's method
for t = 1:101/dt-1
   %compute new active torque in head (driven beam)
   tau_D(1) = cpg(t);
 
   %compute new segment displacements
%    y = (eye(n,n) + (Kb*dt/xi)*(A'*A))\(y - (Kb*dt/xi)*A'*(tau_D));
    y = (eye(n,n) + (Kb*dt)/(xi*delX^4)*(A'*A))\(y-(Kb*dt)/(xi*delX^2)*A'*(tau_D));

   %plot displacements
   plot(y); 
   ylim([0, 1]); 
    ylim([-0.05, 0.05]);
    pause(0.1); 
end

