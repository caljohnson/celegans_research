%Simple_model_with_discrete_beam_mechanics
%Neural oscillators coupled mechanically via a discrete beam equation

%Mechanical parameters
n = 100; %number of discrete segments in beam
delX = 1/n; %segment length
dt = 0.1;   %time step for discrete PDE backward Euler solve
Kb = 10; %bending stiffness
xi = 1; %viscosity
char_scale = 1000; %char scale Kb/xi

%discrete beam PDE matrix
e = ones(n,1);
A = spdiags([e -2*e e], [0 1 2], n-2, n);

%Neural parameters
%assume neural dynamics occur on a fast timescale compared to
%muscle+mechanics s.t. always at steady state - 0 or 1
t_m = 10; %muscular activity timescale
% I_AVB = 0.5; %driving AVB current
eps_h = 0.2; %hysteresis window
m = 1; %proprioceptive signal distance
k_SR = 100; %stretch receptor weight

%guassian kernel for proprioception
grid = -m:m;
sig = 1; %std dev of kernel
G = exp(-grid.^2./(2*sig^2)); %Gaussian kernel of width 2m
G = 2*G/sum(G); %normalize so half of G has area 1    

%inhomogeneous driving current
I_AVB = 0.5*ones(n-2,1);
I_AVB2 = I_AVB;
% I_AVB2 = 0.5*ones(n-2,1); %breaks dorsal-ventral symmetry
% I_AVB(1) = 10e10;

%Neural functions
% I = @(K) I_AVB + k_SR*K;
S = @(I,s) 1.*(I>=0.5+eps_h*(0.5-s)) + 0.*(I<=0.5+eps_h*(0.5-s));

%driving torque- depends pw-linearly on muscle activities AD, AV
pw_lin = @(A) 0*(A<=0) + A.*(0<A & A<=1) + 1.*(A>=1);
% m0 = @(AD, AV) 100*(pw_lin(AD) - pw_lin(AV));
m0 = @(AD, AV) 100*(AD - AV);

%initial states
y = zeros(n,1);
K = zeros(n-2,1);
AD = rand(n-2,1);
AV = rand(n-2,1);
SD = zeros(n-2,1);
SV = zeros(n-2,1);



save = zeros(n-2, 100/dt-1);
for t = 1:100/dt-1
    %Mechanics!
    %compute new active moment in each body segment
    active_moment = m0(AD, AV);
    
    %compute new segment displacements via backward Euler's method=
    y = (eye(n,n) + char_scale*dt/(delX^4)*(A'*A))\(y-char_scale*dt/(delX^2)*A'*active_moment);

    %compute new curvatures
    K = A*y;
    
    %plot displacements
    figure(1);
    plot(y(1:end));
    text(3,0.04, strcat('t= ',num2str(t*dt)));
    ylim([-1, 1]);  
    pause(0.0001); 
    
    %Neural activity!
    
    %proprioceptive integration
    P = zeros(n-2,1);
    %for neurons at the HEAD, less proprioceptive coupling since fewer
    %anterior neurons
    for ii= 1:m
       %guassian kernel for proprioception - need new one for each size
       grid2 = -ii:ii;
       G2 = exp(-grid2.^2./(2*sig^2)); %Gaussian kernel of width 2m
       G2 = 2*G2/sum(G2); %normalize so half of G has area 1   
       for jj=0:ii-1
            P(ii) = P(ii) + K(ii-jj)*G2(ii-jj);
            %signal filtered with Gaussian kernel
       end
    end
    
    %for neurons m posterior from the head, m anterior neurons yield max
    %proprioceptive coupling
    for ii = 1+m:n-2
        for jj = 0:m-1
            P(ii) = P(ii) + K(ii-jj)*G(m-jj);
            %signal filtered with Gaussian kernel
        end
    end
        
    %compute neural states
    SD = S(I_AVB + k_SR*P, SD);
    SV = S(I_AVB2 - k_SR*P, SV);
    
    %integrate muscle activity
    AD = (SD - SV) + (-(SD - SV)  + AD)*exp(-dt/t_m);
    AV = (SV - SD) + (-(SV - SD)  + AV)*exp(-dt/t_m);
   
    %save neural activity
    save(:,t) = AD;
end

figure(2);
surf(save');
view(2); shading flat;

