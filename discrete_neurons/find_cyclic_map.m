%find cyclic map

%case I>0, I<eps
eps = 1;
I = 0.01;

%strength parameter
c = 1;

%time constant
tau = 0.1;

%set thresholds
K_D_OFF = I + eps/2;
K_V_ON = -I+eps/2;
K_D_ON = I-eps/2;
K_V_OFF = -I-eps/2;

%initial a0
a0 = -1;

for i =1:1000
%piece 1 - K(0)=K_V_OFF, K'(0)=a0
%find constants from IC
B1 = (K_V_OFF + a0 - c)/(1-1/tau);
A1 = K_V_OFF - B1 - c;
%set solution form
K = @(t) A1*exp(-t) + B1*exp(-t/tau) + c;

%find t1 s.t. K(t1) = K_V_ON
t1 = fzero(@(t) K(t)-K_V_ON,1);

%set a1 = K'(t1)
a1 = -A1*exp(-t1) - (B1/tau)*exp(-t1/tau);

%Piece 2 - K(0)=K_V_ON, K'(0)=a1
%find constants from IC
B2 = (K_V_ON + a1)/(1-1/tau);
A2 = K_V_ON - B2;
%set solution form
K = @(t) A2*exp(-t) + B2*exp(-t/tau);

%find t2 s.t. K(t2) = K_D_OFF
t2 = fzero(@(t) K(t)-K_D_OFF,0.1);

%set a2 = K'(t2)
a2 = -A2*exp(-t2) - (B2/tau)*exp(-t2/tau);

%--------------------------------
%piece 3 - K(0)=K_D_OFF, K'(0)=a2
%find constants from IC
B3 = (K_D_OFF + a2 + c)/(1-1/tau);
A3 = K_D_OFF - B3 + c;
%set solution form
K = @(t) A3*exp(-t) + B3*exp(-t/tau) - c;

%find t3 s.t. K(t3) = K_D_ON
t3 = fzero(@(t) K(t)-K_D_ON,1);

%set a1 = K'(t3)
a3 = -A3*exp(-t3) - (B3/tau)*exp(-t3/tau);

%-----------------------
%Piece 4 - K(0)=K_D_ON, K'(0)=a3
%find constants from IC
B4 = (K_D_ON + a3)/(1-1/tau);
A4 = K_D_ON - B4;
%set solution form
K = @(t) A4*exp(-t) + B4*exp(-t/tau);

%find t4 s.t. K(t4) = K_V_OFF
t4 = fzero(@(t) K(t)-K_V_OFF,0.1);

%set a0 = K'(t4)
a0 = -A4*exp(-t4) - (B4/tau)*exp(-t4/tau);

end

