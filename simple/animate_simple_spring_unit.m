function [  ] = animate_simple_spring_unit( y, tspan )
%ANIMATE_SIMPLE_SPRING_UNIT animates the model spring unit
%   takes input y = [m, AD, AV] to determine spring parameters at each time
%   step in tspan

pw_lin = @(A) 0*(A<=0) + A*(0<A & A<=1) + 1*(A>=1);

m = y(1,:);
AD = y(2,:);
AV = y(3,:);

L_D0_init = 1;
L_V0_init = 1;
L_D0 = L_D0_init;
L_V0 = L_V0_init;

%draw starting plot
figure(5);xlim([-0.1,1.1]); ylim([-1.1,1.1])
line([1/2-L_D0/2, 1/2 + L_D0/2], [1, 1]); %dorsal spring, starting
line([1/2 - L_V0/2, 1/2+L_V0/2], [-1, -1]); %ventral spring, starting
line([0,1], [0,0]); %center rod, immutable
%side rods
line([0, 1/2-L_D0/2], [0,1]) %left middle-dorsal connection
line([1, 1/2+L_D0/2], [0,1]) %right middle-dorsal connection
line([0, 1/2-L_V0/2], [0,-1]) %left middle-ventral connection
line([1, 1/2+L_V0/2], [0,-1]) %right middle-ventral connection

for t=1:50:size(tspan,2)-1
%     L_D0 = 1.5 - pw_lin(AD(t));
%     L_V0 = 1.5 - pw_lin(AV(t));
    %evaluate current spring lengths using torque m
    L_D = 1 + m(t)/2;
    L_V = 1 - m(t)/2;
   
    %plot
    figure(5); clf; xlim([-0.1,1.1]); ylim([-1.1,1.1])
    line([1/2-L_D/2, 1/2 + L_D/2], [1, 1]); %dorsal spring, starting
    line([1/2 - L_V/2, 1/2+L_V/2], [-1, -1]); %ventral spring, starting
    line([0,1], [0,0]); %center rod, immutable
    %side rods
    line([0, 1/2-L_D/2], [0,1]) %left middle-dorsal connection
    line([1, 1/2+L_D/2], [0,1]) %right middle-dorsal connection
    line([0, 1/2-L_V/2], [0,-1]) %left middle-ventral connection
    line([1, 1/2+L_V/2], [0,-1]) %right middle-ventral connection
    txt = strcat('t= ',int2str(t));
    text(0,-1,txt)
    pause(0.001)
end


end

