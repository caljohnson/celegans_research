%Discrete Beam Eqn model

%parameters
n = 100; %number of discrete segments in beam
delX = 1/n; %segment length

dt = 0.1;
Kb = 10; %bending stiffness
xi = 0.1; %viscosity

%char scale Kb/xi
char_scale = 10; 

%discrete beam PDE matrix
e = ones(n,1);
A = spdiags([e -2*e e], [0 1 2], n-2, n);

%active_moment
active_moment = zeros(n-2,1);

%initial segment displacements
y = zeros(n,1);
  
%solve via backward Euler's method=
for t = 1:1000/dt-1
%     %compute new active torque in each body segment (traveling wave)
%     wave = @(k) 5.*sin((t+k).*dt);
%     active_moment = wave(2:n-1)';
%     
%     %compute new active torques (standing wave)
    active_moment = ones(n-2,1);
    wave = 5*sin(t*dt);
    active_moment = wave*active_moment;

    %compute new segment displacements
    y = (eye(n,n) + char_scale*dt/(delX^4)*(A'*A))\(y-char_scale*dt/(delX^2)*A'*active_moment);

    %plot displacements
    plot(y(1:end));
    text(3,0.04, strcat('t= ',num2str(t*dt)));
    ylim([-1, 1]);  
    pause(0.1); 
end
