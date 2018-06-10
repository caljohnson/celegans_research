
%c_vs_I_oscillations
%runs alpha_star map to find if oscillations exist for various c,I
%and plots c vs I vs alpha_star

%fixed parameters
%params
eps = 2;
% I = 0.01; %must be >0 for now



%iterate over parameters c,tau
cees = 1.1:0.2:5;
eyes = 0.001:0.01:0.1;
tau  = 5; 

alpha_star = zeros(size(cees,2), size(eyes,1));
period_LC = zeros(size(cees,2), size(eyes,1));
amp_K = zeros(size(cees,2), size(eyes,1));

for j = 1:size(cees,2)
 for k = 1:size(eyes,2)
     
     I = eyes(k);
     %results in
     K_V_ON = eps/2-I;
     K_V_OFF = -eps/2-I;
     c = cees(j);
     
     
%MAPS 1 and 2
alpha0andHalf = @(a0,ah) -(K_V_ON -c + tau.*ah)./(a0.*tau + K_V_OFF - c) + ...
    ((K_V_ON-c+ah)./(K_V_OFF-c+a0)).^tau;

alphaHalfand1 = @(ah, a1) -(-K_V_OFF + tau.*a1)./(ah.*tau + K_V_ON) + ...
    ((-K_V_OFF+a1)./(K_V_ON+ah)).^tau;

%times 1/2 and 1
t_half = @(a0,ah) -tau*reallog((K_V_ON - c + ah)/(K_V_OFF-c+a0));
t_one = @(ah,a1) -tau*reallog((-K_V_OFF + a1)/(K_V_ON+ah));

N = 100;
a0s = linspace(-10,0,N)';
ahs = zeros(N,1);
ap1s = zeros(N,1);
a3hs = zeros(N,1);
ap2s = zeros(N,1);

options = optimset('Display','off');

for i=1:N
    if tau>1
        upper_bound = (c-K_V_ON)/tau;
    else
        upper_bound = c-K_V_ON;
    end
    if sign(alpha0andHalf(a0s(i), 0))==sign(alpha0andHalf(a0s(i), upper_bound)) || ...
            ~isfinite(alpha0andHalf(a0s(i), 0)) || ~isfinite(alpha0andHalf(a0s(i), upper_bound))
        exitflag=0;
    else
    [ahs(i), ~, exitflag, ~] = ...
        fzero(@(ah) alpha0andHalf(a0s(i), ah), [0, upper_bound], options);
    if exitflag==1 && isfinite(alphaHalfand1(ahs(i),ahs(i)))
        [ap1s(i), ~, exitflag, ~] = ...
            fzero(@(a1) alphaHalfand1(ahs(i), a1), [ahs(i)], options);
        if exitflag == 1 && ...
                sign(alpha0andHalf(-ap1s(i), 0))~=sign(alpha0andHalf(-ap1s(i), upper_bound)) ...
                && isfinite(alpha0andHalf(-ap1s(i), 0)) && isfinite(alpha0andHalf(-ap1s(i), upper_bound))
            [x, ~, exitflag, ~] = ...
                fzero(@(a3h) alpha0andHalf(-ap1s(i), a3h),[0, upper_bound], options);
            a3hs(i) = -x;
            if exitflag==1 && isfinite(alphaHalfand1(-a3hs(i),-a3hs(i)))
               [x, ~, exitflag, ~]...
                   = fzero(@(a2) alphaHalfand1(-a3hs(i), a2), [-a3hs(i)], options); 
               ap2s(i) = -x;
            else
                exitflag=0;
            end
        else
            exitflag=0;
        end
    else
        exitflag=0;
    end
    end
    if exitflag~=1
        ap2s(i) = nan;
    end
    
    
end

[~, ii] = min(abs(a0s - ap2s));
alpha_star(j,k) = ap2s(ii);
%get LC period for valid c,I
if ~isnan(ap2s(ii))
    period_LC(j,k) = 2*t_half(a0s(ii), ahs(ii)) + 2*t_one(ahs(ii), ap1s(ii));
else
     period_LC(j,k) = nan;
end

%%%%%compute amplitude of LC
kappa_dot = @(t,kappa, A) (-k/mu)*(kappa - (A(1) -A(2))/(2*a*delX));
state_v_1 = discrete_neural_state_init(0, K_V_OFF, K_V_OFF, K_V_ON, 1);
state_d_1 = discrete_neural_state_init(1, K_V_OFF, K_V_OFF, K_V_ON, 0);
muscle_activity = @(t, K,A) [-A(1) + (state_d_1(K) - state_v_1(K)); ...
                -A(2) + (state_v_1(K) - state_d_1(K));];
%integrate system
ode_rhss = @(t,X) [kappa_dot(t,X(1), X(2:3)); muscle_activity(t,X(1),X(2:3));];
init_cond = [K_V_OFF; 0;0;];
[t,y] = ode113(ode_rhss,[0,1e2], init_cond);
%get amplitude of LC from time series
amp_K = max(y(end-2*period_LC:end,1));

end
end

%plot period
[eyes,cees] = meshgrid(eyes(:,1)', cees(1,:));
surf(cees, eyes, period_LC);
h = colorbar; set(h, 'ylim', [0 8])
xlabel('c'); ylabel('I');
title(strcat('tau= ', num2str(tau)));
h = colorbar;
ylabel(h, 'Period')
xlim([1.1,5])

%plot amp
surf(cees, eyes, amp_K);
h = colorbar; set(h, 'ylim', [0 8])
xlabel('c'); ylabel('I');
title(strcat('tau= ', num2str(tau)));
h = colorbar;
ylabel(h, 'Amplitude (\Kappa)')
xlim([1.1,5])

