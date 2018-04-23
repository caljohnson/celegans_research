%c_vs_tau_oscillations
%runs alpha_star map to find if oscillations exist for various c,tau
%and plots c vs tau vs alpha_star

%fixed parameters
%params
eps = 2;
I = 0.1; %must be >0 for now

%results in
K_V_ON = eps/2-I;
K_V_OFF = -eps/2-I;

%iterate over parameters c,tau
cees = 1.1:0.2:10;
% taus = [0.01; 0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8; 0.9; 1.1; 1.5; 5; 10; 50;];
taus = logspace(-7,7,40)';  

alpha_star = zeros(size(cees,2), size(taus,1));

for j = 1:size(cees,2)
 for k = 1:size(taus,1)
     c = cees(j);
     tau = taus(k);
%MAPS 1 and 2
alpha0andHalf = @(a0,ah) -(K_V_ON -c + tau.*ah)./(a0.*tau + K_V_OFF - c) + ...
    ((K_V_ON-c+ah)./(K_V_OFF-c+a0)).^tau;

alphaHalfand1 = @(ah, a1) -(-K_V_OFF + tau.*a1)./(ah.*tau + K_V_ON) + ...
    ((-K_V_OFF+a1)./(K_V_ON+ah)).^tau;

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

% figure(4); clf;
% plot(a0s, ahs, 'o'); hold on
% plot(a0s, ap1s, 'o'); hold on
% plot(a0s, a3hs, 'o'); hold on
% plot(a0s, ap2s, '.-'); hold on
% plot(a0s, a0s, '--k');hold off
% legend('\alpha_{1/2}', '\alpha_{1}', 'y=x');
% xlabel('\alpha_0'); ylabel('\alpha_2');

[~, ii] = min(abs(a0s - ap2s));
alpha_star(j,k) = ap2s(ii);
% if ~isnan(alpha_star(j,k))
%     alpha_star(j,k)=1;
% end
 end
end

[taus,cees] = meshgrid(taus(:,1)', cees(1,:));
surf(cees, log(taus), alpha_star);
xlabel('c'); ylabel('log(\tau)');
title(strcat('I=', num2str(I)));
xlim([1.1,10])

