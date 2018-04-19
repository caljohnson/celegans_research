%finds alpha_0 and alpha_1/2 that yield oscillations

%fixed parameters
%params
eps = 0.5;
c = 10;
I = 0.01;
tau = 0.4;

%results in
K_V_ON = eps/2-I;
K_V_OFF = -eps/2-I;

%MAPS 1 and 2
alpha0andHalf = @(a0,ah) -(K_V_ON -c + tau.*ah)./(a0.*tau + K_V_OFF - c) + ...
    ((K_V_ON-c+ah)./(K_V_OFF-c+a0)).^tau;

alphaHalfand1 = @(ah, a1) -(-K_V_OFF + tau.*a1)./(ah.*tau + K_V_ON) + ...
    ((-K_V_OFF+a1)./(K_V_ON+ah)).^tau;

N = 30;
a0s = linspace(-10,10,N);
ahs = zeros(N,1);
ap1s = zeros(N,1);
a3hs = zeros(N,1);
ap2s = zeros(N,1);

for i=1:N
    ahs(i) = fzero(@(ah) alpha0andHalf(a0s(i), ah), [0, abs(K_V_ON - c)]);
    ap1s(i) = fzero(@(a1) alphaHalfand1(ahs(i), a1), [0,ahs(i)]);
    a3hs(i) = -fzero(@(a3h) alpha0andHalf(-ap1s(i), a3h),[0, abs(K_V_ON - c)]);
    ap2s(i) = -fzero(@(a2) alphaHalfand1(-a3hs(i), a2), [0, -a3hs(i)]);
end

figure(4); clf;
% plot(a0s, ahs, 'o'); hold on
% plot(a0s, ap1s, 'o'); hold on
% plot(a0s, a3hs, 'o'); hold on
plot(a0s, ap2s, 'o-'); hold on
plot(a0s, a0s, '--k');hold off
% legend('\alpha_{1/2}', '\alpha_{1}', 'y=x');
xlabel('\alpha_0'); ylabel('\alpha_2');

