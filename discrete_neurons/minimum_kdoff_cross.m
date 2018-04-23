%minimum K_D_OFF cross
%find minimum solution that crosses the threshold K_D_OFF

%params
eps = 2;
c = 2;
I = 0.01;

%results in
K_V_ON = eps/2-I;
K_V_OFF = -eps/2-I;

tau = [0.01; 0.1; 0.2; 0.3; 0.4; 0.5; 0.9; 1.1; 10;];
% tau = 1.1;

for ii = 1:size(tau,1);
   figure(ii); clf;    
   title_plot = strcat('tau = ', num2str(tau(ii)));
   f = @(ah) -tau(ii)*reallog((-K_V_OFF)./(K_V_ON+ah));
   g = @(ah) -reallog((-K_V_OFF)./(K_V_ON+ah.*tau(ii)));
   ah = 0:0.01:c-K_V_ON;
   plot(ah, f(ah)-g(ah));
   title(title_plot);
   xlabel('\alpha_{1/2}'); ylabel('t_1^f - t_1^g');
   hold off
    
    
end