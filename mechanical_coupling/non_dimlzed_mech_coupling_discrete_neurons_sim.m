%non_dimlzed_Mech_coupling_Discrete_neurons_Sim
%simulates mechanical coupling on discrete neuron-system pair
%with "nondimensionalized" parameters

%simulation runtime
TF = 5e1;

%mechanical params
Gamma = 10; %parameter sweep this - mech coupling strength
tau_f = 5; %mu/k
c_MA = 10; %muscle activity mechanical feedback strength
%timescale of muscle activation
tau_m = 1;

%simulate single-oscillator to get limit cycle and pick initial phases
[ T, cycle, tees ] = get_period_mech_coupling( tau_f, c_MA, tau_m );

% % Mechanically coupled system 
LHS_matrix = tau_f.*[(Gamma/tau_f)+1, (Gamma/tau_f)+1; 1, -1;];
RHS_matrix = [-1, -1; -1, +1;];
%R*[theta1'; theta2;] = L[theta1; theta2;] - c_MA*[M1+M2; M1-M2;]
kappa_dot = @(t,kappa, A) LHS_matrix\(RHS_matrix*[kappa(1); kappa(2);]...
    -c_MA.*[(A(2)-A(1))+(A(4)-A(3));(A(2)-A(1))-(A(4)-A(3));]); 
%GAMMA = 0:
% kappa_dot = @(t,kappa,A) (1/tau_f).*[-kappa(1) - c_MA*(A(2)-A(1));...
%     -kappa(2) - c_MA*(A(4)-A(3));];
   
%neural params
eps = 2;   %%-] These two determine thresholds together
I = 0.01;  %%-]

%results in the thresholds
K_V_ON = eps/2-I;   %K_D_ON is negative of this
K_V_OFF = -eps/2-I; %K_D_OFF is negative of this


%pick initial conditions
T_i = size(tees,1);  %get number of indices in one cycle
% %first at phase 0
% K(1) = cycle(1,1);
% A(1) = cycle(1,2);
% A(2) = cycle(1,3);
% %first at phase 1/2
% K(1) = cycle(round(T_i/2),1);
% A(1) = cycle(round(T_i/2),2);
% A(2) = cycle(round(T_i/2),3);

phases_1 = 1; %round([1; T_i/100; T_i/50; T_i/20; T_i/10; T_i/5;...
%     T_i/4; 3*T_i/8; T_i/2; 5*T_i/8; 3*T_i/4; 7*T_i/8;]);
phases_2 = round([1; T_i/100; T_i/50; T_i/20; T_i/10; T_i/5;...
    T_i/4; 3*T_i/8; T_i/2; 5*T_i/8; 3*T_i/4; 7*T_i/8;]);

max_step = 1e-2;
t0 = 0:max_step:TF;
series = zeros(size(phases_2,1),size(t0,2), 6);
for kk = 1:size(phases_1,1)
for jj = 1:size(phases_2,1)
    phase_1 = phases_1(kk);
    phase_2 = phases_2(jj);
    K(1) = cycle(phase_1,1);
    A(1) = cycle(phase_1,2);
    A(2) = cycle(phase_1,3); 
    K(2) = cycle(phase_2,1);
    A(3) = cycle(phase_2,2);
    A(4) = cycle(phase_2,3); 

%neural functions
state_d_1 = discrete_neural_state_init_coupling(K(1), K_V_OFF, K_V_ON, 0);
state_v_1 = discrete_neural_state_init_coupling(K(1), K_V_OFF, K_V_ON, 1);
state_d_2 = discrete_neural_state_init_coupling2(K(2), K_V_OFF, K_V_ON, 0);
state_v_2 = discrete_neural_state_init_coupling2(K(2), K_V_OFF, K_V_ON, 1);

%muscle eqns:
muscle_activity = @(t, K,A) (1/tau_m).*[-A(1) + (state_d_1(K(1),K(2)) - state_v_1(K(1),K(2))); ...
                -A(2) + (state_v_1(K(1),K(2)) - state_d_1(K(1),K(2)));...
                -A(3) + (state_d_2(K(2), K(1)) - state_v_2(K(2), K(1))); ...
                -A(4) + (state_v_2(K(2), K(1)) - state_d_2(K(2), K(1)));];
            
            

ode_rhss = @(t,X) [kappa_dot(t,X(1:2),X(3:6)); muscle_activity(t,X(1:2),X(3:6));];
init_cond = [K(1); K(2); A(1);A(2);A(3);A(4);];

options = odeset('RelTol',1e-8,'AbsTol',1e-10,  'MaxStep', max_step);
[t,y] = ode23(ode_rhss,[0,TF], init_cond, options);
    
%sample cycle at even intervals
y = interp1(t,y,t0);
t = t0;

series(jj,:,:) = y;
end
end


%plot kappa1 vs kappa2 to see antiphase or in-phase
figure(2);clf;
for jj=1:size(phases_2,1)
    subplot(3,4,jj);  plot(series(jj,:,1),series(jj,:,2), '-'); 
    xlabel('\kappa_1'); ylabel('\kappa_2'); 
    title(strcat('initial phase difference = 1/',num2str(T_i/phase_2)));
% title(strcat('initial phase difference = 1/2 - 1/',num2str(T_i/phase_2)));
end

% phase differences!
% phase_diffs = zeros(size(phases_2,1), size(phases_1,1));
figure(3); clf;
for jj=1:size(phases_2,1)
    y = reshape(series(jj, :,:),TF/max_step+1,6);
    [phase_1,osc_flag1] = map_phase(cycle, y(:, [1 3 4]));
    [phase_2,osc_flag2] = map_phase(cycle, y(:, [2 5 6]));
    if osc_flag1 == 1 && osc_flag2 == 1
        phase_diff = phase_1 - phase_2;
        phase_diff(phase_diff<0) = phase_diff(phase_diff<0) + 1;
    else
        phase_diff = nan*phase_1;
    end

    subplot(3,4,jj);  plot(phase_diff, '-'); 
    ylim([-0.1 1.1]);  %xlim([t(track_inds) t(end)]);
    xlabel('t'); ylabel('\alpha_1 - \alpha_2'); 
     title(strcat('initial phase difference = 1/',num2str(T_i/phase_2)));
% title(strcat('initial phase difference = 1/2 - 1/',num2str(T_i/phase_2)));

    % phase_diffs(jj,kk) = decide_antiphase(y(:,1), y(:,2));
end


% nan_diffs = zeros(size(phase_diffs));
% for jj = 1:size(phase_diffs,1)
%     for kk = 1:size(phase_diffs,2)
%         if isnan(phase_diffs(jj,kk))
%             nan_diffs(jj,kk) = 0;
%         else
%             nan_diffs(jj,kk) = nan;
%         end
%     end
% end
% 
% figure(3); clf;
% % surf(phases_1, phases_2, phase_diffs);
% plot(phases_2./T_i, phase_diffs, 'o');
% hold on; plot(phases_2./T_i, nan_diffs, 'or');
% legend('final phase difference', 'oscillations killed');