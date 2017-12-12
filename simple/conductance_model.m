%simple conductance model for understanding
C_mem = 1;
G_mem = 500;
G_AVB = 150;
V_act = -60;
V_rest = -72;
V_AVB = -20.05;

%activation function - induces hysteresis in model
k_act = 500;
G_act = 20;
I_act = @(V) G_act./(1+exp(-k_act*(V-V_act)));

%leak function
I_leak = @(V) G_mem*(V-V_rest);

%AVB driving function
I_AVB = @(V) G_AVB*(V_AVB - V);

%solve for V_abv
V_avb = @(V) V+(1/G_AVB)*(G_mem*(V-V_rest) - I_act(V));

V = linspace(-60.1, -59.9, 100);
y_avb = V_avb(V);
% y_act = I_act(V);
% y_else = -I_leak(V) + I_AVB(V);

figure(1)
% plot(V,y_act); hold on
% plot(V, -y_else);
plot(y_avb,V);
ylabel('V_{eq}');
xlabel('V_{avb}');