%V neuron
is_V = 1;
K_V_ON = 1;
K_V_OFF = -1;
SV = discrete_neural_state_init(0, -2, K_V_OFF, K_V_ON, is_V);

s1 = zeros(100,1);
s2 = zeros(100,1);

k1 = linspace(-2,2,100);
k2 = fliplr(k1);

for i =1:100
    s1(i) = SV(k1(i));
end
for i =1:100
    s2(i) = SV(k2(i));
end

figure(1); clf;
plot(k1, s1, '->'); hold on
plot(k2, s2, '<-');
xlabel('K');
ylabel('S');
legend('increasing K', 'decreasing K');

%%D neuron
is_V = 0;
SD = discrete_neural_state_init(0, -2, K_V_OFF, K_V_ON, is_V);

s1 = zeros(100,1);
s2 = zeros(100,1);

k1 = linspace(-2,2,100);
k2 = fliplr(k1);

for i =1:100
    s1(i) = SD(k1(i));
end
for i =1:100
    s2(i) = SD(k2(i));
end

figure(2); clf
plot(k1, s1, '->'); hold on
plot(k2, s2, '<-');
xlabel('K');
ylabel('S');
legend('increasing K', 'decreasing K');


