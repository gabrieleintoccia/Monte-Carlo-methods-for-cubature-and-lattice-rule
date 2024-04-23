omega = @(x) x.^2 - x + 1/6;
n = 4001;
s_max = 100;
gamma=ones(s_max,1);
beta=ones(s_max,1);
[z, e2] = fastrank1(n, s_max, omega, gamma, beta)
t=1:s_max;
plot(t,e2(t),'o','DisplayName','n=4001')
hold on
n1=514229;
[z, e2_1] = fastrank1(n1, s_max, omega, gamma, beta)
plot(t,e2_1(t),'o','DisplayName','n=514229')
legend()