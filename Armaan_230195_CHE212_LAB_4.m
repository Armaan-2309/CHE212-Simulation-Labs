clc
clearvars
L = 0.1;
d = 0.003;
N = 5;
h = 50;
k = 50;
T0 = 200;
w = 0.01;
H = 0.01;
delx = L/(N - 1);
Tinf = 40;
g = 0;
p = pi*d;
A1 = p*delx/2;
A2 = H*delx;
Ac = pi/4 * d^2;
A = zeros(N,N);
B = zeros(N,1);
A(1,1)=k/delx;
B(1)=k*(T0)/delx - h*(A1 + Ac)*(Tinf - T0);
for i=2:N - 1
 A(i,i - 1) = k*Ac/delx;
 A(i,i) = -2*(k*Ac/delx)-2*h*(A1);
 A(i,i + 1) = k*Ac/delx;
 B(i) = -2*h*(A1)*Tinf;
end
A(N,N-1) = k*Ac/delx;
A(N,N) = -h*(A1+Ac) - k*Ac/delx;
B(N) = -h*(A1+Ac)*Tinf;
T = A\B;
x = linspace(0,L,N);
figure
plot(x,T,'b');
hold on
%% TRANSIENT CONDUCTION
clc
clearvars
L = 0.1;
d = 0.003;
N = 5;
h = 50;
k = 50;
T0 = 200;
delx = L/(N - 1);
dt = 10;
C = 470;
rho = 7800;
Tinf = 40;
p = pi*d;
A1 = p*delx;
Ac = pi/4 * d^2;
T_old = zeros(1,N);
x = linspace(0,L,N);
for i = 1:N
 T_old(i) = T0;
end
T_new = zeros(1,N);
for i = 1:10
 T_new(1) = 200;
 T_new(N) = (k*Ac/delx*(T_old(N - 1) - T_old(N)) + h*(A1/2 + Ac)*(Tinf -
T_old(N)))*2*dt/(rho*Ac*delx*C) + T_old(N);
 for j = 2:N - 1
 T_new(j) = (k*Ac/delx*(T_old(j - 1) + T_old(j + 1) -2*T_old(j)) + h*A1*(Tinf -
T_old(j)))*dt/(rho*Ac*delx*C) + T_old(j);
 end
 plot(x,T_new)
 % if (i == 1)
 % plot(x,T_new)
 % end
 % if (i == 5)
 % plot(x,T_new)
 % end
 % if (i == 10)
 % plot(x,T_new)
 % end
 % T_old = T_new;

end
% legend('steady state','1st iter','5th iter','10th iter')
legend('steady state','t=10s','t=20s','t=30s','t=40s','t=50s','t=60s','t=70s','t=80s','t=90s','t=100s')
hold off
% plot(x,T_new,'*')
figure
%% Changing dt
clc
clearvars
L = 0.1;
d = 0.003;
N = 5;
h = 50;
k = 50;
T0 = 200;
delx = L/(N - 1);
dt_values = [0.1,0.5,1,3,5,10];
C = 470;
rho = 7800;
Tinf = 40;
p = pi*d;
A1 = p*delx;
Ac = pi/4 * d^2;
T_old = zeros(1,N);
x = linspace(0,L,N);
for i = 1:N
 T_old(i) = T0;
end
T_new = zeros(1,N);
for i = 1:6
 dt = dt_values(i);
 T_new(1) = 200;
 T_new(N) = (k*Ac/delx*(T_old(N - 1) - T_old(N)) + h*(A1/2 + Ac)*(Tinf -
T_old(N)))*2*dt/(rho*Ac*delx*C) + T_old(N);
 for j = 2:N - 1
 T_new(j) = (k*Ac/delx*(T_old(j - 1) + T_old(j + 1) -2*T_old(j)) + h*A1*(Tinf -
T_old(j)))*dt/(rho*Ac*delx*C) + T_old(j);
 end
 plot(x,T_new)
 hold on
 T_old = T_new;

end
legend('dt = 0.1s','dt = 0.5s','dt = 1s','dt = 3s','dt = 5s','dt = 10s')