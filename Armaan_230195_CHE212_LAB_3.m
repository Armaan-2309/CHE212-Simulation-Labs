clc
clearvars


L = 0.005;
N = 100;
h = 4.5;
k = 28;
T0 = 100;
w = 0.01;
H = 0.01;
delx = L/(N - 1);
Ti = 30;
g = 0;
A1 = w*delx;
A2 = H*delx;
Ac = H*w;
p = 2*(w + H);
A = zeros(N,N);
B = zeros(N,1);
A(1,1)=k/delx;
B(1)=k*(T0)/delx - h*(A1 + A2 + Ac)*(Ti - T0);
for i=2:N - 1
A(i,i - 1) = k*Ac/delx;
A(i,i) = -2*(k*Ac/delx)-2*h*(A1 + A2);
A(i,i + 1) = k*Ac/delx;
B(i) = -2*h*(A1 + A2)*Ti;
end
A(N,N-1) = k/delx;
A(N,N) = -h*(A1 + A2) - k/delx;
B(N) = -h*(A1 + A2)*Ti;
T = A\B;
x = linspace(0,L,N);
figure
plot(x,T,'o');
hold on
m = sqrt((h*p)/(k*Ac));
Analy_func=@(x)(Ti+(T0-Ti)*(exp(m*(L - x)) + exp(m*(x - L)))/(exp(m*L) + exp(-m*L)));
plot(x,Analy_func(x),LineWidth=1.5);
title('Graph for h = 4.5 and h = 450');
hold on
Q1 = abs(k*Ac*(T(2)-T(1))/(delx))
Q2 = 0;
for i=2:N-1
    Q2 = Q2 + h*p*delx*(T(i) - Ti);
end
Q2 = Q2
result = 1 - Q1/Q2


%% part 2


L = 0.04;
N = 100;
h = 450;
k = 28;
T0 = 100;
w = 0.01;
H = 0.01;
delx = L/(N - 1);
Ti = 30;
g = 0;
A1 = w*delx;
A2 = H*delx;
Ac = H*w;
p = 2*(w + H);
A = zeros(N,N);
B = zeros(N,1);
A(1,1)=k/delx;
B(1)=k*(T0)/delx - h*(A1 + A2 + Ac)*(Ti - T0);
for i=2:N - 1
A(i,i - 1) = k*Ac/delx;
A(i,i) = -2*(k*Ac/delx)-2*h*(A1 + A2);
A(i,i + 1) = k*Ac/delx;
B(i) = -2*h*(A1 + A2)*Ti;
end
A(N,N-1) = k/delx;
A(N,N) = -h*(A1 + A2) - k/delx;
B(N) = -h*(A1 + A2)*Ti;
T = A\B;
x = linspace(0,L,N);
plot(x,T,'o');
hold on
m = sqrt((h*p)/(k*Ac));
Analy_func=@(x)(Ti+(T0-Ti)*(exp(m*(L - x)) + exp(m*(x - L)))/(exp(m*L) + exp(-m*L)));
plot(x,Analy_func(x),LineWidth=1.5);
legend('h = 4.5','h = 4.5(analyt)','h = 450','h = 450(analyt)')
hold off
Q1 = abs(k*Ac*(T(2)-T(1))/(delx));
Q2 = 0;
for i=2:N-1
    Q2 = Q2 + h*p*delx*(T(i) - Ti);
end
Q2 = Q2;
result = 1 - Q1/Q2;

%% calculating h
clc
clear
clearvars

L = 0.1;
Tb = 100 + 273;
T_inf = 30 + 273;
k = 28;
t = 0.01;
k_air = 0.026;
Pr = 0.7;
v = 1.5e-5;
B = 1/T_inf;
As = L*t;
P = 2*(L + t);
Lc = As/P;
g = 9.8;
Ts = 100 + 273;       % initial assumption
Ra = g*B*(Ts - T_inf)*(Lc^3)*Pr/(v^2);

if 1e2 < Ra && Ra <= 1e7
    Nu = 0.54*Ra^(1/4);
end
if Ra > 1e7 && Ra<=1e11
    Nu = 0.15*Ra^(1/3);
end
h = k_air*Nu/Lc



%% comparision of h values
h_list = [4.5,15.6647,45,450];
figure
for j = 1:4
    L = 0.1;
    N = 50;
    h = h_list(j);
    k = 28;
    T0 = 100;
    w = 0.01;
    H = 0.01;
    delx = L/(N - 1);
    Ti = 30;
    g = 0;
    A1 = w*delx;
    A2 = H*delx;
    Ac = H*w;
    p = 2*(w + H);
    A = zeros(N,N);
    B = zeros(N,1);
    A(1,1)=k/delx;
    B(1)=k*(T0)/delx - h*(A1 + A2 + Ac)*(Ti - T0);
    for i=2:N - 1
    A(i,i - 1) = k*Ac/delx;
    A(i,i) = -2*(k*Ac/delx)-2*h*(A1 + A2);
    A(i,i + 1) = k*Ac/delx;
    B(i) = -2*h*(A1 + A2)*Ti;
    end
    A(N,N-1) = k/delx;
    A(N,N) = -h*(A1 + A2) - k/delx;
    B(N) = -h*(A1 + A2)*Ti;
    T = A\B;
    x = linspace(0,L,N);
    plot(x,T);
    hold on
    m = sqrt((h*p)/(k*Ac));
%     Analy_func=@(x)(Ti+(T0-Ti)*(exp(m*(L - x)) + exp(m*(x - L)))/(exp(m*L) + exp(-m*L)));
%     plot(x,Analy_func(x),LineWidth=1.5);
    title('Num-Analyt Comparison for Temp profile');
    hold on
    Q1 = abs(k*Ac*(T(2)-T(1))/(delx));
    Q2 = 0;
    for i=2:N-1
        Q2 = Q2 + h*p*delx*(T(i) - Ti);
    end
    Q2 = Q2;
    result = 1 - Q1/Q2;
end
legend('4.5','15.6647(obtained)','45','450')
