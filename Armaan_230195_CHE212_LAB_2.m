clc
clearvars
L = 3;
N = 50;
h = 45;
k = 28;
T0 = 100;
w = 0.2;
H = 0.3;
delx = L/(N - 1);
Ti = 30;
g = 0;
A1 = w*delx;
A2 = H*delx;
A3 = H*w;
p = 2*(w + H);
A = zeros(N,N);
B = zeros(N,1);
A(1,1)=k/delx;
B(1)=k*(T0)/delx - h*(A1 + A2 + A3)*(Ti - T0);
for i=2:N - 1
A(i,i - 1) = k*A3/delx;
A(i,i) = -2*(k*A3/delx)-2*h*(A1 + A2);
A(i,i + 1) = k*A3/delx;
B(i) = -2*h*(A1 + A2)*Ti;
end
A(N,N-1) = k/delx;
A(N,N) = -h*(A1 + A2 + A3) - k/delx;
B(N) = -h*(A1 + A2 + A3)*Ti;
T = A\B;
x = linspace(0,L,N);
figure
plot(x,T,'o');
hold on
Analy_func=@(x)(Ti+(T0-Ti)*exp(-x*sqrt((h*p)/(k*A3))));
plot(x,Analy_func(x));
legend('Numerical','Analytical');
title('Num-Analyt Comparison for Temp profile');
hold off
Q1=abs(k*A3*(T0-T(1))/(delx))
sum_temp=T0;
for i=1:N
sum_temp=sum_temp+T(i);
end
Q2=abs(2*h*delx*(H+w+L)*(sum_temp-N*Ti))
