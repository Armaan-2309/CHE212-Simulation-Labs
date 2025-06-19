clearvars
clear all
L = 0.04;
k = 28;
g = 5e6;
N = 5;
delx = L/(N-1);
h = 45;
T_0 = 0;
T_inf = 30;
x = linspace(0,L,N);
A = zeros(N,N);
B = zeros(N,1);
A(1,1) = 1;
B(1) = 0;
figure
for i = 2:N-1
A(i,i-1) = 1;
A(i,i) = -2;
A(i,i+1) = 1;
B(i) = -g*delx.^2/k;
end
A(N,N-1) = -1;
A(N,N) = 1 + h*delx/k;
B(N) = h*delx*T_inf/k + g*delx^2/(2*k);
T = A\B;
analy_func = @(x)(0.5*g*h*L^2/k + g*L + T_inf*h).*x/(h*L + k) - g*(x.^2)/(2*k);
plot(x,T)
hold on
plot(x,analy_func(x),'o')
title('Num-Analyt Comparison for N=5')
legend('Numerical', 'Analytical')
hold off
clearvars
clear all
L = 0.04;
k = 28;
g = 5e6;
N = 50;
delx = L/(N-1);
h = 45;
T_0 = 0;
T_inf = 30;
x = linspace(0,L,N);
A = zeros(N,N);
B = zeros(N,1);
A(1,1) = 1;
B(1) = 0;
figure
for i = 2:N-1
A(i,i-1) = 1;
A(i,i) = -2;
A(i,i+1) = 1;
B(i) = -g*delx.^2/k;
end
A(N,N-1) = -1;
A(N,N) = 1 + h*delx/k;
B(N) = h*delx*T_inf/k + g*delx^2/(2*k);
T = A\B;
analy_func = @(x)(0.5*g*h*L^2/k + g*L + T_inf*h).*x/(h*L + k) - g*(x.^2)/(2*k);
plot(x,T)
hold on
plot(x,analy_func(x),'o')
title('Num-Analyt Comparison for N=50')
legend('Numerical', 'Analytical')
hold off
clearvars
clear all
L = 0.04;
k = 28;
g = [5e6,3e6,1e6,9e5,7e5];
N = 50;
delx = L/(N-1);
h = 45;
T_0 = 0;
T_inf = 30;
figure
for j = 1:5
x = linspace(0,L,N);
A = zeros(N,N);
B = zeros(N,1);
A(1,1) = 1;
B(1) = 0;
for i = 2:N-1
A(i,i-1) = 1;
A(i,i) = -2;
A(i,i+1) = 1;
B(i) = -g(j)*delx.^2/k;
end
A(N,N-1) = -1;
A(N,N) = 1 + h*delx/k;
B(N) = h*delx*T_inf/k + g(j)*delx^2/(2*k);
T = A\B;
plot(x,T,'o')
title("Effect of decreasing g")
hold on
end
legend('5e6','3e6','1e6','9e5','7e5')
hold off
clearvars
clear all
L = 0.04;
k = [28,24,20,16,12];
g = 5e6;
N = 50;
delx = L/(N-1);
h = 45;
T_0 = 0;
T_inf = 30;
figure
for j = 1:5
x = linspace(0,L,N);
A = zeros(N,N);
B = zeros(N,1);
A(1,1) = 1;
B(1) = 0;
for i = 2:N-1
A(i,i-1) = 1;
A(i,i) = -2;
A(i,i+1) = 1;
B(i) = -g*delx.^2/k(j);
end
A(N,N-1) = -1;
A(N,N) = 1 + h*delx/k(j);
B(N) = h*delx*T_inf/k(j) + g*delx^2/(2*k(j));
T = A\B;
plot(x,T,'p')
title("Effect of decreasing k(thermal conductivity")
hold on
end
legend('28','24','20','16','12')