%% Bifurcation diagram for h
clear all, close all, clc

p = @(r) r.*exp(-r);
r1 = linspace(0,1,1000);
r2 = linspace(1,9,1000);
figure;
plot(p(r1),r1,'LineWidth',1.5, 'Color','#60B866','LineStyle','--'); 
hold on;
plot(p(r2),r2,'LineWidth',1.5, 'Color','#60B866'); 
legend('Unstable branch','Stable branch','FontSize',14,'Location','northeast','FontSize',16)
ylabel('\Delta','FontSize',20)
xlabel('h','FontSize',20)
title('h = \psi (\Delta)','FontSize',15)

%% Turing-like bifurcation 
s = 1.5; m = 10; th = 0.5; n = 2^10; L = 10*pi; Nn = 100;
T = 100; tspan = [0 T]; 
W = @(r) 1/sqrt(pi)*exp(-r.^2)- 1/(sqrt(pi)*s)*exp(-(r/s).^2);
F = @(u) 1*(1./(1+exp(-m*u+th)) - 1./(1+exp(th)));

%Creating a grid
dx = 2*L/n; x = -L+[0:n-1]'*dx;

% Form matrix (ring geometry)
M = zeros(n,n);
y = W(x)*dx;
iRows = 1:n;
iShift = -n/2:n/2-1;
for i = 1:n
  M(iRows(i),:) = circshift(y, iShift(i));
end

%setting up an initial condition as random noise in interval (a,b)
a = -1;
b = 1;
u0 = a + (b-a).*rand(size(x));

figure;
A = 4;
[t,U] = ode45(@(t,u) -u + A*M*F(u), tspan, u0);
plot(x,U(end,:),'LineWidth',1.5,'Color','#F59924'); hold on

A = 1;
[t,U] = ode45(@(t,u) -u + A*M*F(u), tspan, u0);
plot(x,U(end,:),'LineWidth',1.5, 'Color','#3F96BF'); hold off
legend('A=4','A=1','FontSize',16)
ylabel('brain activity','FontSize',20)
xlabel('space','FontSize',20)
title('Steady states for different A','FontSize',15)
xlim([-25 25])

%% One bump vs no bump steady state 
m = 10; T = 100; tspan = [0 T]; 
W = @(r) (1-abs(r)).*exp(-abs(r));

%Creating a grid
dx = 2*L/n; x = -L+[0:n-1]'*dx;

% Form matrix (ring geometry)
M = zeros(n,n);
y = W(x)*dx;
iRows = 1:n;
iShift = -n/2:n/2-1;
for i = 1:n
  M(iRows(i),:) = circshift(y, iShift(i));
end

%setting up an initial condition as random noise in interval (a,b)
u0 = 0.9./cosh(0.5*x).^2;

figure;
h = 0.3;
F  = @(u) 1./(1+exp(-m*(u-h)));
[t,U] = ode45(@(t,u) -u + M*F(u), tspan, u0);
plot(x,U(end,:),'LineWidth',1.5,'Color','#F59924'); hold on

h = 0.5;
F  = @(u) 1./(1+exp(-m*(u-h)));
[t,U] = ode45(@(t,u) -u + M*F(u), tspan, u0);
plot(x,U(end,:),'LineWidth',1.5, 'Color','#3F96BF'); hold off
legend('h=0.3','h=0.5','FontSize',16)
ylabel('brain activity','FontSize',20)
xlabel('space','FontSize',20)
title('Steady states for different h','FontSize',15)
xlim([-25 25])