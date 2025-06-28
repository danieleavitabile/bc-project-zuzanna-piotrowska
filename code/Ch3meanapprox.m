%% Single bump mean and variance approximation 
clear all, close all, clc

% rng(14)
% rng(15)
rng(1010)

%setting parameters and variables 
m = 10; n = 2^10; L = 20; NSC = 20; N = 100; T = 100; tspan = [0 T]; 

c = 0.27; d = 0.5; %uniform interval for h with bump or no bump
% c = 0.27; d = 0.33; %uniform interval for h with one bump always

Emc100 = 0; E2mc100 = 0; Emc20 = 0; E2mc20 = 0;
Esc20 = 0; E2sc20 = 0; Esc100 = 0; E2sc100 = 0;

MCcount = 0; SCcount = 0; Hcount = 0; H1count = 0;
MCpeaksval = 0; SCpeaksval = 0; SCpeaknumb = 0;

%setting a mexican hat kernel
W = @(r) (1-abs(r)).*exp(-abs(r));

%Creating a grid (cortex)
dx = 2*L/n; x = -L+[0:n-1]'*dx;

% Form matrix (ring geometry)
M = zeros(n,n);
y = W(x)*dx;
iRows = 1:n;
iShift = -n/2:n/2-1;
for i = 1:n
  M(iRows(i),:) = circshift(y, iShift(i));
end

%Setting up initial conditions
u0 = 0.9./cosh(0.5*x).^2;

%Plotting the solutions and mean approximations
uTfinal = zeros(N,n);
for i = 1:N
    h = unifrnd(c,d);
    F  = @(u) 1./(1+exp(-m*(u-h)));
    [t,U] = ode45(@(t,u) -u + M*F(u), tspan, u0);
    uTfinal(i,:) = U(end,:);
    plot(x,U(end,:),'k','LineWidth',0.001, 'Color',[0.5 0.5 0.5 0.2],'HandleVisibility','off')
    xlim([-L L])
    hold on
    drawnow
    Emc100 = Emc100 + (1/N)*U(end,:);
    E2mc100 = E2mc100 + (1/N)*(U(end,:).^2);
    if i < 21
        Emc20 = Emc20 + (1/20)*U(end,:);
        E2mc20 = E2mc20 + (1/20)*(U(end,:).^2);
    end
end
plot(x,Emc20,'LineWidth',1,'DisplayName','E_{MC}^{20} estimator','Color','#2AC60B'); hold on;
plot(x,Emc100,'r','LineWidth',1,'DisplayName','E_{MC}^{100} estimator'); hold on;

[nodes,weights]=lgwt(NSC, c, d);   
for i = 1:NSC
    h = nodes(i);
    F = @(u) 1./(1+exp(-m*(u-h)));
    [t,U] = ode45(@(t,u) -u + M*F(u), tspan, u0);
    Esc20 = Esc20 + U(end,:).*(1/(d-c)).*weights(i);
    E2sc20 = E2sc20 + ((U(end,:)).^2).*(1/(d-c)).*weights(i);
end

[nodes,weights]=lgwt(100, c, d);   
for i = 1:100
    h = nodes(i);
    F = @(u) 1./(1+exp(-m*(u-h)));
    [t,U] = ode45(@(t,u) -u + M*F(u), tspan, u0);
    Esc100 = Esc100 + U(end,:).*(1/(d-c)).*weights(i);
    E2sc100 = E2sc100 + ((U(end,:)).^2).*(1/(d-c)).*weights(i);
end
plot(x,Esc20,'b','LineWidth',1,'DisplayName','E_{SC}^{20} estimator');
xlabel("space x",'FontSize',20); ylabel("brain activity",'FontSize',20)
legend('show','FontSize',16)
title('Solutions of the neural field equation with h~U[0.27,0.5]','FontSize',15)
plot(x,Esc100,'LineWidth',1.5,'DisplayName','E_{SC}^{100} estimator','Color','#F79210','LineStyle','--');
hold off;
% xlim([-2 2])
% ylim([0.06 0.25])

% plotting variance 
figure; hold on
plot(x, E2mc20 - (Emc20).^2,'LineWidth',1.5,'Color','#2AC60B','DisplayName','Var_{MC}^{20} estimator')
plot(x, E2mc100 - (Emc100).^2,'LineWidth',1.5,'Color','r','DisplayName','Var_{MC}^{100} estimator')
plot(x, E2sc20 - (Esc20).^2,'LineWidth',1.5,'Color','b','DisplayName','Var_{SC}^{20} estimator')
plot(x, E2sc100 - (Esc100).^2,'LineWidth',1.5,'LineStyle','--','Color','#F79210','DisplayName','Var_{SC}^{100} estimator')
xlabel("space x",'FontSize',20); ylabel("variance",'FontSize',20)
title('Variance of the solutions to the neural field equation with h ~ U[0.27 0.5]','FontSize',15)
legend('show','Location','best','FontSize',16)
hold off
% xlim([-2 2])
% ylim([0.08 0.098])

%% 3d histogram of the data 
uVals = linspace(-0.5,1,80);
Z = zeros(n,length(uVals)-1);
for i = 1:n
   counts=histcounts(uTfinal(:,i),uVals);
   Z(i,:) = counts;
end
uvals = linspace(-0.9,1.8,79);
[X, UVALS] = meshgrid(x,uvals);

figure;
surf(UVALS, X, Z')
shading interp
colorbar
xlabel('values of u')
ylabel('x grid - cortex')
view(-45,70)

%% Approximate the height of a bump
Esc100 = 0;
m = 100;
[nodes,weights]=lgwt(100, 0.27, 0.36);   
for i = 1:100
    h = nodes(i);
    F = @(u) 1./(1+exp(-m*(u-h)));
    [t,U] = ode45(@(t,u) -u + M*F(u), tspan, u0);
    Esc100 = Esc100 + U(end,:).*(1/(0.36-0.27)).*weights(i);
    E2sc100 = E2sc100 + ((U(end,:)).^2).*(1/(d-c)).*weights(i);
end
max(Esc100)

%% SC vs MC on the same set 
m = 10; n = 2^10; L = 20; NSC = 20;
c = 0.27; d = 0.5; %uniform interval for h
% c = 0.27; d = 0.33;
T = 100; tspan = [0 T]; 
Emc100 = 0; Esc20 = 0;

% Form matrix (ring geometry)
M = zeros(n,n);
y = W(x)*dx;
iRows = 1:n;
iShift = -n/2:n/2-1;
for i = 1:n
  M(iRows(i),:) = circshift(y, iShift(i));
end

%Setting up the main ODE+initial conditions
rhs = @(t,u) -u + M*F(u);
u0 = 0.9./cosh(0.5*x).^2;

[nodes,weights]=lgwt(NSC, c, d);   
for i = 1:NSC
    h = nodes(i);
    F = @(u) 1./(1+exp(-m*(u-h)));
    [t,U] = ode45(@(t,u) -u + M*F(u), tspan, u0);
    Esc20 = Esc20 + U(end,:).*(1/(d-c)).*weights(i);
    Emc100 = Emc100 + (1/NSC)*U(end,:);
end

figure;
plot(x,Emc100,'r','LineWidth',1,'DisplayName','E_{MC}^{100} estimator'); hold on;
plot(x,Esc20,'b','LineWidth',1,'DisplayName','E_{SC}^{20} estimator'); hold off;