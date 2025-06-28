%% Single bump steady state - error convergence
clear all, close all, clc

rng(1010)

%setting parameters and variables 
m = 10; n = 2^10; L = 20; NSC = 100; T = 100; tspan = [0 T]; 
c = 0.27; d = 0.5; %one vs no bump
% c = 0.27; d = 0.33; %one bump only 
Esc = 0; Nconv = 1:1:20;
Escconv = zeros(length(Nconv),n); Escconvmax = zeros(length(Nconv),1);
Emc = zeros(length(Nconv),n); Emcconvmax = zeros(length(Nconv),1);

%setting tolerance 
opts = odeset('RelTol',5e-10,'AbsTol',5e-10);

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

%initial conditions
u0 = 0.9./cosh(0.5*x).^2;
   
[nodes,weights]=lgwt(NSC, c, d);
for i = 1:NSC
    h = nodes(i);
    f  = @(u) 1./(1+exp(-m*(u-h)));
    [t,U] = ode45(@(t,u) -u + M*f(u), tspan, u0, opts);
    Esc = Esc + U(end,:).*(1/(d-c)).*weights(i);
end

%SC and MC convergence 
for i = 1:length(Nconv)
    [nodes,weights]=lgwt(Nconv(i), c, d);
    for j = 1:Nconv(i)
        h = nodes(j);
        f  = @(u) 1./(1+exp(-m*(u-h)));
        [t,U] = ode45(@(t,u) -u + M*f(u), tspan, u0, opts);
        Escconv(i,:) = Escconv(i,:) + U(end,:).*(1/(d-c)).*weights(j);

        h = unifrnd(c,d);
        f  = @(u) 1./(1+exp(-m*(u-h)));
        [t,U] = ode45(@(t,u) -u + M*f(u), tspan, u0);
        Emc(i,:) = Emc(i,:) + (1/Nconv(i))*U(end,:);
    end
    Escconvmax(i) = max(abs(Esc-Escconv(i,:)));
    Emcconvmax(i) = max(abs(Esc-Emc(i,:)));
end

figure;
plot(Nconv,Emcconvmax,'r.-','LineWidth',1); hold on;
plot(Nconv,0.8.*Nconv.^(-1/2),'LineWidth',1,'Color',"#F9D90A"); hold on;
plot(Nconv,Escconvmax,'b.-','LineWidth',1); hold on;
plot(Nconv,0.4.*Nconv.^(-1),'LineWidth',1,'Color',"#9417E5"); hold off;
title('Errors for different approximations in logscale','FontSize',15)
xlabel('iterations','FontSize',20)
ylabel('error','FontSize',20)
% legend("Error E_{MC}","Error E_{SC}", 'Location', 'southwest','FontSize',16)
yscale log
% xscale log 
xlim([1 20])
ylim([0.00005 1])
legend("Error E_{MC}",'Error bound for E_{MC}',"Error E_{SC}",'Error bound for E_{SC}', 'Location', 'best','FontSize',16)

%% Heatmap
clear all, close all, clc

%setting parameters
heatmap = linspace(0.27,0.5,30);

%setting parameters
m = 10; n = 2^10; L = 20; N = 30;
T = 100; tspan = [0 T]; 

% opts = odeset('RelTol',5e-14,'AbsTol',5e-14);

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

uTfinal = zeros(N,n);
for i = 1:length(heatmap)
    h = heatmap(i);   
    f  = @(u) 1./(1+exp(-m*(u-h)));
    [t,U] = ode45(@(t,u) -u + M*f(u), tspan, u0);
    uTfinal(i,:) = U(end,:);
end

imagesc(uTfinal); hold on;
% line([1,1024], [17,17], 'Color', 'r'); hold off;
colorbar
xticklabels = -L:2:L;
xticks = linspace(1, size(uTfinal, 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
ylabel('h','FontSize',20)
xlabel('space x','FontSize',20)

yticklabels = 0.5:-0.018:0.21; %0.27:0.018:0.5;
yticks = linspace(1, size(uTfinal, 1), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels(:)))
title('Solutions of NF with equally spaced h','FontSize',15)

%% Comparison of SC 
clear all, close all, clc

%setting parameters
m = 10; n = 2^10; L = 20; NSC = 100;
T = 100; tspan = [0 T]; 
c = 0.27; d = 0.5; %one vs no bump
% c = 0.27; d = 0.33; %one bump
Esc = 0;

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

figure;
[nodes,weights]=lgwt(NSC, c, d);
for i = 1:NSC
    h = nodes(i);
    f  = @(u) 1./(1+exp(-m*(u-h)));
    [t,U] = ode45(@(t,u) -u + M*f(u), tspan, u0);
    Esc = Esc + U(end,:).*(1/(d-c)).*weights(i);
end
plot(x, Esc, 'DisplayName','E_{SC}^{100} truth','LineWidth',4,'Color',[1, 0, 0 0.4]); hold on;

N = 1:15:91;
Escdiff = zeros(length(N),n);
for i = 1:length(N)
    [nodes,weights]=lgwt(N(i), c, d);
    for j = 1:N(i)
        h = nodes(j);
        f  = @(u) 1./(1+exp(-m*(u-h)));
        [t,U] = ode45(@(t,u) -u + M*f(u), tspan, u0);
        Escdiff(i,:) = Escdiff(i,:) + U(end,:).*(1/(d-c)).*weights(j);
    end
    txt = ['E_{SC}^{m} with m=',num2str(N(i))];
    plot(x,Escdiff(i,:),'LineWidth',1,'DisplayName',txt)
    xlim([-L L])
    hold on
    drawnow
end

fontsize(16,'points')
legend show
title('Stochastic collocation approximations','FontSize',15)
xlabel('space x','FontSize',20)
ylabel('brain activity','FontSize',20)
% xlim([-4 4])
% ylim([0.22 0.36])