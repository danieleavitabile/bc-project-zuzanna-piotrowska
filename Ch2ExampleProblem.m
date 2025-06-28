%% Random parameter in a simple case - mean and variance
clear all, close all, clc

rng(1) 
% rng(4) %better MC approximation

T = 3; %final time
c = 1; d = 2; %interval for uniform distribution 
NM = 20; %number of solutions for MC approximation
NS = 20; %number of nodes in SC 
u0 = 1; %initial condition 


%setting up needed variables 
tspan = 0:0.05:T;

Emc = 0; E2mc = 0;
Errorms = zeros(1,NM); Errorms2 = zeros(1,NM);

Esc = 0; E2sc = 0;
Errorsc = zeros(1,NM); Errorsc2 = zeros(1,NM);

%Pen and paper results 
ExpHand = @(t) (-exp(-d*t) + exp(-c*t))./((d-c)*t);
VarHand = @(t) (-exp(-2*d*t)+exp(-2*c*t))./(2*t*(d-c)) - ExpHand(t).^2 ;

%plotting solutions and calculating monte carlo approximation for every number 1,...n of nodes 
figure, hold;
for i = 1:NM
    A = unifrnd(c,d);
    opts = odeset('RelTol',5e-6,'AbsTol',5e-6); %6,10,14
    [t,u] = ode45(@(t,u) -A*u, tspan, u0, opts);
    plot(t,u,'k','LineWidth',0.001, 'Color',[0.5 0.5 0.5 0.25], 'HandleVisibility','off');  
    drawnow;
    
    Emc = Emc + (1/NM)*u;
    E2mc = E2mc + (1/NM)*(u.^2);

    for j = i:NM
        Errorms(j) = Errorms(j)+(1/j)*u(end);
        Errorms2(j) = Errorms2(j)+(1/j)*(u(end)).^2;
    end    
end

%calculating stochastic collocation approximation
[nodes,weights]=lgwt(NS, c, d);
for i = 1:NS
    A = nodes(i);
    [t,u] = ode45(@(t,u) -A*u, tspan, u0);
    Esc = Esc + u*(1/(d-c))*weights(i);
    E2sc = E2sc + u.^2*(1/(d-c))*weights(i);
end 

%calculating stochastic collocation approximation for every number 1,...n of nodes
tspansc = [0 T];
for i = 1:NM
    [nodes,weights]=lgwt(i,c,d);
    for j = 1:i
        A = nodes(j);
        opts = odeset('RelTol',5e-6,'AbsTol',5e-6); %6,10,14
        [tsc,u] = ode45(@(t,u) -A*u, tspansc, u0, opts);
        Errorsc(i)=Errorsc(i)+u(end)*(1/(d-c))*weights(j);
        Errorsc2(i)=Errorsc2(i)+(u(end))^2*(1/(d-c))*weights(j);
    end
end 

plot(t,ExpHand(t),'Color','#F79210','LineWidth',4,'LineStyle','--');
plot(t,Emc,'r','LineWidth',1.5);
plot(t,Esc,'Color',[0 0.4470 0.7410],'LineWidth',1.5);
xlabel("t",'FontSize',20)
ylabel("v",'FontSize',20)
legend("True Expectation","E_{MC}","E_{SC}",'FontSize',16)
title('Solutions of the example model with A ~ U[1,2]','FontSize',15)
hold off;

%plotting variance
figure; hold on
plot(t,VarHand(t),'LineWidth',4,'Color','#F79210','LineStyle','--')
plot(t, E2mc - (Emc).^2,'r','LineWidth',2)
plot(t, E2sc - (Esc).^2,'Color',[0 0.4470 0.7410],'LineWidth',2)
legend("True Variance","Var_{MC}", "Var_{SC}",'FontSize',16)
xlabel("t",'FontSize',20); ylabel("variance",'FontSize',20)
title('Variance of solutions of the example model with A ~ U[1,2]','FontSize',15)
hold off

%showing the number of steps used to solve ODE with different tolerances
disp(length(tsc))

%plotting errors with time stepped solutions
m = 1:1:NM;
figure; 
hold on;
ExpHandT = ExpHand(T)*ones(1,NM);
plot(m,abs(ExpHandT-Errorms),'r.-','LineWidth',1)
plot(m,0.05*m.^(-1/2),'LineWidth',1,'Color',"#F9D90A")
plot(m,abs(ExpHandT-Errorsc),'b.-','LineWidth',1);
plot(m,10^3*exp(-m*5),'LineWidth',1,'Color',"#9417E5")
title('Errors for different approximations in loglogscale','FontSize',15)
yscale log
% xscale log
xlabel('number of points','FontSize',20); ylabel('error','FontSize',20);
xlim([1 15]); ylim([10^(-19) 10]);
legend("Error E_{MC}",'Error bound for E_{MC}',"Error E_{SC}",'Error bound for E_{SC}', 'Location', 'southwest','FontSize',16)
hold off;


%% Plotting theoretical errors with solution known 

NM = 20; c = 1; d = 2; T = 3; m = 1:1:NM;
ExpHand = @(t) (-exp(-d*t) + exp(-c*t))./((d-c)*t);

Errorms = zeros(1,NM); Errorms2 = zeros(1,NM);
Errorsc = zeros(1,NM); Errorsc2 = zeros(1,NM);

for i = 1:NM
    [nodes,weights]=lgwt(i,c,d);
    for j = 1:i
        A = nodes(j);
        Errorsc(i)=Errorsc(i)+exp(-A*T)*(1/(d-c))*weights(j);
        Errorsc2(i)=Errorsc2(i)+exp(-2*A*T)*(1/(d-c))*weights(j);
    end
end 

for i = 1:NM
    A = unifrnd(c,d);
    for j = i:NM
        Errorms(j) = Errorms(j)+(1/j)*exp(-A*T);
        Errorms2(j) = Errorms2(j)+(1/j)*exp(-2*A*T);
    end    
end

figure; 
hold on;
plot(m,abs(ExpHandT-Errorms),'r.-','LineWidth',1)
plot(m,0.05*m.^(-1/2),'LineWidth',1,'Color',"#F9D90A")
plot(m,abs(ExpHandT-Errorsc),'b.-','LineWidth',1);
plot(m,10^5*exp(-m*10),'LineWidth',1,'Color',"#9417E5")
title('Errors for different approximations in logscale','FontSize',15)
yscale log
% xscale log
xlim([0.8 15])
ylim([10^(-19) 1])
legend("Error E_{MC}",'Error bound for E_{MC}',"Error E_{SC}",'Error bound for E_{SC}', 'Location', 'east','FontSize',16)
xlabel('number of points','FontSize',20)
ylabel('error','FontSize',20)
hold off;

%% Heatmap
clear all, close all, clc

%setting parameters
heatmap = linspace(1,2,50);

T = 3; tspan = 0:0.05:T; 
u0 = 1;
ufinal = zeros(length(heatmap),length(tspan));

for i = 1:length(heatmap)
    A = heatmap(i);   
    [t,u] = ode45(@(t,u) -A*u, tspan, u0);
    ufinal(i,:) = u(:);
end

imagesc(ufinal);
colorbar
xticklabels = 0:0.5:T;
xticks = linspace(1, size(ufinal, 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
xlabel('time','FontSize',20)
ylabel('values of A','FontSize',20)
yticklabels = 1:0.05:2; %0.27:0.018:0.5;
yticks = linspace(1, size(ufinal, 1), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels(:)))
title('Solutions with equally spaced h','FontSize',15)

%% SC vs MC on the same set 

rng(1)

T = 3; %final time
c = 1; d = 2; %interval for uniform distribution 
NS = 20; %number of nodes
u0 = 1; %initial condition 


%setting up needed variables 
tspan = 0:0.05:T;

Emc = 0; Esc = 0;

%Pen and paper results 
ExpHand = @(t) (-exp(-d*t) + exp(-c*t))./((d-c)*t);

%calculating stochastic collocation approximation
[nodes,weights]=lgwt(NS, c, d);
figure; hold;
for i = 1:NS
    A = nodes(i);
    [t,u] = ode45(@(t,u) -A*u, tspan, u0);
    plot(t,u,'k','LineWidth',0.001, 'Color',[0.5 0.5 0.5 0.25], 'HandleVisibility','off');  
    drawnow;
    Esc = Esc + u*(1/(d-c))*weights(i);
    Emc = Emc + (1/NS)*u;
end 
plot(t,ExpHand(t),'Color','#F79210','LineWidth',4,'LineStyle','--'); hold on;
plot(t,Emc,'r','LineWidth',1.5); hold on;
plot(t,Esc,'Color',[0 0.4470 0.7410],'LineWidth',1.5); hold off;
xlabel("t",'FontSize',20)
ylabel("v",'FontSize',20)
legend("True Expectation","E_{MC}","E_{SC}",'FontSize',16)
title('Comparison of the methods on the same set of values','FontSize',15)