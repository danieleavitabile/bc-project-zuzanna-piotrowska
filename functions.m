%% Plotting functions used in NF model
clear all, close all, clc

m = 10; h = 0.4;

F  = @(u) 1./(1+exp(-m*(u-h))); 
W = @(r) (1-abs(r)).*exp(-abs(r));

%Creating a grid
dx = 2*L/n; x = -L+[0:n-1]'*dx;

figure;
plot(x,F(x),'LineWidth',2,'Color','#29719A')
% title('Firing rate')
xlabel('voltage','FontSize',20)
ylabel('firing rate','FontSize',20)
ylim([-0.05 1.05])
xlim([-2 2])

figure;
plot(x,W(x),'LineWidth',2,'Color','#C15E93')
% title('Synaptic kernel')
xlabel('distance','FontSize',20)
ylabel('connectivity','FontSize',20)
ylim([-0.22 1.05])
xlim([-15 15])