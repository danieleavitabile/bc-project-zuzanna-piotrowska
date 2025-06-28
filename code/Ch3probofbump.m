%% error in probability approximation 
% rng(10)
rng(1010)
% proberror = 20:30:170;
c = 0.27; d = 0.5; m = 100;
proberror = 20:30:200;
Emcerror = zeros(length(proberror),height(x));
Escerror = zeros(length(proberror),height(x));
MCmeancounterror = zeros(1,length(proberror));
SCmeancounterror = zeros(1,length(proberror));

for i = 1:length(proberror)
    [nodes,weights]=lgwt(proberror(i), c, d);
    for j = 1:proberror(i)
        h = unifrnd(c,d);
        F  = @(u) 1./(1+exp(-m*(u-h)));
        [t,U] = ode45(@(t,u) -u + M*F(u), tspan, u0);
        Emcerror(i,:) = Emcerror(i,:) + (1/proberror(i))*U(end,:);

        h = nodes(j);
        F  = @(u) 1./(1+exp(-m*(u-h)));
        [t,U] = ode45(@(t,u) -u + M*F(u), tspan, u0);
        Escerror(i,:) = Escerror(i,:) + U(end,:).*(1/(d-c)).*weights(j);
    end
    MCmeancounterror(i) = max(Emcerror(i,:))./0.7175;
    SCmeancounterror(i) = max(Escerror(i,:))./0.7175;

end

% for i = 1:length(proberror)
%     [nodes,weights]=lgwt(proberror(i), c, d);
%     for j = 1:proberror(i)
%         h = unifrnd(c,d);
%         F  = @(u) 1./(1+exp(-m*(u-h)));
%         [t,U] = ode45(@(t,u) -u + M*F(u), tspan, u0);
%         MCmeancounterror(i) = MCmeancounterror(i) + (1/proberror(i))*max(U(end,:));
% 
%         h = nodes(j);
%         F  = @(u) 1./(1+exp(-m*(u-h)));
%         [t,U] = ode45(@(t,u) -u + M*F(u), tspan, u0);
%         SCmeancounterror(i) = SCmeancounterror(i) + max(U(end,:)).*(1/(d-c)).*weights(j);
%     end
%     MCmeancounterror(i) = MCmeancounterror(i)./0.7175;
%     SCmeancounterror(i) = SCmeancounterror(i)./0.7175;
%  end

theoreticalprob = (1/exp(1)-c)/(d-c);
theory = ones(1,length(proberror)) .* theoreticalprob;

figure;
plot(proberror, theory,'Color','#F79210','LineWidth',1.5); hold on;
plot(proberror, MCmeancounterror,'r.-'); hold on;
plot(proberror, SCmeancounterror,'b.-'); hold off;
legend('True probability','MC approximation','SC approximation','FontSize',16)
title('Probability of seeing a bump for h~U[0.27,0.5]','FontSize',15)
xlabel("number of points",'FontSize',20)
ylabel("probability",'FontSize',20)

%% error in probability approximation for MC only
rng(1010)
c = 0.27; d = 0.5; m = 100;
proberror = 100:100:1000;
Emcerror = zeros(length(proberror),height(x));
MCmeancounterror = zeros(1,length(proberror));

for i = 1:length(proberror)
    for j = 1:proberror(i)
        h = unifrnd(c,d);
        F  = @(u) 1./(1+exp(-m*(u-h)));
        [t,U] = ode45(@(t,u) -u + M*F(u), tspan, u0);
        Emcerror(i,:) = Emcerror(i,:) + (1/proberror(i))*U(end,:);
    end
    MCmeancounterror(i) = max(Emcerror(i,:))./0.7175;
end

proberror = 100:100:1000;
theory = ones(1,length(proberror)) .* theoreticalprob;

figure;
plot(proberror, theory,'Color','#F79210','LineWidth',1.5); hold on;
plot(proberror, MCmeancounterror,'r.-'); hold off;
legend('True probability','MC approximation','FontSize',16)
title('Probability of seeing a bump with h~U[0.27,0.5]','FontSize',15)
xlabel("number of points",'FontSize',20)
ylabel("probability",'FontSize',20)

