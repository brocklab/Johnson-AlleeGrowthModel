% Analyze Allee Effect
close all; clear all; clc
% This script ensures that the Allee effect model (odefunAllee) works for a
% set of input conditions, and analyzes the resutling effects of growth
% with various parameters and initial conditions
%%
t = 1:2:100;
tbig = repmat(t, 5, 1);
tlong = reshape(tbig, [size(tbig,1)*size(tbig,2), 1]);
A = 5;
carcap = 2e2;
g = 0.01;
N0 = [ 2 10 12 20 ];
for i = 1:length(N0)
[tout,Nmodel(:,i)] = ode45(@(t,Nmodel) odefunAllee(t, Nmodel,g, A, carcap ),t, N0(i));
for j = 1:length(Nmodel(:,i))
    if Nmodel(j,i) <0
        Nmodel(j,i) = 0;
    end
end
end

noise = 5*(1-2*randn(length(tout),length(N0)));

Nfake = Nmodel + noise;
figure;
subplot(1,2,1)
hold off
for j = 1:length(N0)
    hold on
    semilogy(tout, Nmodel(:,j), 'LineWidth', 3)
    text(tout(end-27), Nmodel(end-27,j), [ 'N_{0} = ', num2str(N0(j)) ],'HorizontalAlignment','left','VerticalAlignment','bottom','color','k')
   % plot(tout, Nfake(:,j), 'o')
   %ylim([ 0 200])
end
xlabel('time')
ylabel('N')
title(['Expected growth dynamics, A= ', num2str(A), ', K= ', num2str(carcap)])
%legend('N_{0} = 2', 'N_{0} = 5', 'N_{0} = 10', 'N_{0} = 12', 'N_{0} = 20')
legend boxoff

for j = 1:length(N0)
    percapitag(:,j) = Nmodel(:,j)/N0(j);
end

subplot(1,2,2)
hold on
for j = 1:length(N0)
    plot(tout, log(percapitag(:,j)),'LineWidth', 3);
    text(tout(18), log(percapitag(18,j)), [ 'N_{0} = ', num2str(N0(j)) ],'HorizontalAlignment','left','VerticalAlignment','bottom','color','k')
end
%ylim([-1 2])
xlabel('time')
ylabel('log(N/N_{o})')
title(['Expected per capita growth rate, A= ', num2str(A), ', K= ', num2str(carcap)])
%% Repeat and add in diversity as divider

t = 1:2:100;
A = 16;
carcap = 1e8;
g = 0.02;
q = 2;
N0 = [ 2 5 10  12 20];

[tout,Nmodelq] = ode45(@(t,Nmodelq) odefunAlleeqdiv(t, Nmodelq,g, A, carcap, q ),t, N0);

for i = 1:size(Nmodelq,2)
for j = 1:length(Nmodelq(:,i))
    if Nmodelq(j,i) <0
        Nmodelq(j,i) = 0;
    end
end
end

noise = 5*(1-2*randn(length(tout),length(N0)));
Nfake = Nmodelq + noise;
figure;
hold off
for j = 1:length(N0)
    hold on
    plot(tout, Nmodelq(:,j),'LineWidth', 2)
    text(tout(5), Nmodelq(5,j), [ 'N_{0} = ', num2str(N0(j)) ],'HorizontalAlignment','left','VerticalAlignment','bottom','color','k')
   % plot(tout, Nfake(:,j), 'o')
end
xlabel('time')
ylabel('N')
title(['Theoretical Allee Effect with diversity, A= ', num2str(A), ' q= ', num2str(q) ])


%% Repeat and add in diversity as subtractor

t = 1:2:100;
A = 10;
carcap = 1e8;
g = 0.02;
q = 2;
N0 = [ 2 5 10  16 20];

[tout,Nmodelqs] = ode45(@(t,Nmodelqs) odefunAlleeqsub(t, Nmodelqs,g, A, carcap, q ),t, N0);
noise = 5*(1-2*randn(length(tout),length(N0)));
for i = 1:size(Nmodelqs,2)
for j = 1:length(Nmodelqs(:,i))
    if Nmodelqs(j,i) <0
        Nmodelqs(j,i) = 0;
    end
end
end
Nfake = Nmodelqs + noise;
figure;
hold off
for j = 1:length(N0)
    hold on
    plot(tout, Nmodelqs(:,j),'LineWidth', 2)
    text(tout(5), Nmodelqs(5,j), [ 'N_{0} = ', num2str(N0(j)) ],'HorizontalAlignment','left','VerticalAlignment','bottom','color','k')
   % plot(tout, Nfake(:,j), 'o')
end
xlabel('time')
ylabel('N')
title(['Theoretical Allee Effect with diversity as subtractor, A= ', num2str(A), ' q= ', num2str(q) ])
% Make dummy data for fitting
span = length(Nmodel)*length(N0);
noise = 1-2.*randn(span,1);
N_fake = reshape(Nmodel, [span,1]) + 5*noise;
t_fake = reshape(repmat(tout,1, length(N0)), [span,1]);
%% Make phase diagram of simplest Allee effect equation
g = 0.03;
A = 2;
%%
syms g A N(t)
sol = dsolve( diff(N,t)== g*N*(1-(A/N)))

syms g N(t)
A = 10;
solA = dsolve(diff(N,t) == g*N*(1-A/N))

%%
figure;
hold off
plot(tout, Nmodel)
hold on
plot(t_fake, N_fake, 'o')
xlabel('time')
ylabel('N')

%% 
% Make one large vector of time, Ng, Nr, and Ntot for fitting

for j = 1:length(mix)
    t_all(:,j) = mix(j).time;
    Ng_all(:,j) = mix(j).Ngreen;
    Nr_all(:,j) = mix(j).Nred;
    Ntot_all(:,j) = mix(j).Ntot;
    Ng0_all(j,1) = mix(j).Ng0;
    Nr0_all(j,1)= mix(j).Nr0;
    Ntot0_all(j,1)=mix(j).Ntot0;
end


tspan = length(mix)*length(t_all);
t_long = reshape(t_all, [tspan, 1]);
Ng_long = reshape(Ng_all, [tspan,1]);
Nr_long = reshape(Nr_all, [tspan,1]);
Ntot_long = reshape(Ntot_all, [tspan,1]);

%% 
% Allee effect fitting
Nmeas = {Ng_long, Nr_long, Ntot_long};
LB = zeros(3,1);  % Lower Bounds
UB = [Inf Inf Inf]; % Upper Bounds
params0 = [.02; 100];% Initial Guess...
options = optimset('TolFun',1e-12,'Tolx',1e-12,'MaxIter',1000,'Display','off','FinDiffRelStep',1e-3);
carcap = 1e7; % set K, don't worry about it
num_samps = length(mix);

paramsAllee= lsqnonlin(@fit_Allee, params0, LB, UB, options, Nr_long, carcap, Nr0_all, t_long, num_samps);
