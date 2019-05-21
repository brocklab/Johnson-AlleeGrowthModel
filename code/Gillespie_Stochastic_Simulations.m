% This code is first attempt at running stochastic simulations of small
% cell number growth curves to investigate the role of the
% Allee effect in tumor initiation in an in vitro model system.

% ODE form: logistic growth with Allee threshold
%dN/dt = gN(1-A/N)

% for our purposes, we will ignore the carrying capacity term (assume wells
% don't reach carrying capacity/cut off analysis before they do)

% Gillespie Algorithm: 
% P(birth) = dt*[bN - (b-d)A]
% P(death) = dt*[dN ]

close all; clear all; clc

S = load('../out/greenf832.mat');
green= S.green;
S = load('../out/greensum832.mat');
greensum= S.greensum;
%% First try: simple birth death model

% Pull growth rate of bulk population from experiment of 32 cell wells
% g = 0.0188 cells/hr
% variance = 5.56e-5
g = 0.0188;
N_tf= greensum(4).model_gtkoff(end);
netg = g
varg = 5.56e-10;
N0 = 32;
t = 165;
bplusd = ((netg*varg)./N0)*(1./((exp(2*netg*t))-(exp(netg*t))));
d = 0.5*(bplusd-netg);
b = bplusd-d;
%Note this isn't possible... something must be wrong here...
%%
% back-calculate corresponding b/d (see birth-and-death-proc explained for
% derivation)

% assume birth rate = 0.0238, death rate = 0.005
% Set up time and N counter
%Ninit = [1 2 4 8 16 32 48 64 128];
Ninit = [1 2 3 4 5 6 7 8 9];
for m = 1: length(Ninit)
%state(1, :) = Ninit(4); % at time 0, number of cells =N
%tjump(1, :) = 0; % start at time 0
b = 0.0238; % birth rate
d = 0.005; % death rate
A = 0; % Allee Threshold
%birth_n = b*N-(b-d)*A; % birth hazard function
%death_n = d*N; % death hazard function
num_samps = 1000;
num_iters = 100;
take_offs = 0;
state = zeros(num_iters,num_samps);
tstate = zeros(num_iters,num_samps);
state(1, :) = Ninit(m); % at time 0, number of cells =N
tjump(1, :) = 0; % start at time 0
for j = 1:num_samps
    %N = round(normrnd(mu,sigma));
    N=Ninit(m);
    N0 = N;
    time(1)= 0;
for k = 2:num_iters
    birth_n = b*N-(b-d)*A; % birth 
    death_n = d*N; % death 
    if N==0
        N=0;
    else
        r = rand;
        if r< (birth_n)/(birth_n+death_n)
        N = N+1;
        end
        if r>= (birth_n)/(birth_n+death_n)
        N = N-1;
        end
    end
    state(k, j) = N;
    % set time step to be proportional
    tstep =-log(r)/(birth_n + death_n);
    if tstep == Inf
        % make tstep same as previous?
        tstep = 1;
    end
    time = time + tstep;
    tstate(k,j) = time;
    
    % If N goes below 0, cells go extinct, N=0 throughout
    if N <= 0
       state(k:end,j) = 0;
    end
   
end

% fit a g for each sample (1000) at each initial m
    LB = -Inf ;  % Lower Bounds
    UB = Inf; % Upper Bounds
    params0 = 0.02;% Initial Guess...
    options = optimset('TolFun',1e-12,'Tolx',1e-12,'MaxIter',1000,'Display','off','FinDiffRelStep',1e-3);

    [gind]= lsqnonlin(@fit_single_exp_sim, params0, LB, UB, options, state(:,j), tstate(:,j));
    gmat(j, m)= gind;
    
    thres = 0;
     if state(end,j) > thres
        take_offs= take_offs +1;
     end 
    ind(j) = state(end,j)>thres;

end



take_off_pct(m) = 100*(take_offs/num_samps);
% pretend time is constant for now
t = 1:1:100;
% Fit all cell trajectories to a single growth rate (avg growth rate)

LB = -Inf ;  % Lower Bounds
UB = Inf; % Upper Bounds
params0 = 0.01;% Initial Guess...
options = optimset('TolFun',1e-12,'Tolx',1e-12,'MaxIter',1000,'Display','off','FinDiffRelStep',1e-3);
Ntkoff = state(:,ind);
ttkoff = tstate(:,ind);
[gtkoff, resnorm, reslsq]= lsqnonlin(@fit_single_exp_sim, params0, LB, UB, options, Ntkoff, ttkoff);
[gtot, resnorm, reslsq]= lsqnonlin(@fit_single_exp_sim, params0, LB, UB, options, state, tstate);

g_tkoffs(m) = gtkoff;
g_all(m) = gtot;
tmax = max(max(tstate));
tsim = 1:1:0.7*tmax;
N_avg_tkoff = singleexpmodel(gtkoff, N0,tsim);
N_tot = singleexpmodel(gtot, N0, tsim);


subplot(3, 3,m)
%figure;
hold off
for j = 1:num_samps
plot(tstate(:,j), state(:,j))
hold on
end
plot(tsim, N_avg_tkoff, 'k*', 'LineWidth', 2')
plot(tsim, N_tot, 'r*', 'LineWidth',2)
text(tsim(end), N_tot(end-4), ['g_{tot}=',num2str(gtot)]) 
%text(tsim(end), N_avg_tkoff(end-2), ['g_{tkoff}=',num2str(gtkoff)]) 
xlabel ('time')
ylabel('Number of Cells')
%ylim([0, N_avg_tkoff(end)])
title ([ 'N_{seed} =', num2str(Ninit(m)),', N_{sims}= ', num2str(num_samps),', Pct_{take-off}= ',num2str(round(take_off_pct(m),0)), '%'])
end

%% Plot 1 & 2 cell figures and save gmat
save('../out/gmat_sim.mat', 'gmat')
gmat = struct2cell(gmat);
gmat = cell2mat(gmat);
%%
for m = 1:2
subplot(1, 2,m)
%figure;
hold off
for j = 1:num_samps
plot(tstate(:,j), state(:,j))
hold on
end
%plot(tsim, N_avg_tkoff, 'k*', 'LineWidth', 2')
%plot(tsim, N_tot, 'r*', 'LineWidth',2)
%text(tsim(end), N_tot(end-4), ['g_{tot}=',num2str(gtot)]) 
%text(tsim(end), N_avg_tkoff(end-2), ['g_{tkoff}=',num2str(gtkoff)]) 
xlabel ('time')
ylabel('Number of Cells')
ylim([0, 200])
xlim([0, 168])
title ([ 'N_{seed} =', num2str(Ninit(m)),', N_{sims}= ', num2str(num_samps),', Pct_{take-off}= ',num2str(round(take_off_pct(m),0)), '%'])
end
figure;
for j = 1:2
    sigma(j) = std(gmat(:,j));
    mean_g(j)= mean(gmat(:,j));
    med_g(j)= median(gmat(:,j));
    subplot(1,2,j)
    hold off
    histogram(gmat(:,j),[-.02:.001:.05], 'facecolor', 'r')
    xlabel ('time')
ylabel('Frequency')
xlabel('g_{ind}')
ylim([ 0 200])
xlim ([ -0.02 0.05])
title ([ 'N_{seed} =', num2str(Ninit(j)), ', \sigma=', num2str(round(sigma(j),3)),', med=',num2str(round(med_g(j), 4))])
end

%% Distribution of Growth rates for each Ninit
hold off
for j = 1:9
    sigma(j) = std(gmat(:,j));
    mean_g(j)= mean(gmat(:,j));
    med_g(j)= median(gmat(:,j));
    subplot(3,3,j)
    hold off
    histogram(gmat(:,j),[-.02:.001:.05])
    xlabel ('time')
ylabel('Frequency')
xlabel('g_{ind}')
ylim([ 0 200])
xlim ([ -0.02 0.05])
title ([ 'N_{seed} =', num2str(Ninit(j)), ', \sigma=', num2str(round(sigma(j),3)),', med=',num2str(round(med_g(j), 4))])
end

figure;
hold off
bar(sigma)
xlabel('N_{seed}')
ylabel('\sigma of g_{distribution}')
title('Variance in growth rates vs. N_{seed} for stochastic growth model')

%% Look at how percent-take off and growth rate scale with N0
t= 1:1:100;
figure;
subplot(1,4,1)
hold off
plot( Ninit, g_all, '*')
xlabel ('Initial Cell Number')
ylabel ('Average growth rate')
title ('growth rate vs. N_{0}')

subplot(1,4,2)
hold off
bar(sigma)
xlabel('Initial Cell Number')
ylabel('\sigma of g_{distribution}')
title('Variance in g')

subplot(1,4,3)
for m=1:length(Ninit)
    Nmodel(:,m) = singleexpmodel(g_all(m), Ninit(m), t);
    plot(t, Nmodel(:,m), 'LineWidth', 1.5)
    text(t(end-13),Nmodel(end,m), ['N_{0}=', num2str(Ninit(m))])
    hold on
end
xlabel('time (hours)')
ylabel('Number of Cells')
title('Average growth curves')

subplot(1,4,4)
for m = 1:length(Ninit)
plot(t, log(Nmodel(:,m)/Ninit(m)), 'LineWidth', 1.5)
text(t(end-4),log(Nmodel(end,m)/Ninit(m)), ['N_{0}=', num2str(Ninit(m))])
hold on
end
xlabel('time (hours)')
ylabel('log(N/N_{0})')
title('Average per capita growth rate')
%% Plot Percent of Surviving Wells vs. Expected % Surviving Wells if Poisson distributed
PrCSC= take_off_pct(1)/100;
for m = 1:length(Ninit)
expected_prob(m) = 1-((1-PrCSC)^Ninit(m));
expected_pct(m)=100*expected_prob(m);
if expected_pct(m)>100
    expected_pct(m)=100;
end
end
figure;
hold off
plot(Ninit, take_off_pct, 'o-', 'LineWidth',1.5)
hold on
plot(Ninit, expected_pct, 'ro-', 'LineWidth', 1.5)
xlabel('Initial Cell Number')
ylabel('Percent of Surviving Wells')
title('Percent Surviving Wells vs. N_{0}')
legend('Simulated % Surviving Wells', 'Expected % Survivng Wells if Poisson Distributed')
legend boxoff
xlim([ 0 8])