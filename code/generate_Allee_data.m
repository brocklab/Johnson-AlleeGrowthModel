% Moment- Approach fitting of stochastic birth-death-Allee process for small cell
% number data

% This code is first attempt to use moment-approach (so mean and
% variance)of N(t) trajectories to fit small cell number well data.

% We will start by simulating stochastic cell trajectories using a set b
% and d, and A , and see if we can use the following equations for <n(t)> and V(t)
% to perform parameter estimation in Bayesian framework
% For random birth death process we have:

%d<n>/dt=(b-d)<n>(1-A/<n>)


close all; clear all; clc

% Start by generating 20 n0=1 trajectories
% assume birth rate = 0.0238, death rate = 0.005
% Set up time and N counter
Ninit = [2 3 4 5];
%Ninit = 5;
for i = 1:length(Ninit)
b = 0.0233 + .0005; % birth rate
d = 0.0045 + .0005; % death rate
delta(i)= b+d;
A = 2;
%birth_n = b*N; % birth hazard function
%death_n = d*N; % death hazard function
num_samps = 2000;
num_iters = 150;
take_offs = 0;
state = zeros(num_iters,num_samps);
tstate = zeros(num_iters,num_samps);
state(1,:) = Ninit(i); % at time 0, number of cells =N
tjump(1, :) = 0; % start at time 0
ct_extinct = 0;
for j = 1:num_samps

    N=Ninit(i);
    N0 = N;
    time(1)= 0;
for k = 2:num_iters
    birth_n = (b*N)-(b-d)*A; % birth 
    if birth_n <0
        birth_n = 0;
    end
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
    r2=rand;
    tstep =-log(r2)/(birth_n + death_n);
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
    
    thres = 0;
     if state(end,j) > thres
        take_offs= take_offs +1;
     end 
    ind(j) = state(end,j)>thres;


end
  P_takeoff(i)= take_offs/num_samps;
  P_tkoff_theor(i) = 1-((d/b).^Ninit(i));
  end

% find minimum time "measured"
tmin = min(tstate(end, :));

% add constant technical noise to data
sigmaT = 0.1; % guess that on average count is off by one quarter of a cell?
state = state+ round(normrnd(0, sigmaT,size(state)));
% eliminate those below 0
for j=1:size(state,1)
    for i = 1:size(state,2)
        if state(j,i)<0
            state(j,i)=0;
        end
    end
end

%% Plot simulated cell number trajectories out to minimum time reached
figure;
hold off
for j = 1:num_samps
plot(tstate(:,j), state(:,j))
hold on
end
xlim([0, tmin])
xlabel('time (hours)')
ylabel('Number of cells')
title(['Simulated N(t) trajectories for b=', num2str(b), ' & d=', num2str(d)])


%% UNIFORM SAMPLING from stochastic trajectories 
%Want to smooth by sampling from N at tsamp
tstart=0;
tint=2;

tsamp = tstart:tint:90+tstart;
for j = 1:num_samps
    tstoch = tstate(:,j);
    Nstoch = state(:,j);

for i = 1:length(tsamp)
    % find nearest tstate that is less that tsamp
    ind =find(tstoch<=tsamp(i),1,'last');
    tfind = tstoch(ind);
    Nsamp(i,j)=Nstoch(ind);
end

end
mu_data = mean(Nsamp,2);

figure;
hold on 
for j = 1:num_samps
plot(tsamp,Nsamp(:,j), 'r.');
hold on
plot(tstate(:,j), state(:,j),'b.');
end
plot(tsamp, mu_data, 'k-', 'LineWidth',3)
xlabel ('time (hours)')
ylabel('<n> expected')
xlim([0, tsamp(end)])
title('Mean value of N in time for N0=5, A = 2')