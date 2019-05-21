% This script generates simulated data for multiple initial conditions of
% N0 for the b-d-A model. 
% It generates the state matrix, performs uniform sampling, finds the mean
% and variance for each initial cell number, and generates a "musimdatavec"
% "varsimdatavec" and a "timesimdatavec". These will be loaded into the fit
% data stoch all script to perform parameter estimation and model
% selection.
% Note, to simulate data for no Allee effect, just set A =0
close all; clear all; clc;
%%

Ninit = [1 2 3 4 5 6 7];
%Ninit = 7;
musimdatavec = [];
varsimdatavec = [];
timesimdatavec = [];
N0vec = [];
for i = 1:length(Ninit)
b = 0.0233 + .0005; % birth rate
d = 0.0045 + .0005; % death rate
delta= b+d;
A = 0; % change to make this Allee effect or not
%birth_n = b*N; % birth hazard function
%death_n = d*N; % death hazard function
num_samps = 75;
num_iters = 2000;
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


tmin(i) = min(tstate(end, :));

tstart=0;
tint=4;

tsamp = tstart:tint:200+tstart;
tsamp = tsamp';
    for m = 1:num_samps
        tstoch = tstate(:,m);
        Nstoch = state(:,m);

    for n = 1:length(tsamp)
        % find nearest tstate that is less that tsamp
        ind =find(tstoch<=tsamp(n),1,'last');
        tfind = tstoch(ind);
        Nsamp(n,m)=Nstoch(ind);
    end

    end
    mu_data = mean(Nsamp,2);
    noise = 0.*(1-2*randn(length(mu_data),1)); % remove noise
    mu_data = mu_data + noise;
    n_2_data = mean((Nsamp.^2),2) + noise.^2;
    var_data = n_2_data - ((mu_data).^2);
    
    for l=1:length(mu_data)
        if mu_data(l)<0
            mu_data(l)=0;
        end
        if var_data(l)<0
            var_data(l)=0;
        end
    end
    N0list = repmat(Ninit(i),length(mu_data),1);

    musimdatavec = vertcat(musimdatavec, mu_data(2:end));
    varsimdatavec = vertcat(varsimdatavec, var_data(2:end));
    timesimdatavec = vertcat(timesimdatavec, tsamp(2:end));
    N0vec = vertcat(N0vec,N0list(2:end));
  end
%% Plot the simulated data versus the model
figure;
p = [b,d,A];
V0=0;
for i =1:length(Ninit)
plot(tsamp,mu_fxnA(p, tsamp, Ninit(i)),'r-')
hold on
end
plot(timesimdatavec, musimdatavec,'b*')
title(['Allee model & mean of simulated data, N=', num2str(num_samps)])
xlabel('time(hours)')
ylabel('mean cell number')


figure;
for i = 1:length(Ninit)
plot(tsamp,V_fxnA(p, tsamp, Ninit(i), V0),'r-')
hold on
end
plot(timesimdatavec, varsimdatavec,'b*')
title('variance')
title(['Allee model & variance of simulated data, N=', num2str(num_samps)])
xlabel('time(hours)')
ylabel('variance in cell number')

%% Save the mean, variance, and, time vectors if running Allee model
save('../out/musimdatavecA.mat', 'musimdatavec')
save('../out/varsimdatavecA.mat', 'varsimdatavec')
save('../out/timesimdatavecA.mat', 'timesimdatavec')
save('../out/N0vecA.mat', 'N0vec')
save('../out/num_sampsA.mat', 'num_samps')
%% Save the mean, variance, and time vectors for b-d model
save('../out/musimdatavec.mat', 'musimdatavec')
save('../out/varsimdatavec.mat', 'varsimdatavec')
save('../out/timesimdatavec.mat', 'timesimdatavec')
save('../out/N0vec.mat', 'N0vec')
save('../out/num_samps.mat', 'num_samps')

