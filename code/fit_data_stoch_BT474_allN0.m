 %This script fits the Nmat for all initial cell numbers to both the
 %birth-death model and the Allee model

clear all; close all; clc
mudatavec = [];
vardatavec = [];
%% load this chunk to fit actual data 
S = load('../out/BTfit.mat');
BT= S.BT;

S = load('../out/BTsumfit.mat');
BTsum= S.BTsum;

% Need to make one big matrix of mu and vardata with corresponding time vecs in order
mudatavec = [];
vardatavec = [];
timevec = [];
N0vec = [];
numsamps = 0;

for j = 1:length(BTsum)
    
    mudatavec = vertcat(mudatavec, BTsum(j).mu_t(2:end));
    vardatavec = vertcat(vardatavec, BTsum(j).var_t(2:end));
    N0vec = vertcat(N0vec, repmat(BTsum(j).N0, length(BTsum(j).timevec)-1, 1));
    timevec = vertcat(timevec, BTsum(j).timevec(2:end)');
    numsamps = BTsum(j).num + numsamps;

end
%% Run this chunk to load in simulated data from Allee model
mudatavec = load('../out/musimdatavecA.mat');
mudatavec = struct2cell(mudatavec);
mudatavec = cell2mat(mudatavec);
vardatavec = load('../out/varsimdatavecA.mat');
vardatavec = struct2cell(vardatavec);
vardatavec = cell2mat(vardatavec);
timevec = load('../out/timesimdatavecA.mat');
timevec = struct2cell(timevec);
timevec = cell2mat(timevec);
numsamps = load('../out/num_sampsA.mat');
numsamps = struct2cell(numsamps);
numsamps = cell2mat(numsamps);
N0vec = load('../out/N0vecA.mat');
N0vec = struct2cell(N0vec);
N0vec = cell2mat(N0vec);


numsamps = numsamps*length(unique(N0vec));
%% Run this chunk to load in simulated data from b-d model
mudatavec = load('../out/musimdatavec.mat');
mudatavec = struct2cell(mudatavec);
mudatavec = cell2mat(mudatavec);
vardatavec = load('../out/varsimdatavec.mat');
vardatavec = struct2cell(vardatavec);
vardatavec = cell2mat(vardatavec);
timevec = load('../out/timesimdatavec.mat');
timevec = struct2cell(timevec);
timevec = cell2mat(timevec);
numsamps = load('../out/num_samps.mat');
numsamps = struct2cell(numsamps);
numsamps = cell2mat(numsamps);
N0vec = load('../out/N0vec.mat');
N0vec = struct2cell(N0vec);
N0vec = cell2mat(N0vec);

numsamps = numsamps*length(unique(N0vec));

%% Plot to check

figure;
plot(timevec, mudatavec,'*')
xlabel('time')
ylabel('mean cell number')
title('mean cell number all data')

figure;
plot(timevec, vardatavec,'*')
xlabel('time')
ylabel('variance in cell number')
title('variance in cell number all data')

%% DOWNSAMPLE!
% run this chunk to downsample the data from the time vec and the
% corresponding mudatavec and vardatavec

%% Now need to figure out how to generate concatenated model trajectories
% SINGLE EXPONENTIAL MODEL FIT
% rewrote the functions so that they take a time vecto and a N0 vector that
% correspond to the times sampled and the initial condition for many
% samples concatenated into a single vector

p=[0.0238, 0.005];

V0=0;
N = numsamps;
modelfun_mu= @(p)mu_fxn_all(p,timevec,N0vec);
modelfun_V=@(p)V_fxn_all(p, timevec,N0vec, V0);
var_in_mean =  @(p)(1/N).*(V_fxn_all(p, timevec, N0vec, V0)); % vertical
var_in_var = @(p)(1/N).*(V4_fxn_ODE_all(p, N0vec, V0, timevec)-(((N-3)./(N-1)).*(V_fxn_all(p, timevec, N0vec, V0).^2)));
%% Plot some example figures for arbitrary p against the data


figure;
plot(timevec,modelfun_mu(p),'r*')
hold on
plot(timevec, mudatavec,'b*')
title('mean')
legend('model', 'data')

figure;
plot(timevec,modelfun_V(p),'r*')
hold on
plot(timevec, vardatavec,'b*')
title('variance')
legend('model', 'data')

figure;
plot(timevec,var_in_mean(p),'r*')
title('variance in mean')

figure;
plot(timevec, var_in_var(p),'b*')
title('variance invariance')
%% PERFORM FITTING OF ALL DATA TO SINGLE EXPONENTIAL MODEL
pfxform = @(pval)[1 1].*log(pval); %'forward' parameter transform into Reals
pbxform = @(phat)[1 1].*exp(phat); %'backward' parameter transform into model space
yfxform = @(y)log(y); % 'forward' transform for data and model output
ybxform = @(yhat)exp(yhat); % 'inverse' transform for data and model output

%% Two observations to fit on, mu_data and var_data
% both use tsamp
pt1= length(mudatavec)-20;
pt2= length(mudatavec)-5;
% Initial guess
bguess = ((yfxform(mudatavec(pt2)))-(yfxform(mudatavec(pt1))))/(timevec(pt2)-timevec(pt1)) + 0.001; 
dguess = 0.001;
theta = [bguess,dguess];
%%
% Compute the negative log likelihood not using log transforms!

J = @(phat) (sum(((mudatavec-modelfun_mu(pbxform(phat))).^2)./(2.*sqrt(var_in_mean(pbxform(phat)))) + log(sqrt(var_in_mean(pbxform(phat)))) + 0.5*log(2*pi))+...
    sum(((vardatavec-modelfun_V(pbxform(phat))).^2)./(2.*sqrt(var_in_var(pbxform(phat)))) + log(sqrt(var_in_var(pbxform(phat)))) + 0.5*log(2*pi)));

objfun_J = @(phat)(J(phat));
lb = [0 0 ];
ub = [ Inf Inf ];
options = optimset('TolX', 1e-8);
[phatbest_J,fval,exitflag] = fminsearch(objfun_J, pfxform(theta), options);% need a better way to search parameter space,i.e. montecarlo
negLLguess= objfun_J(pfxform(theta))
negLLfit= objfun_J(phatbest_J)
params_best_J= pbxform(phatbest_J)

k = 2;
n = length(timevec);% or should this be the number of trajectories???
%n=numsamps;
AICbd = 2*J(pfxform(params_best_J)) + 2*k
BICbd = 2*J(pfxform(params_best_J)) +log(n)*k
%% Plot the model versus the data for the mean and variance
uniqN0 = unique(N0vec);
tplot = linspace(1,max(timevec),50);

figure;
for i =1:length(uniqN0)
plot(tplot,mu_fxn(params_best_J, tplot, uniqN0(i)),'r-')
hold on
end
plot(timevec, mudatavec,'b*')
title('B-d model fit to mean of the data')
xlabel('time(hours)')
ylabel('mean cell number')


figure;
for i = 1:length(uniqN0)
plot(tplot,V_fxn(params_best_J, tplot, uniqN0(i), V0),'r-')
hold on
end
plot(timevec, vardatavec,'b*')
title('variance')
ylim([0 2500])
title('B-d model fit to variance of the data')
xlabel('time(hours)')
ylabel('variance in cell number')
%% Find the profile likelihood of both parameters from fminsearch
profile = [];
params = params_best_J
fval = negLLfit
threshold = chi2inv(0.95,length(params))/2 + fval;
factor = 0.1; %percent range for profile to run across
numpoints = 10;


profindex = 1; % PROFILE b PARAMETER
    % Profile

    profrangeDown = linspace((params(profindex)), (params(profindex)*(1-factor)),numpoints)'; 
    profrangeUp = linspace((params(profindex)), (params(profindex)*(1+factor)),numpoints)';
    % split into up and down so we can use last fitted value as starting value for next run
    profrange = [profrangeDown;profrangeUp];
    profrange = sort(profrange)
    currfvals = [];
    currparams = [];
    currflags = [];
    paramstemp = [];
    guess = params(2);
    profile = [];
    for m = 1:length(profrange)
        [m] %track progress
        currb = profrange(m);
        guess = profrange(m)-params(1)-params(2);
        objfun_donly = @(d) objfun_J(pfxform([currb,d]));
        [paramstemp, fvaltemp, flagtemp] = fminsearch(objfun_donly, guess);
        currfvals = [currfvals; fvaltemp];
        currflags = [currflags; flagtemp];
        currparams = [currparams; [profrange(m),paramstemp]]; %storing the profiled value too, so the output parameter values are easy to run the model with
    end
    
    profile = horzcat(profrange, real(currfvals));
 profilelikeb=profile
ilo=find(profile(:,2)<threshold);


% 1st column is the parameter values that are "profiled"
% 2nd column is the negLL corresponding to each "profiled" parameter

figure;
plot(profile(:,1), profile(:,2))
hold on
plot(params(1),fval,'r*','LineWidth',2)
plot([profile(1,1) profile(end,1)],[threshold threshold],'r--')
xlabel('b-values')
ylabel('Cost Function Value')
title('Profile Likelihood of Parameter b')

CIonb=[profrange(ilo(1),1), profrange(ilo(end)+1),1]


profindex = 2; % PROFILE d PARAMETER
    % Profile
    factor = 1;
    profrangeDown = linspace((params(profindex)), (params(profindex)*(1-factor)),numpoints)'; 
    profrangeUp = linspace((params(profindex)), (params(profindex)*(1+factor)),numpoints)';
    % split into up and down so we can use last fitted value as starting value for next run
    profrange = [profrangeDown;profrangeUp];
    profrange = sort(profrange)
    currfvals = [];
    currparams = [];
    currflags = [];
    paramstemp = [];
    guess = params(1);
    profile = [];
    for m = 1:length(profrange)
        [m] %track progress
        currd = profrange(m);
        guess =params(1) - params(2)+profrange(m);
        objfun_bonly = @(b) objfun_J(pfxform([b,currd]));
        [paramstemp, fvaltemp, flagtemp] = fminsearch(objfun_bonly, guess);
        currfvals = [currfvals; fvaltemp];
        currflags = [currflags; flagtemp];
        currparams = [currparams; [paramstemp,profrange(m)]]; %storing the profiled value too, so the output parameter values are easy to run the model with
    end
    
    profile = horzcat(profrange, real(currfvals));
 
ilo=find(profile(:,2)<threshold);



% 1st column is the parameter values that are "profiled"
% 2nd column is the negLL corresponding to each "profiled" parameter
profliked= profile;

figure;
plot(profile(:,1), profile(:,2))
hold on
plot(params(2),fval,'r*','LineWidth',2)
plot([profile(1,1) profile(end,1)],[threshold threshold],'r--')
xlabel('d-values')
ylabel('Cost Function Value')
title('Profile Likelihood of Parameter d')

CIond=[profrange(ilo(1),1), profrange(ilo(end)+1),1]




%% MCMC Search algorithm
% 1. Uniformly sample domain and calculate negLL at each point
% 2. Take lowest negLL as initial guess 
% 3. Random walk to search for better neg LL

% set domain of b, d, and A
bvec = [ 0: 0.005: 0.05];
dvec = [ 0: 0.001: 0.01];

% make a 3D mesh
[B,D] = meshgrid(bvec, dvec); % 3D grid of parameters
Bflat = reshape(B,1,[]); % vector of bs
Dflat = reshape(D,1,[]); % vector of ds

% run loop through vector of parameters and calculate negLL

    for i = 1:length(Bflat)
        % set params
        pguess(i,1:2) = horzcat(Bflat(i),Dflat(i));
        negLL(i)= objfun_J(pfxform(pguess(i,1:2)));
        pguess(i,3) = negLL(i);   
    end

NEGLL = reshape(negLL, size(B));
[Jinit,imin] = min(pguess(:,3));
pinit = pguess(imin,1:2)

figure;
hold off;
surf(bvec,dvec,NEGLL(:,:));
hold on
plot3(pinit(1), pinit(2), Jinit, 'r*', 'LineWidth',8)
xlabel( 'b')
ylabel('d')
zlabel('negLL')
title('Initial parameter space search')
%%
J_curr = Jinit; % set current negLL value to inital value to initialize
count_accepted = 0; % start with no changes accepted

% Initialize a vector to store b and d params in
nruns = 10000;
store_params = zeros(nruns, 3); % store J, b,  d,for each run

% Outside the loop: initialize temperature and iterations
T0= 20;
k = 1:1:nruns;
T= T0*exp(-4*k./(nruns)); % sets up gradual cooling


% check that temperature annealing is working
figure;
hold off
set(gca,'LineWidth',1.5,'FontSize',12);
plot(k, T, 'LineWidth', 3)
xlabel('Markov step')
ylabel('Temperature')
title('Temperature Annealing')

bstep = 0.001; % step size for searching birth param
dstep = 0.001; % step size for searchign death param

% set initial guess as b1 and d1
b1=pinit(1);
d1=pinit(2);


store_acc_params = [];
store_acc_params(1,1) = Jinit;
store_acc_params(1,2) = b1;
store_acc_params(1,3) = d1;
theta1 = pinit;
    for k = 1:nruns
        r = rand;
        % set up to only change one parameter at a time, use rand to choose the
        % parameter to update
        if r<0.5
            % update b
            b2 = b1 + bstep*(2*rand-1);
            d2 = d1;
        end
        if r>0.5
            % update d
            d2 = d1 + dstep*(2*rand-1);
            b2 = b1;
        end


    % Constrain search region to domain
        if b2<0 
            b2=0;
        end
        if b2>.05
            b2 = .05;
        end
        if d2<0
            d2=0;
        end
        if d2>.01
            d2 = .01;
        end


        % find the neg LL of the new params
        theta2 = horzcat(b2,d2);
        J_new = objfun_J( pfxform(theta2));
        % store the  negLL and searched parameters
        store_params(k,1) = J_new;
        store_params(k,2) =b2;
        store_params(k,3) = d2;


        prob(k) = exp((J_curr-J_new)./T(k));
        % if Jcurr > Jnew, Jnew is better--the numerator in the exponent is positive &
        % prob>1--> change will always be accepted
        % if Jcurr < Jnew, numerator is negative, & prob <1 --> chance change
        % will be accepted

        if rand < prob(k) % if true, accept the change!
            b1 = b2;
            d1 = d2;

            theta1 = theta2;
            J_curr = J_new;
            count_accepted = count_accepted +1;
            % decrease search step size
            if r<0.5
                % update b step
                bstep = 0.999*bstep;
            end
            if r>0.5
                % update d step
                 dstep = 0.999*dstep;
            end

            params(1,1) = J_curr;
            params(1,2)= b1;
            params(1,3)= d1;
            store_acc_params= vertcat(store_acc_params, params);
        else
            % increase search step size
            if r<0.5
                % update b step
                bstep = 1.001*bstep;
            end
            if r>0.5
                % update d step
                 dstep = 1.001*dstep;
            end
        end
      end
%%
params_best_MC = store_acc_params(end,2:3)
[lowest_LL,ind] = min(store_acc_params(:,1));
lowestLLparams = store_acc_params(ind, 2:3)
negLLMC = lowest_LL;

likelihoodval = exp(-store_acc_params(:,1));
likelihoodvalall = exp(-store_params(:,1));
pointsize = 30;
figure;
scatter(store_acc_params(:,2), store_acc_params(:,3), pointsize, likelihoodval,'filled');
colorbar
xlabel('b')
ylabel('d')
title(' b versus d colored by likelihood for accepted parameters')


%%



% Get credible interval just from accepted parameter values

CrInt(:,1) = prctile(store_acc_params(:,2),[2.5, 97.5])
CrInt(:,2) = prctile(store_acc_params(:,3),[2.5, 97.5])


figure;
plot(timevec, mudatavec, 'b*')
hold on
plot(timevec, modelfun_mu(lowestLLparams),'r*')
legend('data','model')
title(['fit to mean, negLL = ', num2str(lowest_LL)])

figure;
plot(timevec,vardatavec,'b*')
hold on
plot(timevec, modelfun_V(lowestLLparams),'r*')
legend('data','model')
ylim([0 2500])
title(['fit to variance, negLL = ', num2str(lowest_LL)])
%% Make a table of model selection statistics for the b-d model
% BIC = ln(n)k-2ln(L)
% AIC = 2k-2ln(L)
k = 2;
n = N;% or should this be the number of trajectories???
%n=numsamps;
AICbd = 2*J(pfxform(lowestLLparams)) + 2*k
BICbd = 2*J(pfxform(lowestLLparams)) +log(n)*k

%% PERFORM FITTING OF ALL DATA TO B-D-A MODEL
p = [0.0238, 0.005, 2];
V0=0;
N = numsamps;
pfxform = @(pval)[1 1 1].*log(pval); %'forward' parameter transform into Reals
pbxform = @(phat)[1 1 1].*exp(phat); %'backward' parameter transform into model space
yfxform = @(y)log(y); % 'forward' transform for data and model output
ybxform = @(yhat)exp(yhat); % 'inverse' transform for data and model output

% Two observations to fit on, mu_data and var_data
% both use tsamp
% Redefine these
modelfun_mu= @(p)mu_fxnA_all(p,timevec, N0vec); % <n> for params p = b,d, & A
modelfun_V= @(p)V_fxnA_all(p,timevec, N0vec, V0); % variance as a function of parameters
var_in_mean =  @(p)(1/N).*(V_fxnA_all(p, timevec, N0vec, V0)); % variance in mean
var_in_var = @(p)(1/N).*(V4_fxnA_all(p,timevec, N0vec, V0)-(((N-3)./(N-1)).*(V_fxnA_all(p, timevec, N0vec, V0).^2)));

J = @(phat) (sum(((mudatavec-modelfun_mu(pbxform(phat))).^2)./(2.*sqrt(var_in_mean(pbxform(phat)))) + log(sqrt(var_in_mean(pbxform(phat)))) + 0.5*log(2*pi))+...
    sum(((vardatavec-modelfun_V(pbxform(phat))).^2)./(2.*sqrt(var_in_var(pbxform(phat)))) + log(sqrt(var_in_var(pbxform(phat)))) + 0.5*log(2*pi)));


%% Test moment functions with figures to compare to data
test2 = modelfun_mu(p);

figure;
plot(timevec, var_in_var(p), 'g*')

figure;
plot(timevec,modelfun_mu(p),'r*')
hold on
plot(timevec, mudatavec,'b*')
title('mean')
legend('model', 'data')
test = modelfun_V(p);
figure;
plot(timevec,modelfun_V(p),'r*')
hold on
plot(timevec, vardatavec,'b*')
title('variance')
legend('model', 'data')

figure;
plot(timevec,var_in_mean(p),'r*')
title('variance in mean')

figure;
plot(timevec, var_in_var(p),'b*')
title('variance invariance')

%% FITTNG FOR B-D-A model

Aguess = 1;
theta = [bguess,dguess,Aguess];

% test objective function

objfun_JA = @(phat)(J(phat));
testLL = objfun_JA(pfxform(theta))
%%
lb = [0 0 0];
ub = [ Inf Inf  Inf];
[phatbest_JA,fval,exitflag] = fminsearch(objfun_JA, pfxform(theta), options);% need a better way to search parameter space,i.e. montecarlo
negLLguess= objfun_JA(pfxform(theta));
negLLfitA= objfun_JA(phatbest_JA)
params_best_JA= pbxform(phatbest_JA)

%% Plot the model versus the data for the mean and variance

% Make functions that plot mu and variance with mode parameters and N0


figure;
for i =1:length(uniqN0)
plot(tplot,mu_fxnA(params_best_JA, tplot, uniqN0(i)),'r-')
hold on
end
plot(timevec, mudatavec,'b*')
title('Allee model fit to mean of the data')
xlabel('time(hours)')
ylabel('mean cell number')


figure;
for i = 1:length(uniqN0)
plot(tplot,V_fxnA(params_best_JA, tplot, uniqN0(i), V0),'r-')
hold on
end
plot(timevec, vardatavec,'b*')
title('variance')
ylim([0 2500])
title('Allee model fit to variance of the data')
xlabel('time(hours)')
ylabel('variance in cell number')
%%
k = 3;
n = numsamps;% or should this be the number of trajectories???
%n=numsamps;
AICbdA = 2*J(pfxform(params_best_JA)) + 2*k
BICbdA = 2*J(pfxform(params_best_JA)) +log(n)*k

%% MCMC Search algorithm
% 1. Uniformly sample domain and calculate negLL at each point
% 2. Take lowest negLL as initial guess 
% 3. Random walk to search for better neg LL

% set domain of b, d, and A
bvec = [ 0: 0.005: 0.05];
dvec = [ 0: 0.001: 0.01];
Avec = [ 0: 1:10];

% make a 3D mesh
[B,D,A] = meshgrid(bvec, dvec, Avec); % 3D grid of parameters
Bflat = reshape(B,1,[]); % vector of bs
Dflat = reshape(D,1,[]); % vector of ds
Aflat = reshape(A,1, []); % vector of As

% run loop through vector of parameters and calculate negLL

for j = 1:length(Bflat)
    % set params
    pguess(j,1:3) = horzcat(Bflat(j),Dflat(j),Aflat(j));
    negLL(j)= objfun_JA(pfxform(pguess(j,1:3)));
    pguess(j,4) = negLL(j);   
end

NEGLL = reshape(negLL, size(A));
[Jinit,imin] = min(pguess(:,4));
pinit = pguess(imin,1:3)

figure;
hold off;
surf(bvec,dvec,NEGLL(:,:,1));
hold on
plot3(pinit(1), pinit(2), Jinit, 'r*', 'LineWidth',8)
xlabel( 'b')
ylabel('d')
zlabel('negLL')
title('Initial parameter space search')
%%

J_curr = Jinit; % set current negLL value to inital value to initialize
count_accepted = 0; % start with no changes accepted

% Initialize a vector to store b and d params in
nruns = 1000;
store_params = zeros(nruns, 5); % store J, b,  d, sigu, and sig vfor each run

% Outside the loop: initialize temperature and iterations
T0= 50;
k = 1:1:nruns;
T= T0*exp(-1*k./(nruns)); % sets up gradual cooling


% check that temperature annealing is working
figure;
hold off
set(gca,'LineWidth',1.5,'FontSize',12);
plot(k, T, 'LineWidth', 3)
xlabel('Markov step')
ylabel('Temperature')
title('Temperature Annealing')

bstep = 0.001; % step size for searching birth param
dstep = 0.0002; % step size for searchign death param
Astep = 0.2;
% sigustep = 0.005;
% sigvstep = 0.005;
% set initial guess as b1 and d1
b1=pinit(1);
d1=pinit(2);
A1= pinit(3);
% sigu1= theta5(4);
% sigv1 = theta5(5);
store_acc_params = [];
store_acc_params(1,1) = Jinit;
store_acc_params(1,2) = b1;
store_acc_params(1,3) = d1;
store_acc_params(1,4) = A1;
%%
theta1 = pinit;
for k = 1:nruns
    r = rand;
    % set up to only change one parameter at a time, use rand to choose the
    % parameter to update
    if r<0.33
        % update b
        b2 = b1 + bstep*(2*rand-1);
        d2 = d1;
        A2 = A1;
    end
    if r>0.33 && r<0.67
        % update d
        d2 = d1 + dstep*(2*rand-1);
        b2 = b1;
        A2 = A1;
    end
    if r>0.67
        % update A
        A2 = A1 + Astep*(2*rand-1);
        b2 = b1;
        d2 = d1;
    end
    
%     b2 = b1 + bstep*(2*rand-1);
%     d2 = d1 + dstep*(2*rand-1);
%     A2 = A1 + Astep*(2*rand-1);

    
% Constrain search region to domain
    if b2<0 
        b2=0;
    end
    if b2>.05
        b2 = .05;
    end
    if d2<0
        d2=0;
    end
    if d2>.01
        d2 = .01;
    end
    if A2<0
        A2=0;
    end
    if A2>10
        A2 = 10;
    end
    
    % find the neg LL of the new params
    theta2 = horzcat(b2,d2,A2);
    J_new = objfun_JA( pfxform(theta2));
    % store the  negLL and searched parameters
    store_params(k,1) = J_new;
    store_params(k,2) =b2;
    store_params(k,3) = d2;
    store_params(k,4)=A2;
%     store_params(k,5) = sigu2;
%     store_params(k,6) = sigv2;
    
    prob(k) = exp((J_curr-J_new)./T(k));
    % if Jcurr > Jnew, Jnew is better--the numerator in the exponent is positive &
    % prob>1--> change will always be accepted
    % if Jcurr < Jnew, numerator is negative, & prob <1 --> chance change
    % will be accepted
    
    if rand < prob(k) % if true, accept the change!
        b1 = b2;
        d1 = d2;
        A1 = A2;
%         sigu1= sigu2;
%         sigv1 = sigv2;
        theta1 = theta2;
        J_curr = J_new;
        count_accepted = count_accepted +1;
        % decrease search step size
        if r<0.33
            % update b step
            bstep = 0.999*bstep;
        end
        if r>0.33 && r<0.67
            % update d step
             dstep = 0.999*dstep;
        end
        if r>0.67
            % update A
              Astep = 0.999*Astep;
        end

        params(1,1) = J_curr;
        params(1,2)= b1;
        params(1,3)= d1;
        params(1,4)= A1;
        store_acc_params= vertcat(store_acc_params, params);
    else
        % increase search step size
        if r<0.33
            % update b step
            bstep = 1.001*bstep;
        end
        if r>0.33 && r<0.67
            % update d step
             dstep = 1.001*dstep;
        end
        if r>0.67
            % update A
              Astep = 1.001*Astep;
        end

    end
       
end
%%
params_best_MCA = store_acc_params(end,2:4)
[lowest_LLA,ind] = min(store_acc_params(:,1));
lowestLLparamsA = store_acc_params(ind, 2:4)
k=3
AICbdA = 2*J(pfxform(params_best_MCA)) + 2*k
BICbdA = 2*J(pfxform(params_best_MCA)) +log(n)*k

CrIntA(:,1) = prctile(store_acc_params(:,2),[2.5, 97.5])
CrIntA(:,2) = prctile(store_acc_params(:,3),[2.5, 97.5])
CrIntA(:,3) = prctile(store_acc_params(:,4),[2.5,97.5])
%% Likelihood ratio test for model comparison
LLalt = -negLLfitA;
LLnull = -negLLfit;
logLR = 2*(LLalt-LLnull)

figure;
x=0:1:50;
plot(x, chi2pdf(x,1))
P= 1- chi2cdf(logLR,1);

%% Do this from memory from ppt for data
LLalt = -1.896e5;
LLnull = -1.898e5;
logLR = 2*(LLalt-LLnull)
P= 1- chi2cdf(logLR,1)
