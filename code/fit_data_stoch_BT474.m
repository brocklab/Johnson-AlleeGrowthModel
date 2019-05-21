% This script fits the Nmat from each initial cell number using the
% stochastic parameter estimation methods for the birth death and birth
% death Allee model

clear all; close all; clc

S = load('../out/BTfit.mat');
BT= S.BT;

S = load('../out/BTsumfit.mat');
BTsum= S.BTsum;

%% Perform Bayesian parameter estimation on birth-death model

for j =2% 1:length(BTsum)
    
% For now, start by assuming sigma = 1 and fit for b and d
    N0 = BTsum(j).N0;
    V0 = 0;
    mudata = [];
    vardata = [];
    t = [];
    negLLguesst = [];
    negLLfitt = [];
    mu_fitJst = [];
    V_fitJst = [];
    t=BTsum(j).timevec(2:end); % only want to fit data after initial condition
    N= length(t);
    mudata = BTsum(j).mu_t(2:end);
    vardata = BTsum(j).var_t(2:end);
    
% Want to fit parameters b, d, and sigma
% only want to fit data after initial condition
    N= BTsum(j).num;
% need to make N0 and V0 at time tsamp(2)=t0

% For now, start by assuming sigma = 1 and fit for b and d

pfxform = @(pval)[1 1].*log(pval); %'forward' parameter transform into Reals
pbxform = @(phat)[1 1].*exp(phat); %'backward' parameter transform into model space
yfxform = @(y)log(y); % 'forward' transform for data and model output
ybxform = @(yhat)exp(yhat); % 'inverse' transform for data and model output

% Two observations to fit on, mu_data and var_data
% both use tsamp
modelfun_mu = @(p)mu_fxn(p, t, N0); % <n> for params (p) = b & d, vertical
modelfun_V=@(p)V_fxn(p, t, N0, V0); % vertical
var_in_mean =  @(p)(1/N).*(V_fxn(p, t, N0, V0)); % vertical


var_in_var = @(p)(1/N).*(V4_fxn_ODE(p,N0, V0, t)-(((N-3)./(N-1)).*(V_fxn(p, t, N0, V0).^2)));

% Initial guess
bguess = (yfxform(mudata(end-5))-yfxform(mudata(end-20)))/(t(end-5)-t(end-20)) + 0.005; 
dguess = 0.005;
theta = [bguess,dguess];

% Compute the negative log likelihood not using log transforms!

J = @(phat) (sum(((mudata-modelfun_mu(pbxform(phat))).^2)./(2.*sqrt(var_in_mean(pbxform(phat)))) + log(sqrt(var_in_mean(pbxform(phat)))) + 0.5*log(2*pi))+...
    sum(((vardata-modelfun_V(pbxform(phat))).^2)./(2.*sqrt(var_in_var(pbxform(phat)))) + log(sqrt(var_in_var(pbxform(phat)))) + 0.5*log(2*pi)));
J_t = @(phat) ((((mudata-modelfun_mu(pbxform(phat))).^2)./(2.*sqrt(var_in_mean(pbxform(phat)))) + log(sqrt(var_in_mean(pbxform(phat)))) + 0.5*log(2*pi))+...
    (((vardata-modelfun_V(pbxform(phat))).^2)./(2.*sqrt(var_in_var(pbxform(phat)))) + log(sqrt(var_in_var(pbxform(phat)))) + 0.5*log(2*pi)));

objfun_J = @(phat)(J(phat));
lb = [0 0 ];
ub = [ Inf Inf ];
[phatbest_J,fval,exitflag] = fminsearchbnd(objfun_J, pfxform(theta), pfxform(lb), pfxform(ub));% need a better way to search parameter space,i.e. montecarlo
negLLguess= objfun_J(pfxform(theta));
negLLfit= objfun_J(phatbest_J)
params_best_J= pbxform(phatbest_J)
BTsum(j).negLLML= negLLfit;
BTsum(j).paramsML =params_best_J;

% look at it in time
negLLguesst=J_t(pfxform(theta));
negLLfitt = J_t(phatbest_J);

figure;
hold on
plot(t, negLLguesst, 'LineWidth',2)
plot(t, negLLfitt, 'LineWidth', 2)
xlabel ('time')
ylabel('contribution to negative LL')
title('Negative LL in time')
legend( 'initial guess', 'fit param')

mu_fitJst = mu_fxn(params_best_J,t, N0);
V_fitJst= V_fxn(params_best_J,t,N0, V0);
BTsum(j).mu_fitML = mu_fitJst;
BTsum(j).var_fitML = V_fitJst;


% MCMC Search algorithm
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

params_best_MC = store_acc_params(end,2:3)
[lowest_LL,ind] = min(store_acc_params(:,1));
lowestLLparams = store_acc_params(ind, 2:3)
negLLMC = lowest_LL;
negLLMCt = J_t(pfxform(params_best_MC));
BTsum(j).negLLMCt = negLLMCt;
BTsum(j).negLLMC = lowest_LL;
BTsum(j).params_MC= lowestLLparams;
BTsum(j).acc_params = store_acc_params;

mu_fitMC = mu_fxn(params_best_MC,t, N0);
V_fitMC= V_fxn(params_best_MC,t,N0, V0);
BTsum(j).mu_fitMC = mu_fitMC;
BTsum(j).var_fitMC = V_fitMC;


% Get credible interval just from accepted parameter values

BTsum(j).CrInt(:,1) = prctile(store_acc_params(:,2),[2.5, 97.5]);
BTsum(j).CrInt(:,2) = prctile(store_acc_params(:,3),[2.5, 97.5]);
end
%%
figure;
subplot(1,2,1)
plot(t, mudata, 'r*')
hold on
plot(t, mu_fitJst, 'b--', 'LineWidth',3)
plot(t, mu_fitMC, 'k-', 'LineWidth',2)
xlabel('time (hours)')
ylabel('<n>')
title('Bayes fit mean')
legend('data', 'MLfit', 'MCMC fit')
legend boxoff


subplot(1,2,2)
plot(t, vardata, 'g*')
hold on
plot(t, V_fitJst, 'b--', 'LineWidth',3)
plot(t, V_fitMC, 'k-', 'LineWidth',2)
xlabel('time (hours)')
ylabel('Variance')
title('Bayes fit variance')
legend('data', 'MLfit', 'MCMC fit')
legend boxoff

figure;
hold on
plot(t, negLLfitt, 'LineWidth',2)
plot(t, negLLMC, 'LineWidth', 2)
xlabel ('time')
ylabel('contribution to negative LL')
title('Negative LL in time')
legend( 'fminsearch', 'MC fit Params')

figure;
plot(0:1:count_accepted, store_acc_params(:,1))
xlabel('iterations')
ylabel('negLL')
title('negLL')

 figure;
 plot(0:1:count_accepted, store_acc_params(:,2), 'g.')
 hold on
 plot(0:1:count_accepted, store_acc_params(:,3), 'b.')
 legend( 'birth rates', 'death rates')
 xlabel('iterations')
 ylabel('rate')
 title('Accepted birth & death parameters')
 
 
 figure;
 plot(1:1:nruns, store_params(:,2), 'g.')
 hold on
 plot(1:1:nruns, store_params(:,3), 'b.')
 legend( 'birth rates', 'death rates')
 xlabel('iterations')
 ylabel('rate')
 title('Explored birth, death, and A parameters')
 
 figure;
 plot(1:1:nruns, prob, 'r.')
 ylim([0 1])

figure;
subplot(1,2,1)
hist(store_acc_params(:,2))
ylabel('frequency')
xlabel('b')
subplot(1,2,2)
hist(store_acc_params(:,3))
ylabel('frequency')
xlabel('d')


% Plots of two parameters colored by likelihood
likelihoodval = exp(-store_acc_params(:,1));
likelihoodvalall = exp(-store_params(:,1));
pointsize = 30;
figure;
scatter(store_acc_params(:,2), store_acc_params(:,3), pointsize, likelihoodval,'filled');
colorbar
xlabel('b')
ylabel('d')
title(' b versus d colored by likelihood for accepted parameters')
%%  Compute the profile likelihoods for each parameter for each initial cell number
% Use these to determine parameter certainty and find the 95 % CI on the
% b and d parameter estimates

% basic idea: hold 1 parameter constant at a certain value, allow search
% algorithm to find the best possible likelihood with that parameter
% constant, and report likelihood value.

% For first pass, we will just use fminsearch for each step through
% parameters

for j = 2%1:length(BTsum)  
profile = [];
params = BTsum(j).paramsML;
fval = BTsum(j).negLLML;
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
 
ilo=find(profile(:,2)<threshold);

BTsum(j).proflikeb = profile;

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
BTsum(j).CIonb=CIonb;

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
BTsum(j).profliked= profile;

figure;
plot(profile(:,1), profile(:,2))
hold on
plot(params(2),fval,'r*','LineWidth',2)
plot([profile(1,1) profile(end,1)],[threshold threshold],'r--')
xlabel('d-values')
ylabel('Cost Function Value')
title('Profile Likelihood of Parameter d')

CIond=[profrange(ilo(1),1), profrange(ilo(end)+1),1]
BTsum(j).CIond=CIond;

end
%%
paramnames = {'b','d'};
for i = 1:length(paramests)
figure(10+i)
    set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial')
    hold on
    plot(profiles(:,1,i),profiles(:,2,i),'k','LineWidth',2)
    plot(paramests(i),fval,'r*','LineWidth',2)
    plot([profiles(1,1,i) profiles(end,1,i)],[threshold threshold],'r--')
    xlabel(paramnames{i})
    ylabel('Cost Function Value')
%plot parameter relationships
    figure(20+i)
    set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial')
    hold on
    plot(profiles(:,1,i),profiles(:,4:end,i),'LineWidth',2)
    plot(paramests(i),paramests,'r*')
    xlabel(paramnames{i})
    ylabel('Estimated Parameter Value')
    legend(paramnames)
end
%% Perform Bayesian parameter estimation on Allee Model
for j = 3% :length(BTsum)
% start by clearing all variables since size will change each loop
mudata = [];
vardata = [];
t = [];
negLLguesst = [];
negLLfitt = [];
mu_fitJst = [];
V_fitJst = [];
    
    
% Want to fit parameters b, d, A, and  both sigmas
t=BTsum(j).timevec(2:end); % only want to fit data after initial condition
N= length(t);
mudata = BTsum(j).mu_t(2:end);
vardata = BTsum(j).var_t(2:end);
% For now, start by assuming sigma = 1 and fit for b and d
N0 = BTsum(j).N0;
V0 = 0;
pfxform = @(pval)[1 1 1].*log(pval); %'forward' parameter transform into Reals
pbxform = @(phat)[1 1 1].*exp(phat); %'backward' parameter transform into model space
yfxform = @(y)log(y); % 'forward' transform for data and model output
ybxform = @(yhat)exp(yhat); % 'inverse' transform for data and model output

% Two observations to fit on, mu_data and var_data
% both use tsamp
modelfun_mu= @(p)mu_fxnA(p,t, N0); % <n> for params p = b,d, & A
modelfun_V= @(p)V_fxnA(p,t, N0, V0); % variance as a function of parameters
var_in_mean =  @(p)(1/N).*(V_fxnA(p, t, N0, V0)); % variance in mean
var_in_var = @(p)(1/N).*(V4_fxnA(p,t, N0, V0)-(((N-3)./(N-1)).*(V_fxnA(p, t, N0, V0).^2)));
% Initial guess
bguess = (yfxform(mudata(end-5))-yfxform(mudata(end-20)))/(t(end-5)-t(end-20)) + 0.005; 
dguess = 0.005;
Aguess = 1;

theta = [bguess,dguess, Aguess];

% Redo log likelihood function not using log transforms!
Jnew = @(phat) (sum(((mudata-modelfun_mu(pbxform(phat))).^2)./(2.*sqrt(var_in_mean(pbxform(phat)))) + log(sqrt(var_in_mean(pbxform(phat)))) + 0.5*log(2*pi))+...
    sum(((vardata-modelfun_V(pbxform(phat))).^2)./(2.*sqrt(var_in_var(pbxform(phat)))) + log(sqrt(var_in_var(pbxform(phat)))) + 0.5*log(2*pi)));
J_t = @(phat) ((((mudata-modelfun_mu(pbxform(phat))).^2)./(2.*sqrt(var_in_mean(pbxform(phat)))) + log(sqrt(var_in_mean(pbxform(phat)))) + 0.5*log(2*pi))+...
    (((vardata-modelfun_V(pbxform(phat))).^2)./(2.*sqrt(var_in_var(pbxform(phat)))) + log(sqrt(var_in_var(pbxform(phat)))) + 0.5*log(2*pi)));

objfun_J = @(phat)(Jnew(phat));

% need to use MCMC search algorithm to find lowest LL in uniform region of
% parameter space

bvec = [ 0: 0.005: 0.05];
dvec = [ 0: 0.001: 0.01];
Avec = [ 0: 1:10];

% make a 3D mesh
[B,D,A] = meshgrid(bvec, dvec, Avec); % 3D grid of parameters
Bflat = reshape(B,1,[]); % vector of bs
Dflat = reshape(D,1,[]); % vector of ds
Aflat = reshape(A,1, []); % vector of As

% run loop through vector of parameters and calculate negLL

    for i = 1:length(Bflat)
        % set params
        pguess(i,1:3) = horzcat(Bflat(i),Dflat(i),Aflat(i));
        negLL(i)= objfun_J(pfxform(pguess(i,1:3)));
        pguess(i,4) = negLL(i);   
    end

NEGLL = reshape(negLL, size(A));
[Jinit,imin] = min(pguess(:,4));
pinit = pguess(imin,1:3)


figure;
hold off;
surf(bvec,dvec,NEGLL(:,:,pinit(3)+1));
hold on
plot3(pinit(1), pinit(2), Jinit, 'r*', 'LineWidth',8)
xlabel( 'b')
ylabel('d')
zlabel('negLL')
title('Initial parameter space search')


J_curr = Jinit; % set current negLL value to inital value to initialize
count_accepted = 0; % start with no changes accepted

% Initialize a vector to store b and d params in
nruns = 10000;
store_params = zeros(nruns, 4); % store J, b,  d, A each run

% Outside the loop: initialize temperature and iterations
T0= 20;
k = 1:1:nruns;
T= T0*exp(-4*k./(nruns)); % sets up gradual cooling
bstep = 0.001; % step size for searching birth param
dstep = 0.0002; % step size for searchign death param
Astep = 0.2;
b1=pinit(1);
d1=pinit(2);
A1= pinit(3);
store_acc_params = [];
store_acc_params(1,1) = Jinit;
store_acc_params(1,2) = b1;
store_acc_params(1,3) = d1;
store_acc_params(1,4) = A1;
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
        J_new = objfun_J( pfxform(theta2));
        % store the  negLL and searched parameters
        store_params(k,1) = J_new;
        store_params(k,2) =b2;
        store_params(k,3) = d2;
        store_params(k,4)=A2;

        prob(k) = exp((J_curr-J_new)./T(k));
        % if Jcurr > Jnew, Jnew is better--the numerator in the exponent is positive &
        % prob>1--> change will always be accepted
        % if Jcurr < Jnew, numerator is negative, & prob <1 --> chance change
        % will be accepted

        if rand < prob(k) % if true, accept the change!
            b1 = b2;
            d1 = d2;
            A1 = A2;
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
    params_best_J = store_acc_params(end,2:4);
    lowestLLparams = store_acc_params(end, 1);


% look at model fit
mu_fitJst = mu_fxnA(params_best_J,t, N0);
V_fitJst= V_fxnA(params_best_J,t,N0, V0);

% same some things
BTsum(j).paramsA = params_best_J;
BTsum(j).mufit = mu_fitJst;
BTsum(j).varfit = V_fitJst;
end

%% Plot fit versus data

for j = 3%:length(BTsum)
    subplot(1,2,1)
    plot(BTsum(j).timevec(2:end), BTsum(j).mufit, 'r-', 'LineWidth',2)
    hold on
    plot(BTsum(j).timevec, BTsum(j).mu_t,'b*-')
    xlabel('time (hours)')
    ylabel('mean cell number')
    title('fit of mean cell number Allee model')
    
    subplot(1,2,2)
    plot(BTsum(j).timevec(2:end), BTsum(j).varfit, 'g-', 'LineWidth',2)
    hold on
    plot(BTsum(j).timevec, BTsum(j).var_t,'y*-')
    xlabel('time (hours)')
    ylabel('variance in cell number')
    title('fit of variance in cell number Allee model')
end