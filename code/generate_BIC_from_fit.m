% Allee Effect model selection & validation
% In this script, we will generate simulated data from a single model
% structure, and attempt to calibrate the data to three different models.
% We will demonstrate that the BIC value is able to properly identify the
% correct model structure, and demonstrate the ability to provide reliable
% parameter estimates

% One thing we are interested in investigating here is what resolution of
% sampling allows us to identify the correct model structure
close all; clear all; clc

% In order to reuse functions, use model code system
%MODEL CODES:
%1: simple b-d model
%2: strong Allee on birth
%3: strong Allee on death
%4: strong Allee on birth and death equally
%5: weak/strong Allee on birth
%6: weak/strong Allee on death
%7: weak/strong Allee on birth & death equally

%% Step 1: Generate simulated data for a few N0 and high resolution time points
% We will first do this from the Allee model, where we assume Allee is on
% both birth and death
tints = [ 4 8 12 24 48];
for k = 1:5
% Set up simulation
Nsim = 8000;
tsamp = [ 0:tints(k):336];
b = 0.0092;
d = 0.001;
A = 2;
tau = 3;
Ninit = [ 3; 5; 10]; 
paramsA= [b,d,A];
paramsw = [b,d,-A, tau];
params=[b,d];
N0vec = [];
for i =1:length(Ninit)
    N0vec = vertcat(N0vec, repmat(Ninit(i), length(tsamp)-1, 1));
end

% Run strong Allee model, Allee effect on birth
[ NsampAwb,NstatAwb, CstatAwb, musimdatavecAwb, varsimdatavecAwb, timesimdatavec] = run_strAllbmodel(paramsA, tsamp, Ninit, Nsim);
%[ NsampAbd,NstatAbd, CstatAbd, musimdatavecAbd, varsimdatavecAbd, timesimdatavec] = run_strAllbdmodel(paramsA, tsamp, Ninit, Nsim);
%[ Nsampbd,Nstatbd, Cstatbd, musimdatavecbd, varsimdatavecbd, timesimdatavec] = run_bdmodel(params, tsamp, Ninit, Nsim);
sim_mu_clean = musimdatavecAwb;
sim_var_clean = varsimdatavecAwb;
% take mean and variance of data and add noise
% Note, we expect measurement noise will not be constant, but will be
% proportional to the correspnding measurement, so we simulate accordingly
eta = 0.0;

noise_mu = (eta.*(sim_mu_clean)).*(randn(length(sim_mu_clean),1));
sim_mean = sim_mu_clean + noise_mu;

noise_var = (eta.*(sim_var_clean)).*(randn(length(sim_var_clean),1));
sim_var = sim_var_clean + noise_var;




% Perform fitting: first fit simulated mean and variance data to b-d model
% Want to calculate BIC and AIC value from the fitting
% Need to make one big matrix of mu and vardata with corresponding time vecs in order
mudatavec = sim_mean;
vardatavec = sim_var;
timevec = timesimdatavec;
numsamps = Nsim;

N = numsamps;

% For each fitting, we want to write functions in terms of the functions
% that run the model we are fitting to forward. The gen_model_mu runs the
% system of ODEs for the parameters and generates an expected model mean
% and variance for the input p parameters
% information

modelcode = 1;
% Initial guess
bguess = (log(mudatavec(end-2))-log(mudatavec(end-10)))/(timevec(end-2)-timevec(end-10)) + 0.002; 
dguess = 0.002;
theta = [bguess,dguess];

[pbestbd, BICbd, AICbd,negLLfit] = ML_fit(theta, tsamp, Ninit, N,mudatavec, vardatavec, modelcode)

% need to define th inline functions that run forward the mean and the
% variance
modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);


% Fit simulated mean and variance data to the strong Allee model on birth
modelcode = 2;
bguess = 0.008;
%bguess = (log(mudatavec(end-2))-log(mudatavec(end-10)))/(timevec(end-2)-timevec(end-10)) + 0.002; 
dguess = 0.002;
Aguess = 2.5;
theta = [bguess,dguess, Aguess];
[pbeststrAb, BICstrAb, AICstrAb,negLLfitstrAb] = ML_fit(theta, tsamp, Ninit, N,mudatavec, vardatavec, modelcode);

modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);


% Now fit simulated mean and variance data to the strong Allee model on death
modelcode = 3;
%*** Look back at fitting***
% Initial guess
%bguess = (log(mudatavec(end-2))-log(mudatavec(end-10)))/(timevec(end-2)-timevec(end-10)) + 0.002; 
dguess = 0.001;
Aguess = 1.5;
theta = [bguess,dguess, Aguess];

[pbeststrAd, BICstrAd, AICstrAd,negLLfitstrAd] = ML_fit(theta, tsamp, Ninit, N,mudatavec, vardatavec, modelcode);

% need to define th inline functions that run forward the mean and the
% variance
modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);


% Now fit simulated mean and variance data to the strong Allee model on birth & death
modelcode = 4;
%*** Look back at fitting***
% Initial guess
%bguess = (log(mudatavec(end-2))-log(mudatavec(end-10)))/(timevec(end-2)-timevec(end-10)) + 0.002; 
dguess = 0.001;
Aguess = 3;
theta = [bguess,dguess, Aguess];

[pbeststrAbd, BICstrAbd, AICstrAbd,negLLfitstrAbd] = ML_fit(theta, tsamp, Ninit, N,mudatavec, vardatavec, modelcode);

% need to define th inline functions that run forward the mean and the
% variance
modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);

% Fit data to weak Allee model on birth
modelcode = 5;
% Initial guess
%bguess = (log(mudatavec(end-2))-log(mudatavec(end-10)))/(timevec(end-2)-timevec(end-10)) + 0.002; 
dguess = 0.001;
Aguess = -2;
tauguess = 3;
theta = paramsw;
theta = [paramsA, 0]

[pbestAwb, BICAwb, AICAwb,negLLfitAwb] = ML_fit(theta, tsamp, Ninit, N,mudatavec, vardatavec, modelcode);

% need to define th inline functions that run forward the mean and the
% variance

% Fit data to weak Allee model on death
modelcode = 6;
% Initial guess
%bguess = (log(mudatavec(end-2))-log(mudatavec(end-10)))/(timevec(end-2)-timevec(end-10)) + 0.002; 
dguess = 0.001;
Aguess = -2;
tauguess = 3;
theta = paramsw;


[pbestAwd, BICAwd, AICAwd,negLLfitAwd] = ML_fit(theta, tsamp, Ninit, N,mudatavec, vardatavec, modelcode);

% need to define th inline functions that run forward the mean and the
% variance

modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);


% Fit data to weak Allee model on birth & death
modelcode = 7;
% Initial guess
%bguess = (log(mudatavec(end-2))-log(mudatavec(end-10)))/(timevec(end-2)-timevec(end-10)) + 0.002; 
dguess = 0.001;
Aguess = -2;
tauguess = 3;
theta = paramsw;


[pbestAwbd, BICAwbd, AICAwbd,negLLfitAwbd] = ML_fit(theta, tsamp, Ninit, N,mudatavec, vardatavec, modelcode);

% need to define th inline functions that run forward the mean and the
% variance

modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);


% Make vector of BIC values, negLL values, and parameter numbers
num_params = [ 2, 3, 3, 3, 4, 4, 4];
modellist = [1, 2, 3, 4, 5, 6, 7];
BICvals = [ BICbd, BICstrAb, BICstrAd, BICstrAbd, BICAwb, BICAwd, BICAwbd];
negLLvals = [negLLfit, negLLfitstrAb, negLLfitstrAd, negLLfitstrAbd, negLLfitAwb, negLLfitAwd, negLLfitAwbd];

figure;
subplot(3,1,1)
semilogy(modellist, BICvals, 'ko', 'LineWidth', 4)
xlabel('Model')
ylabel('BIC value')
xlim([0 8])
ylim([ min(BICvals)-10 max(BICvals)+5])
title(['BIC sampling every ', num2str(tints(k)),' hours'])


subplot(3,1,2)
semilogy(modellist, negLLvals, 'ro', 'LineWidth', 4)
xlabel('Model')
ylabel(' best negLL')
xlim([0 8])
ylim([ min(negLLvals)-5, max(negLLvals)+5])

subplot(3,1,3)
plot(modellist, num_params, 'bo', 'LineWidth', 4)
xlabel('Model')
ylabel(' number of parameters')
xlim([0 8])
ylim([ 1 5 ])


negLLvalsall(:,k)= negLLvals;
BICvalsall(:,k) = BICvals;

end
save(['../out/negLLvalsbd.mat'], 'negLLvals')
save('../out/BICval4bd.mat', 'BICvals')
%%
for k =1:5
figure;
subplot(3,1,1)
semilogy(modellist, BICvals, 'ko', 'LineWidth', 4)
xlabel('Model')
ylabel('BIC value')
xlim([0 8])
ylim([ min(BICvals)-10 max(BICvals)+5])
title(['BIC sampling every ', num2str(tints(k)),' hours'])


subplot(3,1,2)
semilogy(modellist, negLLvals, 'ro', 'LineWidth', 4)
xlabel('Model')
ylabel(' best negLL')
xlim([0 8])
ylim([ min(negLLvals)-5, max(negLLvals)+5])

subplot(3,1,3)
plot(modellist, num_params, 'bo', 'LineWidth', 4)
xlabel('Model')
ylabel(' number of parameters')
xlim([0 8])
ylim([ 1 5 ])


negLLvalsall(:,k)= negLLvals;
BICvalsall(:,k) = BICvals;

end

%% Load in all negLLvals and BIC vals
negLLvals = [];
BICvals = [];
negLLvals(:,1) = cell2mat(struct2cell(load('../out/negLLvals4bd.mat')));

negLLvals(:,2) = cell2mat(struct2cell(load('../out/negLLvals12bd.mat')));
negLLvals(:,3) = cell2mat(struct2cell(load('../out/negLLvals24bd.mat')));
negLLvals(:,4) = cell2mat(struct2cell(load('../out/negLLvals48bd.mat')));

BICvals(:,1) = cell2mat(struct2cell(load('../out/BICvals4bd.mat')));
BICvals(:,2) = cell2mat(struct2cell(load('../out/BICvals12bd.mat')));
BICvals(:,3) = cell2mat(struct2cell(load('../out/BICvals24bd.mat')));
BICvals(:,4) = cell2mat(struct2cell(load('../out/BICvals48bd.mat')));
%%
timeres = [ 4 12 24 48];
figure;
hold on
for i = 2%1:4
subplot(3,1,1)
hold on
semilogy(modellist, BICvals(:,i), 'o-', 'LineWidth', 4)
%text(4, BICvals(end, i),['Sampling time = every', num2str(timeres(i)),'hours'])
xlabel('Model')
ylabel('BIC value')
xlim([.5 4.55])
ylim([ 200 1.5e4])
legend('4 hours', '12 hours', '24 hours', '48 hours')
legend boxoff

subplot(3,1,2)
hold on
semilogy(modellist, negLLvals(:,i), 'o-', 'LineWidth', 4)
%text(4, negLLvals(end, i),['Sampling time = every', num2str(timeres(i)),'hours'])
xlabel('Model')
ylabel(' best negLL')
xlim([.4 4.55])
ylim([ 200 1e4])
legend('4 hours', '12 hours', '24 hours', '48 hours')
legend boxoff
end
subplot(3,1,3)
plot(modellist, num_params, 'bo', 'LineWidth', 4)
xlabel('Model')
ylabel(' number of parameters')
xlim([0 5])
ylim([ 1 4])

%% Now fit simulated mean and variance data to the weak Allee model Allee on birth
modelcode = 5;
% test out the functions

N=Nsim;
% Initial guess
bguess = 0.009;
dguess = 0.001;
Aguess = -1;
tauguess = 2;
theta = paramsw;
lb = [0 0 -10 0];
ub = [ 1 1 10 10];

pfxform = @(pval)[1 1 1 1].*(pval-lb)./(ub-lb); %'forward' parameter transform (normalized)
pbxform = @(phat)[1 1 1 1].*((ub-lb).*(phat))+lb;
thetanorm = pfxform(theta)
thetaback = pbxform(thetanorm)

%%

[pbestwkAb, BICwkAb, AICwkAb,negLLfitwkAb] = ML_fit(theta, tsamp, Ninit, N,mudatavec, vardatavec, modelcode);
%%
% need to define th inline functions that run forward the mean and the
% variance
modelcode =1
modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);
modelfun_V4= @(p)gen_model_v4(p, tsamp, Ninit, modelcode);
var_in_mean =  @(p)(1/N).*(modelfun_V(p)); % vertical
var_in_var = @(p)(1/N).*(modelfun_V4(p)-(((N-3)./(N-1)).*(modelfun_V(p).^2)));


figure;
plot(timesimdatavec, var_in_var(theta), 'k.')
title('var in var for weak model')

figure;
plot(timesimdatavec, var_in_mean(theta), 'k.')
title('var in mean for weak model')
%%


figure;
hold on
plot(timesimdatavec, sim_mean, 'r*')
plot(timesimdatavec, modelfun_mu(theta), 'k.')
for i = 1:length(Ninit)
    text(tsamp(end-10), CstatAb(end-10,1,i),['N_{0}=', num2str(Ninit(i))])
end
xlabel('time(hours)')
ylabel('mean cell number')
title('Simulated data versus best fit to strong Allee model on b & d')
legend('simulated data mean', 'best fit mean')

figure;
hold on
plot(timesimdatavec, sim_var, 'g*')
plot(timesimdatavec, modelfun_V(theta), 'k.')
for i = 1:length(Ninit)
    text(tsamp(end-10), CstatAb(end-10,3,i),['N_{0}=', num2str(Ninit(i))])
end
xlabel('time(hours)')
ylabel('variance in cell number')
title('Simulated data versus best fit to strong Allee model on b & d')
legend('simulated data mean', 'best fit mean')
