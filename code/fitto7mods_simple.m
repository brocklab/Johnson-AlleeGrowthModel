%% Load in selected mudatavec
% Want to calculate BIC and AIC value from the fitting
% Need to make one big matrix of mu and vardata with corresponding time vecs in order
% Load in raw data, not down sampled
clear all; close all; clc

mudatavec = [];
vardatavec = [];
timevec = [];
mudatavecall = [];
vardatavecall = [];
timevecall = [];
N0vec = [];
numsamps = 0;

S = load('../out/BTsumfit.mat');
BTsum= S.BTsum;
N0list = [ 2 4 10];


for j = 1:length(N0list)
    for i = 1:length(BTsum)
    if BTsum(i).N0==N0list(j)
        mudata = BTsum(i).mu_t(2:end);
        ind = 1:1:length(mudata);
        mudatasmall= mudata(ind);
    mudatavec = vertcat(mudatavec, mudatasmall);
        vardata = BTsum(i).var_t(2:end);
        vardatasmall = vardata(ind);
    vardatavec = vertcat(vardatavec, vardatasmall);
        time = BTsum(i).timevec(2:end)';
        timesmall = time(ind);
    timevec = vertcat(timevec, timesmall);
  
   
    vardatavecall = vertcat(vardatavecall, BTsum(i).var_t(2:end));
    mudatavecall = vertcat(mudatavecall, BTsum(i).mu_t(2:end));
    timevecall = vertcat(timevecall, BTsum(i).timevec(2:end)');
    numsamps = BTsum(i).num + numsamps;

end
    end
end

%% confirm it's same
figure;
for j = 1:length(N0list)
plot(timevecall, mudatavecall, 'r.')
hold on
plot(timevec, mudatavec, 'b*')
end
xlabel('time (hours)')
ylabel('cell number')
title('Down sampled data to every 24 hours')

%% Perform fitting: fit all data to a single b & d
% For each fitting, we want to write functions in terms of the functions
% that run the model we are fitting to forward. The gen_model_mu runs the
% system of ODEs for the parameters and generates an expected model mean
% and variance for the input p parameters
% information

modelcode = 1;
tsamp = 0:4:332; % set sample time (if data is down sampled change accordingly
tsampfit = tsamp(2:end);
tsampnew = tsampfit(ind);
tsamp = [0 tsampnew];
% Initial guess
bguess = (log(mudatavec(end-2))-log(mudatavec(end-4)))/(timevec(end-2)-timevec(end-4)) + 0.002; 
dguess = 0.002;
theta = [bguess,dguess];
Ninit = N0list;
N= numsamps;
[pbestbd, BICbd, AICbd,negLLfitbd] = ML_fit(theta, tsamp, Ninit, N,mudatavec, vardatavec, modelcode)

% need to define th inline functions that run forward the mean and the
% variance
modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);



%% Perform fitting: fit all data to strong Allee model on birth
% For each fitting, we want to write functions in terms of the functions
% that run the model we are fitting to forward. The gen_model_mu runs the
% system of ODEs for the parameters and generates an expected model mean
% and variance for the input p parameters
% information

modelcode = 2;

% Initial guess
bguess = (log(mudatavec(end-2))-log(mudatavec(end-4)))/(timevec(end-2)-timevec(end-4)) + 0.002; 
dguess = 0.0002;
Aguess = 1;
theta = [bguess,dguess, Aguess];
Ninit = N0list;
N= numsamps;
[pbeststrAb, BICstrAb, AICstrAb,negLLfitstrAb] = ML_fit(theta, tsamp, Ninit, N,mudatavec, vardatavec, modelcode);

% need to define th inline functions that run forward the mean and the
% variance
modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);

%% Plot the data compared to best fit model for strong Allee on birth
modelcode = 2;

figure;
subplot(1,2,1)
hold on
plot(timevec, mudatavec, 'r*')
for j = 1:length(N0list)
    modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit(j), modelcode);
plot(tsamp(2:end), modelfun_mu(pbeststrAb), 'b-', 'LineWidth', 2) % test model expected
end
xlabel('time (hours)')
ylabel('<n> time')
legend('data', 'best fit strong Allee')
legend boxoff
title(['Mean  in data fit to strong Allee on birth'])
xlim([ 0 332])
subplot(1,2,2)
hold on
plot(timevec, vardatavec, 'g*')
for j = 1:length(N0list)
    modelfun_V =@(p)gen_model_var(p,tsamp, Ninit(j), modelcode);
plot(tsamp(2:end), modelfun_V(pbeststrAb), 'b-', 'LineWidth', 2) % test model expected
end
xlabel('time (hours)')
ylabel(' Variance in time')
legend('data', 'best fit strong Allee')
legend boxoff
title(['Variance in data fit to strong Allee on birth'])
xlim([ 0 332])


%% Perform fitting: fit all data to strong Allee model on death
% For each fitting, we want to write functions in terms of the functions
% that run the model we are fitting to forward. The gen_model_mu runs the
% system of ODEs for the parameters and generates an expected model mean
% and variance for the input p parameters
% information

modelcode = 3;

% Initial guess
bguess = (log(mudatavec(end-2))-log(mudatavec(end-4)))/(timevec(end-2)-timevec(end-4)) + 0.002; 
dguess = 0.002;
Aguess = 1;
theta = [bguess,dguess, Aguess];
Ninit = N0list;
N= numsamps;
[pbeststrAd, BICstrAd, AICstrAd,negLLfitstrAd] = ML_fit(theta, tsamp, Ninit, N,mudatavec, vardatavec, modelcode);

% need to define th inline functions that run forward the mean and the
% variance
modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);
%%

% Perform fitting: fit all data to strong Allee model on  birth & death
% For each fitting, we want to write functions in terms of the functions
% that run the model we are fitting to forward. The gen_model_mu runs the
% system of ODEs for the parameters and generates an expected model mean
% and variance for the input p parameters
% information

modelcode = 4;
% Initial guess
bguess = (log(mudatavec(end-2))-log(mudatavec(end-4)))/(timevec(end-2)-timevec(end-4)) + 0.002; 
dguess = 0.002;
Aguess = 1;
theta = [bguess,dguess, Aguess];
Ninit = N0list;
N= numsamps;
[pbeststrAbd, BICstrAbd, AICstrAbd,negLLfitstrAbd] = ML_fit(theta, tsamp, Ninit, N,mudatavec, vardatavec, modelcode)

% need to define th inline functions that run forward the mean and the
% variance
modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);



%% Perform fitting: fit all data to weak/strong Allee on birth
% For each fitting, we want to write functions in terms of the functions
% that run the model we are fitting to forward. The gen_model_mu runs the
% system of ODEs for the parameters and generates an expected model mean
% and variance for the input p parameters
% information

modelcode = 5;

% Initial guess
bguess = (log(mudatavec(end-2))-log(mudatavec(end-4)))/(timevec(end-2)-timevec(end-4)) + 0.002; 
dguess = 0.0002;
Aguess = -2;
tauguess = 3;
theta = [bguess,dguess, Aguess, tauguess];
Ninit = N0list;
N= numsamps;
[pbestwkAb, BICwkAb, AICwkAb,negLLfitwkAb] = ML_fit(theta, tsamp, Ninit, N,mudatavec, vardatavec, modelcode)

% need to define th inline functions that run forward the mean and the
% variance
modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);

%% Plot the data compared to best fit model
figure;
subplot(1,2,1)
hold on
plot(timevec, mudatavec, 'r*')
for j = 1:length(N0list)
    modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit(j), modelcode);
plot(tsamp(2:end), modelfun_mu(pbestwkAb), 'b-', 'LineWidth', 2) % test model expected
end
xlabel('time (hours)', 'FontSize',16)
ylabel('Mean cell number (<N(t)>)', 'FontSize', 16)
legend('mean in data', 'best fit weak Allee', 'FontSize',12, 'Location', 'NorthWest')
legend boxoff
%title(['Mean  in data fit to weak Allee on birth'])
xlim([ 0 332])
subplot(1,2,2)
hold on
plot(timevec, vardatavec, 'g*')
for j = 1:length(N0list)
    modelfun_V =@(p)gen_model_var(p,tsamp, Ninit(j), modelcode);
plot(tsamp(2:end), modelfun_V(pbestwkAb), 'b-', 'LineWidth', 2) % test model expected
end
xlabel('time (hours)', 'FontSize',16)
ylabel('Variance in cell number (\Sigma(t))', 'FontSize',16)
legend('variance in data', 'best fit weak Allee', 'FontSize',12, 'Location', 'NorthWest')
legend boxoff
%title(['Variance in data fit to weak Allee on birth'])
xlim([ 0 332])


%% Perform fitting: fit all data to weak/strong Allee on death
% For each fitting, we want to write functions in terms of the functions
% that run the model we are fitting to forward. The gen_model_mu runs the
% system of ODEs for the parameters and generates an expected model mean
% and variance for the input p parameters
% information

modelcode = 6;
tsamp = 0:4:332; % set sample time (if data is down sampled change accordingly
tsampfit = tsamp(2:end);
tsampnew = tsampfit(ind);
tsamp = [0 tsampnew];
% Initial guess
bguess = (log(mudatavec(end-2))-log(mudatavec(end-10)))/(timevec(end-2)-timevec(end-10)) + 0.002; 
dguess = 0.002;
Aguess = -1;
tauguess = 2;
theta = [bguess,dguess, Aguess, tauguess];
Ninit = N0list;
N= numsamps;
[pbestwkAd, BICwkAd, AICwkAd,negLLfitwkAd] = ML_fit(theta, tsamp, Ninit, N,mudatavec, vardatavec, modelcode)

% need to define th inline functions that run forward the mean and the
% variance
modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);


% Perform fitting: fit all data to weak/strong Allee on birth & death
% For each fitting, we want to write functions in terms of the functions
% that run the model we are fitting to forward. The gen_model_mu runs the
% system of ODEs for the parameters and generates an expected model mean
% and variance for the input p parameters
% information

modelcode = 7;

% Initial guess
bguess = (log(mudatavec(end-2))-log(mudatavec(end-4)))/(timevec(end-2)-timevec(end-4)) + 0.002; 
dguess = 0.002;
Aguess = -1;
tauguess = 2;
theta = [bguess,dguess, Aguess, tauguess];
Ninit = N0list;
N= numsamps;
[pbestwkAbd, BICwkAbd, AICwkAbd,negLLfitwkAbd] = ML_fit(theta, tsamp, Ninit, N,mudatavec, vardatavec, modelcode)

% need to define th inline functions that run forward the mean and the
% variance
modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);



%% Make vector of BIC values, negLL values, and parameter numbers
num_params = [ 2, 3, 3, 3, 4, 4, 4];
modellist = [1, 2, 3, 4, 5, 6, 7];
BICvals = [ BICbd, BICstrAb, BICstrAd, BICstrAbd, BICwkAb, BICwkAd, BICwkAbd];
negLLvals = [negLLfitbd, negLLfitstrAb, negLLfitstrAd, negLLfitstrAbd, negLLfitwkAb, negLLfitwkAd, negLLfitwkAbd];
modelnames = {'  ' 'b-d', 'strAb', 'strAd', 'strAbd', 'wkAb', 'wkAd', 'wkAbd', '   '};
shortnames = { 'strAb', 'wkAb'};
figure;
subplot(2,1,1)
semilogy(modellist, BICvals, 'ro', 'LineWidth', 4)
xlabel('Model', 'FontSize',16)
ylabel('BIC value', 'FontSize',16)
xlim([0 8])
ylim([ min(BICvals)-300 max(BICvals)+300])
title('/DeltaBIC values for each model')
set(gca,'xticklabel',modelnames)


% subplot(3,1,2)
% semilogy(modellist, negLLvals, 'ro', 'LineWidth', 4)
% xlabel('Model')
% ylabel(' best negLL')
% xlim([0 8])
% ylim([ min(negLLvals)-5, max(negLLvals)+5])
% title('negLL values for each model (sampling every 4 hours)')
% set(gca,'xticklabel',modelnames)

subplot(2,1,2)
plot(modellist, num_params, 'bo', 'LineWidth', 4)
xlabel('Model')
ylabel(' number of parameters')
xlim([0 8])
ylim([ 1 5 ])
title ('number of parameters')
set(gca,'xticklabel',modelnames)
%% FInd weights

deltaBIC = BICvals- min(BICvals)
sumexps = sum(exp(-0.5*deltaBIC))

for i = 1:length(deltaBIC)
weights(i) = exp(-0.5*deltaBIC(i))./sumexps
end

PstrAb = weights(2)./sum(weights)
PwkAb = weights(5)./sum(weights)
factor = weights(2)./weights(5)
%%
figure;
subplot(3,1,1)
plot(modellist, deltaBIC, 'ro', 'LineWidth', 4)
xlabel('Model', 'FontSize',16)
ylabel('\DeltaBIC from data', 'FontSize',16)
xlim([0 8])
title('\DeltaBIC values for each model fit to experimental data', 'FontSize',14)
set(gca,'xticklabel',modelnames, 'FontSize',14)


subplot(3,1,2)
plot(modellist, weights, 'go', 'LineWidth', 4)
xlabel('Model', 'FontSize',16)
ylabel('BIC weight','FontSize',16)
xlim([0 8])
title('BIC weights','FontSize',16)
set(gca,'xticklabel',modelnames, 'FontSize',14)

subplot(3,1,3)
plot(modellist, num_params, 'b*', 'LineWidth', 4)
xlabel('Model','FontSize',16)
ylabel(' number of parameters','FontSize',16)
xlim([0 8])
ylim([ 1 5 ])
title ('number of parameters','FontSize',16)
set(gca,'xticklabel',modelnames, 'FontSize',14)

%% What if we just compare those two
deltaBIC2 = [0, BICvals(5)-min(BICvals)];
sumexps2 = sum(exp(-0.5*deltaBIC2))

for i = 1:length(deltaBIC2)
weights2(i) = exp(-0.5*deltaBIC2(i))./sumexps2
end

PstrAb2 = weights2(1)./sum(weights2)
PwkAb2 = weights2(2)./sum(weights2)
factor2 = weights2(1)./weights2(2)

%%
figure;
plot([1 2], [(BICvals(2)) (BICvals(5))], 'ro', 'LineWidth',4)
xlabel('Model')
% xlim([ 0 3])
ylim([ 1800 1860])
ylabel('BIC value')
title('BIC model comparison')
set(gca,'xticklabel',shortnames)

save('../out/BICvals24hrs.mat', 'BICvals')
save('../out/weights24hrs.mat', 'weights')
%% Load in data sets
BICvals4hrs = load('../out/BICvals4hrs.mat');
BICvals4hrs = cell2mat(struct2cell(BICvals4hrs));
BICvals12hrs = load('../out/BICvals12hrs.mat');
BICvals12hrs = cell2mat(struct2cell(BICvals12hrs));
BICvals24hrs = load('../out/BICvals24hrs.mat');
BICvals24hrs = cell2mat(struct2cell(BICvals24hrs));

BICall = vertcat(sort(BICvals4hrs), sort(BICvals12hrs), sort(BICvals24hrs))

modelnames = {'  ' 'b-d', 'strAb', 'strAd', 'strAbd', 'wkAb', 'wkAd', 'wkAbd', '   '};
timeres = {'every 4 hours', 'every 12 hours', 'every 24 hours'}
figure;
bar(BICall)
set(gca,'xticklabel',timeres)
ylabel('BIC value')
legend('wkAb', 'strAb', 'wkAbd', 'strAbd', 'wkAd', 'strAd', 'b-d')
title('BIC values with decreasing time resolution')
