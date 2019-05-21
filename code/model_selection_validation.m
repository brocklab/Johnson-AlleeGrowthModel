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

% Set up simulation
Nsim = 5000;
tsamp = [ 0:4:160];
b = 0.0238;
d = 0.005;
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
[ NsampAb,NstatAb, CstatAb, musimdatavecAb, varsimdatavecAb, timesimdatavec] = run_strAllbmodel(paramsA, tsamp, Ninit, Nsim);
%[ NsampAbd,NstatAbd, CstatAbd, musimdatavecAbd, varsimdatavecAbd, timesimdatavec] = run_strAllbdmodel(paramsA, tsamp, Ninit, Nsim);
%[ Nsampbd,Nstatbd, Cstatbd, musimdatavecbd, varsimdatavecbd, timesimdatavec] = run_bdmodel(params, tsamp, Ninit, Nsim);
sim_mu_clean = musimdatavecAb;
sim_var_clean = varsimdatavecAb;
% take mean and variance of data and add noise
% Note, we expect measurement noise will not be constant, but will be
% proportional to the correspnding measurement, so we simulate accordingly
%% Plot the stochastic trajectories


%%
eta = 0.05;

noise_mu = (eta.*(sim_mu_clean)).*(randn(length(sim_mu_clean),1));
sim_mean = sim_mu_clean + noise_mu;

noise_var = (eta.*(sim_var_clean)).*(randn(length(sim_var_clean),1));
sim_var = sim_var_clean + noise_var;

%% Plot the current simulated data for this run
modelcode=2;
N=Nsim;
V0=0;
modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);
modelfun_V4= @(p)gen_model_v4(p, tsamp, Ninit, modelcode);
var_in_mean =  @(p)(1/N).*(modelfun_V(p)); % vertical
var_in_var = @(p)(1/N).*(modelfun_V4(p)-(((N-3)./(N-1)).*(modelfun_V(p).^2)));

%% 4th order variance is fine
figure;
hold on
plot(timesimdatavec, modelfun_V4(paramsA), 'b*') % test model expected
for i = 1:length(Ninit)
    plot(tsamp, CstatAb(:,6,i), 'c.')
    text(tsamp(end-10), CstatAb(end-10,6,i),['N_{0}=', num2str(Ninit(i))])
end

legend(' v4 from function', 'expected v4')
%% variance is fine
figure;
hold on
plot(timesimdatavec, sim_var, 'r*')
plot(timesimdatavec, modelfun_V(paramsA), 'b*') % test model expected
for i = 1:length(Ninit)
    plot(tsamp, CstatAb(:,3,i), 'c.')
    text(tsamp(end-10), CstatAb(end-10,3,i),['N_{0}=', num2str(Ninit(i))])
end
legend(' v from data', 'vfrom function', 'expected v')
%% mean is fine
figure;
hold on
plot(timesimdatavec, modelfun_mu(paramsA), 'b*') % test model expected
plot(timesimdatavec, sim_mean, 'r*')
hold on
for i = 1:length(Ninit)
    plot(tsamp, CstatAb(:,1,i), 'k.')
    text(tsamp(end-10), CstatAb(end-10,1,i),['N_{0}=', num2str(Ninit(i))])
end
legend(' mu from function', 'mu data', 'expected mu')
%% Plot data against model for this run
figure;
hold on
plot(timesimdatavec, sim_mean, 'r*')
%plot(timesimdatavec, modelfun_mu(paramsA), 'b*') % test model expected
%variance
for i = 1:length(Ninit)
    plot(tsamp, CstatAb(:,1,i), 'k.')
    text(tsamp(end-10), CstatAb(end-10,1,i),['N_{0}=', num2str(Ninit(i))])
end
xlabel('time(hours)')
ylabel('mean cell number')
title(['Simulated mean from Allee model, b =', num2str(round(b,4)),', d=', num2str(round(d,4)),', A=', num2str(A),', \eta=', num2str(eta),', N_{sim}=', num2str(Nsim)])
legend('simulated data mean', 'expected mean')

figure;
hold on
plot(timesimdatavec, sim_var, 'g*')
hold on
%plot(timesimdatavec, modelfun_V(paramsA), 'b*') % test model expected
%variance
for i = 1:length(Ninit)
    plot(tsamp, CstatAb(:,3,i), 'k.')
    text(tsamp(end-10), CstatAb(end-10,3,i),['N_{0}=', num2str(Ninit(i))])
end
xlabel('time(hours)')
ylabel('variance in cell number')
title(['Simulated variance from Allee model, b =', num2str(round(b,4)),', d=', num2str(round(d,4)),', A=', num2str(A),', \eta=', num2str(eta),', N_{sim}=', num2str(Nsim)])
legend('simulated data variance', 'expected variance')
%% Perform fitting: first fit simulated mean and variance data to b-d model
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

%% Plot the data compared to the best fit model
figure;
hold on
plot(timesimdatavec, sim_mean, 'r*')
plot(timesimdatavec, modelfun_mu(pbestbd), 'k.')
for i = 1:length(Ninit)
    text(tsamp(end-10), CstatAb(end-10,1,i),['N_{0}=', num2str(Ninit(i))])
end
xlabel('time(hours)')
ylabel('mean cell number')
title('Simulated data versus best fit of mean')
legend('simulated data mean', 'best fit mean')

figure;
hold on
plot(timesimdatavec, sim_var, 'g*')
plot(timesimdatavec, modelfun_V(pbestbd), 'k.')
for i = 1:length(Ninit)
    text(tsamp(end-10), CstatAb(end-10,3,i),['N_{0}=', num2str(Ninit(i))])
end
xlabel('time(hours)')
ylabel('variance in cell number')
title('Simulated data versus best fit of variance')
legend('simulated data variance', 'best fit variance')
%% Fit simulated mean and variance data to the strong Allee model on birth
modelcode = 2;
bguess = 0.0238;
%bguess = (log(mudatavec(end-2))-log(mudatavec(end-10)))/(timevec(end-2)-timevec(end-10)) + 0.002; 
dguess = 0.005;
Aguess = 2;
theta = [bguess,dguess, Aguess];
[pbeststrAb, BICstrAb, AICstrAb,negLLfitstrAb] = ML_fit(theta, tsamp, Ninit, N,mudatavec, vardatavec, modelcode);

modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);

%% Plot the data compared to the best fit model
modelcode = 2;
modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);

figure;
hold on
plot(timesimdatavec, sim_mean, 'r*')
for i = 1:length(Ninit)
modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit(i), modelcode);  
plot(tsamp(2:end), modelfun_mu(pbeststrAb), 'k-', 'LineWidth',2)
end
for i = 1:length(Ninit)
    text(tsamp(end-10), CstatAb(end-10,1,i),['N_{0}=', num2str(Ninit(i))], 'FontSize',14)
end
xlabel('time (hours)', 'FontSize', 16)
ylabel('Mean cell number (<n(t)>)', 'FontSize',16)
%title('Simulated data versus best fit to strong Allee model')
legend('simulated data mean', 'best fit mean', 'Location', 'NorthWest', 'FontSize',14)
legend boxoff


figure;
hold on
plot(timesimdatavec, sim_var, 'g*')
for i = 1:length(Ninit)
modelfun_V =@(p)gen_model_var(p,tsamp, Ninit(i), modelcode);  
plot(tsamp(2:end), modelfun_V(pbeststrAb), 'k-', 'LineWidth',2)
end
for i = 1:length(Ninit)
    text(tsamp(end-10), CstatAb(end-10,3,i),['N_{0}=', num2str(Ninit(i))], 'FontSize',14)
end
xlabel('time (hours)', 'FontSize', 16)
ylabel('Variance in cell number (\Sigma(t))', 'FontSize',16)
%title('Simulated data versus best fit to strong Allee model')
legend('simulated data variance', 'best fit variance','Location', 'NorthWest', 'FontSize',14)
legend boxoff

%% Profile likelihood of strong Allee on birth model
% Use these to determine parameter certainty and find the 95 % CI on the
% b and d parameter estimate from all data

num_params = 3;

 % Check fxn

        fbest = negLLfitstrAb;
        factor = 0.02;
        numpoints = 10;
        params = pbeststrAb;
        [profilesstrAb] = profile_likelihood(params, tsamp, Ninit, numsamps ,mudatavec, vardatavec, modelcode, factor, numpoints);
        
        threshold = chi2inv(0.95,length(params))/2 + 1.05*fbest;
        ilo1=find(profilesstrAb(:,2,1)<threshold);
        if ilo1(end) ==numpoints*2 ||ilo1(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIstrAb(1,:) = [NaN, NaN];
        else
             CIstrAb(1,:) = [profilesstrAb(ilo1(1),1,1), profilesstrAb(ilo1(end),1,1)];
        end
        ilo2 = find(profilesstrAb(:,2,2)<threshold);
        if ilo2(end) ==numpoints*2 ||ilo2(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIstrAb(2,:)= [NaN NaN];
        else
             CIstrAb(2,:) = [profilesstrAb(ilo2(1),1,2), profilesstrAb(ilo2(end),1,2)];
        end
        ilo3 = find(profilesstrAb(:,2,3)<threshold);
        if ilo3(end) ==numpoints*2 ||ilo3(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIstrAb(3,:)= [NaN NaN];
        else
             CIstrAb(3,:) = [profilesstrAb(ilo3(1),1,3), profilesstrAb(ilo3(end),1,3)];
        end
    % this output a profiles 3D matrix whose first z dimension contains the
    % profiled parameters and the corresponding negLL
    % the CI contains the upper and lower CI on first the b parameter than
    % the d parameter
%% Plot profiles of strong Allee on birth to view CI
% Think about how to visualize this better
 figure;
         subplot(1, 3,1)
        plot(profilesstrAb(:,1,1), profilesstrAb(:,2,1), 'LineWidth', 2)
        hold on
        plot(pbeststrAb(1),negLLfitstrAb,'r*','LineWidth',2)
        plot([paramsA(1), paramsA(1)], [negLLfitstrAb-10, negLLfitstrAb+10], 'g-', 'LineWidth', 2)
        plot([profilesstrAb(1,1,1) profilesstrAb(end,1,1)],[threshold threshold],'r--')
        xlabel(' profiled b values')
        xlim([profilesstrAb(1,1,1) profilesstrAb(end,1,1)])
        %legend ('profiles', 'MLE', 'true parameter','thershold')
        %legend boxoff
        ylabel('Cost Function Value')
        title('b', 'FontSize', 16)
        %title(['b CI = [',num2str(CIstrAb(1,:)),']'])
        %xlim ([ 0.023, 0.0246])

        subplot(1,3,2)
        plot(profilesstrAb(:,1,2), profilesstrAb(:,2,2), 'LineWidth',2)
        hold on
        plot(pbeststrAb(2),negLLfitstrAb,'r*','LineWidth',2)
        plot([paramsA(2), paramsA(2)], [negLLfitstrAb-50, negLLfitstrAb+50], 'g-', 'LineWidth', 2)
        plot([profilesstrAb(1,1,2) profilesstrAb(end,1,2)],[threshold threshold],'r--')
        xlabel(' profiles d values')
        xlim ([profilesstrAb(1,1,2) profilesstrAb(end,1,2)])
        %legend ('profiles', 'MLE', 'true parameter','thershold')
        %legend boxoff
        ylabel('Cost Function Value')
        title('d', 'FontSize', 16)
        %title(['d CI = [',num2str(CIstrAb(2,:)),']'])
        
        subplot(1,3,3)
        plot(profilesstrAb(:,1,3), profilesstrAb(:,2,3), 'LineWidth',2)
        hold on
        plot(pbeststrAb(3),negLLfitstrAb,'r*','LineWidth',2)
        plot([paramsA(3), paramsA(3)], [negLLfitstrAb-100, negLLfitstrAb+100], 'g-', 'LineWidth', 2)
        plot([profilesstrAb(1,1,3) profilesstrAb(end,1,3)],[threshold threshold],'r--')
        xlabel(' profiled A values')
        xlim([profilesstrAb(1,1,3) profilesstrAb(end,1,3)])
        %legend ('profiles', 'MLE', 'true parameter','thershold')
        %legend boxoff
        ylabel('Cost Function Value')
        title('A', 'FontSize', 16)
        xlim ([1.7 2.2])
        %title(['A CI = [',num2str(CIstrAb(3,:)),']'])
        



%% Now fit simulated mean and variance data to the strong Allee model on death
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

%% Plot the data compared to the best fit model
figure;
hold on
plot(timesimdatavec, sim_mean, 'r*')
plot(timesimdatavec, modelfun_mu(pbeststrAd), 'k.')
for i = 1:length(Ninit)
    text(tsamp(end-10), CstatAb(end-10,1,i),['N_{0}=', num2str(Ninit(i))])
end
xlabel('time(hours)')
ylabel('mean cell number')
title('Simulated data versus best fit to strong Allee model Allee on death')
legend('simulated data mean', 'best fit mean')

figure;
hold on
plot(timesimdatavec, sim_var, 'g*')
plot(timesimdatavec, modelfun_V(pbeststrAd), 'k.')
for i = 1:length(Ninit)
    text(tsamp(end-10), CstatAb(end-10,3,i),['N_{0}=', num2str(Ninit(i))])
end
xlabel('time(hours)')
ylabel('variance in cell number')
title('Simulated data versus best fit to strong Allee model Allee on death')
legend('simulated data variance', 'best fit variance')

%% Now fit simulated mean and variance data to the strong Allee model on birth & death
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
%% Plot the data compared to the best fit model
figure;
hold on
plot(timesimdatavec, sim_mean, 'r*')
plot(timesimdatavec, modelfun_mu(pbeststrAbd), 'k.')
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
plot(timesimdatavec, modelfun_V(pbeststrAbd), 'k.')
for i = 1:length(Ninit)
    text(tsamp(end-10), CstatAb(end-10,3,i),['N_{0}=', num2str(Ninit(i))])
end
xlabel('time(hours)')
ylabel('variance in cell number')
title('Simulated data versus best fit to strong Allee model on b & d')
legend('simulated data mean', 'best fit mean')

%% Fit data to weak Allee model on birth
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
%%
modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);

% test

figure; 
plot(timesimdatavec, modelfun_mu(theta))
lb = [0 0 -10 0];
ub = [ 2 2 10 10];
       
pfxform = @(pval)[1 1 1 1].*(pval-lb)./(ub-lb); %'forward' parameter transform (normalized)
pbxform = @(phat) [1 1 1 1].*((ub-lb).*(phat))+lb;
theta
norm_theta = pfxform(theta)
thetaback = pbxform(norm_theta)
%% Plot the data compared to the best fit model
figure;
hold on
plot(timesimdatavec, sim_mean, 'r*')
plot(timesimdatavec, modelfun_mu(pbestAwb), 'k.')
%plot(timesimdatavec, modelfun_mu(paramsw), 'c*')
for i = 1:length(Ninit)
    text(tsamp(end-4), CstatAb(end-4,1,i),['N_{0}=', num2str(Ninit(i))])
end
xlabel('time(hours)')
ylabel('mean cell number')
title('Simulated data versus best fit to weak Allee model on b ')
legend('simulated data mean', 'best fit mean', 'true model')

figure;
hold on
plot(timesimdatavec, sim_var, 'g*')
plot(timesimdatavec, modelfun_V(pbestAwb), 'k.')
%plot(timesimdatavec, modelfun_V(paramsw), 'c*')
for i = 1:length(Ninit)
    text(tsamp(end-4), CstatAb(end-4,3,i),['N_{0}=', num2str(Ninit(i))])
end
xlabel('time(hours)')
ylabel('variance in cell number')  
title('Simulated data versus best fit to weak Allee model on b')
legend('simulated data variance', 'best fit variance', 'true model')

%% Fit data to weak Allee model on death
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

%% Plot the data compared to the best fit model
figure;
hold on
plot(timesimdatavec, sim_mean, 'r*')
plot(timesimdatavec, modelfun_mu(pbestAwd), 'k.')
%plot(timesimdatavec, modelfun_mu(paramsw), 'c*')
for i = 1:length(Ninit)
    text(tsamp(end-4), CstatAb(end-4,1,i),['N_{0}=', num2str(Ninit(i))])
end
xlabel('time(hours)')
ylabel('mean cell number')
title('Simulated data versus best fit to weak Allee model on d ')
legend('simulated data mean', 'best fit mean', 'true model')

figure;
hold on
plot(timesimdatavec, sim_var, 'g*')
plot(timesimdatavec, modelfun_V(pbestAwd), 'k.')
%plot(timesimdatavec, modelfun_V(paramsw), 'c*')
for i = 1:length(Ninit)
    text(tsamp(end-4), CstatAb(end-4,3,i),['N_{0}=', num2str(Ninit(i))])
end
xlabel('time(hours)')
ylabel('variance in cell number')  
title('Simulated data versus best fit to weak Allee model on d')
legend('simulated data variance', 'best fit variance', 'true model')
%% Fit data to weak Allee model on birth & death
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

%% Plot the data compared to the best fit model
figure;
hold on
plot(timesimdatavec, sim_mean, 'r*')
plot(timesimdatavec, modelfun_mu(pbestAwbd), 'k.')
%plot(timesimdatavec, modelfun_mu(paramsw), 'c*')
for i = 1:length(Ninit)
    text(tsamp(end-4), CstatAb(end-4,1,i),['N_{0}=', num2str(Ninit(i))])
end
xlabel('time(hours)')
ylabel('mean cell number')
title('Simulated data versus best fit to weak Allee model on b & d ')
legend('simulated data mean', 'best fit mean', 'true model')

figure;
hold on
plot(timesimdatavec, sim_var, 'g*')
plot(timesimdatavec, modelfun_V(pbestAwbd), 'k.')
%plot(timesimdatavec, modelfun_V(paramsw), 'c*')
for i = 1:length(Ninit)
    text(tsamp(end-4), CstatAb(end-4,3,i),['N_{0}=', num2str(Ninit(i))])
end
xlabel('time(hours)')
ylabel('variance in cell number')  
title('Simulated data versus best fit to weak Allee model on b & d')
legend('simulated data variance', 'best fit variance', 'true model')
%% Make vector of BIC values, negLL values, and parameter numbers
num_params = [ 2, 3, 3, 3, 4, 4, 4];
modellist = [1, 2, 3, 4, 5, 6, 7];
BICvals = [ BICbd, BICstrAb, BICstrAd, BICstrAbd, BICAwb, BICAwd, BICAwbd];
negLLvals = [negLLfit, negLLfitstrAb, negLLfitstrAd, negLLfitstrAbd, negLLfitAwb, negLLfitAwd, negLLfitAwbd];
modelnames = {'b-d', 'strAb', 'strAd', 'strAbd', 'wkAb', 'wkAd', 'wkAbd', };
deltaBIC = BICvals- min(BICvals);

sumexps = sum(exp(-0.5*deltaBIC));

for i = 1:length(deltaBIC)
weights(i) = exp(-0.5*deltaBIC(i))./sumexps;
end

figure;
subplot(3,1,1)
plot(modellist, deltaBIC, 'ro', 'LineWidth', 4)
xlabel('Model')
ylabel('\DeltaBIC value')
xlim([ 0.5 7.5])
title('\DeltaBIC values for each model')
set(gca,'xticklabel',modelnames)

subplot(3,1,2)
plot(modellist, weights, 'go', 'LineWidth', 4)
xlabel('Model')
ylabel('BIC weight')
xlim([ 0.5 7.5])
title('BIC weights for each model')
set(gca,'xticklabel',modelnames)

subplot(3,1,3)
plot(modellist, num_params, 'bo', 'LineWidth', 4)
xlabel('Model')
xlim([ 0.5 7.5])
ylim([ 1 5])
ylabel('number of parameters')
title ('Model Complexity')
set(gca,'xticklabel',modelnames)




save('../out/negLLvals4bd.mat', 'negLLvals')
save('../out/BICval4bd.mat', 'BICvals')


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
