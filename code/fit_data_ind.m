% function [pure] = fitData( mix, avg_sd)
% function fits the purepopulation mixtures to a model of exponential
% growth using the red cell count and green cell count individual
% subpopulation growth dynamic data from the incucyte
close all; clear all; clc
S = load('../out/pure.mat');
pure= S.pure;
%%
% single exponential
pfxform = @(pval)[1].*log(pval); %'forward' parameter transform into Reals
pbxform = @(phat)[1].*exp(phat); %'backward' parameter transform into model space
yfxform = @(y)log(y); % 'forward' transform for data and model output
ybxform = @(yhat)exp(yhat); % 'inverse' transform for data and model output

sigma = 0.1; % * avg_sd... need to find measurement error in Incucyte collection

for j = 1:length(pure)
    time = pure(j).time;
    ind0 = find(time == 0);
    N = pure(j).post_sg; % change depend on type of processing we choose
    N0 = pure(j).Nplate; % start by using what we think FACS plated

    
    % pick indices in the middle roughly for guessing
    l = length(time); % length of each replicate
    % picked 0.6 because seemed to be exponential growth for all curves...
    lup = 120; % indices for guessing g
    llow = N0; 

    
    % Going to give the model N0
    modelfungrowth = @(p)singleexpmodel(p, N0,time); % p is single parameter, N0g is initial Vol

    
    % INITIAL GUESSES BASED ON DATA
    gguess= (yfxform(N(lup))-yfxform(N(llow)))/(time(lup)-time(llow));
 
    % Initial guess matrices
    theta = gguess; % g only

    

    
    % loglikelihood calculation
    %cell growth parameters
    loglikelihood = @(phat)sum(log(normpdf(yfxform(N),yfxform(modelfungrowth(pbxform(phat))), sigma)));
    LLcontrib= @(phat)(log(normpdf(yfxform(N),yfxform(modelfungrowth(pbxform(phat))), sigma)));
    
    objfun = @(phat)-loglikelihood(phat);

    phatbest = fminsearch(objfun, pfxform(theta)); % find best fitting parameters

    
    % save parameters and model into mix structure
    pure(j).phat = pbxform(phatbest);

    
    pure(j).tseries = 1:1:time(end); % make model time only go from 1:end
    % make model cell number based on average of N0g inputs

    pure(j).Nmodel = singleexpmodel(pbxform(phatbest),N0, time); % model fit for plotting
    if any(isnan( pure(j).Nmodel)) 
        pure(j).Nmodel = zeros([length(pure(j).Nmodel),1]);
    end
    pure(j).LL = LLcontrib(phatbest);
    
    pure(j).residuals = pure(j).Nmodel - pure(j).post_sg;
    pure(j).Nbar = mean(pure(j).post_sg); % average value of volumes over all measured times
    pure(j).Rsq = 1- (sum((pure(j).residuals).^2)./(sum((pure(j).Nbar-pure(j).post_sg).^2)));
    
    
    % AIC
    % AIC for single exponential model
    num_params1 = 1;
    n = length(pure(j).time);
    AIC = -2*loglikelihood(phatbest) + 2*num_params1;
    pure(j).AIC= AIC;
    
    % Now try using lsqnonlin
    
    LB = -Inf ;  % Lower Bounds
    UB = Inf; % Upper Bounds
    params0 = 0;% Initial Guess...
    options = optimset('TolFun',1e-12,'Tolx',1e-12,'MaxIter',1000,'Display','off','FinDiffRelStep',1e-3);
    [lsq_g, resnorm, reslsq]= lsqnonlin(@fit_singleexp, params0, LB, UB, options, N, N0, time);
    pure(j).lsq_g = lsq_g;
    pure(j).residuals_lsq = reslsq;
    pure(j).Rsqlsq = 1- (sum((reslsq.^2))./(sum((pure(j).Nbar-pure(j).post_sg).^2)));
    pure(j).AIClsq = n*log(sum(reslsq.^2)./n) +2*num_params1;
    pure(j).Nmodellsq = singleexpmodel(lsq_g,N0, time); % model fit for plotting



    
end
%% Fit to the Allee effect

% Allee effect fitting
Nmeas = [];
Nmeas0 = [];
tlong = [];
for i = 1:length(pure)
    Nmeas = vertcat(Nmeas, pure(i).post_sg);
    tlong = vertcat(tlong, pure(i).time);
    Nmeas0 = vertcat(Nmeas0, pure(i).Nplate);
end
LB = zeros(3,1);  % Lower Bounds
UB = [Inf Inf Inf]; % Upper Bounds
params0 = [.02; 8; 1e7];% Initial Guess... g, A, carcap
options = optimset('TolFun',1e-12,'Tolx',1e-12,'MaxIter',1000,'Display','off','FinDiffRelStep',1e-3);

num_samps = length(pure);
%%
paramsAllee= lsqnonlin(@fit_Allee, params0, LB, UB, options, Nmeas, Nmeas0, tlong, num_samps);
%%
save('../out/puref.mat', 'pure')