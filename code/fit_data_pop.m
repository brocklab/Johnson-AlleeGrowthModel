% function [pure] = fitDataPop( mix, avg_sd)
% function fits the purepopulation mixtures to a model of exponential
% growth using the red cell count and green cell count individual
% subpopulation growth dynamic data from the incucyte
close all; clear all; clc
S = load('../out/pure.mat');
pure= S.pure;
%% First put data into one big structure
% Pop_structure should contain 
% Vector of length  times measured x number of individuals)
% Vector of length of  cell nums measured x number of individuals 
% Vector of length number of individuals that corresponds to the N0 seeded in
% sorting
% Vector of length N that corresponds to the measured N0 (in case we go
% back to this)
pop.timevec = [];
pop.Nvec = [];
pop.Nplate = [];
pop.N0meas =[];
for j = 98:103%1:length(pure)
    pop.timevec= vertcat(pop.timevec, pure(j).time);
    pop.Nvec = vertcat(pop.Nvec, pure(j).cellnum);
    pop.Nplate = vertcat(pop.Nplate, pure(j).Nplate);
    pop.N0meas = vertcat(pop.N0meas, pure(j).N0avg);
end
%% Next need to set up an model system
% This model needs to find all of the times that are 0, and at that point
% needs to produce the Ns starting from the next Nplate.... test this with
% dummy data...

% Going to give the model N0

% test
tbig = repmat(t, 5, 1);
tlong = reshape(tbig', [size(tbig,1)*size(tbig,2), 1]);
A = 16;
carcap = 1e8;
g = 0.02;
N0 = [ 2 5 10 12 20 ];
P = [g, A, carcap];
   
%modelfunAllee = @(p)Alleemodel(p, time, N0); % p is single parameter, N0g is initial Vol
 modelfunAllee = @(p)Alleemodel(p, tlong, N0); % p is single parameter, N0g is initial Vol
 Nmodellong = modelfunAllee(P);
 
 noise = 0.5*(1-2*randn(length(tlong),length(N0)));

Nfake = Nmodellong + noise;
 figure (1)
 hold off
 plot(tlong, Nmodellong, '.')
 hold on
 plot(tlong, Nfake, '*')
 xlabel ('time')
 ylabel('Allee model cell number')
 title('Test Allee Model')
 legend ('model A = 16, g = 0.02', 'model + noise')
 
    
%% Set up fitting 

% Transform parameters into log space (g, A, k)
pfxform = @(pval)[1 1 1].*log(pval); %'forward' parameter transform into Reals
pbxform = @(phat)[1 1 1].*exp(phat); %'backward' parameter transform into model space
yfxform = @(y)log(y); % 'forward' transform for data and model output
ybxform = @(yhat)exp(yhat); % 'inverse' transform for data and model output


% should this be the standard deviation of the individual growth rates??
sigma = 1; % * avg_sd... need to find measurement error in Incucyte collection

time = round(pop.timevec, 0);
N0 = pop.Nplate;
N = pop.Nvec;

% give the model N0 and long time vector
 modelfunAllee = @(p)Alleemodel(p, time, N0); % p is set of params A
 
 % Model guesses ( should be data driven... but right now 
 gguess = 0.115;
 Aguess = 2;
 carcapguess = 1e4;
 %test
 Nmodel = modelfunAllee([gguess, Aguess, carcapguess]); % What does this error mean?
 
 %%

 
 theta = [gguess, Aguess, carcapguess];
 
 % loglikelihood calculation
 loglikelihood = @(phat)sum(log(normpdf(yfxform(N),yfxform(modelfunAllee(pbxform(phat))), sigma)));
 LLcontrib= @(phat)(log(normpdf(yfxform(N),yfxform(modelfunAllee(pbxform(phat))), sigma)));
    
    objfun = @(phat)-loglikelihood(phat);

    phatbest = fminsearch(objfun, pfxform(theta)); % find best fitting parameters
%%
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