% Preliminary data fit to BT-474s of all inital cell numbers to a single
% exponential model

% Metrics we want to compute for this data set:
% for each individual well:
    % growth rate
    % R-squared value
    % AIC value
    % model fit
    
% for groups of wells:
    % percent of take-off within each Nseed
    % percapita growth rate 
    % variance of per capita growth rate

close all; clear all; clc

%% Run this chunk to load in 2 and 5 cell per well data sets
S = load('../out/BTall.mat');
BT= S.BT;
%% DONT RUN Run this chunk to load pilot study data
S = load('../out/BTpsdata.mat');
BT= S.BT;
%%
for j = 1:length(BT)
    if BT(j).badfit ==0
    % want to fit on all but t=0 (we give it this)
    time = [];
    N = [];
    
    time = BT(j).time;
    n = length(time);
    num_params = 1;
    ind0 = find(time == 0);
    % again fit on all but t=0

    N = BT(j).cellnum; % change depend on type of processing we choose
    N0 = BT(j).N0; % start by using what we think FACS plated
    BT(j).Nbar = mean(N); % average value of volumes over all measured times
    
    % First fit each curve to single exponential over entire time course
    LB = -Inf ;  % Lower Bounds
    UB = Inf; % Upper Bounds
    params0 = 0.1;% Initial Guess...
    options = optimset('TolFun',1e-12,'Tolx',1e-12,'MaxIter',1000,'Display','off','FinDiffRelStep',1e-3);
    [g, resnorm, reslsq]= lsqnonlin(@fit_singleexp, params0, LB, UB, options, N, N0, time);
    BT(j).g=g;
    BT(j).residuals_g = reslsq;
    BT(j).Rsq = 1- (sum((reslsq.^2))./(sum((BT(j).Nbar-N).^2)));
    BT(j).AICg = n*log(sum(reslsq.^2)./n) +2*num_params;
    BT(j).Nmodel = singleexpmodel(g,N0, BT(j).time); % model fit for plotting (include t=0)
    
    
    BT(j).maxN = max(N);
    BT(j).last = mean(BT(j).cellnum(end-5:end));
    thres = 3*BT(j).Nseed;
  
    end
end
%% Go through all of the trajectories and compare to fits
% Use this to go through the actual counts and write down the numbers of
% the wells that should be "excluded" from the count by eye
for i =1:length(BT)
    hold off
    plot(BT(i).time, BT(i).cellnum, '*')
    hold  on
    plot(BT(i).time, BT(i).Nmodel, 'r-')
    xlabel('time(hours)')
    ylabel('Number of cells')
    title(['N_{0}=', num2str(BT(i).N0),', R-sq=', num2str(BT(i).Rsq), ', badfit=',num2str(BT(i).badfit), ', dieoff=', num2str(BT(i).dieoff), ', persist=', num2str(BT(i).persist)])%, ', well=', BT(i).well,' date=', BT(i).date])
    
   pause
end