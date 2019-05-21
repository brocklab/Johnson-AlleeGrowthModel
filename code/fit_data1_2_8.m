% This script fits the wells of 1, 1&2, and 2 from the green12 data set
% (12-20 expt) to a single exponential model.
% It also plots the raw data, truncates the data, at 2 weeks, and smooths
% the data


close all; clear all; clc
S = load('../out/green128.mat');
green= S.green;
%% Truncate data to two weeks, add 2 days, set N0
for j = 1:length(green)
    t_tot = green(j).time;
    ind = t_tot < 168 ; % make this 2 weeks exactly
    new_t = t_tot(ind);
    green(j).t = new_t;
    green(j).N = green(j).cellnum(ind);
    post_sg = smooth(green(j).N,30,'sgolay');
    green(j).Nsmooth = post_sg;
end
%% Fit data to a single exponential
for j = 1:length(green)
    time = green(j).t;
    n = length(time);
    num_params1 = 1;
    ind0 = find(time == 0);
    N = green(j).N; % change depend on type of processing we choose
    N0 = green(j).N(1); % start by using what we think FACS plated
    green(j).Nbar = mean(green(j).Nsmooth); % average value of volumes over all measured times
    
    % First fit each curve to single exponential over entire time course
    LB = -Inf ;  % Lower Bounds
    UB = Inf; % Upper Bounds
    params0 = 0;% Initial Guess...
    options = optimset('TolFun',1e-12,'Tolx',1e-12,'MaxIter',1000,'Display','off','FinDiffRelStep',1e-3);
    [g, resnorm, reslsq]= lsqnonlin(@fit_singleexp, params0, LB, UB, options, N, N0, time);
    green(j).g = g;
    green(j).residuals_g = reslsq;
    green(j).Rsqlsq = 1- ((sum(reslsq.^2))./(sum((green(j).Nbar-green(j).Nsmooth).^2)));
    green(j).AICg = n*log(sum(reslsq.^2)./n) +2*num_params1;
    green(j).Nmodelg = singleexpmodel(g,N0, time); % model fit for plotting   
end

%% Plot raw data for an example trajectory

figure;
j = 64;
hold off
plot(green(j).t, green(j).N, '*','color',green(j).color)
hold on
plot(green(j).t, green(j).Nmodelg, 'm','LineWidth',1.5)
xlabel('time (hours)');
ylabel('cell count');
    % xlim ([ 0, (mix(j).time(end)+10)])
if (max(green(j).cellnum))~=0
ylim ([ 0, 20])
end

title([' Raw data: N_{initial}=' num2str(green(j).N0actual),', N_{col}=' num2str(green(j).num_colonies')])

%% Cycle through all data
for j = 59:length(green)
    subplot(1,2,1)
    hold off
        pause
        hold off
        plot(green(j).time, green(j).cellnum, '*','color',green(j).color)
        xlabel('time (hours)');
        ylabel('cell count');
            % xlim ([ 0, (mix(j).time(end)+10)])
        if (max(green(j).cellnum))~=0
        ylim ([ 0, (max(green(j).cellnum))])
        end
        
       title([' Raw data: N_{initial}=' num2str(green(j).N0actual),', N_{col}=' num2str(green(j).num_colonies')])
 subplot(1,2,2)
        hold off
        plot(green(j).t, green(j).N, '.','color',green(j).color)
        hold on
        plot(green(j).t, green(j).Nsmooth, 'k')
        plot(green(j).t, green(j).Nmodelg, 'm','LineWidth',1.5)
        legend('raw data', 'smoothed data', ['model, R-sq=', num2str(round(green(j).Rsqlsq, 2))])
        xlabel('time (hours)');
        ylabel('cell count');
            % xlim ([ 0, (mix(j).time(end)+10)])
        if (max(green(j).N))~=0
        ylim ([ 0, (max(green(j).cellnum))])
        end
       title([' Processed Raw data: N_{initial}=' num2str(green(j).N0actual),', N_{col}=' num2str(green(j).num_colonies')])

 end
       
xlabel('time (hours)');
ylabel('cell count');
% xlim ([ 0, (mix(j).time(end)+10)])
% ylim ([ 0, (max(max(mix(j).Ngreen, mix(j).Nred)+10))])
title(['N_{initial}=' num2str(green(j).N0actual),', N_{colonies}=', num2str(green(j).num_colonies')])
hold off
%% Find take-off pct from Num_colonies alone
for i = 1:length(green)
    N0_list(i) = green(i).N0actual;
end
uniqueN0= unique(N0_list);
ct = zeros(length(uniqueN0), 4);
% ct will contained number of cells starting at N0=0 to (here 3)
% number of wells with colonies in that bit
% percent of wells that formed colonies
% expected percent of wells that formed colonies based on Poisson
% distribution of chance of finding 1 CSC

for j = 1:length(uniqueN0)
    for i = 1:length(green)
    if green(i).N0actual == uniqueN0(j)
        ct(j,1) = ct(j,1)+1;
       if green(i).num_colonies >0
       ct(j,2) = ct(j,2)+1;
       end
    end
    end
end
ct(:,3) = ct(:,2)./ct(:,1)*100; % percent-take-offs experimental
PrnCSC = 1-(ct(2,3)./100);
% 
for j = 2:length(uniqueN0)
    n0 = j-1;
    ct(j,4) = (1-(PrnCSC)^n0)*100;
end

% 1 cell well only
P_est= ct(2,3)
Pext=100-P_est
%% Plot Percent Take-offs actual vs. Expected
figure;
plot(1:1:7, ct(2:end, 3), '-bo', 'LineWidth', 2)
hold on
plot(1:1:7, ct(2:end, 4), '-ro', 'LineWidth', 2)
legend('Observed take-off percents', 'Poisson expected take-off percents')
legend boxoff
xlabel ('N_{0}')
ylabel('Percent take-offs by colony counting')
% do this for just 1 and 2 cell wells
figure;
plot(1:1:2, ct(2:3,3), '-bo', 'LineWidth', 2)
hold on
plot(1:1:2, ct(2:3,4), '-ro', 'LineWidth', 2)
legend('Observed take-off percents', 'Poisson expected take-off percents')
legend boxoff
xlabel ('N_{0}')
ylabel('Percent take-offs by colony counting')

figure;
bar(1:1:7, ct(2:end, 3))

%% Plot raw data by initial cell number
figure (1)
hold off;
for j = 1:length(green)
    plot(green(j).t, green(j).Nsmooth, 'color', green(j).color, 'LineWidth', 1)
    hold on
end
xlim([0, green(j).t(end)])
xlabel('Time (hours)')
ylabel('N')
title('All wells (N = 210) raw data')
%%

for i= 1:7
subplot(2,4,i)
hold off
for j = 1:length(green)
    if green(j).N0actual ==i
    plot(green(j).t, green(j).Nsmooth, 'color', green(j).color, 'LineWidth', 1)
    end
    hold on
end
xlim([0, green(j).t(end)])
ylim([0, 600])
xlabel('Time (hours)')
ylabel('N')
title(['N_{0} = ', num2str(i), ', N_{wells}= ', num2str(ct(i+1, 1))])
end

% Find take-off pct from Num_colonies alone
for i = 1:length(green)
    N0_list(i) = green(i).N0actual;
end
uniqueN0= unique(N0_list);
ct = zeros(length(uniqueN0), 4);
% ct will contained number of cells starting at N0=0 to (here 3)
% number of wells with colonies in that bit
% percent of wells that formed colonies
% expected percent of wells that formed colonies based on Poisson
% distribution of chance of finding 1 CSC

for j = 1:length(uniqueN0)
    for i = 1:length(green)
    if green(i).N0actual == uniqueN0(j)
        ct(j,1) = ct(j,1)+1;
       if green(i).num_colonies >0
       ct(j,2) = ct(j,2)+1;
       end
    end
    end
end
ct(:,3) = ct(:,2)./ct(:,1)*100; % percent-take-offs experimental
PrnCSC = 1-(ct(2,3)./100);
% 
for j = 2:length(uniqueN0)
    n0 = j-1;
    ct(j,4) = (1-(PrnCSC)^n0)*100;
end
%% Filter out bad R-squared values and where num_cols = 0
ct_clean = 0;
ct_death = 0;

for j = 1:length(green)
    if green(j).num_colonies ~=0 && green(j).Rsqlsq > 0.2
        ct_clean = ct_clean +1;
        green(j).clean = 1;
        green(j).death = 0;
    else
        green(j).clean = 0;
        green(j).death = 0;
    end

end
for j = 1:length(green)
    if green(j).num_colonies ==0 && green(j).N0actual ~= 0
        ct_death = ct_death +1;
        green(j).death = 1;
    end
end

% Want to combine clean fits with wells that die to come up with
% distribution
for j = 1:length(green)
    ind(j) = green(j).death == 1 || green(j).clean == 1;
end

greenan = green(ind);
for j = 1:length(greenan)
    if greenan(j).death == 1
        greenan(j). g = 0;
    end
end

%% Plot greenan

% Plot raw data by initial cell number
figure (1)
hold off;
for j = 1:length(greenan)
    plot(greenan(j).t, greenan(j).Nsmooth, 'color', green(j).color, 'LineWidth', 1)
    hold on
end
xlim([0, green(j).t(end)])
ylim([0, 200])
xlabel('Time (hours)')
ylabel('N')
title('Clean Wells (N = 124) raw data')

subplot(1,2,1)
hold off
for j = 1:length(greenan)
    if greenan(j).N0actual ==1
        plot(greenan(j).t, greenan(j).Nsmooth, 'color', greenan(j).color, 'LineWidth', 1)
        hold on
    end
end
xlim([0, greenan(j).t(end)])
ylim([0, 200])
xlabel('Time (hours)')
ylabel('N')
title('N_{0} = 1, N_{wells}=84')

subplot(1,2,2)
hold off
for j = 1:length(greenan)
    if greenan(j).N0actual ==2
        plot(greenan(j).t, greenan(j).Nsmooth, 'color', greenan(j).color, 'LineWidth', 1)
        hold on
    end
end
xlim([0, greenan(j).t(end)])
ylim([0, 200])
xlabel('Time (hours)')
ylabel('N')
title('N_{0} = 2, N_{wells} = 23')


%%
% Make greenan into greensum
N0actual = [1 2 3 4 5 6 7];
for i = 1:length(N0actual)
    greensum(i).gtot = [];
    greensum(i).ct_death= [];
end

for i = 1:length(N0actual) % number of unique seed numbers
    greensum(i).Noactual = N0actual(i);
    for j = 1:length(greenan)
        N0 = green(j).N0actual;
        time = green(j).t;
        if greenan(j).N0actual == N0actual(i)
            greensum(i).gtot = vertcat(greensum(i).gtot, greenan(j).g);
            greensum(i).ct_death = vertcat(greensum(i).ct_death, greenan(j).death);
        end
    end 
end

for j = 1:length(greensum)
    greensum(j).sigma = std(greensum(j).gtot);
    greensum(j).mean = mean(greensum(j).gtot);
    greensum(j).media = mean(greensum(j).gtot);
end
%% Plot growth rate distributions for 1 & 2 cell wells

for j = 1:2
    subplot(1,2,j)
    hold off
    histogram(greensum(j).gtot,[-.02:.001:.05])
    xlabel ('time')
ylabel('Frequency')
xlabel('g_{ind}')
%ylim([ 0 2])
xlim ([ -0.02 0.05])
title ([ 'N_{0} =', num2str(j), ', \sigma=', num2str(round(greensum(j).sigma,3)),', med=',num2str(round(greensum(j).media, 4))])
end

save('../out/greensum18.mat', 'greensum')
save('../out/greenan18.mat', 'greenan')
%% Perform KS test of g distributions
S = load('../out/greensum18.mat');
greensum= S.greensum;

gmat = load('../out/gmat_sim.mat');
gmat = struct2cell(gmat);
gmat = cell2mat(gmat);
%%
for i = 1:2
[h(i),p(i)] = kstest2(greensum(i).gtot, gmat(:,i))
end
figure;
for j = 1:2
    subplot(1,2,j)
    hold off
    histogram(gmat(:,j), [-.02:.001:.05], 'facecolor', 'r')
    hold on
    histogram(greensum(j).gtot,[-.02:.001:.05], 'facecolor', 'b')
    xlabel ('time')
ylabel('Frequency')
xlabel('g_{ind}')
%ylim([ 0 2])
xlim ([ -0.02 0.05])
legend ('stoch sim g_{distrib}','observed g_{distrib}')
legend boxoff
%title ([ 'N_{0} =', num2str(j), ', \sigma=', num2str(round(greensum(j).sigma,3)),', med=',num2str(round(greensum(j).media, 4))])
end

%% Fit all data to Allee and single exponential model
 % Will need to make long t vectors and long models for input from
 % structure greenan
 close all; clear all; clc
 S = load('../out/greenan18.mat');
 greenan= S.greenan;
 S=load('../out/greensum18.mat');
 greensum = S.greensum;
 %% Make data and t long
 tlong = [];
 Ndatalong = [];
 N0 = [];
 for j = 1:length(greenan)
     tlong = vertcat(tlong,greenan(j).t);
     Ndatalong = vertcat(Ndatalong, greenan(j).N);
     N0 = vertcat(N0, greenan(j).N0actual);
 end
 
 %% Fit with lsqnonlin
 
  % First fit each curve to single exponential over entire time course
    LB = -Inf;  % Lower Bounds
    UB = Inf; % Upper Bounds
    params0 = 0.0006;% Initial Guess...
    LBA = [-Inf, 0];
    UBA = [Inf, Inf];
    p0A= [ 0.0006, 1];

    options = optimset('TolFun',1e-12,'Tolx',1e-12,'MaxIter',1000,'Display','off','FinDiffRelStep',1e-3);
    [g, resnorm, resexp]= lsqnonlin(@fit_simmodelsingexplong, params0, LB, UB, options, tlong, Ndatalong, N0);
    [params_A, resnorm, resAllee] = lsqnonlin(@fit_simmodelAlleelong, p0A, LBA, UBA, options, tlong, Ndatalong, N0);
    %% 
    %g, tlong,Ndata, igood, N0
    num_params1 = 1;
    num_paramsA = 2;
    n = length(Ndatalong);
    Nbar = mean(Ndatalong);
    Rsqexp = 1- ((sum(resexp.^2))./(sum((Nbar-Ndatalong).^2)));
    RsqAll= 1- ((sum(resAllee.^2))./(sum((Nbar-Ndatalong).^2)));
    AICexp = n*log(sum(resexp.^2)./n) +2*num_params1;
    AICAll = n*log(sum(resAllee.^2)./n) +2*num_paramsA;
  %%  
    N0sim = [1 2 3 4 5 6 7];
    num_runs = length(N0sim);
    t = 0:1:tlong(end);
    tbig = repmat(t, num_runs, 1);
    num_meas = length(t);
    vert_length = num_meas*num_runs;
    tlongsim = reshape(tbig', vert_length,1);
    
    isim = find(tlongsim>=0);
    %%
    Nmodelexp = simmodelsingexplong(g, tlongsim, isim, N0sim); % model fit for plotting 
    NmodelAllee = simmodelAlleelong(params_A, tlongsim, isim, N0sim); % model fit for plotting
    %% Plot model outputs
    figure;
    subplot(1,2,1)
    plot(tlong, Ndatalong, 'k.', 'LineWidth', .2)
    hold on
    plot(tlongsim, NmodelAllee, 'g.', 'Linewidth', 2)
    xlabel ('time(hours)')
    ylabel('Cell Number')
    xlim([0 tlongsim(end)])
    title ('Allee Model Fit AIC = 56627')
    subplot(1,2,2)
    plot(tlong, Ndatalong, 'k.', 'LineWidth', .2)
    hold on
    plot(tlongsim, Nmodelexp, 'r.', 'LineWidth', 2)
    xlabel ('time(hours)')
    ylabel('Cell Number')
    xlim([0 tlongsim(end)])
    title ('Exponential Model Fit AIC = 56625')
 
 %% Fit Allee effect with Bayesian

% Allee effect equation params g and A
pfxform = @(pval)[1 1].*log(pval); %'forward' parameter transform into Reals
pbxform = @(phat)[1 1].*exp(phat); %'backward' parameter transform into model space
pfxformg = @(pval)[1].*log(pval);
pbxformg = @(phat)[1].*exp(phat);
yfxform = @(y)log(y); % 'forward' transform for data and model output
ybxform = @(yhat)exp(yhat); % 'inverse' transform for data and model output

% Find 0s and mark those indices. Will use censoring on 0 measurements with
% assumption of error being both additive and proportional (i.e. 0
% measurement has a SD of error of 1 cell, so we will sum over the interval
% from 0 to 1 for those data points)
LLOQ = 1;
icens = find(Ndatalong <LLOQ);
igood = find(Ndatalong>= LLOQ);

modelfungood= @(p)simmodelAlleelong(p, tlong, igood, N0);
modelfuncens= @(p)simmodelAlleelong(p, tlong, icens, N0);
modelfunsinggood = @(p)simmodelsingexplong(p, tlong, igood, N0);
modelfunsingcens = @(p)simmodelsingexplong(p, tlong, icens, N0);

gguess = 0.0182;
Aguess = 1.46e-12;
sigma = 0.075;
theta = [gguess, Aguess]; % g and A
thetag = gguess;
Ndatalonghigh = ones(length(icens),1);
Ndatalonglow = zeros(length(icens),1);

loglikelihood = @(phat)(sum(log(normpdf(yfxform(Ndatalong(igood)),yfxform(modelfungood(pbxform(phat))), sigma)))+...
    sum(log(normcdf(yfxform(Ndatalonghigh), yfxform(modelfuncens(pbxform(phat))),1))));
loglikelihoodg = @(phat)(sum(log(normpdf(yfxform(Ndatalong(igood)),yfxform(modelfunsinggood(pbxformg(phat))), sigma)))+...
    sum(log(normcdf(yfxform(Ndatalonghigh), yfxform(modelfunsingcens(pbxformg(phat))),1))));

% Probably need to add normcdf term for upper interval - normcdf of lower
% interval
% minimize objective function for each structural model
objfun = @(phat)-loglikelihood(phat);
objfunsing = @(phat)-loglikelihoodg(phat);
options = optimset('MaxFunEvals',1e5, 'MaxIter', 1e5);
phatbest = fminsearch(objfun, pfxform(theta), options); % find best fitting parameters
phatbestsing = fminsearch(objfunsing, pfxformg(thetag), options);

params_Bayes = pbxform(phatbest);
params_sing = pbxformg(phatbestsing);

iall = find(tlong>=0);

Nfitlong= simmodelAlleelong(params_Bayes, tlong, iall, N0);
Nfitsing = simmodelsingexplong(params_sing, tlong, iall, N0);
% Now find chi-squared, R-squared, and AIC value for Allee and single
% exponential model
%
%residuals (model- measured)
resAllee = Nfitlong-Ndatalong;
resexp = Nfitsing-Ndatalong;
% average N
Nbar = (mean(Ndatalong));
%R-squared
RsqAllee = 1- (sum((resAllee).^2)./(sum((Nbar-Ndatalong).^2)));
Rsqexp = 1- (sum((resexp).^2)./(sum((Nbar-Ndatalong).^2)));
num_p_Allee = 2; % g and A
num_p_exp = 1; %g
AICAllee = -2*loglikelihood(phatbest) + 2*num_p_Allee;
AICexp = -2*loglikelihoodg(phatbestsing) + 2*num_p_exp;
delta_AIC = AICAllee- AICexp; % if value is negative then Allee is better model
% will run a loop through different values of A and eta and find the
% corresponding delta_AIC at each point. This will be plotted on a heat map
% to be used for parameter identifiability

%%
save('../out/greenf.mat', 'green')
save('../out/greensum.mat', 'greensum')