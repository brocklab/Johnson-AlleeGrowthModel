% This script fits the wells of 1, 1&2, and 2 from the green12 data set
% (12-20 expt) to a single exponential model.
% It also plots the raw data, truncates the data, at 2 weeks, and smooths
% the data


close all; clear all; clc
S = load('../out/green12.mat');
green= S.green;
%% Truncate data to two weeks, add 2 days, set N0
for j = 1:length(green)
    t_tot = green(j).time;
    ind = t_tot < 120; % make this match green8 last time
    new_t = t_tot(ind)
    green(j).t(1) = 0;
    green(j).t = vertcat(green(j).t, new_t);
    green(j).N = vertcat(green(j).N0actual, green(j).cellnum(ind));
    post_sg = smooth(green(j).N(2:end),30,'sgolay');
    green(j).Nsmooth = vertcat(green(j).N0actual, post_sg);
end
%% Fit data to a single exponential
for j = 1:length(green)
    time = green(j).t;
    n = length(time);
    num_params1 = 1;
    ind0 = find(time == 0);
    N = green(j).Nsmooth; % change depend on type of processing we choose
    N0 = green(j).N0actual; % start by using what we think FACS plated
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

%% Plot raw data


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
%% Get distribution of gs for those that have colonies and reasonable R-
for j = 1:3
greensum.gdistrib(j) = 
end
%%
for i = 1:length(uniqueN0-1)
for j = 1:length(green)
    if green(j).N0actual == uniqueN0(i+1) && green(j).num_colonies ~=0 && green(j).Rsqlsq>0.3
        greensum.g(i)=vertcat(greensum.g(i), green(j).g);
    end
    
end
end

    
%%
for j = 1:length(green)
    time = green(j).time;
    n = length(time);
    num_params1 = 1;
    ind0 = find(time == 0);
    N = green(j).cellnum; % change depend on type of processing we choose
    N0 = green(j).Nseed; % start by using what we think FACS plated
    green(j).Nbar = mean(green(j).cellnum); % average value of volumes over all measured times
    
    % First fit each curve to single exponential over entire time course
    LB = -Inf ;  % Lower Bounds
    UB = Inf; % Upper Bounds
    params0 = 0;% Initial Guess...
    options = optimset('TolFun',1e-12,'Tolx',1e-12,'MaxIter',1000,'Display','off','FinDiffRelStep',1e-3);
    [gtot, resnorm, reslsq]= lsqnonlin(@fit_singleexp, params0, LB, UB, options, N, N0, time);
    green(j).gtot = gtot;
    green(j).residuals_gtot = reslsq;
    green(j).Rsqlsq = 1- (sum((reslsq.^2))./(sum((green(j).Nbar-green(j).cellnum).^2)));
    green(j).AIClsq = n*log(sum(reslsq.^2)./n) +2*num_params1;
    green(j).Nmodellsq = singleexpmodel(gtot,N0, time); % model fit for plotting
    
    green(j).maxN = max(green(j).cellnum);
    green(j).last = mean(green(j).cellnum(end-5:end));
    thres = 8*green(j).Nseed;
    if green(j).last > thres
        green(j).takeoff = 1;
    else
        green(j).takeoff = 0;
    end
    
   
        
end

ct = zeros(4,3);

for j = 1:length(green)
    if green(j).Nseed == 8
        ct(1,2) = ct(1,2)+1;
        if green(j).takeoff == 1
            ct(1,1) = ct(1,1)+1;
        end
    end
    
    if green(j).Nseed == 16
        ct(2,2) = ct(2,2)+1;
        if green(j).takeoff == 1
            ct(2,1) = ct(2,1)+1;
        end
    end
    
    
    if green(j).Nseed == 24
        ct(3,2) = ct(3,2)+1;
        if green(j).takeoff == 1
            ct(3,1) = ct(3,1)+1;
        end
    end
    if green(j).Nseed == 32
        ct(4,2) = ct(4,2)+1;
        if green(j).takeoff == 1
            ct(4,1) = ct(4,1)+1;
        end
    end
    
    
end

% gives number that took off, number total, and % take-off
ct(:,3) = 100.*(ct(:,1))./(ct(:,2));

% Make green sum to group by Nseed
Nseed = [8 16 24 32];
for i = 1:length(Nseed)
    greensum(i).gtot = [];
    greensum(i).tkoffs = [];
    greensum(i).gtkoff = [];
end


greensum(1).cc = [1 0 0];
greensum(2).cc = [0 1 0];
greensum(3).cc = [0 1 1];
greensum(4).cc = [0 0 1];


for i = 1:length(Nseed) % number of unique seed numbers
    greensum(i).Nseed = Nseed(i);
    for j = 1:length(green)
        N0 = green(j).Nseed;
        time = green(j).time;
        if green(j).Nseed == Nseed(i)
            greensum(i).gtot = vertcat(greensum(i).gtot, green(j).gtot);
            greensum(i).tkoffs = vertcat(greensum(i).tkoffs, green(j).takeoff);
            if green(j).takeoff ==1
            greensum(i).gtkoff= vertcat(greensum(i).gtkoff, green(j).gtot);% model fit for plotting
            end
        end
    end
    
    greensum(i).pct_takeoff = ct(i,3); 
end
%

for i= 1:length(Nseed)
    
    greensum(i).var_g = (std(greensum(i).gtot))^2;
    greensum(i).var_gtkoff = (std(greensum(i).gtkoff))^2;
    greensum(i).norm_var = greensum(i).var_g/greensum(i).Nseed;
    greensum(i).mean_g = mean(greensum(i).gtot);
    greensum(i).mean_gtkoff = mean(greensum(i).gtkoff);  
    greensum(i).model_gtkoff= singleexpmodel(greensum(i).mean_gtkoff,greensum(i).Nseed, time);% model fit for plotting
            
end
%%
save('../out/greenf.mat', 'green')
save('../out/greensum.mat', 'greensum')