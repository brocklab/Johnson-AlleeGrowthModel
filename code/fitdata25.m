% This script fits the wells of 2 & 5 from combined data sets and finds growth rate and other metrics.

% Metrics we want to compute for this data set:
% for each individual well:
    % time to take-off
    % growth rate after take off
% for groups of wells:
    % percent of take-off within each Nseed
    % percapita growth rate 
    % variance of per capita growth rate

close all; clear all; clc
S = load('../out/greenc.mat');
greenc= S.greenc;
%%
for j = 1:length(greenc)
    time = greenc(j).time;
    n = length(time);
    num_params1 = 1;
    ind0 = find(time == 0);
    N = greenc(j).cellnum; % change depend on type of processing we choose
    N0 = greenc(j).Nseed; % start by using what we think FACS plated
    greenc(j).Nbar = mean(greenc(j).cellnum); % average value of volumes over all measured times
    
    % First fit each curve to single exponential over entire time course
    LB = -Inf ;  % Lower Bounds
    UB = Inf; % Upper Bounds
    params0 = 0;% Initial Guess...
    options = optimset('TolFun',1e-12,'Tolx',1e-12,'MaxIter',1000,'Display','off','FinDiffRelStep',1e-3);
    [gtot, resnorm, reslsq]= lsqnonlin(@fit_singleexp, params0, LB, UB, options, N, N0, time);
    greenc(j).gtot = gtot;
    greenc(j).residuals_gtot = reslsq;
    greenc(j).Rsqlsq = 1- (sum((reslsq.^2))./(sum((greenc(j).Nbar-greenc(j).cellnum).^2)));
    greenc(j).AIClsq = n*log(sum(reslsq.^2)./n) +2*num_params1;
    greenc(j).Nmodellsq = singleexpmodel(gtot,N0, time); % model fit for plotting
    
    greenc(j).maxN = max(greenc(j).cellnum);
    greenc(j).last = mean(greenc(j).cellnum(end-5:end));
    thres = 10*greenc(j).Nseed;
    if greenc(j).last > thres
        greenc(j).takeoff = 1;
    else
        greenc(j).takeoff = 0;
    end
    
   
        
end
%%
ct = zeros(2,3);

for j = 1:length(greenc)
    if greenc(j).Nseed == 2
        ct(1,2) = ct(1,2)+1;
        if greenc(j).takeoff == 1
            ct(1,1) = ct(1,1)+1;
        end
    end
    
    if greenc(j).Nseed == 5
        ct(2,2) = ct(2,2)+1;
        if greenc(j).takeoff == 1
            ct(2,1) = ct(2,1)+1;
        end
    end
    
    
    
end

% gives number that took off, number total, and % take-off
ct(:,3) = 100.*(ct(:,1))./(ct(:,2));
%%
Nseed = [2 5];
greencsum(1).gtot = [];
greencsum(1).tkoffs = [];
greencsum(1).gtkoff = [];
greencsum(2).gtot = [];
greencsum(2).tkoffs = [];
greencsum(2).gtkoff = [];



greencsum(1).cc = [1 0 0];
greencsum(2).cc = [0 1 0];



for i = 1:2 % number of unique seed numbers
    greencsum(i).Nseed = Nseed(i);
    for j = 1:length(greenc)
        N0 = greenc(j).Nseed;
        time = greenc(j).time;
        if greenc(j).Nseed == Nseed(i)
            greencsum(i).gtot = vertcat(greencsum(i).gtot, greenc(j).gtot);
            greencsum(i).tkoffs = vertcat(greencsum(i).tkoffs, greenc(j).takeoff);
            if greenc(j).takeoff ==1
            greencsum(i).gtkoff= vertcat(greencsum(i).gtkoff, greenc(j).gtot);% model fit for plotting
            end
        end
    end
    
    greencsum(i).pct_takeoff = ct(i,3); 
end
%%

for i= 1:2
    
    greencsum(i).var_g = (std(greencsum(i).gtot))^2;
    greencsum(i).var_gtkoff = (std(greencsum(i).gtkoff))^2;
    greencsum(i).norm_var = greencsum(i).var_g/greencsum(i).Nseed;
    greencsum(i).mean_g = mean(greencsum(i).gtot);
    greencsum(i).mean_gtkoff = mean(greencsum(i).gtkoff);  
    greencsum(i).model_gtkoff= singleexpmodel(greencsum(i).mean_gtkoff,greencsum(i).Nseed, time);% model fit for plotting
    greencsum(i).time = time;        
end
%%
save('../out/greencf.mat', 'greenc')
save('../out/greencsum.mat', 'greencsum')