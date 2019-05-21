% This script fits the wells of 8, 16, 24, and 32 cells from the green data set
% (9-28 and 12-1 expt) to a single exponential model.

% Metrics we want to compute for this data set:
% for each individual well:
    % time to take-off
    % growth rate
% for groups of wells:
    % percent of take-off within each Nseed
    % percapita growth rate 
    % variance of per capita growth rate

close all; clear all; clc
S = load('../out/green832.mat');
green= S.green;


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
save('../out/greenf832.mat', 'green')
save('../out/greensum832.mat', 'greensum')