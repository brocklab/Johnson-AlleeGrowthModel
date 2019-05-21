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
%% Run this chunk to load pilot study data
S = load('../out/BTpsdata.mat');
BT= S.BT;
%%
for j = 1:length(BT)
    % want to fit on all but t=0 (we give it this)
    time = [];
    N = [];
    time = BT(j).time(2:end,:);
    n = length(time);
    num_params = 1;
    ind0 = find(time == 0);
    % again fit on all but t=0
    N = BT(j).cellnum(2:end,:); % change depend on type of processing we choose
    N0 = BT(j).N0; % start by using what we think FACS plated
    BT(j).Nbar = mean(BT(j).cellnum); % average value of volumes over all measured times
    
    % First fit each curve to single exponential over entire time course
    LB = -Inf ;  % Lower Bounds
    UB = Inf; % Upper Bounds
    params0 = 0.1;% Initial Guess...
    options = optimset('TolFun',1e-12,'Tolx',1e-12,'MaxIter',1000,'Display','off','FinDiffRelStep',1e-3);
    [g, resnorm, reslsq]= lsqnonlin(@fit_singleexp, params0, LB, UB, options, N, N0, time);
    BT(j).g=g;
    BT(j).residuals_g = reslsq;
    BT(j).Rsq = 1- (sum((reslsq.^2))./(sum((BT(j).Nbar-BT(j).cellnum).^2)));
    BT(j).AICg = n*log(sum(reslsq.^2)./n) +2*num_params;
    BT(j).Nmodel = singleexpmodel(g,N0, BT(j).time); % model fit for plotting (include t=0)
    % Next fit to single exponential death model
    [d, resnorm, reslsqd]= lsqnonlin(@fit_singleexpd, params0, LB, UB, options, N, N0, time);
    BT(j).d = d;
    BT(j).residuals_d = reslsqd;
    BT(j).Rsqd = 1- (sum((reslsqd.^2))./(sum((BT(j).Nbar-BT(j).cellnum).^2)));
    BT(j).AICd = n*log(sum(reslsqd.^2)./n) +2*num_params;
    
    BT(j).maxN = max(BT(j).cellnum);
    BT(j).last = mean(BT(j).cellnum(end-5:end));
    thres = 3*BT(j).Nseed;
    if BT(j).last > thres
        BT(j).takeoff = 1;
    else
        BT(j).takeoff = 0;
    end
    
    if BT(j).Rsq <0.6 || BT(j).Rsqd< 0.6
        BT(j).badfit =1;
    else 
        BT(j).badfit = 0;
    end
    
    if BT(j).AICd < BT(j).AICg
        BT(j).dieout = 1;
    else
        BT(j).dieout = 0;
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
    title(['N_{0}=', num2str(BT(i).N0),', R-sq=', num2str(BT(i).Rsq)])%, ', well=', BT(i).well,' date=', BT(i).date])
    pause
end

%% Plot some things

% growth rate versus No
figure;
for i = 1:length(BT)
    if BT(i).badfit == 0
    semilogx((BT(i).N0), BT(i).g, 'b.')
    hold on
    end
end
xlabel('N_{0}')
ylabel('growth rate')
title('Growth rate versus initial cell number')

% plot the cell trajectories colored by No
% first do all on one plot
figure;
subplot(1,2,1)
for i = 1:length(BT)
    if BT(i). badfit == 0
        plot(BT(i).time, BT(i).cellnum, 'color', BT(i).color)
        hold on
    end
end
xlim ([ 0 336])
xlabel ('time (hours)')
ylabel('number of cells')
title ('Cell number trajectories colored by N0')
subplot(1,2,2)
for i = 1:length(BT)
    if BT(i). badfit == 0
        plot(BT(i).time, log(BT(i).cellnum./BT(i).N0), 'color', BT(i).color)
        hold on
    end
end
xlabel ('time (hours)')
ylabel('log(N/N0)')
xlim ([ 0 336])
title ('Normalized cell number vs. N0')

%% Separate by N0
uniqN0=[4 8 16];
figure;
for j = 1:3
subplot(1,3,j)
for i = 1:length(BT)
    if  BT(i).badfit == 0 && BT(i).Nseed == uniqN0(j)
        plot(BT(i).time, BT(i).cellnum, 'color', BT(i).color)
        hold on
    end
end

xlim ([ 0 250])
ylim([ 0 300])
xlabel ('time (hours)')
ylabel('number of cells')
title ('Cell number trajectories colored by N0')
end

%% Plot example trajectories
figure;
for j = 1:length(BT)
    hold off
    plot(BT(j).time, BT(j).cellnum,'*')
    hold on
    plot(BT(j).time, BT(j).Nmodel, '-')
    xlabel('time')
    ylabel('N')
    title(['N_{0}= ',num2str(BT(j).N0),', R-squared =', num2str(BT(j).Rsq), ', g=', num2str(BT(j).g)])
    pause

end


%% Make a new structure which summarizes for each initial cell number
    % percent of take-off within each Nseed
    % percapita growth rate 
    % variance of per capita growth rate
 
 % find groups by N0
 for i = 1:length(BT)
     N0list(i,1) = BT(i).N0;
 end
 
 uniqN0= unique(N0list); % get list of unique N0s
 
 for i = 1:length(uniqN0)
    BTsum(i).N0 = uniqN0(i);
    BTsum(i).gtot = []; % this will be a vector of all the gs
    BTsum(i).tkoffs = 0; % this will be a count of take offs
    BTsum(i).num = 0; % this will be the number of trajectories for that cell number
    BTsum(i).tmin = [];
 end

 for i = 1:length(uniqN0) % number of unique seed numbers
    for j = 1:length(BT)
        N0 = BT(j).N0;
        time = BT(j).time;
        if BT(j).N0 == uniqN0(i)
            BTsum(i).color = BT(j).color;
             if BT(j). badfit == 0
            BTsum(i).tkoffs = BTsum(i).tkoffs + BT(j).takeoff;
            BTsum(i).num = BTsum(i).num +1;
            BTsum(i).tmin = vertcat(BTsum(i).tmin, BT(j).time(end));
            if BT(j).takeoff ==1 
            BTsum(i).gtot = vertcat(BTsum(i).gtot, BT(j).g);
            end
             end
        end
    end
 end
 
 % need to make a uniform time vector for each initial cell number
 % this will likely change for each cell number ( higher N, less time)
 for i = 1:length(uniqN0)
     BTsum(i).timevec =0:4:((round(min(BTsum(i).tmin))./4)*4); % makes a uniform time vector
     BTsum(i).Nmat = [];
 end

 for i= 1:length(uniqN0)
    BTsum(i).pct_takeoff = 100*(BTsum(i).tkoffs./BTsum(i).num);
    BTsum(i).var_g = (std(BTsum(i).gtot))^2;
    BTsum(i).norm_var = BTsum(i).var_g/BTsum(i).N0;
    BTsum(i).mean_g = mean(BTsum(i).gtot);
    BTsum(i).model_g= singleexpmodel(BTsum(i).mean_g,BTsum(i).N0, BTsum(i).timevec);% model fit for plotting           
 end
 
% This loop makes an a matrix of each N trajectory sampled uniformly in
% time which will be used for the stochastic parameter estimation
%%
for j = 1:length(BTsum)
    tsamp = BTsum(j).timevec;
 for k = 1:length(BT)
     % note change depending on whether or not we want to include badfits
     if BT(k).N0 == BTsum(j).N0 && BT(k).badfit ==0
         tstoch = BT(k).time;
         Nstoch = BT(k).cellnum;
         Nsamp = [];
        for i = 1:length(tsamp)
        % find nearest tstate that is less that tsamp
        ind =find(tstoch<=tsamp(i),1,'last');
        tfind = tstoch(ind);
        Nsamp(i,1)=Nstoch(ind);
        end
        BTsum(j).Nmat= horzcat(BTsum(j).Nmat, Nsamp);
     end
end

end
%%
% Now find the measured mean and variance at each uniformly sampled time
% point
for j = 1:length(BTsum)
    Nsamp = [];
    mu_t = [];
    n_2_t = [];
    var_t = [];
    Nsamp = BTsum(j).Nmat;
    mu_t = mean(Nsamp,2);
    n_2_t = mean((Nsamp.^2),2);
    var_t = n_2_t - ((mu_t).^2);
    BTsum(j).mu_t = mu_t;
    BTsum(j).var_t = var_t;
end
%% Plot some more things

% average growth rate versus cell number
figure;
for i =1:length(uniqN0)
    errorbar((BTsum(i).N0), BTsum(i).mean_g, 1.96*sqrt(BTsum(i).var_g), 'bo')
    set(gca,'XScale','log')
    hold on
end
xlabel ('N0')
ylabel ('mean growth rate')
%xlim ([ 0, 1000])
title ('Average growth rate vs. N_{0} w/ exclusions')

%% mean growth rate model trajectories
figure;
subplot(1,2,1)
for i = 1:length(BTsum)
        plot(BTsum(i).timevec, BTsum(i).model_g, 'color', BTsum(i).color)
        hold on
        text(BTsum(i).timevec(end-10), BTsum(i).model_g(end-10),['N_{0}=', num2str(BTsum(i).N0)], 'FontSize',8)
end
xlim ([ 0 336])
xlabel ('time (hours)')
ylabel('number of cells')
title ('Mean growth rate model by N0')
subplot(1,2,2)
for i = 1:length(BTsum)
        plot(BTsum(i).timevec, log(BTsum(i).model_g./BTsum(i).N0), 'color', BTsum(i).color)
        hold on
        text(BTsum(i).timevec(end-10), log(BTsum(i).model_g(end-10)./BTsum(i).N0),['N_{0}=', num2str(BTsum(i).N0)], 'FontSize',8)

end
xlabel ('time (hours)')
ylabel('log(N/N0)')
xlim ([ 0 336])
title ('Normalized cell number vs. N0')

%% plot the uniform time sampled data for N0=2

figure;
subplot (1,2,1)
for j = 1:BTsum(2).num
plot(BTsum(2).timevec, BTsum(2).Nmat(:,j),'b.')
hold on
end
plot(BTsum(2).timevec, BTsum(2).mu_t,'r-', 'LineWidth', 3)
xlabel ('time')
xlim([0 336])
ylabel('number of cells')
title('Cell number trajectories for N0=2 with mean')
subplot (1,2,2)
for j = 1:BTsum(2).num
plot(BTsum(2).timevec, BTsum(2).Nmat(:,j),'b.')
hold on
end
plot(BTsum(2).timevec, BTsum(2).var_t,'g-', 'LineWidth', 3)
xlim([ 0 336])
xlabel ('time')
ylabel('number of cells')
title('Cell number trajectories for N0=2 with variance')

%% Plot uniform time sample for all initial cell numbers
figure;
for i = 1:length(BTsum)-1
subplot(1,length(BTsum),i)
    for j = 1:BTsum(i).num
    plot(BTsum(i).timevec, BTsum(i).Nmat(:,j),'b.')
    hold on
    end
 plot(BTsum(i).timevec, BTsum(i).mu_t,'r-', 'LineWidth', 3)
xlim([0 336])
xlabel ('time')
ylabel('number of cells')
title(['N0= ', num2str(BTsum(i).N0),' mean'])
end
%%
figure;
for i = 1:length(BTsum)-1
subplot(1,length(BTsum),i)
for j = 1:BTsum(i).num
plot(BTsum(i).timevec, BTsum(i).Nmat(:,j),'b.')
hold on
end
plot(BTsum(i).timevec, BTsum(i).var_t,'g-', 'LineWidth', 3)
xlabel ('time')
xlim([0 336])
ylabel('number of cells')
title(['N0= ', num2str(BTsum(i).N0),' variance'])
end
%% Save the two fitted data structures

% this saves the raw data structure for all imported BT-474 data sets
save('../out/BTfit.mat', 'BT')
save('../out/BTsumfit.mat', 'BTsum')