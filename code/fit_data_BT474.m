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

%% Run this chunk to load in form all seeding numbers well data sets
S = load('../out/BTall.mat');
BT= S.BT;
%% Run this chunk to load pilot study data
S = load('../out/BTpsdata.mat');
BT= S.BT;
%%
for j = 1:length(BT)
    if BT(j).badfit == 0 && BT(j).persist == 0
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
    BT(j).lastN = (BT(j).cellnum(end));
    end
    if BT(j).persist == 1
        BT(j).g=0;
        BT(j).cellnum = repmat(BT(j).N0,length(BT(j).time),1);
        BT(j).Nmodel = repmat(BT(j).N0,length(BT(j).time),1);
    end
    
   
    if BT(j).dieoff==1
      [d, resnorm, reslsq]= lsqnonlin(@fit_singleexpd, params0, LB, UB, options, N, N0, time);
       BT(j).g=-d; 
    end
end
%%
for j = 1:length(BT)
    if ~isempty(BT(j).Rsq)
    if BT(j).Rsq < 0.6 && BT(j). badfit == 0 && BT(j).dieoff == 0 && BT(j).persist == 0
        BT(j).badfit =1;
    end
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
%subplot(1,2,1)
for i = 1:length(BT)
    if BT(i). badfit == 0
        plot(BT(i).time, BT(i).cellnum, 'color', BT(i).color)
        hold on
    end
end
xlim ([ 0 328])

xlabel ('time (hours)')
ylabel('number of cells')
title ('Cell number trajectories colored by N0')
%subplot(1,2,2)
figure;
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
%% Plot example of raw data
figure;
for i = 1:length(BT)
    if BT(i). badfit == 0 && BT(i).N0==4
        plot(BT(i).time, BT(i).cellnum, 'color', BT(i).color)
        hold on
    end
%     if BT(i). badfit == 0 && BT(i).N0==5
%         plot(BT(i).time, BT(i).cellnum, 'color', BT(i).color)
%         hold on
%     end
    if BT(i). badfit == 0 && BT(i).N0==12
        plot(BT(i).time, BT(i).cellnum, 'color', BT(i).color)
        hold on
    end
end
xlim ([ 0 328])
xlabel ('time (hours)')
ylabel('number of cells')
title ('Cell number trajectories colored by N0')

%% Separate by N0
uniqN0 = [];
N0list = [];
for j = 1:length(BT)
    if ~isnan(BT(j).N0)
    N0list =vertcat(N0list, BT(j).N0);
    end
    
end
uniqN0 = unique(N0list);
%%
figure;
for j = 1:length(uniqN0)
subplot(1,length(uniqN0),j)
for i = 1:length(BT)
    if  BT(i).badfit == 0 && BT(i).N0 == uniqN0(j)
        plot(BT(i).time, BT(i).cellnum, 'color', BT(i).color)
        hold on
    end
end

xlim ([ 0 250])
ylim([ 0 300])
xlabel ('time (hours)')
ylabel('N')
title (['N0=', num2str(uniqN0(j))])
end

%% Plot example trajectories
figure;
for j = 1:length(BT)
    if BT(j).badfit == 0 && BT(j).persist ==0
    hold off
    plot(BT(j).time, BT(j).cellnum,'*')
    hold on
    plot(BT(j).time, BT(j).Nmodel, '-')
    xlabel('time')
    ylabel('N')
    title(['N_{0}= ',num2str(BT(j).N0),', R-squared =', num2str(BT(j).Rsq), ', g=', num2str(BT(j).g)])
    end
    pause

end


%% Make a new structure which summarizes for each initial cell number
    % percent of take-off within each Nseed
    % percapita growth rate 
    % variance of per capita growth rate
 
 % find groups by N0

 
 for i = 1:length(uniqN0)
    BTsum(i).N0 = uniqN0(i);
    BTsum(i).gtot = []; % this will be a vector of all the gs
    BTsum(i).dieoffs = 0; % this will be a count of take offs
    BTsum(i).badfits = 0;
    BTsum(i).persisters = 0;
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
            BTsum(i).dieoffs = BTsum(i).dieoffs + BT(j).dieoff;
            BTsum(i).persisters = BTsum(i).persisters +BT(j).persist;
            BTsum(i).badfits = BTsum(i).badfits + BT(j).badfit;
            BTsum(i).num = BTsum(i).num +1;
            BTsum(i).tmin = vertcat(BTsum(i).tmin, BT(j).time(end));
            if BT(j).badfit ==0
            BTsum(i).gtot = vertcat(BTsum(i).gtot, BT(j).g);
            end
             end
        end
    end
 end
 
 % need to make a uniform time vector for each initial cell number
 % this will likely change for each cell number ( higher N, less time)
 for i = 1:length(uniqN0)
     BTsum(i).timevec =0:4:324; % makes a uniform time vector
     BTsum(i).Nmat = [];
 end

 for i= 1:length(uniqN0)
    BTsum(i).pct_takeoff = 100*(1-BTsum(i).dieoffs-BTsum(i).persisters)./(BTsum(i).num);
    BTsum(i).var_g = (std(BTsum(i).gtot))^2;
    BTsum(i).norm_var = BTsum(i).var_g/BTsum(i).N0;
    BTsum(i).mean_g = mean(BTsum(i).gtot);
    BTsum(i).model_g= singleexpmodel(BTsum(i).mean_g,BTsum(i).N0, BTsum(i).timevec);% model fit for plotting           
 end
 
% This loop makes an a matrix of each N trajectory sampled uniformly in
% time which will be used for the stochastic parameter estimation
%%
for j = 1:length(BTsum)
    tsamp = 0:4:324;
 for k = 1:length(BT)
     % note change depending on whether or not we want to include badfits
     if BT(k).N0 == BTsum(j).N0 && BT(k).badfit ==0
        Nsamp=BT(k).cellnum(1:82);
     end
        BTsum(j).Nmat= horzcat(BTsum(j).Nmat, Nsamp);
        Nsamp = [];
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
title ('Average growth rate vs. N_{0}')

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

%% Same plot but sample a few
figure;
subplot(1,3,1)
timemat = [];
timetotvec2 = [];
Ntotvec2 = [];
for j = 1:length(BTsum)
    if BTsum(j).N0 == 2
        timemat = repmat(BTsum(j).timevec, size(BTsum(j).Nmat,2),1);
        timemat = timemat';

timetotvec1 = reshape(timemat,size(timemat,1)*size(timemat,2), 1);
Ntotvec1 = reshape(BTsum(j).Nmat, size(timemat,1)*size(timemat,2),1);
    end
    
end
hold on
timemat = [];
timetotvec2 = [];
Ntotvec2 = [];
for j = 1:length(BTsum)
    if BTsum(j).N0 == 4
        timemat = repmat(BTsum(j).timevec, size(BTsum(j).Nmat,2),1);
        timemat = timemat';

timetotvec2 = reshape(timemat,size(timemat,1)*size(timemat,2), 1);
Ntotvec2 = reshape(BTsum(j).Nmat, size(timemat,1)*size(timemat,2),1);
    end
end
timemat3 = [];
timetotvec3 = [];
Ntotvec = [];
for j = 1:length(BTsum)
    if BTsum(j).N0 == 10
        timemat = repmat(BTsum(j).timevec, size(BTsum(j).Nmat,2),1);
        timemat = timemat';

timetotvec3 = reshape(timemat,size(timemat,1)*size(timemat,2), 1);
Ntotvec3 = reshape(BTsum(j).Nmat, size(timemat,1)*size(timemat,2),1);
    end

end
plot(timetotvec1, Ntotvec1,'.', 'color', 'r')
    hold on
  plot(timetotvec2, Ntotvec2,'.', 'color', 'g')
    hold on  
plot(timetotvec3, Ntotvec3,'.', 'color', 'b')
    hold on
xlabel('time (hours)')
ylabel ('Number of cells')
title ('Raw data cell number trajectories')
legend ('N_{0}= 2', 'N_{0}= 4', 'N_{0}= 10', 'Location', 'NorthWest')
legend boxoff
xlim ([0 328])


subplot(1,3,2)
for i = 1:length(BTsum)
        if BTsum(i).N0 ==2
        plot(BTsum(i).timevec, BTsum(i).model_g, 'color', 'r', 'Linewidth', 2)
        hold on
        text(BTsum(i).timevec(end-10), BTsum(i).model_g(end-10),['N_{0}=', num2str(BTsum(i).N0)], 'FontSize',8)
        end
        if BTsum(i).N0 ==4
        plot(BTsum(i).timevec, BTsum(i).model_g, 'color', 'g','Linewidth', 2)
        hold on
        text(BTsum(i).timevec(end-10), BTsum(i).model_g(end-10),['N_{0}=', num2str(BTsum(i).N0)], 'FontSize',8)
        end
        if BTsum(i).N0 ==10
        plot(BTsum(i).timevec, BTsum(i).model_g, 'color', 'b', 'Linewidth', 2)
        hold on
        text(BTsum(i).timevec(end-10), BTsum(i).model_g(end-10),['N_{0}=', num2str(BTsum(i).N0)], 'FontSize',8)
        end
end
xlim ([ 0 328])
xlabel ('time (hours)')
ylabel('number of cells')
title ('Mean growth rate model by N0')
subplot(1,3,3)
for i = 1:length(BTsum)
        if BTsum(i).N0 ==2
        plot(BTsum(i).timevec, log(BTsum(i).model_g./BTsum(i).N0), 'color', 'r','Linewidth', 2)
        hold on
        text(BTsum(i).timevec(end-10), log(BTsum(i).model_g(end-10)./BTsum(i).N0),['N_{0}=', num2str(BTsum(i).N0)], 'FontSize',8)
        end
        if BTsum(i).N0 ==4
        plot(BTsum(i).timevec, log(BTsum(i).model_g./BTsum(i).N0), 'color', 'g','Linewidth', 2)
        hold on
        text(BTsum(i).timevec(end-10), log(BTsum(i).model_g(end-10)./BTsum(i).N0),['N_{0}=', num2str(BTsum(i).N0)], 'FontSize',8)
        end
        if BTsum(i).N0 ==10
        plot(BTsum(i).timevec, log(BTsum(i).model_g./BTsum(i).N0), 'color', 'b','Linewidth', 2)
        hold on
        text(BTsum(i).timevec(end-10), log(BTsum(i).model_g(end-10)./BTsum(i).N0),['N_{0}=', num2str(BTsum(i).N0)], 'FontSize',8)
        end
        
end
xlabel ('time (hours)')
ylabel('log(N/N0)')
xlim ([ 0 328])
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
xlim([0 BTsum(2).timevec(end)])
ylabel('number of cells')
title('Cell number trajectories for N0=2 with mean')
subplot (1,2,2)
for j = 1:BTsum(2).num
plot(BTsum(2).timevec, BTsum(2).Nmat(:,j),'b.')
hold on
end
plot(BTsum(2).timevec, BTsum(2).var_t,'g-', 'LineWidth', 3)
xlim([ 0 BTsum(2).timevec(end)])
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