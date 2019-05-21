% This script is going to run the load and fit data codes for the BT-474
% pilot study data at high N0.

% This will then be passed to the "fit to 7 models" code that will fit the
% data to our stochastic models to show that the Allee effect is not
% relevant at high cell densities (hopefully)
close all; clear all; clc

 [N, T] =xlsread('../data/BT-474_Nseed_large.xls');
 Nr = N;

 
%%
%( don't run initially, we already have this separately analyzed)
for i = 1:60 % size matrix, first row, second column
    BT(i).time = round(Nr(:,1));
    BT(i).cellnum = Nr(:,i+1);
%     BT(i).N0 = N2(1, i-sz(1,3)+1);
    BT(i).Nseed = []; % this will have to be updated on each loop
    BT(i).date = '12-4-18';
    BT(i).well = T(1, i+1);
end
for i = 1:60
    cell512 = { 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9', 'B10', 'B11', 'C8', 'C9', 'C10', 'C11',...
        'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7',...
        'D8', 'D9', 'D10', 'D11'};
        cell1024= {'E10', 'E11', 'E6', 'E7', 'E8', 'E9', 'E2', 'E3', 'E4', 'E5','F2', 'F3', 'F4', 'F5'...
       'F6', 'F7', 'F8', 'F9','F10', 'F11', 'G10', 'G11', 'G6', 'G7', 'G8', 'G9', 'G10', 'G2', 'G3', 'G4', 'G5'};
    if contains(BT(i).well, cell512 )
        BT(i).Nseed = 512;
    end
    if contains(BT(i).well, cell1024)
        BT(i).Nseed = 1024;
    end
end
%% Order by N0
Afields = fieldnames(BT);
Acell = struct2cell(BT);
szarray= size(Acell);            % Notice that the this is a 3 dimensional array.
                            % For MxN structure array with P fields, the size
                            % of the converted cell array is PxMxN

% Convert to a matrix
Acell = reshape(Acell, szarray(1), []);      % Px(MxN)

% Make each field a column
Acell = Acell';                         % (MxN)xP

% Sort by 3rd field "N0"
Acell = sortrows(Acell, 3);

% Put back into original cell array format
Acell = reshape(Acell', szarray);

% Convert to Struct
BT = cell2struct(Acell, Afields, 1);
%%  Add color by N0
for j = 1:length(BT)
    N0(j) = BT(j).Nseed;
end
colorsets = varycolor(length(unique(N0)));
uniqN0= unique(N0);

for i = 1:length(BT)
    BT(i).color = [];
    for j = 1:length(uniqN0)
        if BT(i).Nseed==uniqN0(j)
            BT(i).color =colorsets(j,:);
        end
    end
end

%% Fit each trajectory to single exponential growth model
for j = 1:length(BT)
    % want to fit on all but t=0 (we give it this)
    time = [];
    N = [];
    
    time = BT(j).time;
    n = length(time);
    num_params = 1;
    ind0 = find(time == 0);
    % again fit on all but t=0

    N = BT(j).cellnum; % change depend on type of processing we choose
    N0 = BT(j).cellnum(1); % start by using what we think FACS plated
    BT(j).Nbar = mean(N); % average value of volumes over all measured times
    BT(j).N0 = BT(j).cellnum(1);
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

for j = 1:length(BT)
    if ~isempty(BT(j).Rsq)
    if BT(j).Rsq < 0.6
        BT(j).badfit =1;
    else BT(j).badfit =0;
    end
    end
%     if BT(j).N0 > 1.2*BT(j).Nseed ||BT(j).N0 < 0.8*BT(j).Nseed
%         BT(j).badfit =1;
%     end
  
     
end
%% Separate by N0
uniqN0 = [];
N0list = [];
for j = 1:length(BT)
    if ~isnan(BT(j).Nseed)
    N0list =vertcat(N0list, BT(j).Nseed);
    end
    
end
uniqN0=unique(N0list);
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
%%
 for i = 1:length(uniqN0) % number of unique seed numbers
    for j = 1:length(BT)
        N0 = BT(j).Nseed;
        time = BT(j).time;
        if BT(j).Nseed == uniqN0(i)
            BTsum(i).color = BT(j).color;
            if BT(j). badfit == 0
            BTsum(i).num = BTsum(i).num +1;
            BTsum(i).tmin = vertcat(BTsum(i).tmin, BT(j).time(end));
            BTsum(i).gtot = vertcat(BTsum(i).gtot, BT(j).g);
            end
         end
        end
    end

 % need to make a uniform time vector for each initial cell number
 % this will likely change for each cell number ( higher N, less time)
 for i = 1:length(uniqN0)
     BTsum(i).timevec =Nr(:,1); % makes a uniform time vector
     BTsum(i).Nmat = [];
 end

 for i= 1:length(uniqN0)
    BTsum(i).var_g = (std(BTsum(i).gtot))^2;
    BTsum(i).norm_var = BTsum(i).var_g/BTsum(i).N0;
    BTsum(i).mean_g = mean(BTsum(i).gtot);
    
    BTsum(i).model_g= singleexpmodel(BTsum(i).mean_g,BTsum(i).mean_N0, BTsum(i).timevec);% model fit for plotting           
 end
 
% This loop makes an a matrix of each N trajectory sampled uniformly in
% time which will be used for the stochastic parameter estimation
%%
for j = 1:length(BTsum)
    tsamp = Nr(:,1);
 for k = 1:length(BT)
     % note change depending on whether or not we want to include badfits
     if BT(k).Nseed == BTsum(j).N0 && BT(k).badfit ==0
        Nsamp=BT(k).cellnum;
     else
         Nsamp = [];
     end
        BTsum(j).Nmat= horzcat(BTsum(j).Nmat, Nsamp);
        Nsamp = [];
     end
end
for i= 1:length(uniqN0)
    BTsum(i).var_g = (std(BTsum(i).gtot))^2;
    BTsum(i).norm_var = BTsum(i).var_g/BTsum(i).N0;
    BTsum(i).mean_g = mean(BTsum(i).gtot);
    BTsum(i).mean_N0 = mean(BTsum(i).Nmat(1,:));
    BTsum(i).model_g= singleexpmodel(BTsum(i).mean_g,BTsum(i).mean_N0, BTsum(i).timevec);% model fit for plotting           
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

figure;
for j = 1:length(BTsum)
    for i = 1:BTsum(j).num
    plot(BTsum(j).timevec, BTsum(j).Nmat(:,i), 'color', BTsum(j).color)
    plot(BTsum(j).timevec, BTsum(j).model_g,'k-')
    hold on
    end
    text(BTsum(j).timevec(end-20), BTsum(j).Nmat(end-20,1), ['N_{seed}= ', num2str(BTsum(j).N0)])
    text(BTsum(j).timevec(end-10), BTsum(j).Nmat(end-10,1), ['g= ', num2str(BTsum(j).mean_g),' +/- ', num2str(1.96*sqrt(BTsum(j).var_g))])
end

xlim([ 0 Nr(end,1)])
xlabel ('time (hours)')
ylabel('Number of cells')
title('Growth trajectories for normal seeding density wells')

figure;
subplot(1,3,1)
for j = 1:length(BTsum)
    for i = 1:BTsum(j).num
    plot(BTsum(j).timevec, BTsum(j).Nmat(:,i), 'color', BTsum(j).color)
    hold on
    plot(BTsum(j).timevec, BTsum(j).model_g,'k-')
    
    end
    text(BTsum(j).timevec(end-20), BTsum(j).Nmat(end-20,1), ['N_{seed}= ', num2str(BTsum(j).N0)])
    %text(BTsum(j).timevec(end-10), BTsum(j).Nmat(end-10,1), ['g= ', num2str(BTsum(j).mean_g),' +/- ', num2str(1.96*sqrt(BTsum(j).var_g))])
end
xlim([ 0 Nr(end,1)])
xlabel ('time (hours)', 'FontSize', 14)
ylabel('Number of cells', 'FontSize',14)
%title('A         Growth trajectories', 'FontSize', 14)


subplot(1,3,2)
for j = 1:length(BTsum)
    plot(BTsum(j).timevec, log(BTsum(j).mu_t./BTsum(j).mean_N0), 'color', BTsum(j).color, 'LineWidth',3)
    hold on
end
 %text(BTsum(j).timevec(end-20), BTsum(j).Nmat(end-20,1), ['N_seed= ', num2str(BTsum(j).N0)])
 text(BTsum(1).timevec(end-50), log(BTsum(1).mu_t(end-50)./BTsum(1).mean_N0), ['g_{Nseed=512}= ', num2str(round(BTsum(1).mean_g,4)),' +/- ', num2str(round(1.96*sqrt(BTsum(1).var_g),5))])
 text(BTsum(2).timevec(end-30), log(BTsum(2).mu_t(end-30)./BTsum(1).mean_N0), ['g_{Nseed=1024}= ', num2str(round(BTsum(2).mean_g,4)),' +/- ', num2str(round(1.96*sqrt(BTsum(2).var_g),5))])

xlim([ 0 Nr(end,1)])
xlabel ('time (hours)','FontSize', 14)
legend('N_{seed}=512', 'N_{seed}=1024', 'Location', 'NorthWest')
legend boxoff
ylabel('log(N/N_{0})','FontSize', 14)
%title('B         Normalized growth', 'FontSize', 14)

% average growth rate versus cell number
subplot(1,3,3)
for i =1:length(uniqN0)
    errorbar((BTsum(i).N0), BTsum(i).mean_g, 1.96*sqrt(BTsum(i).var_g), 'o', 'color', BTsum(i).color,'LineWidth', 2)
    %set(gca,'XScale','log')
    hold on
end
xlabel ('N_{0}','FontSize', 14)
ylabel ('mean growth rate','FontSize', 14)
xlim ([ 256, 1024+256])
ylim([0.01, 0.013])
%title ('C        Average growth rate vs. N_{0}', 'FontSize',14)
%%

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
%% Find variance in N0

for j = 1:length(uniqN0)
    N0s = BTsum(j).Nmat(1,:);
    stdevN0s = std(N0s);
    BTsum(j).V0 = stdevN0s.^2;
end
figure;
subplot(1,2,1)
for j = 1:length(uniqN0)
    plot(BTsum(j).timevec, BTsum(j).mu_t)
    hold on
end
subplot(1,2,2)
for j = 1:length(uniqN0)
    plot(BTsum(j).timevec, BTsum(j).var_t)
    hold on
end
%% Save the two fitted data structures

% this saves the raw data structure for all imported BT-474 data sets
save('../out/BTfitps.mat', 'BT')
save('../out/BTsumfitps.mat', 'BTsum')
