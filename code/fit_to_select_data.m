% This script selects a subset of the BT-474 data and fits each initial
% condition to the b-d model, and then fits all N0 and data to each of the
% 7 models 

clear all; close all; clc

S = load('../out/BTfit.mat');
BT= S.BT;

S = load('../out/BTsumfit.mat');
BTsum= S.BTsum;
%% Plot some example trajectories from BT sum
figure;
timemat = [];
timetotvec2 = [];
Ntotvec2 = [];
N0list = [2 4 10];
for j = 1:length(BTsum)
    if BTsum(j).N0 == N0list(1)
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
    if BTsum(j).N0 == N0list(2)
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
    if BTsum(j).N0 == N0list(3)
        timemat = repmat(BTsum(j).timevec, size(BTsum(j).Nmat,2),1);
        timemat = timemat';

timetotvec3 = reshape(timemat,size(timemat,1)*size(timemat,2), 1);
Ntotvec3 = reshape(BTsum(j).Nmat, size(timemat,1)*size(timemat,2),1);
    end

end
figure;
subplot(1,3,1)
plot(timetotvec1, Ntotvec1,'.', 'color', 'r')
    hold on
  plot(timetotvec2, Ntotvec2,'.', 'color', 'g')
    hold on  
plot(timetotvec3, Ntotvec3,'.', 'color', 'b')
    hold on
xlabel('time (hours)', 'FontSize', 16)
ylabel ('Cell number (N(t)) ', 'FontSize', 16)
title ('Raw data cell number trajectories', 'FontSize', 14)
legend (['N_{0}= ',num2str(N0list(1))], ['N_{0}= ',num2str(N0list(2))], ['N_{0}= ',num2str(N0list(3))], 'Location', 'NorthWest', 'FontSize', 14)
legend boxoff
xlim ([0 328])


subplot(1,3,2)
plot(BTsum(N0list(1)).timevec, BTsum(N0list(1)).mu_t, 'r-', 'LineWidth',3)
hold on
plot(BTsum(N0list(2)).timevec, BTsum(N0list(2)).mu_t, 'g-', 'LineWidth',3)
plot(BTsum(N0list(3)).timevec, BTsum(N0list(3)).mu_t, 'b-', 'LineWidth',3)
%text(ttot(end-20,i), N1tot(end-20,i),['Ni=', num2str(N1_init(i))]','HorizontalAlignment','left','VerticalAlignment','bottom','color','k' )
legend (['N_{0}= ',num2str(N0list(1))], ['N_{0}= ',num2str(N0list(2))], ['N_{0}= ',num2str(N0list(3))], 'Location', 'NorthWest', 'FontSize', 14)
legend boxoff
xlabel('time (hours)', 'FontSize', 16)
ylabel('Mean cell number <N(t)>', 'FontSize', 16)
xlim([ 0 328])
%title('dN/dt=gN solution')
title('Average growth', 'FontSize', 14)

subplot(1,3,3)
plot(BTsum(N0list(1)).timevec, log(BTsum(N0list(1)).mu_t./N0list(1)), 'r-', 'LineWidth',3)
hold on
plot(BTsum(N0list(2)).timevec, log(BTsum(N0list(2)).mu_t./N0list(2)), 'g-', 'LineWidth',3)
plot(BTsum(N0list(3)).timevec, log(BTsum(N0list(3)).mu_t./N0list(3)), 'b-', 'LineWidth',3)
text(BTsum(N0list(1)).timevec(30), log(BTsum(N0list(1)).mu_t(35))./N0list(1), ['g_{avg N0=2}=', num2str(round(BTsum(N0list(1)).mean_g, 5)),' +/- ' num2str(round(1.96*sqrt(BTsum(N0list(1)).var_g),5))])
text(BTsum(N0list(2)).timevec(35), log(BTsum(N0list(2)).mu_t(40))./N0list(2), ['g_{avg N0=5}=', num2str(round(BTsum(N0list(2)).mean_g, 5)),' +/- ' num2str(round(1.96*sqrt(BTsum(N0list(2)).var_g),5))])
text(BTsum(N0list(3)).timevec(40), log(BTsum(N0list(3)).mu_t(45))./N0list(3), ['g_{avg N0=12}=', num2str(round(BTsum(N0list(3)).mean_g, 5)),' +/- ' num2str(round(1.96*sqrt(BTsum(N0list(3)).var_g),5))])
legend(['N_{0}= ',num2str(N0list(1))], ['N_{0}= ',num2str(N0list(2))], ['N_{0}= ',num2str(N0list(3))], 'Location', 'NorthWest', 'FontSize', 14)
legend boxoff
xlabel('time (hours)', 'FontSize', 16)
ylabel('log(N(t)/N_{0})', 'FontSize', 16)
xlim([ 0 328])
%title('dN/dt=gN solution')
title('Normalized growth rate', 'FontSize', 14)
%% Perform multiparametric t-test
% Make grouped data vector
g2 = BTsum(2).gtot;
g4pr = BTsum(5).gtot;
g4 = vertcat(g4pr, NaN*ones(length(g2)-length(g4pr),1))
g10pr = BTsum(12).gtot;
g10 = vertcat(g10pr, NaN*ones(length(g2)-length(g10pr),1))
y = horzcat(g2, g4, g10);
[h,p,stats] = anova1(y)
%%
figure;
 
plot(timetotvec3, Ntotvec3,'.', 'color', 'b')
hold on
plot(timetotvec2, Ntotvec2,'.', 'color', 'g')
 plot(timetotvec1, Ntotvec1,'.', 'color', 'r')
 
xlabel('time (hours)','FontSize', 18)
ylabel ('Number of cells','FontSize', 18)
title ('Raw data','FontSize', 18)
legend ('N_{0}= 10', 'N_{0}= 4', 'N_{0}= 2', 'Location', 'NorthWest')
legend boxoff
xlim ([0 328])

 
 figure;
 subplot(1,2,1)
plot(BTsum(2).timevec, BTsum(2).mu_t, 'r-', 'LineWidth',3)
hold on
plot(BTsum(4).timevec, BTsum(4).mu_t, 'g-', 'LineWidth',3)
plot(BTsum(10).timevec, BTsum(10).mu_t, 'b-', 'LineWidth',3)
%text(ttot(end-20,i), N1tot(end-20,i),['Ni=', num2str(N1_init(i))]','HorizontalAlignment','left','VerticalAlignment','bottom','color','k' )
% legend ('N_{0}= 2', 'N_{0}= 4', 'N_{0}= 10', 'Location', 'NorthWest')
% legend boxoff
xlabel('time (hours)','FontSize', 18)
ylabel('Mean','FontSize', 18)
xlim([ 0 328])
%title('dN/dt=gN solution')
title('Mean', 'FontSize', 18)

subplot(1,2,2)
plot(BTsum(2).timevec, BTsum(2).var_t, 'r-', 'LineWidth',3)
hold on
plot(BTsum(4).timevec, BTsum(4).var_t, 'g-', 'LineWidth',3)
plot(BTsum(10).timevec, BTsum(10).var_t, 'b-', 'LineWidth',3)
%text(ttot(end-20,i), N1tot(end-20,i),['Ni=', num2str(N1_init(i))]','HorizontalAlignment','left','VerticalAlignment','bottom','color','k' )
% legend ('N_{0}= 2', 'N_{0}= 4', 'N_{0}= 10', 'Location', 'NorthWest')
% legend boxoff
xlabel('time (hours)', 'FontSize', 18)
ylabel('Variance','FontSize', 18)
xlim([ 0 328])
%title('dN/dt=gN solution')
title('Variance','FontSize', 18)
%%
subplot(3,3,3)
i=1
text(ttot(55,i)+5, log(N1tot(55, i)./(N1_init(i))), ['N_{0}=', num2str(N1_init(i)),', ', num2str(N1_init(i+1)),', &', num2str(N1_init(i+2))])
hold on
for i = 1:length(N1_init)
plot(ttot(:,i), log(N1tot(:,i)./(N1_init(i))), 'c-', 'LineWidth',3)
%text(ttot(55,i), log(N1tot(55, i)./(N1_init(i))), ['N_{0}=', num2str(N1_init(i)),',', num2str(N1_init(i+1),'&', num2str(N1_init(i+2))])
hold on
end
ylim([ 0 5])
xlabel('time (hours)')
ylabel('log(N/N0)')
%title('Normalized growth rate dN/dt=g(N)')
title('Normalized growth rate')
% Make a plot of instantaneous growth rate in time
for i = 1:length(N1_init)
    for j = 1:length(tsamp)-1
    percapg(j,i) = (N1tot(j+1,i)-N1tot(j,i))./(N1tot(j+1,i)*tint);
    end
end

subplot(3,3,2)
for i = 2%1:length(N1_init)
    %subplot(1,length(N1_init),i)
   plot(ttot(1:end-1,i), percapg(:,i), 'c-', 'LineWidth',3)
   text(ttot(55,i), percapg(75,i)+.001, ['N_{0}=',num2str(N1_init(i) )])
   hold off
   xlim([ 0 tsamp(end)])
   ylim([ 0 0.05])
   xlabel ('time (hours)')
   ylabel('per capita growth rate')
   title(['Per capita growth rate'])
   %title(['Per capita growth rate for dN/dt =gN, N_{0}=',num2str(N1_init(i)),', A=', num2str(A)])
end



%% Fit to a subset of the data, here we pick 2, 4, and 10-

for j =1:length(N0list)
for i = 1:length(BTsum)
    if BTsum(i).N0==N0list(j)
        g0list(j) = BTsum(i).mean_g;
        varglist(j)= BTsum(i).var_g;
% For now, start by assuming sigma = 1 and fit for b and d
    N0 = BTsum(i).N0;
    V0 = 0;
    t = [];
    negLLguesst = [];
    negLLfitt = [];
    mu_fitJst = [];
    V_fitJst = [];
    tsamp=0:4:324; 
    N= BTsum(i).num;
    mudata(:,j) = BTsum(i).mu_t(2:end);
    vardata(:,j) = BTsum(i).var_t(2:end);
    
    end 
end
end
%% Add in code to downsample data
ind = 1:1:length(mudata);
mudatanew = mudata(ind,:);
vardatanew = vardata(ind, :);
tsampfit = tsamp(2:end);
tsampnew = tsampfit(ind);

%% confirm it's same
figure;
for j = 1:length(N0list)
plot(tsampfit, mudata(:,j), 'r*')
hold on
plot(tsampnew, mudatanew(:,j), 'b.')
end
xlabel('time (hours)')
ylabel('cell numer')
title('down sampled data for fitting')
tsamp = [0 tsampnew];
mudata = mudatanew;
vardata = vardatanew;
%% Add in code to plot trajectories & log normalized trajectories with growth rate

%%
for j =1:length(N0list)
for i = 1:length(BTsum)
    if BTsum(i).N0==N0list(j) 
    N= BTsum(i).num;
    modelcode = 1;
    % Initial guess
    bguess = (log(mudata(end-2,j))-log(mudata(end-4,j)))/(tsamp(end-2)-tsamp(end-4)) + 0.002; 
    dguess = 0.002;
    theta = [bguess,dguess];
    
    % fit to birth-death model individually for each set of N0
    [pbestbd(j,:), BICbd(j), AICbd(j),negLLfit(j)] = ML_fit(theta, tsamp, N0list(j), N,mudata(:,j), vardata(:,j), modelcode);
    
    BTsum(i).pbestbd = pbestbd(j,:);
    BTsum(i).negLLfitbd = negLLfit(j);
end
end
end

%% Plot fitting results

figure;
for j = 1:length(N0list)
    modelfun_mu =@(p)gen_model_mu(p,tsamp, N0list(j), modelcode);
    modelfun_V = @(p)gen_model_var(p, tsamp, N0list(j), modelcode);
    subplot(2,3,j)
    hold on 
    plot(tsamp(2:end), mudata(:,j),'r*')
    plot(tsamp(2:end), modelfun_mu(pbestbd(j,:)), 'k-', 'LineWidth',2)
    xlabel ('time (hours)', 'FontSize', 14)
    ylabel ('Mean cell number (N(t))', 'FontSize',14)
    legend ('mean data', 'fit', 'Location', 'NorthWest', 'FontSize',12)
    legend boxoff
    title(['N_{0}= ',num2str(N0list(j))], 'FontSize',14)%,', b= ', num2str(round(pbestbd(j,1),6)),', d= ', num2str(round(pbestbd(j,2),6))])
    xlim([ 0 328])
    
    subplot(2,3,j+3)
    hold on
    plot(tsamp(2:end), vardata(:,j),'g*')
    plot(tsamp(2:end), modelfun_V(pbestbd(j,:)), 'k-', 'LineWidth',2)
    xlabel ('time (hours)', 'FontSize', 14)
    ylabel ('Variance in cell number (\Sigma(t))', 'FontSize',14)
    legend ('variance in data', 'fit', 'Location', 'NorthWest', 'FontSize',12)
    legend boxoff
    title(['N_{0}= ',num2str(N0list(j))], 'FontSize',14)%,', b= ', num2str(round(pbestbd(j,1),6)),', d= ', num2str(round(pbestbd(j,2),6))])
    xlim([ 0 328])
end
%%
figure;
errorbar(N0list, g0list, sqrt(varglist)*1.96./2, 'ko', 'LineWidth', 3)
xlabel('N_{0}')
ylabel('average growth rate')
title('Average growth rate by N_{0}')
xlim ([ 1 13])
%ylim ([ 0 0.01])
%%
figure;
subplot(1,2,1)
plot(N0list, pbestbd(:,1), 'r*','LineWidth', 3)
hold on
plot(N0list, pbestbd(:,2), 'b*','LineWidth', 3)
xlabel('N_{0}', 'FontSize',18)
ylabel('parameter', 'FontSize',16)
legend ('b', 'd','FontSize',14)
legend boxoff
xlim ([ 0 13])
ylim ([ 0 0.015])
title('b & d parameter estimates by N_{0}', 'FontSize',16)

subplot(1,2,2)
plot(N0list, pbestbd(:,1) - pbestbd(:,2), 'k*', 'LineWidth',3)
xlabel('N_{0}')
ylabel('g= b-d')
xlim ([ 0 13])
ylim ([ 0.006 0.01])
title('g = b-d parameter estimates by N_{0}')

%%  Compute the profile likelihoods for each parameter for each initial cell number
% Use these to determine parameter certainty and find the 95 % CI on the
% b and d parameter estimates

% basic idea: hold 1 parameter constant at a certain value, allow search
% algorithm to find the best possible likelihood with that parameter
% constant, and report likelihood value.

% For first pass, we will just use fminsearch for each step through
% parameters

%WRITE THIS INTO A FUNCTION!
num_params = 2;

 % Check fxn
 for j =1 % first for N0=2 sicne d is so low
for i = 1:length(BTsum)
    if BTsum(i).N0==N0list(j)
        fbest = BTsum(i).negLLfitbd;
        factor =0.01 ;
        numpoints = 10;
        params = BTsum(i).pbestbd;
        N0 = N0list(j);
        N= BTsum(i).num;
        [profiles] = profile_likelihood(params, tsamp, N0, N,mudata(:,j), vardata(:,j), modelcode, factor, numpoints);
        BTsum(i).profiles = profiles;
        threshold = chi2inv(0.95,length(params))/2 + fbest;
        BTsum(i).threshold = threshold;
        ilo1=find(profiles(:,2,1)<threshold);
        if ilo1(end) ==numpoints*2 ||ilo1(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CI(1,:) = [NaN, NaN];
        else
             CI(1,:) = [profiles(ilo1(1),1,1), profiles(ilo1(end),1,1)];
        end
        ilo2 = find(profiles(:,2,2)<threshold);
        if ilo2(end) ==numpoints*2 ||ilo1(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CI(2,:)= [NaN NaN];
        else
             CI(2,:) = [profiles(ilo2(1),1,2), profiles(ilo2(end),1,2)];
        end
        BTsum(i).CI = CI;
    end
end
 end
 %% Then for N0=4,10 since 
 for j =3%:length(N0list)
for i = 1:length(BTsum)
    if BTsum(i).N0==N0list(j)
        fbest = BTsum(i).negLLfitbd;
        factor = 0.03;
        numpoints = 10;
        params = BTsum(i).pbestbd;
        N0 = N0list(j);
        N= BTsum(i).num;
        [profiles] = profile_likelihood(params, tsamp, N0, N,mudata(:,j), vardata(:,j), modelcode, factor, numpoints);
        BTsum(i).profiles = profiles;
        threshold = chi2inv(0.95,length(params))/2 + fbest;
        BTsum(i).threshold = threshold;
        ilo1=find(profiles(:,2,1)<threshold);
        if ilo1(end) ==numpoints*2 ||ilo1(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CI(1,:) = [NaN, NaN];
        else
             CI(1,:) = [profiles(ilo1(1),1,1), profiles(ilo1(end),1,1)];
        end
        ilo2 = find(profiles(:,2,2)<threshold);
        if ilo2(end) ==numpoints*2 ||ilo1(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CI(2,:)= [NaN NaN];
        else
             CI(2,:) = [profiles(ilo2(1),1,2), profiles(ilo2(end),1,2)];
        end
        BTsum(i).CI = CI;
    end
end
 end
 %%
for j = 1:length(N0list)
    for i = 1:length(BTsum)
    if BTsum(i).N0==N0list(j)
    CImat(:,:,j) = BTsum(i).CI;
    end
    end
end
 
 %%
 figure;
 for j = 1:length(N0list)
     for i = 1:length(BTsum)
         if BTsum(i).N0 == N0list(j)
         subplot(2,3,j)
        plot(BTsum(i).profiles(:,1,1), BTsum(i).profiles(:,2,1), 'LineWidth', 2)
        hold on
        plot(BTsum(i).pbestbd(1),BTsum(i).negLLfitbd,'r*','LineWidth',2)
        plot([BTsum(i).profiles(1,1,1) BTsum(i).profiles(end,1,1)],[BTsum(i).threshold BTsum(i).threshold],'r--')
        xlabel(' profiled parameter values')
        ylabel('Negative Log Likelihood')
        %xlim([ BTsum(i).CI(1,1)-0.1*BTsum(i).CI(1,1) BTsum(i).CI(1,2)+0.1*BTsum(i).CI(1,2)])
        title(['N_{0}=', num2str(N0list(j)),' Profile Likelihood of b'])

        subplot(2,3,j+3)
        plot(BTsum(i).profiles(:,1,2), BTsum(i).profiles(:,2,2), 'LineWidth',2)
        hold on
        plot(BTsum(i).pbestbd(2),BTsum(i).negLLfitbd,'r*','LineWidth',2)
        plot([BTsum(i).profiles(1,1,2) BTsum(i).profiles(end,1,2)],[BTsum(i).threshold BTsum(i).threshold],'r--')
        xlabel(' profiled parameter values')
        ylabel('Negative Log Likelihood')
        %xlim([ BTsum(i).CI(2,1)-0.1*BTsum(i).CI(2,1) BTsum(i).CI(2,2)+0.1*BTsum(i).CI(2,2)])
        title(['N_{0}=', num2str(N0list(j)),' Profile Likelihood of d'])
         end
     end
 end
 
 figure;
 for j = 1:length(N0list)
     for i = 1:length(BTsum)
         if BTsum(i).N0 == N0list(j)
         subplot(2,3,j)
        plot(BTsum(i).profiles(:,1,1), exp(-BTsum(i).profiles(:,2,1)), '--','Linewidth',2)
        hold on
        plot(BTsum(i).pbestbd(1),exp(-BTsum(i).negLLfitbd),'r*','LineWidth',2)
        %plot([BTsum(i).profiles(1,1,1) BTsum(i).profiles(end,1,1)],[BTsum(i).threshold BTsum(i).threshold],'r--')
        xlabel(' profiled b value')
        ylabel('Likelihood')
        %xlim([ BTsum(i).CI(1,1)-0.1*BTsum(i).CI(1,1) BTsum(i).CI(1,2)+0.1*BTsum(i).CI(1,2)])
        title(['N_{0}=', num2str(N0list(j)),' Profile Likelihood of b'])
        subplot(2,3,j+3)
        plot(BTsum(i).profiles(:,1,2), exp(-BTsum(i).profiles(:,2,2)), '--','Linewidth',2)
        hold on
        plot(BTsum(i).pbestbd(2),exp(-BTsum(i).negLLfitbd),'r*','LineWidth',2)
        %plot([BTsum(i).profiles(1,1,1) BTsum(i).profiles(end,1,1)],[BTsum(i).threshold BTsum(i).threshold],'r--')
        xlabel(' profiled d value')
        ylabel('Likelihood')
        %xlim([ BTsum(i).CI(1,1)-0.1*BTsum(i).CI(1,1) BTsum(i).CI(1,2)+0.1*BTsum(i).CI(1,2)])
        title(['N_{0}=', num2str(N0list(j)),' Profile Likelihood of d'])
         end
     end
 end
 %%
figure;
for j = 1:length(N0list)
     for i = 1:length(BTsum)
    if BTsum(i).N0==N0list(j)
errorbar(BTsum(i).N0, BTsum(i).pbestbd(1),0.5.*(BTsum(i).CI(1,1)-pbestbd(1)), 0.5*(BTsum(i).CI(1,2)-pbestbd(1)),'r*','LineWidth', 3)
hold on
errorbar(BTsum(i).N0, BTsum(i).pbestbd(2), 0.5.*(BTsum(i).CI(2,1)-pbestbd(2)),0.5*( BTsum(i).CI(2,2)-pbestbd(2)),'b*','LineWidth', 3)
% plot(BTsum(i).N0, BTsum(i).CI(1,1), 'k.')
% plot(BTsum(i).N0, BTsum(i).CI(1,2), 'k.')
% plot(BTsum(i).N0, BTsum(i).pbestbd(2), 'b*','LineWidth', 3)
% plot(BTsum(i).N0, BTsum(i).CI(2,1), 'k.')
% plot(BTsum(i).N0, BTsum(i).CI(2,2), 'k.')

    end
     end
end
xlabel('N_{0}')
ylabel('parameter value')
legend ('birth rate', 'death rate', 'Location', 'NorthWest')
legend boxoff
xlim ([ 0 11])
% ylim ([ 0 0.015])
%title('b & d parameter estimates by N_{0}')

%% Load in selected mudatavec
% Want to calculate BIC and AIC value from the fitting
% Need to make one big matrix of mu and vardata with corresponding time vecs in order
% Load in raw data, not down sampled
clear all; close all; clc
%%
mudatavec = [];
vardatavec = [];
timevec = [];
mudatavecall = [];
vardatavecall = [];
timevecall = [];
N0vec = [];
numsamps = 0;

S = load('../out/BTsumfit.mat');
BTsum= S.BTsum;
N0list = [ 2 4 10];


for j = 1:length(N0list)
    for i = 1:length(BTsum)
    if BTsum(i).N0==N0list(j)
        mudata = BTsum(i).mu_t(2:end);
        ind = 1:1:length(mudata);
        mudatasmall= mudata(ind);
    mudatavec = vertcat(mudatavec, mudatasmall);
        vardata = BTsum(i).var_t(2:end);
        vardatasmall = vardata(ind);
    vardatavec = vertcat(vardatavec, vardatasmall);
        time = BTsum(i).timevec(2:end)';
        timesmall = time(ind);
    timevec = vertcat(timevec, timesmall);
  
   
    vardatavecall = vertcat(vardatavecall, BTsum(i).var_t(2:end));
    mudatavecall = vertcat(mudatavecall, BTsum(i).mu_t(2:end));
    timevecall = vertcat(timevecall, BTsum(i).timevec(2:end)');
    numsamps = BTsum(i).num + numsamps;

end
    end
end

%% confirm it's same
figure;
for j = 1:length(N0list)
plot(timevecall, mudatavecall, 'r*')
hold on
plot(timevec, mudatavec, 'b.')
end
xlabel('time (hours)')
ylabel('cell numer')
title('down sampled data for fitting')

%% Perform fitting: fit all data to a single b & d
% For each fitting, we want to write functions in terms of the functions
% that run the model we are fitting to forward. The gen_model_mu runs the
% system of ODEs for the parameters and generates an expected model mean
% and variance for the input p parameters
% information

modelcode = 1;
tsamp = 0:4:332; % set sample time (if data is down sampled change accordingly
tsampfit = tsamp(2:end);
tsampnew = tsampfit(ind);
tsamp = [0 tsampnew];
%% Initial guess
bguess = (log(mudatavec(end-2))-log(mudatavec(end-4)))/(timevec(end-2)-timevec(end-4)) + 0.002; 
dguess = 0.002;
theta = [bguess,dguess];
Ninit = N0list;
N= numsamps;
[pbestbd, BICbd, AICbd,negLLfitbd] = ML_fit(theta, tsamp, Ninit, N,mudatavec, vardatavec, modelcode)

% need to define th inline functions that run forward the mean and the
% variance
modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);

%% Plot the data compared to best fit model
figure;
subplot(1,2,1)
hold on
plot(timevec, mudatavec, 'r*')
for j = 1:length(N0list)
    modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit(j), modelcode);
plot(tsamp(2:end), modelfun_mu(pbestbd), 'b-', 'LineWidth', 2) % test model expected
end
xlabel('time (hours)')
ylabel('Mean cell number')
legend('data', 'best fit b-d')
legend boxoff
title(['Mean: simple birth-death model, BIC= ', num2str(BICbd)])
xlim([ 0 332])
subplot(1,2,2)
hold on
plot(timevec, vardatavec, 'g*')
for j = 1:length(N0list)
    modelfun_V =@(p)gen_model_var(p,tsamp, Ninit(j), modelcode);
plot(tsamp(2:end), modelfun_V(pbestbd), 'b-', 'LineWidth', 2) % test model expected
end
xlabel('time (hours)')
ylabel(' Variance')
legend('data', 'best fit b-d')
legend boxoff
title(['Variance: simple birth-death model, BIC= ', num2str(BICbd)])
xlim([ 0 332])
%% Profile likelihood of b-d model 
% Use these to determine parameter certainty and find the 95 % CI on the
% b and d parameter estimate from all data

num_params = 2;

 % Check fxn

        fbest = negLLfitbd;
        factor = 0.05;
        numpoints = 10;
        params = pbestbd;
        [profilesbd] = profile_likelihood(params, tsamp, N0list, numsamps ,mudatavec, vardatavec, modelcode, factor, numpoints);
        threshold = chi2inv(0.95,length(params))/2 + fbest;
        ilo1=find(profilesbd(:,2,1)<threshold);
     
        if ilo1(end) ==numpoints*2 ||ilo1(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIbd(1,:) = [NaN, NaN];
        else
             CIbd(1,:) = [profilesbd(ilo1(1),1,1), profilesbd(ilo1(end),1,1)];
        end
        ilo2 = find(profilesbd(:,2,2)<threshold);
        if ilo2(end) ==numpoints*2 ||ilo1(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIbd(2,:)= [NaN NaN];
        else
             CIbd(2,:) = [profilesbd(ilo2(1),1,2), profilesbd(ilo2(end),1,2)];
        end

    % this output a profiles 3D matrix whose first z dimension contains the
    % profiled parameters and the corresponding negLL
    % the CI contains the upper and lower CI on first the b parameter than
    % the d parameter

%% Plot profiles to view CI

 figure;
         subplot(1, 2,1)
        plot(profilesbd(:,1,1), profilesbd(:,2,1), '--', 'LineWidth', 2)
        hold on
        plot(pbestbd(1),(negLLfitbd),'r*','LineWidth',2)
        plot([profilesbd(1,1,1) profilesbd(end,1,1)],[(threshold) (threshold)],'r--')
        xlabel(' profiles param values')
        ylabel('Cost Function Value')
        title(['Profile Likelihood of b CI = [',num2str(CIbd(1,:)),']'])
        
        subplot(1,2,2)
        plot(profilesbd(:,1,2), profilesbd(:,2,2), '--', 'LineWidth', 2)
        hold on
        plot(pbestbd(2),negLLfitbd,'r*','LineWidth',2)
        plot([profilesbd(1,1,2) profilesbd(end,1,2)],[(threshold) (threshold)],'r--')
        xlabel(' profiles param values')
        ylabel('Cost Function Value')
        title([' Profile Likelihood of d CI = [',num2str(CIbd(2,:)),']'])
    

%% Perform fitting: fit all data to strong Allee model on birth
% For each fitting, we want to write functions in terms of the functions
% that run the model we are fitting to forward. The gen_model_mu runs the
% system of ODEs for the parameters and generates an expected model mean
% and variance for the input p parameters
% information

modelcode = 2;

% Initial guess
bguess = (log(mudatavec(end-2))-log(mudatavec(end-4)))/(timevec(end-2)-timevec(end-4)) + 0.0002; 
dguess = 0.0002;
Aguess = 2;
theta = [bguess,dguess, Aguess];
Ninit = N0list;
N= numsamps;
[pbeststrAb, BICstrAb, AICstrAb,negLLfitstrAb] = ML_fit(theta, tsamp, Ninit, N,mudatavec, vardatavec, modelcode);

% need to define th inline functions that run forward the mean and the
% variance
modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);

%% Plot the data compared to best fit model fr strong Allee on birth
modelcode = 2;

figure;
subplot(1,2,1)
hold on
plot(timevec, mudatavec, 'r*')
for j = 1:length(N0list)
    modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit(j), modelcode);
plot(tsamp(2:end), modelfun_mu(pbeststrAb), 'b-', 'LineWidth', 2) % test model expected
end
xlabel('time (hours)')
ylabel('Mean cell number')
legend('data', 'best fit strong Allee (b)')
legend boxoff
title(['Mean: strong Allee on birth model, BIC= ', num2str(BICstrAb)])
xlim([ 0 328])
subplot(1,2,2)
hold on
plot(timevec, vardatavec, 'g*')
for j = 1:length(N0list)
    modelfun_V =@(p)gen_model_var(p,tsamp, Ninit(j), modelcode);
plot(tsamp(2:end), modelfun_V(pbeststrAb), 'b-', 'LineWidth', 2) % test model expected
end
xlabel('time (hours)')
ylabel(' Variance')
legend('data', 'best fit strong Allee (b)')
legend boxoff
title(['Variance: strong Allee on birth model, BIC= ', num2str(BICstrAb)])
xlim([ 0 328])
%% Profile likelihood of strong Allee on birth model
% Use these to determine parameter certainty and find the 95 % CI on the
% b and d parameter estimate from all data

num_params = 3;

 % Check fxn

        fbest = negLLfitstrAb;
        factor = 0.05;
        numpoints = 10;
        params = pbeststrAb;
        [profilesstrAb] = profile_likelihood(params, tsamp, N0list, numsamps ,mudatavec, vardatavec, modelcode, factor, numpoints);
        threshold = chi2inv(0.95,length(params))/2 + fbest;
        ilo1=find(profilesstrAb(:,2,1)<threshold);
        if ilo1(end) ==numpoints*2 ||ilo1(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIstrAb(1,:) = [NaN, NaN];
        else
             CIstrAb(1,:) = [profilesstrAb(ilo1(1),1,1), profilesstrAb(ilo1(end),1,1)];
        end
        ilo2 = find(profilesstrAb(:,2,2)<threshold);
        if ilo2(end) ==numpoints*2 ||ilo2(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIstrAb(2,:)= [NaN NaN];
        else
             CIstrAb(2,:) = [profilesstrAb(ilo2(1),1,2), profilesstrAb(ilo2(end),1,2)];
        end
        ilo3 = find(profilesstrAb(:,2,3)<threshold);
        if ilo3(end) ==numpoints*2 ||ilo3(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIstrAb(3,:)= [NaN NaN];
        else
             CIstrAb(3,:) = [profilesstrAb(ilo3(1),1,3), profilesstrAb(ilo3(end),1,3)];
        end
    % this output a profiles 3D matrix whose first z dimension contains the
    % profiled parameters and the corresponding negLL
    % the CI contains the upper and lower CI on first the b parameter than
    % the d parameter
%% Plot profiles of strong Allee on birth to view CI

 figure;
         subplot(1, 3,1)
        plot(profilesstrAb(:,1,1), profilesstrAb(:,2,1))
        hold on
        plot(pbeststrAb(1),negLLfitstrAb,'r*','LineWidth',2)
        plot([profilesstrAb(1,1,1) profilesstrAb(end,1,1)],[threshold threshold],'r--')
        xlabel(' profiled b values')
        xlim([profilesstrAb(1,1,1) profilesstrAb(end,1,1)])
        ylabel('Cost Function Value')
        title(['Profile Likelihood of b CI = [',num2str(CIstrAb(1,:)),']'])

        subplot(1,3,2)
        plot(profilesstrAb(:,1,2), profilesstrAb(:,2,2))
        hold on
        plot(pbeststrAb(2),negLLfitstrAb,'r*','LineWidth',2)
        plot([profilesstrAb(1,1,2) profilesstrAb(end,1,2)],[threshold threshold],'r--')
        xlabel(' profiles d values')
        xlim ([profilesstrAb(1,1,2) profilesstrAb(end,1,2)])
        ylabel('Cost Function Value')
        title([' Profile Likelihood of d CI = [',num2str(CIstrAb(2,:)),']'])
        
        subplot(1,3,3)
        plot(profilesstrAb(:,1,3), profilesstrAb(:,2,3))
        hold on
        plot(pbeststrAb(3),negLLfitstrAb,'r*','LineWidth',2)
        plot([profilesstrAb(1,1,3) profilesstrAb(end,1,3)],[threshold threshold],'r--')
        xlabel(' profiled A values')
        xlim([profilesstrAb(1,1,3) profilesstrAb(end,1,3)])
        ylabel('Cost Function Value')
        title([' Profile Likelihood of A CI = [',num2str(CIstrAb(3,:)),']'])
        


%% Perform fitting: fit all data to strong Allee model on death
% For each fitting, we want to write functions in terms of the functions
% that run the model we are fitting to forward. The gen_model_mu runs the
% system of ODEs for the parameters and generates an expected model mean
% and variance for the input p parameters
% information

modelcode = 3;

% Initial guess
bguess = (log(mudatavec(end-2))-log(mudatavec(end-4)))/(timevec(end-2)-timevec(end-4)) + 0.002; 
dguess = 0.002;
Aguess = 1;
theta = [bguess,dguess, Aguess];
Ninit = N0list;
N= numsamps;
[pbeststrAd, BICstrAd, AICstrAd,negLLfitstrAd] = ML_fit(theta, tsamp, Ninit, N,mudatavec, vardatavec, modelcode);

% need to define th inline functions that run forward the mean and the
% variance
modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);

%% Plot the data compared to best fit model
figure;
subplot(1,2,1)
hold on
plot(timevec, mudatavec, 'r*')
for j = 1:length(N0list)
    modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit(j), modelcode);
plot(tsamp(2:end), modelfun_mu(pbeststrAd), 'b-', 'LineWidth', 2) % test model expected
end
xlabel('time (hours)')
ylabel('Mean cell number')
legend('data', 'best fit strong Allee (d)')
legend boxoff
title(['Mean: strong Allee on death model, BIC= ', num2str(BICstrAd)])
xlim([ 0 328])
subplot(1,2,2)
hold on
plot(timevec, vardatavec, 'g*')
for j = 1:length(N0list)
    modelfun_V =@(p)gen_model_var(p,tsamp, Ninit(j), modelcode);
plot(tsamp(2:end), modelfun_V(pbeststrAd), 'b-', 'LineWidth', 2) % test model expected
end
xlabel('time (hours)')
ylabel(' Variance')
legend('data', 'best fit strong Allee (d)')
legend boxoff
title(['Variance:strong Allee on death model, BIC= ', num2str(BICstrAd)])
xlim([ 0 328])

%% Profile likelihood of strong Allee on death model
% Use these to determine parameter certainty and find the 95 % CI on the
% b and d parameter estimate from all data

num_params = 3;

 % Check fxn

        fbest = negLLfitstrAd;
        factor = 0.05;
        numpoints = 10;
        params = pbeststrAd;
        [profilesstrAd] = profile_likelihood(params, tsamp, N0list, numsamps ,mudatavec, vardatavec, modelcode, factor, numpoints);
        threshold = chi2inv(0.95,length(params))/2 + fbest;
        ilo1=find(profilesstrAd(:,2,1)<threshold);
        if ilo1(end) ==numpoints*2 ||ilo1(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIstrAd(1,:) = [NaN, NaN];
        else
             CIstrAd(1,:) = [profilesstrAd(ilo1(1),1,1), profilesstrAd(ilo1(end),1,1)];
        end
        ilo2 = find(profilesstrAd(:,2,2)<threshold);
        if ilo2(end) ==numpoints*2 ||ilo2(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIstrAd(2,:)= [NaN NaN];
        else
             CIstrAd(2,:) = [profilesstrAd(ilo2(1),1,2), profilesstrAd(ilo2(end),1,2)];
        end
        ilo3 = find(profilesstrAd(:,2,3)<threshold);
        if ilo3(end) ==numpoints*2 ||ilo3(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIstrAd(3,:)= [NaN NaN];
        else
             CIstrAd(3,:) = [profilesstrAd(ilo3(1),1,3), profilesstrAd(ilo3(end),1,3)];
        end
    % this output a profiles 3D matrix whose first z dimension contains the
    % profiled parameters and the corresponding negLL
    % the CI contains the upper and lower CI on first the b parameter than
    % the d parameter

%% Plot profiles of strong Allee on death to view CI

 figure;
         subplot(1, 3,1)
        plot(profilesstrAd(:,1,1), profilesstrAd(:,2,1))
        hold on
        plot(pbeststrAd(1),negLLfitstrAd,'r*','LineWidth',2)
        plot([profilesstrAd(1,1,1) profilesstrAd(end,1,1)],[threshold threshold],'r--')
        xlabel(' profiled b values')
       % xlim([profilesstrAd(1,1,1) profilesstrAd(end,1,1)])
        ylabel('Cost Function Value')
        title(['Profile Likelihood of b CI = [',num2str(CIstrAd(1,:)),']'])

        subplot(1,3,2)
        plot(profilesstrAd(:,1,2), profilesstrAd(:,2,2))
        hold on
        plot(pbeststrAd(2),negLLfitstrAd,'r*','LineWidth',2)
        plot([profilesstrAd(1,1,2) profilesstrAd(end,1,2)],[threshold threshold],'r--')
        xlabel(' profiles d values')
        xlim ([profilesstrAd(1,1,2) profilesstrAd(end,1,2)])
        ylabel('Cost Function Value')
        title([' Profile Likelihood of d CI = [',num2str(CIstrAd(2,:)),']'])
        
        subplot(1,3,3)
        plot(profilesstrAd(:,1,3), profilesstrAd(:,2,3))
        hold on
        plot(pbeststrAd(3),negLLfitstrAd,'r*','LineWidth',2)
        plot([profilesstrAd(1,1,3) profilesstrAd(end,1,3)],[threshold threshold],'r--')
        xlabel(' profiled A values')
        xlim([profilesstrAd(1,1,3) profilesstrAd(end,1,3)])
        ylabel('Cost Function Value')
        title([' Profile Likelihood of A CI = [',num2str(CIstrAd(3,:)),']'])
        



%% Perform fitting: fit all data to strong Allee model on  birth & death
% For each fitting, we want to write functions in terms of the functions
% that run the model we are fitting to forward. The gen_model_mu runs the
% system of ODEs for the parameters and generates an expected model mean
% and variance for the input p parameters
% information

modelcode = 4;
% Initial guess
bguess = (log(mudatavec(end-2))-log(mudatavec(end-4)))/(timevec(end-2)-timevec(end-4)) + 0.002; 
dguess = 0.002;
Aguess = 1;
theta = [bguess,dguess, Aguess];
Ninit = N0list;
N= numsamps;
[pbeststrAbd, BICstrAbd, AICstrAbd,negLLfitstrAbd] = ML_fit(theta, tsamp, Ninit, N,mudatavec, vardatavec, modelcode)

% need to define th inline functions that run forward the mean and the
% variance
modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);

%% Plot the data compared to best fit model
figure;
subplot(1,2,1)
hold on
plot(timevec, mudatavec, 'r*')
for j = 1:length(N0list)
    modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit(j), modelcode);
plot(tsamp(2:end), modelfun_mu(pbeststrAbd), 'b-', 'LineWidth', 2) % test model expected
end
xlabel('time (hours)')
ylabel('Mean cell number')
legend('data', 'best fit strong Allee (b & d)')
legend boxoff
title(['Mean: strong Allee on birth & death model, BIC= ', num2str(BICstrAbd)])
xlim([ 0 328])
subplot(1,2,2)
hold on
plot(timevec, vardatavec, 'g*')
for j = 1:length(N0list)
    modelfun_V =@(p)gen_model_var(p,tsamp, Ninit(j), modelcode);
plot(tsamp(2:end), modelfun_V(pbeststrAbd), 'b-', 'LineWidth', 2) % test model expected
end
xlabel('time (hours)')
ylabel(' Variance ')
legend('data', 'best fit strong Allee (b & d)')
legend boxoff
title(['Variance: strong Allee on birth & death model, BIC= ', num2str(BICstrAbd)])
xlim([ 0 328])


%% Profile likelihood of strong Allee on birth & death model
% Use these to determine parameter certainty and find the 95 % CI on the
% b and d parameter estimate from all data

num_params = 3;

 % Check fxn

        fbest = negLLfitstrAbd;
        factor = 0.05;
        numpoints = 10;
        params = pbeststrAbd;
        [profilesstrAbd] = profile_likelihood(params, tsamp, N0list, numsamps ,mudatavec, vardatavec, modelcode, factor, numpoints);
        threshold = chi2inv(0.95,length(params))/2 + fbest;
        ilo1=find(profilesstrAbd(:,2,1)<threshold);
        if ilo1(end) ==numpoints*2 ||ilo1(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIstrAbd(1,:) = [NaN, NaN];
        else
             CIstrAbd(1,:) = [profilesstrAbd(ilo1(1),1,1), profilesstrAbd(ilo1(end),1,1)];
        end
        ilo2 = find(profilesstrAbd(:,2,2)<threshold);
        if ilo2(end) ==numpoints*2 ||ilo2(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIstrAbd(2,:)= [NaN NaN];
        else
             CIstrAbd(2,:) = [profilesstrAbd(ilo2(1),1,2), profilesstrAbd(ilo2(end),1,2)];
        end
        ilo3 = find(profilesstrAbd(:,2,3)<threshold);
        if ilo3(end) ==numpoints*2 ||ilo3(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIstrAbd(3,:)= [NaN NaN];
        else
             CIstrAbd(3,:) = [profilesstrAbd(ilo3(1),1,3), profilesstrAbd(ilo3(end),1,3)];
        end
    % this output a profiles 3D matrix whose first z dimension contains the
    % profiled parameters and the corresponding negLL
    % the CI contains the upper and lower CI on first the b parameter than
    % the d parameter
%% Plot the profiles for strong Allee model on birth & death
    
    figure;
         subplot(1, 3,1)
        plot(profilesstrAbd(:,1,1), profilesstrAbd(:,2,1))
        hold on
        plot(pbeststrAbd(1),negLLfitstrAbd,'r*','LineWidth',2)
        plot([profilesstrAbd(1,1,1) profilesstrAbd(end,1,1)],[threshold threshold],'r--')
        xlabel(' profiled b values')
       % xlim([profilesstrAd(1,1,1) profilesstrAd(end,1,1)])
        ylabel('Cost Function Value')
        title(['Profile Likelihood of b CI = [',num2str(CIstrAbd(1,:)),']'])

        subplot(1,3,2)
        plot(profilesstrAbd(:,1,2), profilesstrAbd(:,2,2))
        hold on
        plot(pbeststrAbd(2),negLLfitstrAbd,'r*','LineWidth',2)
        plot([profilesstrAbd(1,1,2) profilesstrAbd(end,1,2)],[threshold threshold],'r--')
        xlabel(' profiles d values')
        xlim ([profilesstrAbd(1,1,2) profilesstrAbd(end,1,2)])
        ylabel('Cost Function Value')
        title([' Profile Likelihood of d CI = [',num2str(CIstrAbd(2,:)),']'])
        
        subplot(1,3,3)
        plot(profilesstrAbd(:,1,3), profilesstrAbd(:,2,3))
        hold on
        plot(pbeststrAbd(3),negLLfitstrAbd,'r*','LineWidth',2)
        plot([profilesstrAbd(1,1,3) profilesstrAbd(end,1,3)],[threshold threshold],'r--')
        xlabel(' profiled A values')
        xlim([profilesstrAbd(1,1,3) profilesstrAbd(end,1,3)])
        ylabel('Cost Function Value')
        title([' Profile Likelihood of A CI = [',num2str(CIstrAbd(3,:)),']'])
        


%% Perform fitting: fit all data to weak/strong Allee on birth
% For each fitting, we want to write functions in terms of the functions
% that run the model we are fitting to forward. The gen_model_mu runs the
% system of ODEs for the parameters and generates an expected model mean
% and variance for the input p parameters
% information

modelcode = 5;

% Initial guess
bguess = (log(mudatavec(end-2))-log(mudatavec(end-4)))/(timevec(end-2)-timevec(end-4)) + 0.0002; 
dguess = 0.0002;
Aguess = -2;
tauguess = 3;
theta = [bguess,dguess, Aguess, tauguess];
Ninit = N0list;
N= numsamps;
[pbestwkAb, BICwkAb, AICwkAb,negLLfitwkAb] = ML_fit(theta, tsamp, Ninit, N,mudatavec, vardatavec, modelcode)

% need to define th inline functions that run forward the mean and the
% variance
modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);

%% Plot the data compared to best fit model
figure;
subplot(1,2,1)
hold on
plot(timevec, mudatavec, 'r*')
for j = 1:length(N0list)
    modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit(j), modelcode);
plot(tsamp(2:end), modelfun_mu(pbestwkAb), 'b-', 'LineWidth', 2) % test model expected
end
for j = 1:length(N0list)
    modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit(j), modelcode);
    y= modelfun_mu(pbestwkAb)
text(tsamp(end-10*j), y(end-7*j),['N_{0}=', num2str(N0list(j))], 'FontSize',14)  % test model expected
end
xlabel('time (hours)', 'FontSize',16)
ylabel('Mean cell number (<n(t)>)', 'FontSize',16)
legend('mean in data', 'best fit weak Allee (b) model', 'Location', 'NorthWest', 'FontSize',14)
legend boxoff
%title(['Mean fit ', num2str(BICwkAb)])
xlim([ 0 328])
subplot(1,2,2)
hold on
plot(timevec, vardatavec, 'g*')
for j = 1:length(N0list)
    modelfun_V =@(p)gen_model_var(p,tsamp, Ninit(j), modelcode);
plot(tsamp(2:end), modelfun_V(pbestwkAb), 'b-', 'LineWidth', 2) % test model expected
end
for j = 1:length(N0list)
    modelfun_V =@(p)gen_model_var(p,tsamp, Ninit(j), modelcode);
    y= modelfun_V(pbestwkAb)
text(tsamp(end-10*j), y(end-7*j),['N_{0}=', num2str(N0list(j))], 'FontSize',14)  % test model expected
end
xlabel('time (hours)', 'FontSize',16)
ylabel('Variance in cell number (\Sigma(t))', 'FontSize',16)
legend('variance in data', 'best fit weak Allee (b) model', 'Location', 'NorthWest', 'FontSize',14)
legend boxoff
%title(['Variance in data fit to weak Allee(b) model, BIC= ', num2str(BICwkAb)])
xlim([ 0 328])

%% Profile likelihood of weak/strong Allee on birth 
% Use these to determine parameter certainty and find the 95 % CI on the
% b and d parameter estimate from all data

num_params = 4;

 % Check fxn

        fbest = negLLfitwkAb; 
        factor = 0.01;
        numpoints = 10;
        params = pbestwkAb;
        [profileswkAb] = profile_likelihood(params, tsamp, N0list, numsamps ,mudatavec, vardatavec, modelcode, factor, numpoints);
        threshold = chi2inv(0.95,length(params))/2 + fbest;
        ilo1=find(profileswkAb(:,2,1)<threshold);
        if ilo1(end) ==numpoints*2 ||ilo1(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIwkAb(1,:) = [NaN, NaN];
        else
             CIwkAb(1,:) = [profileswkAb(ilo1(1),1,1), profileswkAb(ilo1(end),1,1)];
        end
        ilo2 = find(profileswkAb(:,2,2)<threshold);
        if ilo2(end) ==numpoints*2 ||ilo2(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIwkAb(2,:)= [NaN NaN];
        else
             CIwkAb(2,:) = [profileswkAb(ilo2(1),1,2), profileswkAb(ilo2(end),1,2)];
        end
        ilo3 = find(profileswkAb(:,2,3)<threshold);
        if ilo3(end) ==numpoints*2 ||ilo3(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIwkAb(3,:)= [NaN NaN];
        else
             CIwkAb(3,:) = [profileswkAb(ilo3(1),1,3), profileswkAb(ilo3(end),1,3)];
        end
        ilo4 = find(profileswkAb(:,2,4)<threshold);
        if ilo4(end) ==numpoints*2 ||ilo4(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIwkAb(4,:)= [NaN NaN];
        else
             CIwkAb(4,:) = [profileswkAb(ilo4(1),1,4), profileswkAb(ilo4(end),1,4)];
        end
    % this output a profiles 3D matrix whose first z dimension contains the
    % profiled parameters and the corresponding negLL
    % the CI contains the upper and lower CI on first the b parameter than
    % the d parameter
%% Plot profiles for weak/strong Allee on birth
    
    figure;
         subplot(1, 4,1)
        plot(profileswkAb(:,1,1), profileswkAb(:,2,1), 'LineWidth',2)
        hold on
        plot(pbestwkAb(1),negLLfitwkAb,'r*','LineWidth',2)
        plot([profileswkAb(1,1,1) profileswkAb(end,1,1)],[threshold threshold],'r--')
        xlabel(' profiled b values','FontSize',16)
        xlim([0.01005 0.0102])
        ylabel('Negative Log Likelihood','FontSize',16)
        title('b', 'FontSize', 18)
        %title(['b CI = [',num2str(CIwkAb(1,:)),']'])

        subplot(1,4,2)
        plot(profileswkAb(:,1,2), profileswkAb(:,2,2),'LineWidth',2)
        hold on
        plot(pbestwkAb(2),negLLfitwkAb,'r*','LineWidth',2)
        plot([profileswkAb(1,1,2) profileswkAb(end,1,2)],[threshold threshold],'r--')
        xlabel(' profiles d values', 'FontSize',16)
        xlim([-2e-4 3e-4])
        ylabel('Negative Log Likelihood', 'FontSize',16)
        title('d', 'FontSize', 18)
        %title(['d CI = [',num2str(CIwkAb(2,:)),']'])
        
        subplot(1,4,3)
        plot(profileswkAb(:,1,3), profileswkAb(:,2,3),'LineWidth',2)
        hold on
        plot(pbestwkAb(3),negLLfitwkAb,'r*','LineWidth',2)
        plot([profileswkAb(1,1,3) profileswkAb(end,1,3)],[threshold threshold],'r--')
        xlabel(' profiled A values','FontSize',16)
        ylabel('Negative Log Likelihood','FontSize',16)
        title('A', 'FontSize', 18)
        xlim([-4.5, -2])
        %title(['A CI = [',num2str(CIwkAb(3,:)),']'])
        
        subplot(1,4,4)
        plot(profileswkAb(:,1,4), profileswkAb(:,2,4),'LineWidth',2)
        hold on
        plot(pbestwkAb(4),negLLfitwkAb,'r*','LineWidth',2)
        plot([profileswkAb(1,1,4) profileswkAb(end,1,4)],[threshold threshold],'r--')
        xlabel(' profiled \tau values','FontSize',16)
        xlim ([ 6 9.5])
        ylabel('Negative Log Likelihood','FontSize',16)
        title('\tau', 'FontSize', 18)
        %title(['\tau CI = [',num2str(CIwkAb(4,:)),']'])



%% Perform fitting: fit all data to weak/strong Allee on death
% For each fitting, we want to write functions in terms of the functions
% that run the model we are fitting to forward. The gen_model_mu runs the
% system of ODEs for the parameters and generates an expected model mean
% and variance for the input p parameters
% information

modelcode = 6;
tsamp = 0:4:332; % set sample time (if data is down sampled change accordingly
tsampfit = tsamp(2:end);
tsampnew = tsampfit(ind);
tsamp = [0 tsampnew];
% Initial guess
bguess = (log(mudatavec(end-2))-log(mudatavec(end-10)))/(timevec(end-2)-timevec(end-10)) + 0.002; 
dguess = 0.002;
Aguess = -1;
tauguess = 2;
theta = [bguess,dguess, Aguess, tauguess];
Ninit = N0list;
N= numsamps;
[pbestwkAd, BICwkAd, AICwkAd,negLLfitwkAd] = ML_fit(theta, tsamp, Ninit, N,mudatavec, vardatavec, modelcode)

% need to define th inline functions that run forward the mean and the
% variance
modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);

%% Plot the data compared to best fit model
figure;
subplot(1,2,1)
hold on
plot(timevec, mudatavec, 'r*')
for j = 1:length(N0list)
    modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit(j), modelcode);
plot(tsamp(2:end), modelfun_mu(pbestwkAd), 'b-', 'LineWidth', 2) % test model expected
end
xlabel('time (hours)')
ylabel('Mean cell number')
legend('data', 'best fit weak Allee (d)')
legend boxoff
title(['Mean: weak Allee on death model, BIC= ', num2str(BICwkAd)])
xlim([ 0 328])
subplot(1,2,2)
hold on
plot(timevec, vardatavec, 'g*')
for j = 1:length(N0list)
    modelfun_V =@(p)gen_model_var(p,tsamp, Ninit(j), modelcode);
plot(tsamp(2:end), modelfun_V(pbestwkAd), 'b-', 'LineWidth', 2) % test model expected
end
xlabel('time (hours)')
ylabel(' Variance')
legend('data', 'best fit weak Allee (d)')
legend boxoff
title(['Variance: weak Allee on death model, BIC= ', num2str(BICwkAd)])
xlim([ 0 332])
%% Profile likelihood of weak/strong Allee on death
% Use these to determine parameter certainty and find the 95 % CI on the
% b and d parameter estimate from all data

num_params = 4;

 % Check fxn

        fbest = negLLfitwkAd;
        factor = 0.05;
        numpoints = 10;
        params = pbestwkAd;
        [profileswkAd] = profile_likelihood(params, tsamp, N0list, numsamps ,mudatavec, vardatavec, modelcode, factor, numpoints);
        threshold = chi2inv(0.95,length(params))/2 + fbest;
        ilo1=find(profileswkAd(:,2,1)<threshold);
        if ilo1(end) ==numpoints*2 ||ilo1(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIwkAd(1,:) = [NaN, NaN];
        else
             CIwkAd(1,:) = [profileswkAd(ilo1(1),1,1), profileswkAd(ilo1(end),1,1)];
        end
        ilo2 = find(profileswkAd(:,2,2)<threshold);
        if ilo2(end) ==numpoints*2 ||ilo2(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIwkAd(2,:)= [NaN NaN];
        else
             CIwkAd(2,:) = [profileswkAd(ilo2(1),1,2), profileswkAd(ilo2(end),1,2)];
        end
        ilo3 = find(profileswkAd(:,2,3)<threshold);
        if ilo3(end) ==numpoints*2 ||ilo3(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIwkAd(3,:)= [NaN NaN];
        else
             CIwkAd(3,:) = [profileswkAd(ilo3(1),1,3), profileswkAd(ilo3(end),1,3)];
        end
        ilo4 = find(profileswkAd(:,2,4)<threshold);
        if ilo4(end) ==numpoints*2 ||ilo4(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIwkAd(4,:)= [NaN NaN];
        else
             CIwkAd(4,:) = [profileswkAd(ilo4(1),1,4), profileswkAd(ilo4(end),1,4)];
        end
    % this output a profiles 3D matrix whose first z dimension contains the
    % profiled parameters and the corresponding negLL
    % the CI contains the upper and lower CI on first the b parameter than
    % the d parameter
    %%
    save('../out/profilesbd.mat', 'profilesbd')
    save('../out/profilesstrAb.mat', 'profilesstrAb')
    save('../out/profilesstrAd.mat', 'profilesstrAd')
    save('../out/profilesstrAbd.mat', 'profilesstrAbd')
    save('../out/profileswkAb.mat', 'profileswkAb')
    save('../out/profileswkAd.mat', 'profileswkAd')
%% Plot profiles for weak/strong Allee on death
    
    figure;
         subplot(1, 4,1)
        plot(profileswkAd(:,1,1), profileswkAd(:,2,1))
        hold on
        plot(pbestwkAd(1),negLLfitwkAd,'r*','LineWidth',2)
        plot([profileswkAd(1,1,1) profileswkAd(end,1,1)],[threshold threshold],'r--')
        xlabel(' profiled b values')
       % xlim([profilesstrAd(1,1,1) profilesstrAd(end,1,1)])
        ylabel('Cost Function Value')
        title(['Profile Likelihood of b CI = [',num2str(CIwkAd(1,:)),']'])

        subplot(1,4,2)
        plot(profileswkAd(:,1,2), profileswkAd(:,2,2))
        hold on
        plot(pbestwkAd(2),negLLfitwkAd,'r*','LineWidth',2)
        plot([profileswkAd(1,1,2) profileswkAd(end,1,2)],[threshold threshold],'r--')
        xlabel(' profiles d values')
        %xlim ([profileswkAd(1,1,2) profileswkAd(end,1,2)])
        ylabel('Cost Function Value')
        title([' Profile Likelihood of d CI = [',num2str(CIwkAd(2,:)),']'])
        
        subplot(1,4,3)
        plot(profileswkAd(:,1,3), profileswkAd(:,2,3))
        hold on
        plot(pbestwkAd(3),negLLfitwkAd,'r*','LineWidth',2)
        plot([profileswkAd(1,1,3) profileswkAd(end,1,3)],[threshold threshold],'r--')
        xlabel(' profiled A values')
        xlim([profileswkAd(1,1,3) profileswkAd(end,1,3)])
        ylabel('Cost Function Value')
        title([' Profile Likelihood of A CI = [',num2str(CIwkAd(3,:)),']'])
        
        subplot(1,4,4)
        plot(profileswkAd(:,1,4), profileswkAd(:,2,4))
        hold on
        plot(pbestwkAd(4),negLLfitwkAd,'r*','LineWidth',2)
        plot([profileswkAd(1,1,4) profileswkAd(end,1,4)],[threshold threshold],'r--')
        xlabel(' profiled tau values')
        xlim([profileswkAd(1,1,4) profileswkAd(end,1,4)])
        ylabel('Cost Function Value')
        title([' Profile Likelihood of \tau CI = [',num2str(CIwkAd(4,:)),']'])


%% Perform fitting: fit all data to weak/strong Allee on birth & death
% For each fitting, we want to write functions in terms of the functions
% that run the model we are fitting to forward. The gen_model_mu runs the
% system of ODEs for the parameters and generates an expected model mean
% and variance for the input p parameters
% information

modelcode = 7;

% Initial guess
bguess = (log(mudatavec(end-2))-log(mudatavec(end-4)))/(timevec(end-2)-timevec(end-4)) + 0.002; 
dguess = 0.002;
Aguess = -1;
tauguess = 2;
theta = [bguess,dguess, Aguess, tauguess];
Ninit = N0list;
N= numsamps;
[pbestwkAbd, BICwkAbd, AICwkAbd,negLLfitwkAbd] = ML_fit(theta, tsamp, Ninit, N,mudatavec, vardatavec, modelcode)

% need to define th inline functions that run forward the mean and the
% variance
modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);

%% Plot the data compared to best fit model
figure;
subplot(1,2,1)
hold on
plot(timevec, mudatavec, 'r*')
for j = 1:length(N0list)
    modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit(j), modelcode);
plot(tsamp(2:end), modelfun_mu(pbestwkAbd), 'b-', 'LineWidth', 2) % test model expected
end
xlabel('time (hours)')
ylabel('Mean cell number')
legend('data', 'best fit weak Allee (b & d)')
legend boxoff
title(['Mean: weak Allee birth & death model, BIC= ', num2str(BICwkAbd)])
xlim([ 0 328])
subplot(1,2,2)
hold on
plot(timevec, vardatavec, 'g*')
for j = 1:length(N0list)
    modelfun_V =@(p)gen_model_var(p,tsamp, Ninit(j), modelcode);
plot(tsamp(2:end), modelfun_V(pbestwkAbd), 'b-', 'LineWidth', 2) % test model expected
end
xlabel('time (hours)')
ylabel(' Variance')
legend('data', 'best fit weak Allee')
legend boxoff
title(['Variance: weak Allee birth & death model, BIC= ', num2str(BICwkAbd)])
xlim([ 0 332])

%% Profile likelihood of weak/strong Allee on birth & death
% Use these to determine parameter certainty and find the 95 % CI on the
% b and d parameter estimate from all data

num_params = 4;

 % Check fxn

        fbest = negLLfitwkAbd;
        factor = 0.05;
        numpoints = 10;
        params = pbestwkAbd;
        [profileswkAbd] = profile_likelihood(params, tsamp, N0list, numsamps ,mudatavec, vardatavec, modelcode, factor, numpoints);
        threshold = chi2inv(0.95,length(params))/2 + fbest;
        ilo1=find(profileswkAbd(:,2,1)<threshold);
        if ilo1(end) ==numpoints*2 ||ilo1(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIwkAbd(1,:) = [NaN, NaN];
        else
             CIwkAbd(1,:) = [profileswkAbd(ilo1(1),1,1), profileswkAbd(ilo1(end),1,1)];
        end
        ilo2 = find(profileswkAbd(:,2,2)<threshold);
        if ilo2(end) ==numpoints*2 ||ilo2(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIwkAbd(2,:)= [NaN NaN];
        else
             CIwkAbd(2,:) = [profileswkAbd(ilo2(1),1,2), profileswkAbd(ilo2(end),1,2)];
        end
        ilo3 = find(profileswkAbd(:,2,3)<threshold);
        if ilo3(end) ==numpoints*2 ||ilo3(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIwkAbd(3,:)= [NaN NaN];
        else
             CIwkAbd(3,:) = [profileswkAbd(ilo3(1),1,3), profileswkAd(ilo3(end),1,3)];
        end
        ilo4 = find(profileswkAbd(:,2,4)<threshold);
        if ilo4(end) ==numpoints*2 ||ilo4(1)==1
            display('All negLL vals below threshold. Increase factor.')
            CIwkAbd(4,:)= [NaN NaN];
        else
             CIwkAbd(4,:) = [profileswkAbd(ilo4(1),1,4), profileswkAbd(ilo4(end),1,4)];
        end
    % this output a profiles 3D matrix whose first z dimension contains the
    % profiled parameters and the corresponding negLL
    % the CI contains the upper and lower CI on first the b parameter than
    % the d parameter


%% Plot profiles for weak/strong Allee on birth death
    
    figure;
         subplot(1, 4,1)
        plot(profileswkAbd(:,1,1), profileswkAbd(:,2,1))
        hold on
        plot(pbestwkAbd(1),negLLfitwkAbd,'r*','LineWidth',2)
        plot([profileswkAbd(1,1,1) profileswkAbd(end,1,1)],[threshold threshold],'r--')
        xlabel(' profiled b values')
       % xlim([profilesstrAd(1,1,1) profilesstrAd(end,1,1)])
        ylabel('Cost Function Value')
        title(['Profile Likelihood of b CI = [',num2str(CIwkAbd(1,:)),']'])

        subplot(1,4,2)
        plot(profileswkAbd(:,1,2), profileswkAbd(:,2,2))
        hold on
        plot(pbestwkAbd(2),negLLfitwkAbd,'r*','LineWidth',2)
        plot([profileswkAbd(1,1,2) profileswkAbd(end,1,2)],[threshold threshold],'r--')
        xlabel(' profiles d values')
        %xlim ([profileswkAd(1,1,2) profileswkAd(end,1,2)])
        ylabel('Cost Function Value')
        title([' Profile Likelihood of d CI = [',num2str(CIwkAbd(2,:)),']'])
        
        subplot(1,4,3)
        plot(profileswkAbd(:,1,3), profileswkAbd(:,2,3))
        hold on
        plot(pbestwkAbd(3),negLLfitwkAbd,'r*','LineWidth',2)
        plot([profileswkAbd(1,1,3) profileswkAbd(end,1,3)],[threshold threshold],'r--')
        xlabel(' profiled A values')
        xlim([profileswkAbd(1,1,3) profileswkAbd(end,1,3)])
        ylabel('Cost Function Value')
        title([' Profile Likelihood of A CI = [',num2str(CIwkAbd(3,:)),']'])
        
        subplot(1,4,4)
        plot(profileswkAbd(:,1,4), profileswkAbd(:,2,4))
        hold on
        plot(pbestwkAbd(4),negLLfitwkAbd,'r*','LineWidth',2)
        plot([profileswkAbd(1,1,4) profileswkAbd(end,1,4)],[threshold threshold],'r--')
        xlabel(' profiled tau values')
        xlim([profileswkAbd(1,1,4) profileswkAbd(end,1,4)])
        ylabel('Cost Function Value')
        title([' Profile Likelihood of \tau CI = [',num2str(CIwkAbd(4,:)),']'])

   save('../out/profileswkAbd.mat', 'profileswkAbd')

%% Make vector of BIC values, negLL values, and parameter numbers
num_params = [ 2, 3, 3, 3, 4, 4, 4];
modellist = [1, 2, 3, 4, 5, 6, 7];
BICvals = [ BICbd, BICstrAb, BICstrAd, BICstrAbd, BICwkAb, BICwkAd, BICwkAbd];
negLLvals = [negLLfit, negLLfitstrAb, negLLfitstrAd, negLLfitstrAbd, negLLfitwkAb, negLLfitwkAd, negLLfitwkAbd];
modelnames = {'b-d', 'strAb', 'strAd', 'strAbd', 'wkAb', 'wkAd', 'wkAbd'};
shortnames = {'strAb', 'wkAb'};
figure;
subplot(3,1,1)
semilogy(modellist, BICvals, 'ko', 'LineWidth', 4)
xlabel('Model')
ylabel('BIC value')
xlim([0 8])
ylim([ min(BICvals)-10 max(BICvals)+5])
title('BIC values for each model (sampling every 4 hours)')
set(gca,'xticklabel',modelnames)


subplot(3,1,2)
semilogy(modellist, negLLvals, 'ro', 'LineWidth', 4)
xlabel('Model')
ylabel(' best negLL')
xlim([0 8])
ylim([ min(negLLvals)-5, max(negLLvals)+5])
title('negLL values for each model (sampling every 4 hours)')
set(gca,'xticklabel',modelnames)
subplot(3,1,3)
plot(modellist, num_params, 'bo', 'LineWidth', 4)
xlabel('Model')
ylabel(' number of parameters')
xlim([0 8])
ylim([ 1 5 ])
title ('number of parameters')
set(gca,'xticklabel',modelnames)
figure;
bar([1 2], [(BICvals(2)) (BICvals(5))])
xlabel('Model')
ylim([ 7800 7855])
ylabel('BIC value')
title('BIC model comparison')
set(gca,'xticklabel',shortnames)

