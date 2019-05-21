% Analyze Allee Effect with closed form solution
close all; clear all; clc
% This script ensures that the Allee effect works for a
% set of input conditions, and analyzes the resulting effects of growth
% with various parameters and initial conditions

% Growth Dynamics Function:
% dN/dt = gN(N/A-1)
% N(t)= (N0-A)exp(gt/A) + A
%% Run once  for single value of A, g, and eta
t = 0:2:98;
N0 = [8 16 32 ];
num_runs = length(N0);
tbig = repmat(t, num_runs, 1);
num_meas = length(t);
vert_length = num_meas*num_runs;
tlong = reshape(tbig', vert_length,1);
A = 5;
eta = 0.2;
% take g and sigma from 32 cell well growth curves
g = 0.0188;
sigma = 0.0075;

params = [g,A];

% for i = 1:length(N0)
% Nmodel(:,i)= (N0(i)-A).*exp(g*tbig(i,:)) + A;
% for j = 1:length(Nmodel(:,i))
%     if Nmodel(j,i) <0
%         Nmodel(j,i) = 0;
%     end
% end
% end
iall = find(tlong>=0);
%%
N0 = [ 5, 10, 20];
Nmodel = simmodelAllee(params, tbig, N0);
%Nmodel = round(Nmodel,0);
Nmodellong = simmodelAlleelong(params, tlong,iall, N0);
Nmodellong = simmodelAlleelong_distrib(params,sigma, tlong, N0);
%Nmodellong = round(Nmodellong, 0);

% test

%%
eta = 1;
noise = eta*(1-2*randn(length(tbig),length(N0)));
Nfake = Nmodel + noise;
% eliminate negative values in noisy data
for i = 1:length(N0)
for j = 1:length(Nfake(:,i))
  if Nfake(j,i) <0
        Nfake(j,i) = 0;
  end
end
end
Nfake = round(Nfake,0);
Nfakelong = reshape(Nfake, [size(Nfake,1)*size(Nfake,2),1]);
figure;
subplot(1,2,1)
hold off
for j = 1:length(N0)
    hold on
    semilogy(tbig(j,:), Nmodel(:,j), 'LineWidth', 3)
    hold on
    %semilogy(tbig(j,:), Nfake(:,j), 'LineWidth', 2)
    text(tbig(j, end-15), Nmodel(end-15,j), [ 'N_{0} = ', num2str(N0(j)) ],'HorizontalAlignment','left','VerticalAlignment','bottom','color','k')
   % plot(tout, Nfake(:,j), 'o')
   %ylim([ 0 200])
end
hold on
%plot(tlong, Nmodellong,'*')
xlabel('time')
ylabel('N')
title(['Allee effect growth dynamics, A= ', num2str(A)])
%legend('N_{0} = 2', 'N_{0} = 5', 'N_{0} = 10', 'N_{0} = 12', 'N_{0} = 20')
legend boxoff

for j = 1:length(N0)
    percapitag(:,j) = Nmodel(:,j)/N0(j);
end

subplot(1,2,2)
hold on
for j = 1:length(N0)
    plot(tbig(j,:), log(percapitag(:,j)),'LineWidth', 3);
    text(tbig(j,30), log(percapitag(30,j)), [ 'N_{0} = ', num2str(N0(j)) ],'HorizontalAlignment','left','VerticalAlignment','bottom','color','k')
end
%ylim([-1 2])
xlabel('time')
ylabel('log(N/N_{o})')
title(['Allee effect per capita growth rate, A= ', num2str(A)])

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
LLOQ = 2;
icens = find(Nfakelong <LLOQ);
igood = find(Nfakelong>= LLOQ);

modelfungood= @(p)simmodelAlleelong(p, tlong, igood, N0);
modelfuncens= @(p)simmodelAlleelong(p, tlong, icens, N0);
modelfunsinggood = @(p)simmodelsingexplong(p, tlong, igood, N0);
modelfunsingcens = @(p)simmodelsingexplong(p, tlong, icens, N0);

gguess = g + 0.001;
Aguess = A+.5;
sigma = eta;
theta = [gguess, Aguess]; % g and A
thetag = gguess + 0.001;
Nfakelonghigh = ones(length(icens),1);
Nfakelonglow = zeros(length(icens),1);

loglikelihood = @(phat)(sum(log(normpdf(yfxform(Nfakelong(igood)),yfxform(modelfungood(pbxform(phat))), sigma)))+...
    sum(log(normcdf(yfxform(Nfakelonghigh), yfxform(modelfuncens(pbxform(phat))),1))));
loglikelihoodg = @(phat)(sum(log(normpdf(yfxform(Nfakelong(igood)),yfxform(modelfunsinggood(pbxformg(phat))), sigma)))+...
    sum(log(normcdf(yfxform(Nfakelonghigh), yfxform(modelfunsingcens(pbxformg(phat))),1))));

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
Nfit = simmodelAllee(params_Bayes, tbig, N0);
Nfitsing = simmodelsingexplong(params_sing, tlong, iall, N0);
% Now find chi-squared, R-squared, and AIC value for Allee and single
% exponential model
%
%residuals (model- measured)
resAlleesq = Nfit-Nfake;
resAllee = reshape(resAlleesq, vert_length, 1);
resexp = Nfitsing-Nfakelong;
% average N
Nbar = mean(mean(Nfake));
%R-squared
RsqAllee = 1- (sum((resAllee).^2)./(sum((Nbar-Nfakelong).^2)));
Rsqexp = 1- (sum((resexp).^2)./(sum((Nbar-Nfakelong).^2)));
num_p_Allee = 2; % g and A
num_p_exp = 1; %g
AICAllee = -2*loglikelihood(phatbest) + 2*num_p_Allee;
AICexp = -2*loglikelihoodg(phatbestsing) + 2*num_p_exp;
delta_AIC = AICAllee- AICexp; % if value is negative then Allee is better model
% will run a loop through different values of A and eta and find the
% corresponding delta_AIC at each point. This will be plotted on a heat map
% to be used for parameter identifiability


%% Plot fitted model versus "data"
% Compare visually the data to the model predicted by the Bayesian fitted
% parameters



figure;
hold off
for j = 1:length(N0)
    hold on
    semilogy(tbig(j,:), Nfit(:,j), 'LineWidth', 3,'color', 'b')
    hold on
    semilogy(tbig(j,:), Nfake(:,j), 'LineWidth', 2, 'color', 'g')
    text(tbig(j, end-27), Nfit(end-27,j), [ 'N_{0} = ', num2str(N0(j)) ],'HorizontalAlignment','left','VerticalAlignment','bottom','color','k')
   % plot(tout, Nfake(:,j), 'o')
   %ylim([ 0 200])
end
semilogy(tlong, Nfitsing, 'r.')
xlabel('time')
ylabel('N')
ylim([0 100])
%title('Allee vs. single exponential model')
%title(['Simulated Data vs. Fitted Model, A= ', num2str(A), ', A_{fit}=',num2str(params_Bayes(2)), 'vs. Single Exponential Fit, g_{fit}=', num2str(params_sing)])
%legend('N_{0} = 2', 'N_{0} = 5', 'N_{0} = 10', 'N_{0} = 12', 'N_{0} = 20')
legend boxoff

%% Run for all values of A, eta, and N0
% range of A and eta
A = [1:1:20];
% eta = [.1:.025:.575];
%  N0 = [ 2 5 10 12 20 ];
%  [ALLEE, ETA] = meshgrid(A, eta); % grid of parameters
% Aflat = reshape(ALLEE, 1, []);
% Etaflat = reshape(ETA, 1, []);
% run a loop through all combos of A & eta
% Redo with varying sigma
eta = 0.1;
sig = [.0005:.0004:.0084];
[ALLEE, SIG] = meshgrid(A, sig); % grid of parameters
Aflat = reshape(ALLEE, 1, []);
Sigflat = reshape(SIG, 1, []);
%%
for j = 1:length(Aflat)

% [ RsqAllee(j), Rsqexp(j), AICAllee(j), AICexp(j)] = fitgrowthdataAlleeexp(Aflat(j), g, Etaflat(j), N0);
sigma = Sigflat(j);
[ RsqAllee(j), Rsqexp(j), AICAllee(j), AICexp(j)] = fitgrowthdataAlleeexp_distrib(Aflat(j), g,sigma, eta, N0);
end

%% Next use vectors to reconstruct grid and make image
delta_AIC = -(AICAllee- AICexp);

del_AIC = reshape(delta_AIC, size(ALLEE));
% put fitted patient parameter fit in figure
figure; 
minmax = @(x)([min(x) max(x)]);
imagesc(minmax(A),minmax(sig),del_AIC);
%[C,h]=contourf(GAM,PHI,DTTPDK); clabel(C,h); colorbar
caxis(10000*[-1 1]);
cmap = [[1-[0:0.8/31:.8]';zeros(32,1)] [zeros(32,1);[.2:.8/31:1]'] zeros(64,1)]; % green to red colormap
colormap(cmap); colorbar;
xlabel('A (Allee threshold)');
ylabel('\sigma (spread in g_{distribution})');
title('\delta AIC (relative improvement in Allee vs. single exponential model)');
set(gca,'Ydir','normal'); hold on;

