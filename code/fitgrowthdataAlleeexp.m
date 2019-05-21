function [ RsqAllee, Rsqexp, AICAllee, AICexp] = fitgrowthdataAlleeexp(A, g, eta, N0)

t = 1:2:100;
tbig = repmat(t, 5, 1);
num_meas = length(t);
num_runs = 5; % length of N0
vert_length = num_meas*num_runs;
tlong = reshape(tbig', vert_length,1);

params = [g,A];


iall = find(tlong>=0);
Nmodel = simmodelAllee(params, tbig, N0);
%Nmodel = round(Nmodel,0);
Nmodellong = simmodelAlleelong(params, tlong,iall, N0);
%Nmodellong = round(Nmodellong, 0);

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

for j = 1:length(N0)
    percapitag(:,j) = Nmodel(:,j)/N0(j);
end



% Fit Allee effect with Bayesian

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
end