% Start by looking at single cell well data only
close all; clear all; clc
S = load('../out/green128.mat');
green= S.green;

for j = 1:length(green)
    ind(j)= green(j).N0actual ==1;
    raw = green(j).cellnum;
    cellnum_sm = smooth(raw,20,'sgolay');
    cellnum_sm = raw;
    green(j).N= cellnum_sm;
    green(j).t_real=green(j).time;
end

green1=green(ind);

%% Plot raw data
figure;
hold off
for j =1:length(green1)
    plot(green1(j).t_real, green1(j).N, 'color','b')
    hold on
end
xlabel('time (hours)')
ylabel('cell number')
xlim([0,100])
ylim([0,50])
title('Raw data cell number vs. time N_{0}=1')

%% UNIFORM SAMPLING from stochastic trajectories 
%Want to smooth by sampling from N at tsamp
num_samps = length(green1);
tstart =0; % this is time that you started measuring
tint = 5; % this is roughly how much time was spent in between images of data
tsamp = tstart:tint:100+tstart;
for j = 1:length(green1)%
    tstoch = green1(j).t_real;
    Nstoch = green1(j).N;

for i = 1:length(tsamp)
    % find nearest tstate that is less that tsamp
    ind =find(tstoch<=tsamp(i),1,'last');
    tfind = tstoch(ind);
    Nsamp(i,j)=Nstoch(ind);
end
    green1(j).tsamp = tsamp;
    green1(j).Nsamp = Nsamp;
end
%%
for j = 1:size(Nsamp,1)
    for i =1:size(Nsamp,2)
        if Nsamp(j,i)<0
            Nsamp(j,i)=0;
        end
    end
end
Nsamp = round(Nsamp,0);
mu_data = mean(Nsamp,2);
n_2_data = mean((Nsamp.^2),2);
var_data = n_2_data - ((mu_data).^2);

figure;
subplot(1,2,1)
hold on 
for j = 1:num_samps
plot(tsamp,Nsamp(:,j), 'b.');
hold on
end
plot(tsamp, mu_data, 'r-', 'LineWidth',3)
xlabel ('time (hours)')
ylabel('<n> measured')
xlim([0, tsamp(end)])
title('Measured mean value of N in time for N0=1')

subplot(1,2,2)
for j = 1:num_samps
plot(tsamp,Nsamp(:,j), 'b.');
hold on
end
plot(tsamp, var_data, 'g-', 'LineWidth',3)
xlabel ('time (hours)')
ylabel(' variance measured')
xlim([0, tsamp(end)])
title('Measured variance in time for N0=1')



%% Expected values for <n>, <n^2>, and variance using closed form solutions

% data values for mean, variance, and fourth moment
mu_data = mean(Nsamp,2);
n_2_data = mean((Nsamp.^2),2);
var_data = n_2_data - ((mu_data).^2);
sum_sq_diff= zeros(length(tsamp),1);
sum_quad_diff=zeros(length(tsamp),1);
for j = 1:num_samps
    sq_diff=(Nsamp(:,j)-mu_data).^2;
    quad_diff = (Nsamp(:,j)-mu_data).^4;
    sum_sq_diff = sum_sq_diff + sq_diff;
    sum_quad_diff = sum_quad_diff + quad_diff;
end
var_data2= sum_sq_diff/num_samps; % this just confirms that both ways of calculating variance are the same
var_data4 = sum_quad_diff/num_samps;

%% Perform Bayesian parameter estimation
% Want to fit parameters b, d, and sigma
t=tsamp(2:end); % only want to fit data after initial condition
N= num_samps;
% For now, start by assuming sigma = 1 and fit for b and d
sigma = 0.1;
pfxform = @(pval)[1 1].*log(pval); %'forward' parameter transform into Reals
pfxform4 = @(pval)[1 1 1 1].*log(pval); % for fitting with 4 params
pbxform = @(phat)[1 1].*exp(phat); %'backward' parameter transform into model space
pbxform4 = @(phat)[1 1 1 1].*exp(phat); % for fitting with 4 params
yfxform = @(y)log(y); % 'forward' transform for data and model output
ybxform = @(yhat)exp(yhat); % 'inverse' transform for data and model output


% guess what b-d should be from mu data
gguess = (yfxform(mu_data(end))-yfxform(mu_data(end-10)))/(tsamp(end)-tsamp(end-10)); 
%gguess= 0.023;
bguess = gguess + .01;
dguess = 0.01;
theta = [bguess,dguess];
sigmaT = 0.05; % random guess?
sig1 = sigmaT;
sig2 =sigmaT;
theta4 = horzcat(theta, sig1, sig2);
% can use this for now since all same starting value
N0 = mu_data(1);
V0 = var_data(1);
% Two observations to fit on, mu_data and var_data
% both use tsamp
modelfun_mu = @(p)mu_fxn(p, t, N0); % <n> for params (p) = b & d, vertical
modelfun_mu4 = @(p)mu_fxn(p(1:2), t, N0); % for fitting with sigmas
modelfun_n2= @(p)second_mom_fxn(p, t, N0); % Var for params (p) = b & d
modelfun_V=@(p)V_fxn(p, t, N0, V0); % vertical
modelfun_V4= @(p)V_fxn(p(1:2),t,N0, V0);
modelfun_both = @(p)mu_V_fxn(p, t, N0);
var_in_mean =  @(p)(1/N).*(V_fxn(p, t, N0, V0)) + (sigmaT.^2).*ones(length(t),1); % vertical
var_in_mean4 = @(p)(1/N).*(V_fxn(p(1:2), t, N0, V0)) + ((p(3).^2)).*ones(length(t),1);
var_in_mean_data = (1/N).*var_data;


% To find fourth order moment, take parameters, simulate 100 trajectories,
% and find the corresponding v4


var_in_var = @(p)(1/N).*(V4_fxn_ODE(p,N0, V0, t)-(((N-3)./(N-1)).*(V_fxn(p, t, N0, V0).^2))) + 2*sigmaT;
var_in_var4=@(p)(1/N).*(V4_fxn_ODE(p(1:2),N0, V0, t)-(((N-3)./(N-1)).*(V_fxn(p(1:2), t, N0, V0).^2))) + ((p(4).^2).*ones(length(t),1));
var_in_var_data = (1/N).*(var_data4(2:end)-(((N-3)/(N-1)).*(var_data(2:end).^2)));
%var_in_var = @(p)(1/N).*(var_data4(2:end)- (((N-3)/(N-1)).*((V_fxn(p, t, N0, V0))'.^2)));


% Don't fit on first time point because variance will be 0
v_data= var_data(2:end);
m_data = mu_data(2:end);
ydata = horzcat(m_data,v_data);
mudata = m_data;
vardata= v_data;

%% Try fitting all 4 parameters (b,d, sigu, sigv)

J4= @(phat)((0.5*sum((log(2*pi*var_in_mean4(pbxform4(phat)))) + (((yfxform(modelfun_mu4(pbxform4(phat))) - yfxform(mudata))./sqrt(var_in_mean4(pbxform4(phat)))).^2))) +...
 (0.5*sum((log(2*pi*var_in_var4(pbxform4(phat)))) + (((yfxform(modelfun_V4(pbxform4(phat))) - yfxform(vardata))./sqrt(var_in_var4(pbxform4(phat)))).^2))));

objfun_J4 = @(phat)(J4(phat));
lb = [0 0 sigmaT sigmaT];
ub = [ Inf Inf Inf Inf];
%phatbest_J4 = fminsearch(objfun_J4, pfxform4(theta4));
[phatbest_J4,fval,exitflag] = fminsearchbnd(objfun_J4, pfxform4(theta4), pfxform4(lb), pfxform4(ub))% need a better way to search parameter space,i.e. montecarlo
negLLguess4= objfun_J4(pfxform4(theta4));
params_best_J4= pbxform4(phatbest_J4);
mu_fitJ4 = mu_fxn(params_best_J4(1:2),t, N0);
V_fitJ4= V_fxn(params_best_J4(1:2),t,N0, V0);


figure;
subplot(1,2,1)
plot(tsamp, mu_data, 'r*')
hold on
plot(t, mu_fitJ4, 'b--', 'LineWidth',3)
xlabel('time (hours)')
ylabel('<n>')
title('Bayes fit mean')
legend('data', 'fit')
legend boxoff


subplot(1,2,2)
plot(tsamp, var_data, 'g*')
hold on
plot(t, V_fitJ4, 'b--', 'LineWidth',3)

xlabel('time (hours)')
ylabel('Variance')
title('Bayes fit variance')
legend('data', 'fit')
legend boxoff

%% Doesn't work well-- appears not to change from initial guess

%Take a look at J without summing...
J_t= @(phat)((0.5*((log(2*pi*var_in_mean4(pbxform4(phat)))) + (((yfxform(modelfun_mu4(pbxform4(phat))) - yfxform(mudata))./sqrt(var_in_mean4(pbxform4(phat)))).^2))) +...
 (0.5*((log(2*pi*var_in_var4(pbxform4(phat)))) + (((yfxform(modelfun_V4(pbxform4(phat))) - yfxform(vardata))./sqrt(var_in_var4(pbxform4(phat)))).^2))));
% using 4th order variance in data
% J_t= @(phat)((0.5*((log(2*pi*var_in_mean(pbxform(phat)))) + (((yfxform(modelfun_mu(pbxform(phat))) - yfxform(mudata))./sqrt(var_in_mean(pbxform(phat)))).^2))) +...
%  (0.5*((log(2*pi*var_in_var_data)) + (((yfxform(modelfun_V(pbxform(phat))) - yfxform(vardata))./sqrt(var_in_var_data)).^2))));
negLLguesst=J_t(pfxform4(theta4));
negLLfind = J_t(phatbest_J4);
figure;
hold on
plot(t, negLLguesst, 'LineWidth',2)
plot(t, negLLfind, 'LineWidth', 2)
xlabel ('time')
ylabel('contribution to negative LL')
title('Negative LL in time')
legend('initial guess', 'fit param')

%%  Write Markov-Chain Monte Carlo Search Algorithm

% write a function that finds the negativeLL for the mudata, vardata, and
% parameters


J_init = negLLfxn4( pfxform4(theta4), mudata, vardata, N, t, N0, V0 );
J_curr = J_init; % set current negLL value to inital value to initialize
count_accepted = 0; % start with no changes accepted

% Initialize a vector to store b and d params in
nruns = 5000;
store_params = zeros(nruns, 5); % store J, b,  d, sigu, and sig vfor each run
store_acc_params = zeros(nruns, 5 ); % store the converged parameters and J, b, d, sigu, sigv,

% Outside the loop: initialize temperature and iterations
T0= 3;
k = 1:1:nruns;
T= T0*exp(-1*k./(nruns)); % sets up gradual cooling

% prob_test= exp((J_curr-negLLtrue4)./T(1))

% check that temperature annealing is working
figure;
hold off
set(gca,'LineWidth',1.5,'FontSize',12);
plot(k, T, 'LineWidth', 3)
xlabel('Markov step')
ylabel('Temperature')
title('Temperature Annealing')

bstep = 0.005; % step size for searching birth param
dstep = 0.005; % step size for searchign death param
sigustep = 0.001;
sigvstep = 0.001;
% set initial guess as b1 and d1
b1=theta4(1);
d1=theta4(2);
sigu1= theta4(3);
sigv1 = theta4(4);


%%
for k = 1:nruns
    b2 = b1 + bstep*(2*rand-1);
    d2 = d1 + dstep*(2*rand-1);
    sigu2 = sigu1 + sigustep*(2*rand-1);
    sigv2 = sigv1 + sigvstep*(2*rand-1);
    
    % Constrain search region
    if b2<0 
        b2=0;
    end
    if b2>theta4(1)+.05
        b2 = theta4(1)+.05;
    end
    if d2<0
        d2=0;
    end
    if d2>theta4(2)+.05
        d2 = theta4(2)+.05;
    end
    if sigu2<sigmaT
        sigu2 = sigmaT;
    end
    if sigu2>2*sigmaT
        sigu2= 2*sigmaT;
    end
    if sigv2<0.9*sigmaT
        sigv2 = 0.9*sigmaT;
    end
    if sigv2>2*sigmaT
        sigv2= 2*sigmaT;
    end

    % find the neg LL of the new params
    theta2 = horzcat(b2,d2, sigu2, sigv2);
    J_new = negLLfxn4( pfxform4(theta2), mudata, vardata, N, t, N0, V0 );
    % store the  negLL and searched parameters
    store_params(k,1) = J_new;
    store_params(k,2) =b2;
    store_params(k,3) = d2;
    store_params(k,4) = sigu2;
    store_params(k,5) = sigv2;
    
    prob(k) = exp((J_curr-J_new)./T(k));
    % if Jcurr > Jnew, Jnew is better--the numerator in the exponent is positive &
    % prob>1--> change will always be accepted
    % if Jcurr < Jnew, numerator is negative, & prob <1 --> chance change
    % will be accepted
    
    if rand < prob(k) % if true, accept the change!
        b1 = b2;
        d1 = d2;
        sigu1= sigu2;
        sigv1 = sigv2;
        theta1 = theta2;
        J_curr = J_new;
        count_accepted = count_accepted +1;
        % decrease search step size
        bstep = 0.999*bstep;
        dstep = 0.999*dstep;
        sigustep = 0.99*sigustep;
        sigvstep = 0.99*sigvstep;
    else
        % increase search step size
        bstep = 1.001*bstep;
        dstep = 1.001*dstep;
        sigustep = 1.00001*sigustep;
        sigvstep = 1.00001*sigvstep;
    end
    store_acc_params(k,1) = J_curr;
    store_acc_params(k,2)= b1;
    store_acc_params(k,3)= d1;
    store_acc_params(k,4)= sigu1;
    store_acc_params(k,5)= sigv1;
       
end
%%
params_best_MC = store_acc_params(end,2:end)

negLLfitfmin=J_t((phatbest_J4));
negLLMC = J_t(pfxform4(params_best_MC));
figure;
hold on
plot(t, negLLfitfmin, 'LineWidth',2)
plot(t, negLLMC, 'LineWidth', 2)
xlabel ('time')
ylabel('contribution to negative LL')
title('Negative LL in time')
legend( 'fminsearch', 'MC fit Params')
%%

figure;
plot(1:1:nruns, store_acc_params(:,1))
xlabel('iterations')
ylabel('negLL')
title('negLL')

 figure;
 plot(1:1:nruns, store_acc_params(:,2), 'g.')
 hold on
 plot(1:1:nruns, store_acc_params(:,3), 'b.')
%  plot(1:1:nruns, store_acc_params(:,4),'r.')
%   plot(1:1:nruns, store_acc_params(:,5),'m.')
 legend( 'birth rates', 'death rates', '\sigma\mu', '\sigmavar')
 xlabel('iterations')
 ylabel('rate')
 title('Accepted birth and death parameters')
 
 figure;
 plot(1:1:nruns, store_params(:,2), 'g.')
 hold on
 plot(1:1:nruns, store_params(:,3), 'b.')
%  plot(1:1:nruns, store_params(:,4), 'r.')
%   plot(1:1:nruns, store_params(:,5), 'm.')
 legend( 'birth rates', 'death rates', '\sigma\mu', '\sigmavar')
 xlabel('iterations')
 ylabel('rate')
 title('Explored birth and death parameters')
 
 figure;
 plot(1:1:nruns, prob, 'r.')
 ylim([0 1])
%%
figure;
subplot(1,4,1)
hist(store_acc_params(:,2))
ylabel('frequency')
xlabel('b')
subplot(1,4,2)
hist(store_acc_params(:,3))
ylabel('frequency')
xlabel('d')
subplot(1,4,3)
hist(store_acc_params(:,4))
ylabel('frequency')
xlabel('\sigma\mu')
subplot(1,4,4)
hist(store_acc_params(:,5))
ylabel('frequency')
xlabel('\sigmavar')

% Plots of two parameters colored by likelihood
likelihoodval = exp(-store_acc_params(:,1));
likelihoodvalall = exp(-store_params(:,1));
pointsize = 30;
figure;
scatter(store_acc_params(:,2), store_acc_params(:,3), pointsize, likelihoodval,'filled');
colorbar
xlabel('b')
%xlim([ b-.005 b + .005])
%ylim([0 d+.005])
ylabel('d')
title(' b versus d colored by likelihood for accepted parameters')

figure;
scatter(store_acc_params(:,3), store_acc_params(:,4), pointsize, likelihoodval,'filled');
colorbar
xlabel('\sigma\mu')
%xlim([ b-.005 b + .005])
%ylim([0 d+.005])
ylabel('\sigmavar')
title(' \sigmamu vs \sigmavar colored by likelihood for accepted parameters')

figure;
scatter(store_acc_params(:,3), store_acc_params(:,4), pointsize, likelihoodval,'filled');
colorbar
xlabel('d')
%xlim([ b-.005 b + .005])
%ylim([0 d+.005])
ylabel('\sigma\mu')
title(' \sigma\mu vs d colored by likelihood for accepted parameters')
%%

std_dev = std(store_acc_params(:, 2:end))

CI(1,:)= params_best_MC - 1.96*std_dev;
CI(2,:) = params_best_MC + 1.96*std_dev;
figure;
subplot(2,2,1)
plot(store_acc_params(:,2), likelihoodval,'ro')
hold on
xlabel ('b')
ylabel('likelihood')
title (['b=',num2str(round(params_best_MC(1), 4)),' +/- ', num2str(round(1.96*std_dev(1),4))])
subplot(2,2,2)
plot(store_acc_params(:,3), likelihoodval,'bo')
xlabel ('d')
ylabel('likelihood')
title (['d=',num2str(round(params_best_MC(2), 4)),' +/- ', num2str(round(1.96*std_dev(2),4))])
subplot(2,2,3)
plot(store_acc_params(:,4), likelihoodval,'go')
xlabel ('\sigma_{\mu}')
ylabel('likelihood')
title (['\sigma_{\mu}=',num2str(round(params_best_MC(3), 4)),' +/- ', num2str(round(1.96*std_dev(3),4))])
subplot(2,2,4)
plot(store_acc_params(:,5), likelihoodval,'mo')
xlabel ('\sigma_{var}')
ylabel('likelihood')
title (['\sigma_{var}=',num2str(round(params_best_MC(4), 4)),' +/- ', num2str(round(1.96*std_dev(4),4))])

