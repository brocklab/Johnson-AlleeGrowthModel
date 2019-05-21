% Moment- Approach fitting of stochastic birth-death-Allee process for small cell
% number data

% This code is first attempt to use moment-approach (so mean and
% variance)of N(t) trajectories to fit small cell number well data.

% We will start by simulating stochastic cell trajectories using a set b
% and d, and A , and see if we can use the following equations for <n(t)> and V(t)
% to perform parameter estimation in Bayesian framework
% For random birth death process we have:

%d<n>/dt=(b-d)<n>(1-A/<n>)


close all; clear all; clc

% Start by generating 20 n0=1 trajectories
% assume birth rate = 0.0238, death rate = 0.005
% Set up time and N counter
Ninit = [1 2 3 4 5 6 7];
Ninit = 5;
for i = 1:length(Ninit)
b = 0.0233 + .0005; % birth rate
d = 0.0045 + .0005; % death rate
delta= b+d;
A = 2;
%birth_n = b*N; % birth hazard function
%death_n = d*N; % death hazard function
num_samps = 500;
num_iters = 150;
take_offs = 0;
state = zeros(num_iters,num_samps);
tstate = zeros(num_iters,num_samps);
state(1,:) = Ninit(i); % at time 0, number of cells =N
tjump(1, :) = 0; % start at time 0
ct_extinct = 0;
for j = 1:num_samps

    N=Ninit(i);
    N0 = N;
    time(1)= 0;
for k = 2:num_iters
    birth_n = (b*N)-(b-d)*A; % birth 
    if birth_n <0
        birth_n = 0;
    end
    death_n = d*N; % death 
    if N==0
        N=0;
    else
        r = rand;
        if r< (birth_n)/(birth_n+death_n)
        N = N+1;
        end
        if r>= (birth_n)/(birth_n+death_n)
        N = N-1;
        end
    end
    state(k, j) = N;
    % set time step to be proportional
    r2=rand;
    tstep =-log(r2)/(birth_n + death_n);
    if tstep == Inf
        % make tstep same as previous?
        tstep = 1;
    end
    time = time + tstep;
    tstate(k,j) = time;
    
    % If N goes below 0, cells go extinct, N=0 throughout
    if N <= 0
       state(k:end,j) = 0;
    end
   
end
    
    thres = 0;
     if state(end,j) > thres
        take_offs= take_offs +1;
     end 
    ind(j) = state(end,j)>thres;


end
  P_takeoff(i)= take_offs/num_samps;
  P_tkoff_theor(i) = 1-((d/b).^Ninit(i));
  end

% find minimum time "measured"
tmin = min(tstate(end, :));

% add constant technical noise to data
sigmaT = 0; % guess that on average count is off by one quarter of a cell?
state = state+ round(normrnd(0, sigmaT,size(state)));
% eliminate those below 0
for j=1:size(state,1)
    for i = 1:size(state,2)
        if state(j,i)<0
            state(j,i)=0;
        end
    end
end
%% Plot P_est with simulated data vs expected value
figure;
plot(delta,P_takeoff, '*', 'LineWidth',2)
hold on
plot(delta, P_tkoff_theor,'-','LineWidth',2)
xlabel('\delta = b+d')
ylabel('P_{establishment}')
title('Simulated vs. Expected P_{establishment} for increasing \delta')
legend('simulated', 'expected')
legend boxoff


P_est_theor=1-((d/b).^Ninit);

figure;
plot(Ninit, P_est_theor,'-', 'LineWidth',2)
hold on
plot(Ninit,P_takeoff, '*', 'LineWidth',2)
xlabel('Initial cell number')
ylabel('P_est')
title('P_est vs. N_{0} for constant \delta and b-d')
legend('expected if b-d process', 'measured from Allee sim')
legend boxoff
ylim([0 1])
%% Plot simulated cell number trajectories out to minimum time reached
figure;
hold off
for j = 1:num_samps
plot(tstate(:,j), state(:,j))
hold on
end
xlim([0, tmin])
xlabel('time (hours)', 'FontSize', 16)
ylabel('Number of cells (N(t))', 'FontSize', 16)
title('Simulated stochastic trajectories for N_{0}=5', 'FontSize', 14)
%title(['Simulated N(t) trajectories for b=', num2str(b), ' & d=', num2str(d)])


%% UNIFORM SAMPLING from stochastic trajectories 
%Want to smooth by sampling from N at tsamp
tstart=0;
tint=2;

tsamp = tstart:tint:100+tstart;
for j = 1:num_samps
    tstoch = tstate(:,j);
    Nstoch = state(:,j);

for i = 1:length(tsamp)
    % find nearest tstate that is less that tsamp
    ind =find(tstoch<=tsamp(i),1,'last');
    tfind = tstoch(ind);
    Nsamp(i,j)=Nstoch(ind);
end

end
mu_data = mean(Nsamp,2);
figure;
hold on 
for j = 1:num_samps
plot(tsamp,Nsamp(:,j), 'r.');
hold on
plot(tstate(:,j), state(:,j),'b.');
end
plot(tsamp, mu_data, 'k-', 'LineWidth',3)
xlabel ('time (hours)')
ylabel('number of cells')
xlim([0, tsamp(end)])
%title('Mean value of N in time for N0=5, A = 2')
%% Expected values for Allee model
% need to change all ODEs
% tsamp =tsamp(5:end);
% mu_data = mean(Nsamp,2);
% n_2_data = mean((Nsamp.^2),2);
% n_3_data = mean((Nsamp.^3),2);
% n_4_data = mean((Nsamp.^4),2);
% var_data = n_2_data - ((mu_data).^2);
% var4_data = n_4_data - ((mu_data).^4);
% V0= var_data(5);
% N0 = mu_data(5);

V0=0;
C_init(1)=N0;
C_init(2)=N0.^2;
C_init(3) = V0;
C_init(4)= N0.^3;
C_init(5)= N0.^4;
C_init(6) = V0;


f = @(t,C) [((b-d)*C(1)-(b-d)*A); % dn/dt
             2*C(2).*(b-d) - (2*C(1).*(b-d)*A) + (C(1).*(b+d))-((b-d)*A);   % dn2/dt
            (2*C(2).*(b-d) - (2*C(1).*(b-d)*A) + (C(1).*(b+d))-((b-d)*A)-(2*C(1).*((b-d)*C(1)-(b-d)*A))); % dV/dt
            (3*C(4).*(b-d)) + (3*C(2).*(b+d)-3*C(2)*(b-d)*A) - (3*C(1).*(b-d)*A)-((b-d)*A);%dn3dt
            (4*C(5).*(b-d)) + (6*C(4).*(b+d))-(4*C(4)*(b-d)*A)+(4*C(2).*(b-d))-(6*C(2).*(b-d)*A)+ (C(1).*(b+d))+...
            (4*C(1).*(b-d)*A)-((b-d)*A); %dn4dt
            (4*C(5).*(b-d)) + (6*C(4).*(b+d))-(4*C(4).*(b-d)*A)+(4*C(2).*(b-d))-(6*C(2).*(b-d)*A)+ (C(1).*(b+d))+...
            (4*C(1).*(b-d)*A)-((b-d)*A)-(4.*((C(1)).^3)*(((b-d).*C(1)-(b-d)*A)))];

options1 = odeset('Refine',1);  
options = odeset(options1,'NonNegative',1:6);
[t,C]=ode45(f, tsamp,C_init, options);
mu_C= C(:,1);
n2_C= C(:,2);
v2_C=C(:,3);
n3_C = C(:,4);
n4_C = C(:,5);
v4_C=C(:,6);

%% Plot Sampled Data from Stochastic Simulations alongside expected values of mean and variance
mu_data = mean(Nsamp,2);
n_2_data = mean((Nsamp.^2),2);
n_3_data = mean((Nsamp.^3),2);
n_4_data = mean((Nsamp.^4),2);
var_data = n_2_data - ((mu_data).^2);
var4_data = n_4_data - ((mu_data).^4);
%var_data = (std(Nsamp,0,2).^2);

figure;
subplot(1,2,1)
plot(tsamp, mu_data(1:end), 'r*')
hold on
plot(tsamp, mu_C, 'k-', 'LineWidth',2)
xlabel('time (hours)')
ylabel('Mean cell number')
%title('Strong Allee on birth mean', 'FontSize', 14)
title('Mean', 'FontSize', 14)
legend('mean simulated data', 'expected mean', 'FontSize',12, 'Location', 'NorthWest')
legend boxoff

% subplot(1,3,2)
% plot(tsamp, n_2_data(1:end), 'g*')
% hold on
% plot(tsamp, n2_C, 'b.', 'LineWidth',2)
% xlabel('time (hours)')
% ylabel('<n2>')
% title('Expected vs. simulated <n2>')
% legend('<n2> in simulated data', 'expected <n2>')
% legend boxoff

subplot(1,2,2)
plot(tsamp, var_data(1:end), 'g*')
hold on
plot(tsamp, v2_C, 'k-', 'LineWidth',2)
xlabel('time (hours)')
ylabel('Variance')
%title('Strong Allee on birth variance', 'FontSize', 14)
title ('Variance', 'FontSize',14)
legend('variance in simulated data', 'expected variance', 'FontSize',12, 'Location', 'NorthWest')
legend boxoff
%%
figure;
subplot(1,3,1)
plot(tsamp, n_3_data(1:end), 'm*')
hold on
plot(tsamp, n3_C, 'b-', 'LineWidth',2)
xlabel('time (hours)')
ylabel('<n3>')
title('Expected vs. simulated <n3>')
legend('<n3> in simulated data', 'expected <n3>')
legend boxoff

subplot(1,3,2)
plot(tsamp, n_4_data(1:end), 'g*')
hold on
plot(tsamp, n4_C, 'b-', 'LineWidth',2)
xlabel('time (hours)')
ylabel('<n4>')
title('Expected vs.  simulated <n4>')
legend('<n4> in simulated data', 'expected <n4>')
legend boxoff

subplot(1,3,3)
plot(tsamp, var4_data(1:end), 'r-')
hold on
plot(tsamp, v4_C, 'b-', 'LineWidth',1)
xlabel('time (hours)')
ylabel('Var4')
title('Expected vs.  simulated Var4')
legend('<n4> in simulated data', 'expected <n4>')
legend boxoff
%% Check your function

p = horzcat(b,d,A);
C= moments_Allee_ODE(p, N0, V0, tsamp);

figure;
plot(tsamp, C(:,6),'b-')
hold on
plot(tsamp, var4_data, 'g*')
xlabel('time (hours)')
ylabel('Variance')
title('Simulated fourth moment')
legend('fourth moment SSA fxn', 'fourth moment simulated data')
legend boxoff
%% Test mufxn
ttry=tsamp(4:end);
figure;
plot(tsamp, mu_data)
hold on
plot(ttry, mu_fxnA(p,ttry, N0))
title ('mu comparison')

figure;
plot(tsamp, var_data)
hold on
plot(ttry, V_fxnA(p, ttry,N0, V0))
title ('variance comparison')
figure;
plot(tsamp, var4_data)
hold on
plot(ttry, V4_fxnA(p,ttry,N0,V0))
title('variance 4 comparison')

%% Perform Bayesian parameter estimation
% Want to fit parameters b, d, A, and  both sigmas
t=tsamp(2:end); % only want to fit data after initial condition
N= length(t);
% For now, start by assuming sigma = 1 and fit for b and d


pfxform5 = @(pval)[1 1 1 1 1].*log(pval); %'forward' parameter transform into Reals
pfxform = @(pval)[1 1 1].*log(pval); %'forward' parameter transform into Reals
pbxform5 = @(phat)[1 1 1 1 1].*exp(phat); %'backward' parameter transform into model space
pbxform = @(phat)[1 1 1].*exp(phat); %'backward' parameter transform into model space
yfxform = @(y)log(y); % 'forward' transform for data and model output
ybxform = @(yhat)exp(yhat); % 'inverse' transform for data and model output

% Two observations to fit on, mu_data and var_data
% both use tsamp
modelfun_mu5= @(p)mu_fxnA(p(1:3),t, N0); % <n> for params p = b,d, & A
modelfun_V5= @(p)V_fxnA(p(1:3),t, N0, V0); % variance as a function of parameters
modelfun_mu= @(p)mu_fxnA(p,t, N0); % <n> for params p = b,d, & A
modelfun_V= @(p)V_fxnA(p,t, N0, V0); % variance as a function of parameters
var_in_mean5 =  @(p)(1/N).*(V_fxnA(p(1:3), t, N0, V0)) + (p(4).^2).*ones(length(t),1); % variance in mean
var_in_var5 = @(p)(1/N).*(V4_fxnA(p(1:3),t, N0, V0)-(((N-3)./(N-1)).*(V_fxnA(p(1:3), t, N0, V0).^2))) + ((p(5).^2).*ones(length(t),1));
var_in_mean =  @(p)(1/N).*(V_fxnA(p, t, N0, V0)); % variance in mean
var_in_var = @(p)(1/N).*(V4_fxnA(p,t, N0, V0)-(((N-3)./(N-1)).*(V_fxnA(p, t, N0, V0).^2)));
% Initial guess
bguess = b+.001;
dguess = d+.001;
Aguess = A+.5;

theta = [bguess,dguess, Aguess];

sig1 = 0.1;
sig2 =0.1;

theta5 = horzcat(theta, sig1, sig2);

% Don't fit on first time point because variance will be 0
v_data= var_data(2:end);
m_data = mu_data(2:end);
ydata = horzcat(m_data,v_data);
mudata = m_data;
vardata= v_data;
%% Redo log likelihood function not using log transforms!
p5 = horzcat(p, sig1, sig2);
%J5= @(phat)((0.5*sum((log(2*pi*var_in_mean5(pbxform5(phat)))) + (((yfxform(modelfun_mu5(pbxform5(phat))) - yfxform(mudata))./sqrt(var_in_mean5(pbxform5(phat)))).^2))) +...
 %(0.5*sum((log(2*pi*var_in_var5(pbxform5(phat)))) + (((yfxform(modelfun_V5(pbxform5(phat))) - yfxform(vardata))./sqrt(var_in_var5(pbxform5(phat)))).^2))));

Jnew = @(phat) (sum(((mudata-modelfun_mu5(pbxform5(phat))).^2)./(2.*sqrt(var_in_mean5(pbxform5(phat)))) + log(sqrt(var_in_mean5(pbxform5(phat)))) + 0.5*log(2*pi))+...
    sum(((vardata-modelfun_V5(pbxform5(phat))).^2)./(2.*sqrt(var_in_var5(pbxform5(phat)))) + log(sqrt(var_in_var5(pbxform5(phat)))) + 0.5*log(2*pi)));
J_t = @(phat) ((((mudata-modelfun_mu5(pbxform5(phat))).^2)./(2.*sqrt(var_in_mean5(pbxform5(phat)))) + log(sqrt(var_in_mean5(pbxform5(phat)))) + 0.5*log(2*pi))+...
    (((vardata-modelfun_V5(pbxform5(phat))).^2)./(2.*sqrt(var_in_var5(pbxform5(phat)))) + log(sqrt(var_in_var5(pbxform5(phat)))) + 0.5*log(2*pi)));

objfun_Jnew = @(phat)(Jnew(phat));
lb = [0 0 0 0.1 0.1];
ub = [ Inf Inf Inf Inf Inf];
[phatbest_Jnew,fval,exitflag] = fminsearchbnd(objfun_Jnew, pfxform5(theta5), pfxform5(lb), pfxform5(ub));% need a better way to search parameter space,i.e. montecarlo
negLLguessnew= objfun_Jnew(pfxform5(theta5))
negLLtruenew= objfun_Jnew(pfxform5(p5))
negLLfit= objfun_Jnew(phatbest_Jnew)
params_best_Jnew= pbxform5(phatbest_Jnew)

% look at it in time
negLLtruet=J_t(pfxform5(p5));
negLLguesst=J_t(pfxform5(theta5));
negLLfitt = J_t(phatbest_Jnew);
figure;
plot(t, negLLtruet,'*', 'LineWidth',2)
hold on
plot(t, negLLguesst, 'LineWidth',2)
plot(t, negLLfitt, 'LineWidth', 2)
xlabel ('time')
ylabel('contribution to negative LL')
title('Negative LL in time')
legend('true params', 'initial guess', 'fit param')

%% Redo without fitting for sigu and sig var


Jnew = @(phat) (sum(((mudata-modelfun_mu(pbxform(phat))).^2)./(2.*sqrt(var_in_mean(pbxform(phat)))) + log(sqrt(var_in_mean(pbxform(phat)))) + 0.5*log(2*pi))+...
    sum(((vardata-modelfun_V(pbxform(phat))).^2)./(2.*sqrt(var_in_var(pbxform(phat)))) + log(sqrt(var_in_var(pbxform(phat)))) + 0.5*log(2*pi)));
J_t = @(phat) ((((mudata-modelfun_mu(pbxform(phat))).^2)./(2.*sqrt(var_in_mean(pbxform(phat)))) + log(sqrt(var_in_mean(pbxform(phat)))) + 0.5*log(2*pi))+...
    (((vardata-modelfun_V(pbxform(phat))).^2)./(2.*sqrt(var_in_var(pbxform(phat)))) + log(sqrt(var_in_var(pbxform(phat)))) + 0.5*log(2*pi)));

objfun_Jnew = @(phat)(Jnew(phat));
lb = [0 0 0 ];
ub = [ Inf Inf Inf ];
[phatbest_Jnew,fval,exitflag] = fminsearchbnd(objfun_Jnew, pfxform(theta), pfxform(lb), pfxform(ub));% need a better way to search parameter space,i.e. montecarlo
negLLguessnew= objfun_Jnew(pfxform(theta))
negLLtruenew= objfun_Jnew(pfxform(p))
negLLfit= objfun_Jnew(phatbest_Jnew)
params_best_Jnew= pbxform(phatbest_Jnew)

% look at it in time
negLLtruet=J_t(pfxform(p));
negLLguesst=J_t(pfxform(theta));
negLLfitt = J_t(phatbest_Jnew);
figure;
plot(t, negLLtruet,'*', 'LineWidth',2)
hold on
plot(t, negLLguesst, 'LineWidth',2)
plot(t, negLLfitt, 'LineWidth', 2)
xlabel ('time')
ylabel('contribution to negative LL')
title('Negative LL in time')
legend('true params', 'initial guess', 'fit param')

mu_fitJst = mu_fxnA(params_best_Jnew,t, N0);
V_fitJst= V_fxnA(params_best_Jnew,t,N0, V0);

figure;
subplot(1,2,1)
plot(tsamp, mu_data, 'r*')
hold on
plot(t, mu_fitJst, 'b--', 'LineWidth',3)
plot(tsamp, mu_C, 'k-', 'LineWidth',2)
xlabel('time (hours)')
ylabel('<n>')
title('Bayes fit mean')
legend('data', 'fit', 'expected')
legend boxoff


subplot(1,2,2)
plot(tsamp, var_data, 'g*')
hold on
plot(t, V_fitJst, 'b--', 'LineWidth',3)
plot(tsamp, v2_C, 'k-', 'LineWidth',2)
xlabel('time (hours)')
ylabel('Variance')
title('Bayes fit variance')
legend('data', 'fit', 'expected')
legend boxoff


%% Troubleshoot by looking at negLL in time
J_t= @(phat)((0.5*((log(2*pi*var_in_mean(pbxform(phat)))) + (((yfxform(modelfun_mu(pbxform(phat))) - yfxform(mudata))./sqrt(var_in_mean(pbxform(phat)))).^2))) +...
 (0.5*((log(2*pi*var_in_var(pbxform(phat)))) + (((yfxform(modelfun_V(pbxform(phat))) - yfxform(vardata))./sqrt(var_in_var(pbxform(phat)))).^2))));
% using 4th order variance in data
% J_t= @(phat)((0.5*((log(2*pi*var_in_mean(pbxform(phat)))) + (((yfxform(modelfun_mu(pbxform(phat))) - yfxform(mudata))./sqrt(var_in_mean(pbxform(phat)))).^2))) +...
%  (0.5*((log(2*pi*var_in_var_data)) + (((yfxform(modelfun_V(pbxform(phat))) - yfxform(vardata))./sqrt(var_in_var_data)).^2))));
negLLtruet=J_t(pfxform(p));
negLLguesst=J_t(pfxform(theta));
negLLfind = J_t(phatbest_Jsigt);
figure;
plot(t, negLLtruet,'*', 'LineWidth',2)
hold on
plot(t, negLLguesst, 'LineWidth',2)
plot(t, negLLfind, 'LineWidth', 2)
xlabel ('time')
ylabel('contribution to negative LL')
title('Negative LL in time')
legend('true params', 'initial guess', 'fit param')
%% Redo MCMC Search algorithm
% 1. Uniformly sample domain and calculate negLL at each point
% 2. Take lowest negLL as initial guess 
% 3. Random walk to search for better neg LL

% set domain of b, d, and A
bvec = [ 0: 0.005: 0.05];
dvec = [ 0: 0.001: 0.01];
Avec = [ 0: 1:10];

% make a 3D mesh
[B,D,A] = meshgrid(bvec, dvec, Avec); % 3D grid of parameters
Bflat = reshape(B,1,[]); % vector of bs
Dflat = reshape(D,1,[]); % vector of ds
Aflat = reshape(A,1, []); % vector of As

% run loop through vector of parameters and calculate negLL

for j = 1:length(Bflat)
    % set params
    pguess(j,1:3) = horzcat(Bflat(j),Dflat(j),Aflat(j));
    negLL(j)= objfun_Jnew(pfxform(pguess(j,1:3)));
    pguess(j,4) = negLL(j);   
end

NEGLL = reshape(negLL, size(A));
[Jinit,imin] = min(pguess(:,4));
pinit = pguess(imin,1:3)

figure;
hold off;
surf(bvec,dvec,NEGLL(:,:,3));
hold on
plot3(pinit(1), pinit(2), Jinit, 'r*', 'LineWidth',8)
xlabel( 'b')
ylabel('d')
zlabel('negLL')
title('Initial parameter space search')



%%  Write Markov-Chain Monte Carlo Search Algorithm


J_true= objfun_Jnew( pfxform(p));
J_curr = Jinit; % set current negLL value to inital value to initialize
count_accepted = 0; % start with no changes accepted

% Initialize a vector to store b and d params in
nruns = 10000;
store_params = zeros(nruns, 5); % store J, b,  d, sigu, and sig vfor each run

% Outside the loop: initialize temperature and iterations
T0= 20;
k = 1:1:nruns;
T= T0*exp(-4*k./(nruns)); % sets up gradual cooling

 prob_test= exp((J_curr-J_true)./T(1))

% check that temperature annealing is working
figure;
hold off
set(gca,'LineWidth',1.5,'FontSize',12);
plot(k, T, 'LineWidth', 3)
xlabel('Markov step')
ylabel('Temperature')
title('Temperature Annealing')

bstep = 0.001; % step size for searching birth param
dstep = 0.0002; % step size for searchign death param
Astep = 0.2;
% sigustep = 0.005;
% sigvstep = 0.005;
% set initial guess as b1 and d1
b1=pinit(1);
d1=pinit(2);
A1= pinit(3);
% sigu1= theta5(4);
% sigv1 = theta5(5);
store_acc_params = [];
store_acc_params(1,1) = Jinit;
store_acc_params(1,2) = b1;
store_acc_params(1,3) = d1;
store_acc_params(1,4) = A1;
%%
theta1 = pinit;
for k = 1:nruns
    r = rand;
    % set up to only change one parameter at a time, use rand to choose the
    % parameter to update
    if r<0.33
        % update b
        b2 = b1 + bstep*(2*rand-1);
        d2 = d1;
        A2 = A1;
    end
    if r>0.33 && r<0.67
        % update d
        d2 = d1 + dstep*(2*rand-1);
        b2 = b1;
        A2 = A1;
    end
    if r>0.67
        % update A
        A2 = A1 + Astep*(2*rand-1);
        b2 = b1;
        d2 = d1;
    end
    
%     b2 = b1 + bstep*(2*rand-1);
%     d2 = d1 + dstep*(2*rand-1);
%     A2 = A1 + Astep*(2*rand-1);

    
% Constrain search region to domain
    if b2<0 
        b2=0;
    end
    if b2>.05
        b2 = .05;
    end
    if d2<0
        d2=0;
    end
    if d2>.01
        d2 = .01;
    end
    if A2<0
        A2=0;
    end
    if A2>10
        A2 = 10;
    end
    
    % find the neg LL of the new params
    theta2 = horzcat(b2,d2,A2);
    J_new = objfun_Jnew( pfxform(theta2));
    % store the  negLL and searched parameters
    store_params(k,1) = J_new;
    store_params(k,2) =b2;
    store_params(k,3) = d2;
    store_params(k,4)=A2;
%     store_params(k,5) = sigu2;
%     store_params(k,6) = sigv2;
    
    prob(k) = exp((J_curr-J_new)./T(k));
    % if Jcurr > Jnew, Jnew is better--the numerator in the exponent is positive &
    % prob>1--> change will always be accepted
    % if Jcurr < Jnew, numerator is negative, & prob <1 --> chance change
    % will be accepted
    
    if rand < prob(k) % if true, accept the change!
        b1 = b2;
        d1 = d2;
        A1 = A2;
%         sigu1= sigu2;
%         sigv1 = sigv2;
        theta1 = theta2;
        J_curr = J_new;
        count_accepted = count_accepted +1;
        % decrease search step size
        if r<0.33
            % update b step
            bstep = 0.999*bstep;
        end
        if r>0.33 && r<0.67
            % update d step
             dstep = 0.999*dstep;
        end
        if r>0.67
            % update A
              Astep = 0.999*Astep;
        end

        params(1,1) = J_curr;
        params(1,2)= b1;
        params(1,3)= d1;
        params(1,4)= A1;
        store_acc_params= vertcat(store_acc_params, params);
    else
        % increase search step size
        if r<0.33
            % update b step
            bstep = 1.001*bstep;
        end
        if r>0.33 && r<0.67
            % update d step
             dstep = 1.001*dstep;
        end
        if r>0.67
            % update A
              Astep = 1.001*Astep;
        end

    end
       
end
%%
params_best_MC = store_acc_params(end,2:4)
[lowest_LL,ind] = min(store_acc_params(:,1));
lowestLLparams = store_acc_params(ind, 2:4)
diff_trueMC= params_best_MC(1:3) - p;
diff_trueML=lowestLLparams-p;
pct_error(:) = 100*(abs(diff_trueMC(1:3))./p)
pct_errorML(:) = 100*(abs(diff_trueML(1:3))./p)
negLLMC = J_t(pfxform(params_best_MC));
figure;
plot(t, negLLtruet, 'LineWidth',2)
hold on
plot(t, negLLfitt, 'LineWidth',2)
plot(t, negLLMC, 'LineWidth', 2)
xlabel ('time')
ylabel('contribution to negative LL')
title('Negative LL in time')
legend('true params', 'fminsearch', 'MC fit Params')
%%

figure;
plot(0:1:count_accepted, store_acc_params(:,1))
xlabel('iterations')
ylabel('negLL')
title('negLL')

 figure;
 plot(0:1:count_accepted, store_acc_params(:,2), 'g.')
 hold on
 plot(0:1:count_accepted, store_acc_params(:,3), 'b.')
 %plot(1:1:nruns, store_acc_params(:,4), 'r.')
%  plot(1:1:nruns, store_acc_params(:,4),'r.')
%   plot(1:1:nruns, store_acc_params(:,5),'m.')
 legend( 'birth rates', 'death rates')
 xlabel('iterations')
 ylabel('rate')
 title('Accepted birth, death, and A parameters')
 
 figure;
 plot(0:1:count_accepted, store_acc_params(:,4), 'r.')
%  plot(1:1:nruns, store_acc_params(:,4),'r.')
%   plot(1:1:nruns, store_acc_params(:,5),'m.')
 legend('A')
 xlabel('iterations')
 ylabel('rate')
 title('Accepted A parameters')
 
 figure;
 plot(1:1:nruns, store_params(:,2), 'g.')
 hold on
 plot(1:1:nruns, store_params(:,3), 'b.')
 %plot(1:1:nruns, store_params(:,4), 'r.')
%  plot(1:1:nruns, store_params(:,4), 'r.')
%   plot(1:1:nruns, store_params(:,5), 'm.')
 legend( 'birth rates', 'death rates', 'A','\sigma\mu', '\sigmavar')
 xlabel('iterations')
 ylabel('rate')
 title('Explored birth, death, and A parameters')
 
 figure;
 plot(1:1:nruns, prob, 'r.')
 ylim([0 1])
%%
figure;
subplot(1,3,1)
hist(store_acc_params(:,2))
ylabel('frequency')
xlabel('b')
subplot(1,3,2)
hist(store_acc_params(:,3))
ylabel('frequency')
xlabel('d')
subplot(1,3,3)
hist(store_acc_params(:,4))
ylabel('frequency')
xlabel('A')
% subplot(1,5,4)
% hist(store_acc_params(:,5))
% ylabel('frequency')
% xlabel('\sigmamu')
% subplot(1,5,5)
% hist(store_acc_params(:,6))
% ylabel('frequency')
% xlabel('\sigmavar')
%%

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
scatter(store_acc_params(:,2), store_acc_params(:,4), pointsize, likelihoodval,'filled');
colorbar
xlabel('b')
%xlim([ b-.005 b + .005])
%ylim([0 d+.005])
ylabel('A')
title(' b versus A colored by likelihood for accepted parameters')

figure;
scatter(store_acc_params(:,3), store_acc_params(:,4), pointsize, likelihoodval,'filled');
colorbar
xlabel('d')
%xlim([ b-.005 b + .005])
%ylim([0 d+.005])
ylabel('A')
title(' d versus A colored by likelihood for accepted parameters')

%% Try to make 3D plot colored by likelihood

figure;
scatter3(store_acc_params(:,2),store_acc_params(:,3), store_acc_params(:,4), pointsize, likelihoodval,'filled')
xlabel ('b')
ylabel('d')
zlabel('A')
figure;
scatter3(store_params(:,2),store_params(:,3), store_params(:,4), pointsize, likelihoodvalall,'filled')
xlabel ('b')
ylabel('d')
zlabel('A')
%%

std_dev = std(store_acc_params(:, 2:end))
credible_interval = prctile(store_acc_params(:, 2:end),[ 2.5, 97.5])

CI(1,:)= params_best_MC - 1.96*std_dev;
CI(2,:) = params_best_MC + 1.96*std_dev;
figure;
subplot(1,3,1)
plot(store_acc_params(:,2), likelihoodval,'ro')
hold on
xlabel ('b')
ylabel('likelihood')
title (['b=',num2str(round(params_best_MC(1), 4)),' +/- ', num2str(round(1.96*std_dev(1),4))])
subplot(1,3,2)
plot(store_acc_params(:,3), likelihoodval,'bo')
xlabel ('d')
ylabel('likelihood')
title (['d=',num2str(round(params_best_MC(2), 4)),' +/- ', num2str(round(1.96*std_dev(2),4))])
subplot(1,3,3)
plot(store_acc_params(:,4), likelihoodval,'go')
xlabel ('A')
ylabel('likelihood')
title (['A=',num2str(round(params_best_MC(3), 4)),' +/- ', num2str(round(1.96*std_dev(3),4))])

%%

std_dev = std(store_acc_params(:, 2:end))
