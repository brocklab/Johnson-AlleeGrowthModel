% Moment- Approach fitting of stochastic birth-death process for small cell
% number data

% This code is first attempt to use moment-approach (so mean and
% variance)of N(t) trajectories to fit small cell number well data.

% We will start by simulating stochastic cell trajectories using a set b
% and d, and see if we can use the following equations for <n(t)> and V(t)
% to perform parameter estimation in Bayesian framework
% For random birth death process we have:

%d<n>/dt=(b-d)<n>
%<n>= n0exp[(b-d)t]
%dV/dt =  (b-d) V + (b+d) <n>
%V= [(b+d)/2(b-d)]*n0*exp((b-d)t)*[exp((b-d)t)-1]

close all; clear all; clc

% Start by generating 20 n0=1 trajectories
% assume birth rate = 0.0238, death rate = 0.005
% Set up time and N counter
Ninit = [1 2 4 8 16 32 48 64 128];
Ninit = 1;
for i = 1:length(Ninit)
b = 0.0233 + .0005; % birth rate .0005*i to get b+d to increase
d = 0.0045 + .0005; % death rate
delta(i)= b+d;
%birth_n = b*N; % birth hazard function
%death_n = d*N; % death hazard function
num_samps = 1000;
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
    birth_n = b*N; % birth 
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
sigmaT = 0.1; % guess that on average count is off by one quarter of a cell?
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
% use this figure when want to see the effect of increased both b and d
% alpha) but holding b-d constant
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
% take this figure when holding b-d constant and just want to see effect on
% P est
% figure;
% plot(Ninit(1:5), P_est_theor(1:5),'-', 'LineWidth',2)
% hold on
% plot(Ninit(1:5), P_takeoff(1:5), '*', 'LineWidth',2)
% xlabel('Initial cell number')
% ylabel('P_est')
% title('P_{est} vs. N_{0} for constant \delta and b-d')
% legend('expected', 'measured from simulation')
% legend boxoff
%% Plot simulated cell number trajectories out to minimum time reached
figure;
hold off
for j = 1:num_samps
plot(tstate(:,j), state(:,j))
hold on
end
xlim([0, tmin])
xlabel('time (hours)')
ylabel('Number of cells')
title(['Simulated N(t) trajectories for b=', num2str(b), ' & d=', num2str(d)])


%% UNIFORM SAMPLING from stochastic trajectories 
%Want to smooth by sampling from N at tsamp
tstart=0;
tint=5;

 tsamp = tstart:tint:100+tstart;
%tsamp = tstart:tint:round(tmin);
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
% Plot Sampled Data from Stochastic Simulations alongside expected values of mean and variance

%% Expected values for <n>, <n^2>, and variance using numerical solver
tsamp =tsamp(1:end);
mu_data = mean(Nsamp,2);
n_2_data = mean((Nsamp.^2),2);
n_3_data = mean((Nsamp.^3),2);
n_4_data = mean((Nsamp.^4),2);
var_data = n_2_data - ((mu_data).^2);
var4_data = n_4_data - ((mu_data).^4);



V0 = 0;
C_init(1)=N0;
C_init(2)=N0.^2;
C_init(3)= N0.^3;
C_init(4)= N0.^4;
C_init(5) = V0;
C_init(6)= V0;

f = @(t,C) [(b-d)*C(1); % dn/dt
            2*C(2)*(b-d) + C(1)*(b+d);   % dn2/dt
            3*C(3)*(b-d) + 3*C(2)*(b+d) + C(1)*(b-d); %dn3/dt
            4*C(4)*(b-d)+ 6*C(3)*(b+d) + 4*C(2)*(b-d) + C(1)*(b+d); % dn4/dt
            2*C(5)*(b-d) + (b+d)*C(1); %dV2dt
            4*C(4)*(b-d)+ 6*C(3)*(b+d) + 4*C(2)*(b-d) + C(1)*(b+d)-4*(C(1).^3)*(b-d)*C(1)];%dV4dt
            %4*C(4)*(b-d)+6*C(3)*(b+d) - 12*C(3)*(b-d)+ 4*C(2)*(b-d) - 12*C(2)*(b+d)-4*C(1)*(b-d)+C(1)*(b+d)]; %dV4dt
options1 = odeset('Refine',1);  
options = odeset(options1,'NonNegative',1:6);
[t,C]=ode45(f, tsamp,C_init, options);
mu_C= C(:,1);
n2_C= C(:,2);
n3_C = C(:,3);
n4_C = C(:,4);
v2_C = C(:,5);
v4_C=C(:,6);
%%
mu_data = mean(Nsamp,2);
n_2_data = mean((Nsamp.^2),2);
n_3_data = mean((Nsamp.^3),2);
n_4_data = mean((Nsamp.^4),2);
var_data = n_2_data - ((mu_data).^2);
var4_data = n_4_data - ((mu_data).^4);
%var_data = (std(Nsamp,0,2).^2);

figure;
plot(tsamp, mu_data(1:end), 'r*')
hold on
plot(tsamp, mu_C, 'k-', 'LineWidth',2)
xlabel('time (hours)')
ylabel('Mean (<n>')
title('Mean (t)')
legend('simulated mean', 'model mean')
legend boxoff

figure;
plot(tsamp, var_data(1:end), 'g*')
hold on
plot(tsamp, v2_C, 'k-', 'LineWidth',2)
xlabel('time (hours)')
ylabel('Variance')
title('Variance (t)')
legend('simulated variance', 'model variance')
legend boxoff

%%


figure;
subplot(1,4,1)
plot(tsamp, mu_data(1:end), 'r*')
hold on
plot(tsamp, mu_C, 'k-', 'LineWidth',2)
xlabel('time (hours)')
ylabel('<n>')
title('Expected vs. simulated mu')
legend('mean n simulated data', 'expected mean n')
legend boxoff

subplot(1,4,2)
plot(tsamp, n_2_data(1:end), 'g*')
hold on
plot(tsamp, n2_C, 'b-', 'LineWidth',2)
xlabel('time (hours)')
ylabel('<n2>')
title('Expected vs. simulated <n2>')
legend('<n2> in simulated data', 'expected <n2>')
legend boxoff

subplot(1,4,3)
plot(tsamp, n_3_data(1:end), 'm*')
hold on
plot(tsamp, n3_C, 'b-', 'LineWidth',2)
xlabel('time (hours)')
ylabel('<n3>')
title('Expected vs. simulated <n3>')
legend('<n3> in simulated data', 'expected <n3>')
legend boxoff

subplot(1,4,4)
plot(tsamp, n_4_data(1:end), 'm*')
hold on
plot(tsamp, n4_C, 'b-', 'LineWidth',2)
xlabel('time (hours)')
ylabel('<n4>')
title('Expected vs.  simulated <n4>')
legend('<n4> in simulated data', 'expected <n3>')
legend boxoff

figure;
subplot(1,2,1)
plot(tsamp, var_data(1:end), 'm*')
hold on
plot(tsamp, v2_C, 'b-', 'LineWidth',2)
xlabel('time (hours)')
ylabel('Variance')
title('Expected vs.  simulated Variance')
legend('Variance in simulated data', 'expected Variance')
legend boxoff

subplot(1,2,2)
plot(tsamp, var4_data(1:end), 'r*')
hold on
plot(tsamp, v4_C, 'b-', 'LineWidth',2)
xlabel('time (hours)')
ylabel('Variance')
title('Expected vs.  simulated 4th order var')
legend('4th order var in simulated data', 'expected 4th order var')
legend boxoff
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
%%
%var_data = (std(Nsamp,0,2).^2);
p = horzcat(b,d);
mu_n = mu_fxn(p,tsamp, N0);
V_n= V_fxn(p,tsamp,N0, V0);

n_2_n= second_mom_fxn(p,tsamp,N0);
%V_n = n_2_n - (mu_n).^2;

figure;
hold on 
for j = 1:num_samps
plot(tsamp,Nsamp(1:end,j), 'r.');
hold on
plot(tstate(:,j), state(:,j),'b.');
end
plot(tsamp, mu_n, 'k-', 'LineWidth',3)
xlabel ('time (hours)')
ylabel('<n> expected')
xlim([0, tsamp(end)])
title('Mean value of N in time for N0=1 closed form')

figure;
plot(tsamp, mu_data, 'r*')
hold on
plot(tsamp, mu_n, 'k-', 'LineWidth',2)
xlabel('time (hours)')
ylabel('<n>')
title(['Expected vs. simulated mean cell number, N_{0}=',num2str(Ninit), ', N_{runs}=',num2str(num_samps)])
legend('mean n simulated data', 'expected mean n')
legend boxoff


figure;
plot(tsamp, var_data(1:end), 'g*')
hold on
plot(tsamp, var_data2(1:end), 'b*')
plot(tsamp, V_n, 'k-', 'LineWidth',2)
xlabel('time (hours)')
ylabel('Variance')
title('Expected vs. simulated variance closed form')
legend('variance in simulated data', 'expected variance')
legend boxoff


%% try your V4 fxn
%[V4] = V4_fxn(p, N0, tstart, tint);
 V4 =V4_fxn_ODE(p, N0, V0, tsamp);
%V4samp = interp1(tsamp,tODE,V4);
% Next need to write 
figure;
plot(tsamp, V4,'g-')
hold on
plot(tsamp, var4_data, 'b*')
xlabel('time (hours)')
ylabel('Variance')
title('Simulated fourth moment')
legend('fourth moment SSA fxn', 'fourth moment simulated data')
legend boxoff

%% Find mean and variance at each time point
mu_n = mu_fxn(p,tsamp, N0);
%V_n = n_2_n - (mu_n).^2;
V0=0;
V_n = V_fxn(p,tsamp(2:end),N0, V0);


mu_data = mean(Nsamp,2);
n_2_data = mean((Nsamp.^2),2);
var_data = n_2_data - ((mu_data).^2);
%var_data = (std(Nsamp,0,2).^2);

figure;

plot(tsamp, mu_data, 'r*')
hold on
plot(tsamp, mu_n, 'k-', 'LineWidth',2)
xlabel('time (hours)')
ylabel('<n>')
title(['Expected vs. simulated mean cell number, N_{0}=',num2str(Ninit), ', N_{runs}=',num2str(num_samps)])
legend('mean n simulated data', 'expected mean n')
legend boxoff

figure;
plot(tsamp, var_data, 'g*')
hold on
plot(tsamp(2:end), V_n, 'k-', 'LineWidth',2)
xlabel('time (hours)')
ylabel('Variance')
title('Expected vs. simulated variance')
legend('variance in simulated data', 'expected variance')
legend boxoff
%% Expected variance as b+d increases
for i = 1:5
    % increase b and d at every iteration
    b_test(i)= b+ .01*i;
    d_test(i)= d + .01*i;
    p_test(i,:) = vertcat(b_test(i), d_test(i));
    var_test(:,i) = V_fxn(p_test(i,:),tsamp,N0, V0);
    delta(i) = (b_test(i)+d_test(i))./2;
    P_takeoff(i) = (1-(d_test(i)./b_test(i)).^N0);
end
%
figure;
for i = 1:5
plot(tsamp, var_test(:,i), 'LineWidth',2)
hold on
text(tsamp(end-4), var_test(end-2,i), ['\delta=', num2str(delta(i)),', P_{est}=' num2str(round(P_takeoff(i),3))])
end
xlabel('time (hours)')
ylabel('Expected Variance')
title(['Variance over time for constant b-d =', num2str(b-d), ', increasing b+d(\delta)'])

figure;
plot(delta, P_takeoff, 'o-', 'LineWidth',2)
xlabel('\delta = b+d')
ylabel('P_{establishment}')
title('Theoretical \delta vs. P_est for constant b-d')
%% Perform Bayesian parameter estimation
% Want to fit parameters b, d, and sigma
t=tsamp(2:end); % only want to fit data after initial condition
N=num_samps;
% need to make N0 and V0 at time tsamp(2)=t0

% For now, start by assuming sigma = 1 and fit for b and d

pfxform = @(pval)[1 1].*log(pval); %'forward' parameter transform into Reals
pfxform4 = @(pval)[1 1 1 1].*log(pval); % for fitting with 4 params
pbxform = @(phat)[1 1].*exp(phat); %'backward' parameter transform into model space
pbxform4 = @(phat)[1 1 1 1].*exp(phat); % for fitting with 4 params
yfxform = @(y)log(y); % 'forward' transform for data and model output
ybxform = @(yhat)exp(yhat); % 'inverse' transform for data and model output

% Two observations to fit on, mu_data and var_data
% both use tsamp
modelfun_mu = @(p)mu_fxn(p, t, N0); % <n> for params (p) = b & d, vertical
modelfun_mu4 = @(p)mu_fxn(p(1:2), t, N0); % for fitting with sigmas
modelfun_V=@(p)V_fxn(p, t, N0, V0); % vertical
modelfun_V4= @(p)V_fxn(p(1:2),t,N0, V0);
var_in_mean =  @(p)(1/N).*(V_fxn(p, t, N0, V0)) + (sigmaT.^2).*ones(length(t),1); % vertical
var_in_mean4 = @(p)(1/N).*(V_fxn(p(1:2), t, N0, V0)) + ((p(3).^2)).*ones(length(t),1);
var_in_mean_data = (1/N).*var_data;


% To find fourth order moment, take parameters, simulate 100 trajectories,
% and find the corresponding v4


var_in_var = @(p)(1/N).*(V4_fxn_ODE(p,N0, V0, t)-(((N-3)./(N-1)).*(V_fxn(p, t, N0, V0).^2))) + 2*sigmaT;
var_in_var4=@(p)(1/N).*(V4_fxn_ODE(p(1:2),N0, V0, t)-(((N-3)./(N-1)).*(V_fxn(p(1:2), t, N0, V0).^2))) + ((p(4).^2).*ones(length(t),1));
var_in_var_data = (1/N).*(var_data4(2:end)-(((N-3)/(N-1)).*(var_data(2:end).^2)));
%var_in_var = @(p)(1/N).*(var_data4(2:end)- (((N-3)/(N-1)).*((V_fxn(p, t, N0, V0))'.^2)));

% Initial guess
bguess = b +.005;
dguess = d + .001;
theta = [bguess,dguess];

sig1 = 0.1;
sig2 = 0.1;

theta4 = horzcat(theta, sig1, sig2);

% Don't fit on first time point because variance will be 0
v_data= var_data(2:end);
m_data = mu_data(2:end);
n2_data = n_2_data;
ydata = horzcat(m_data,v_data);
mudata = m_data;
vardata= v_data;
%% Compute the negative log likelihood not using log transforms!

Jnew = @(phat) (sum(((mudata-modelfun_mu(pbxform(phat))).^2)./(2.*sqrt(var_in_mean(pbxform(phat)))) + log(sqrt(var_in_mean(pbxform(phat)))) + 0.5*log(2*pi))+...
    sum(((vardata-modelfun_V(pbxform(phat))).^2)./(2.*sqrt(var_in_var(pbxform(phat)))) + log(sqrt(var_in_var(pbxform(phat)))) + 0.5*log(2*pi)));
J_t = @(phat) ((((mudata-modelfun_mu(pbxform(phat))).^2)./(2.*sqrt(var_in_mean(pbxform(phat)))) + log(sqrt(var_in_mean(pbxform(phat)))) + 0.5*log(2*pi))+...
    (((vardata-modelfun_V(pbxform(phat))).^2)./(2.*sqrt(var_in_var(pbxform(phat)))) + log(sqrt(var_in_var(pbxform(phat)))) + 0.5*log(2*pi)));

objfun_Jnew = @(phat)(Jnew(phat));
lb = [0 0 ];
ub = [ Inf Inf ];
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

mu_fitJst = mu_fxn(params_best_Jnew,t, N0);
V_fitJst= V_fxn(params_best_Jnew,t,N0, V0);

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

%% Redo MCMC Search algorithm
% 1. Uniformly sample domain and calculate negLL at each point
% 2. Take lowest negLL as initial guess 
% 3. Random walk to search for better neg LL

% set domain of b, d, and A
bvec = [ 0: 0.005: 0.05];
dvec = [ 0: 0.001: 0.01];

% make a 3D mesh
[B,D] = meshgrid(bvec, dvec); % 3D grid of parameters
Bflat = reshape(B,1,[]); % vector of bs
Dflat = reshape(D,1,[]); % vector of ds

% run loop through vector of parameters and calculate negLL

for j = 1:length(Bflat)
    % set params
    pguess(j,1:2) = horzcat(Bflat(j),Dflat(j));
    negLL(j)= objfun_Jnew(pfxform(pguess(j,1:2)));
    pguess(j,3) = negLL(j);   
end

NEGLL = reshape(negLL, size(B));
[Jinit,imin] = min(pguess(:,3));
pinit = pguess(imin,1:2)

figure;
hold off;
surf(bvec,dvec,NEGLL(:,:));
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
store_params = zeros(nruns, 3); % store J, b,  d,for each run

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
dstep = 0.001; % step size for searchign death param

% set initial guess as b1 and d1
b1=pinit(1);
d1=pinit(2);


store_acc_params = [];
store_acc_params(1,1) = Jinit;
store_acc_params(1,2) = b1;
store_acc_params(1,3) = d1;

%%
theta1 = pinit;
for k = 1:nruns
    r = rand;
    % set up to only change one parameter at a time, use rand to choose the
    % parameter to update
    if r<0.5
        % update b
        b2 = b1 + bstep*(2*rand-1);
        d2 = d1;
    end
    if r>0.5
        % update d
        d2 = d1 + dstep*(2*rand-1);
        b2 = b1;
    end

    
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
    
    
    % find the neg LL of the new params
    theta2 = horzcat(b2,d2);
    J_new = objfun_Jnew( pfxform(theta2));
    % store the  negLL and searched parameters
    store_params(k,1) = J_new;
    store_params(k,2) =b2;
    store_params(k,3) = d2;

    
    prob(k) = exp((J_curr-J_new)./T(k));
    % if Jcurr > Jnew, Jnew is better--the numerator in the exponent is positive &
    % prob>1--> change will always be accepted
    % if Jcurr < Jnew, numerator is negative, & prob <1 --> chance change
    % will be accepted
    
    if rand < prob(k) % if true, accept the change!
        b1 = b2;
        d1 = d2;

        theta1 = theta2;
        J_curr = J_new;
        count_accepted = count_accepted +1;
        % decrease search step size
        if r<0.5
            % update b step
            bstep = 0.999*bstep;
        end
        if r>0.5
            % update d step
             dstep = 0.999*dstep;
        end
        
        params(1,1) = J_curr;
        params(1,2)= b1;
        params(1,3)= d1;
        store_acc_params= vertcat(store_acc_params, params);
    else
        % increase search step size
        if r<0.5
            % update b step
            bstep = 1.001*bstep;
        end
        if r>0.5
            % update d step
             dstep = 1.001*dstep;
        end

    end
       
end
%%
params_best_MC = store_acc_params(end,2:3)
[lowest_LL,ind] = min(store_acc_params(:,1));
lowestLLparams = store_acc_params(ind, 2:3)
diff_trueMC= params_best_MC(1:2) - p;
diff_trueML=lowestLLparams-p;
pct_error(:) = 100*(abs(diff_trueMC(1:2))./p)
pct_errorML(:) = 100*(abs(diff_trueML(1:2))./p)
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

figure;
plot(0:1:count_accepted, store_acc_params(:,1))
xlabel('iterations')
ylabel('negLL')
title('negLL')

 figure;
 plot(0:1:count_accepted, store_acc_params(:,2), 'g.')
 hold on
 plot(0:1:count_accepted, store_acc_params(:,3), 'b.')
 legend( 'birth rates', 'death rates')
 xlabel('iterations')
 ylabel('rate')
 title('Accepted birth & death parameters')
 
 
 figure;
 plot(1:1:nruns, store_params(:,2), 'g.')
 hold on
 plot(1:1:nruns, store_params(:,3), 'b.')
 legend( 'birth rates', 'death rates')
 xlabel('iterations')
 ylabel('rate')
 title('Explored birth, death, and A parameters')
 
 figure;
 plot(1:1:nruns, prob, 'r.')
 ylim([0 1])
p%%
figure;
subplot(1,2,1)
hist(store_acc_params(:,2))
ylabel('frequency')
xlabel('b')
subplot(1,2,2)
hist(store_acc_params(:,3))
ylabel('frequency')
xlabel('d')


% Plots of two parameters colored by likelihood
likelihoodval = exp(-store_acc_params(:,1));
likelihoodvalall = exp(-store_params(:,1));
pointsize = 30;
figure;
scatter(store_acc_params(:,2), store_acc_params(:,3), pointsize, likelihoodval,'filled');
colorbar
xlabel('b')
ylabel('d')
title(' b versus d colored by likelihood for accepted parameters')




