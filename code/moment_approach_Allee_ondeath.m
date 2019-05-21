% Allee effect on death--> results in the exact same set of ODEs
% therefore there's no way to distinguish if Allee effect is due to birth
% alone or death alone, as the second moment and variance come out
% identical
%Moment- Approach fitting of stochastic birth-death-Allee process for small cell
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
b = 0.0092;
d = 0.001;
delta= b+d;
A = 2;
%birth_n = b*N; % birth hazard function
%death_n = d*N; % death hazard function
num_samps = 5000;
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
    birth_n = (b*N); % birth 
    if birth_n <0
        birth_n = 0;
    end
    death_n = (d*N+(b-d)*A); % death 
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
xlabel('time (hours)')
ylabel('Number of cells')
title(['Simulated N(t) trajectories for b=', num2str(b), ' & d=', num2str(d)])


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
ylabel('<n> expected')
xlim([0, tsamp(end)])
title('Mean value of N in time for N0=5, A = 2')

%% Compare moment approach approximations
V0=0;
C_init(1)=N0;
C_init(2)=N0.^2;
C_init(3) = V0;
C_init(4)= N0.^3;
C_init(5)= N0.^4;
C_init(6) = V0;


f = @(t,C) [((b-d)*C(1)-(b-d)*A); % dn/dt
             2*C(2).*(b-d) - (2*C(1).*(b-d)*A) + (C(1).*(b+d))+((b-d)*A);   % dn2/dt
            (2*C(2).*(b-d) - (2*C(1).*(b-d)*A) + (C(1).*(b+d))+((b-d)*A)-(2.*C(1).*((b-d)*C(1)-(b-d)*A))); % dV/dt
            (3*C(4).*(b-d)) + (3*C(2).*(b+d))-(3*C(2).*(b-d)*A) + (3*C(1).*(b-d)*A) + ((b-d).*C(1))-((b-d)*A);%dn3dt
            (4*C(5).*(b-d)) + (6*C(4).*(b+d))-(4*C(4)*(b-d)*A)+(4*C(2).*(b-d))+(6*C(2).*(b-d)*A)+ (C(1).*(b+d))+...
            (4*C(1).*(b-d)*A)+((b-d)*A); %dn4dt
            (4*C(5).*(b-d)) + (6*C(4).*(b+d))-(4*C(4).*(b-d)*A)+(4*C(2).*(b-d))+(6*C(2).*(b-d)*A)+ (C(1).*(b+d))+...
            (4*C(1).*(b-d)*A)+((b-d)*A)-(4.*((C(1)).^3)*(((b-d).*C(1)-(b-d)*A)))];

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
plot(tsamp, (n2_C-((mu_C).^2)))
hold on
plot(tsamp, v2_C, 'r.')


figure;
subplot(1,2,1)
plot(tsamp, mu_data(1:end), 'r*')
hold on
plot(tsamp, mu_C, 'k-', 'LineWidth',2)
xlabel('time (hours)')
ylabel('mean cell number')
title('Strong Allee on death mean', 'FontSize', 14)
legend('mean n simulated data', 'expected mean n')
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
title('Strong Allee on death variance', 'FontSize', 14)
legend('Variance in simulated data', 'expected Variance')
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

%% Plot Sampled Data from Stochastic Simulations alongside expected values of mean and variance
mu_data = mean(Nsamp,2);
n_2_data = mean((Nsamp.^2),2);
n_3_data = mean((Nsamp.^3),2);
n_4_data = mean((Nsamp.^4),2);
var_data = n_2_data - ((mu_data).^2);
var4_data = n_4_data - ((mu_data).^4);
%var_data = (std(Nsamp,0,2).^2);

figure;
subplot(1,3,1)
plot(tsamp, mu_data(1:end), 'r*')
hold on
plot(tsamp, mu_C, 'k-', 'LineWidth',2)
xlabel('time (hours)')
ylabel('<n>')
title('Expected vs. simulated <n> Allee model')
legend('mean n simulated data', 'expected mean n')
legend boxoff

subplot(1,3,2)
plot(tsamp, n_2_data(1:end), 'g*')
hold on
plot(tsamp, n2_C, 'b.', 'LineWidth',2)
xlabel('time (hours)')
ylabel('<n2>')
title('Expected vs. simulated <n2>')
legend('<n2> in simulated data', 'expected <n2>')
legend boxoff
%%
subplot(1,3,3)
plot(tsamp, var_data(1:end), 'm*')
hold on
plot(tsamp, v2_C, 'b.', 'LineWidth',2)
xlabel('time (hours)')
ylabel('Variance')
title('Expected vs.  simulated Variance')
legend('Variance in simulated data', 'expected Variance')
legend boxoff