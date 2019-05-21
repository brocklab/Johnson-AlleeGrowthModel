% This script is a master script to run functions that give the necessary
% outputs from each model structure we are investigating. We want to both
% run simulations of the Gillespie algorithm for a given N and parameter
% set. We want to extract the summary statistics from each simulated
% data set, as well as the expected moments derived from the CME. Then, we
% can run all of the functions and check to see how they vary depending on
% both the parameters and the assumptions about the Allee effect on either
% the birth or death for the three deterministic model structures:
% 1. dN/dt = gN "bdmodel"
% 2. dN/dt = g(N-A) * Allee on birth, death, or both: "strAllb", "strAlld",
% "strAllbd"
% 3. dN/dt = gN(1-(A+tau)/(N+tau)) * Allee on birth, death, or both:
% "wkAllb", "wkAlld", "wkAllbd"

% Leads to 7 total models to look at.
% We want to probably make comprehensive sheet of mean and variances for
% each of these models, but in the fitting we need up to the 4th moment and
% 4th order variance.

% Function structure:
% [ Nsamp,Nstat(<n>, <n2>, <Var>, <n3>, <n4>, <V4>), Cstat(same
% order), musimdatavec, varsimdatavec, timesimdatavec] = 
% run_modelname(params, tsamp, N0, Nsim)


% Nsamp is a matrix of numTP(#of time points) rows x Nsim columns,
% x numN0 zs
% tsamp is the corresponding time vector of the Nstat
% Nstat is the moments and variances for each time point (so numTP x num
% moments columns x numN0 zs)
% Cstat is the same
% musimdatavec, varsimedatavec, and timesimdatavec are the concatenated
% versions for all NO of each of these.

% Let's set up some conditions we will set for all of our models
close all; clear all; clc

Nsim = 1000;
tsamp = [ 0:4:336];
b = 0.0092;
d = 0.001;
A = 3;
Ninit = [ 5; 10; 15]; 

%% Birth-death model only
% start by writing a function that generates Nsamp, Nstat, Cstat,
% musimdatavec, varsimdatavec, and timesimdatavec for the simple birth
% death model

paramsbd = [b,d];
%Nstat(<n>, <n2>, <Var>, <n3>, <n4>, <V4>)
[ Nsampbd,Nstatbd, Cstatbd, musimdatavecbd, varsimdatavecbd, timesimdatavecbd] = run_bdmodel(paramsbd, tsamp, Ninit, Nsim);

stats = {'mean', '2nd moment', 'variance', '3rd moment', '4th moment', '4th order variance'};
% Plot the moments for each initial cell number
for i = 1:3
    figure;
    hold off
for j = 1:6
    subplot(2,3,j)
    plot(tsamp, Nstatbd(:,j,i), '*', 'LineWidth', 1)
    hold on
    plot(tsamp, Cstatbd(:,j,i), '.', 'LineWidth',2)
    xlabel ('time (hours)')
    ylabel([stats{j}])
    legend([stats{j},' simulated data'], ['expected ', stats{j}])
    title([stats{j},' N_{0}= ', num2str(Ninit(i))])
    legend boxoff
end
end

%% Plot the mean and variance for all three initial cell numbers
figure;
plot(timesimdatavecbd, musimdatavecbd,'r*', 'LineWidth',2)
hold on
for i = 1:length(Ninit)
    plot(tsamp, Cstatbd(:,1,i), 'k.')
    text(tsamp(end-10), Cstatbd(end-10, 1, i), ['N_{0}=', num2str(Ninit(i))])
end
xlabel('time (hours)')
ylabel('mean number of cells')
title('Expected and simulated mean for different N_{0}')

figure;
plot(timesimdatavecbd, varsimdatavecbd,'g*', 'LineWidth',2)
hold on
for i = 1:length(Ninit)
    plot(tsamp, Cstatbd(:,3,i), 'k.')
    text(tsamp(end-10), Cstatbd(end-10, 3, i), ['N_{0}=', num2str(Ninit(i))])
end
xlabel('time (hours)')
ylabel('variance in number of cells')
title('Expected and simulated variance for different N_{0}')
%% Now repeat function making for the strong Allee model on birth
% need to do so for our three different stochastic structures:
A = 3;
paramsAb = [b,d, A];
%Nstat(<n>, <n2>, <Var>, <n3>, <n4>, <V4>)
[ NsampAb,NstatAb, CstatAb, musimdatavecAb, varsimdatavecAb, timesimdatavecAb] = run_strAllbmodel(paramsAb, tsamp, Ninit, Nsim);

stats = {'mean', '2nd moment', 'variance', '3rd moment', '4th moment', '4th order variance'};
% Plot the moments for each initial cell number
for i = 1:3
    figure;
    hold off
for j = 1:6
    subplot(2,3,j)
    plot(tsamp, NstatAb(:,j,i), '*', 'LineWidth', 1)
    hold on
    plot(tsamp, CstatAb(:,j,i), '.', 'LineWidth',2)
    xlabel ('time (hours)')
    ylabel([stats{j}])
    legend([stats{j},' simulated data'], ['expected ', stats{j}])
    title([stats{j},' N_{0}= ', num2str(Ninit(i))])
    legend boxoff
end
end


%% Compare the mean and variance for b-d and b-d-Ab model
figure;
plot(timesimdatavecbd, musimdatavecbd,'m*', 'LineWidth',2)
hold on
plot(timesimdatavecAb, musimdatavecAb,'r*', 'LineWidth',2)
for i = 1:length(Ninit)
    plot(tsamp, Cstatbd(:,1,i), 'k.')
    text(tsamp(end-10), Cstatbd(end-10, 1, i), ['N_{0}=', num2str(Ninit(i))])
end
hold on
for i = 1:length(Ninit)
    plot(tsamp, CstatAb(:,1,i), 'k.')
    text(tsamp(end-10), CstatAb(end-10, 1, i), ['N_{0}=', num2str(Ninit(i))])
end
legend('b-d model', 'b-d-Ab model')
xlabel('time (hours)')
ylabel('mean number of cells')
title(['Difference in means for A= ', num2str(A),', b & d constant'])

figure;
plot(timesimdatavecbd, varsimdatavecbd,'g*', 'LineWidth',2)
hold on
plot(timesimdatavecAb, varsimdatavecAb,'c*', 'LineWidth',2)
for i = 1:length(Ninit)
    plot(tsamp, Cstatbd(:,3,i), 'k.')
    text(tsamp(end-10), Cstatbd(end-10, 3, i), ['N_{0}=', num2str(Ninit(i))])
end
hold on
for i = 1:length(Ninit)
    plot(tsamp, CstatAb(:,3,i), 'k.')
    text(tsamp(end-10), CstatAb(end-10, 3, i), ['N_{0}=', num2str(Ninit(i))])
end
legend('b-d model', 'b-d-Ab model')
xlabel('time (hours)')
ylabel('variance in number of cells')
title(['Difference in variances for A= ', num2str(A),', b & d constant'])
%% Continue making functions for other two strong model structures

% NOTE: Must still be having issues with the moments here because we see
% that when A>1 the mean and variance don't match...

% Problem is that these should be exactly the same
A = 3;
paramsAd = [b,d, A];
%Nstat(<n>, <n2>, <Var>, <n3>, <n4>, <V4>)
[ NsampAd,NstatAd, CstatAd, musimdatavecAd, varsimdatavecAd, timesimdatavecAd] = run_strAlldmodel(paramsAd, tsamp, Ninit, Nsim);

stats = {'mean', '2nd moment', 'variance', '3rd moment', '4th moment', '4th order variance'};
% Plot the moments for each initial cell number
for i = 1:3
    figure;
    hold off
for j = 1:6
    subplot(2,3,j)
    plot(tsamp, NstatAd(:,j,i), '*', 'LineWidth', 1)
    hold on
    plot(tsamp, CstatAd(:,j,i), '.', 'LineWidth',2)
    xlabel ('time (hours)')
    ylabel([stats{j}])
    legend([stats{j},' simulated data'], ['expected ', stats{j}])
    title([stats{j},' N_{0}= ', num2str(Ninit(i))])
    legend boxoff
end
end
%% Check strong Allee model on  birth and death
[ NsampAbd,NstatAbd, CstatAbd, musimdatavecAbd, varsimdatavecAbd, timesimdatavecAbd] = run_strAllbdmodel(paramsAd, tsamp, Ninit, Nsim);

stats = {'mean', '2nd moment', 'variance', '3rd moment', '4th moment', '4th order variance'};
% Plot the moments for each initial cell number
for i = 1:3
    figure;
    hold off
for j = 1:6
    subplot(2,3,j)
    plot(tsamp, NstatAbd(:,j,i), '*', 'LineWidth', 1)
    hold on
    plot(tsamp, CstatAbd(:,j,i), '.', 'LineWidth',2)
    xlabel ('time (hours)')
    ylabel([stats{j}])
    legend([stats{j},' simulated data'], ['expected ', stats{j}])
    title([stats{j},' N_{0}= ', num2str(Ninit(i))])
    legend boxoff
end
end




%%
figure;
plot(timesimdatavecAd, musimdatavecAd,'r*', 'LineWidth',2)
hold on
for i = 1:length(Ninit)
    plot(tsamp, CstatAd(:,1,i), 'k.')
    text(tsamp(end-10), CstatAd(end-10, 1, i), ['N_{0}=', num2str(Ninit(i))])
end
xlabel('time (hours)')
ylabel('mean number of cells')
title('Expected and simulated mean for different N_{0}')

figure;
plot(timesimdatavecAd, varsimdatavecAd,'g*', 'LineWidth',2)
hold on
for i = 1:length(Ninit)
    plot(tsamp, CstatAd(:,3,i), 'k.')
    text(tsamp(end-10), CstatAd(end-10, 3, i), ['N_{0}=', num2str(Ninit(i))])
end
xlabel('time (hours)')
ylabel('variance in number of cells')
title('Expected and simulated variance for different N_{0}')
%% Compare mean and variance from Allee on birth and Allee on death

figure;
plot(timesimdatavecAb, musimdatavecAb,'r*', 'LineWidth',2)
hold on
plot(timesimdatavecAd, musimdatavecAd,'b*', 'LineWidth',2)
for i = 1:length(Ninit)
    plot(tsamp, CstatAb(:,1,i), 'k.')
    text(tsamp(end-10), CstatAb(end-10, 1, i), ['N_{0}=', num2str(Ninit(i))])
end
hold on
for i = 1:length(Ninit)
    plot(tsamp, CstatAd(:,1,i), 'k.')
    text(tsamp(end-10), CstatAd(end-10, 1, i), ['N_{0}=', num2str(Ninit(i))])
end
legend('b-d-Ab model', 'b-d-Ad model')
xlabel('time (hours)')
ylabel('mean number of cells')
title(['Difference in means for A= ', num2str(A),', A on P(b) or P(d)'])

figure;
plot(timesimdatavecAb, varsimdatavecAb,'g*', 'LineWidth',2)
hold on
plot(timesimdatavecAd, varsimdatavecAd,'y*', 'LineWidth',2)
for i = 1:length(Ninit)
    plot(tsamp, CstatAb(:,3,i), 'k.')
    text(tsamp(end-10), CstatAb(end-10, 3, i), ['N_{0}=', num2str(Ninit(i))])
end
hold on
for i = 1:length(Ninit)
    plot(tsamp, CstatAd(:,3,i), 'k.')
    text(tsamp(end-10), CstatAd(end-10, 3, i), ['N_{0}=', num2str(Ninit(i))])
end
legend('b-d-Ab model', 'b-d-Ad model')
xlabel('time (hours)')
ylabel('variance in number of cells')
title(['Difference in variances for A= ', num2str(A),', A on P(b) or P(d)'])

%% Make functions for the weak/strong model structure
% Problem is that these should be exactly the same
A = -1;
tau = 2;
paramsAw = [b,d, A, tau];
Ninit = [ 10; 15; 20; 30]; 
%Nstat(<n>, <n2>, <Var>, <n3>, <n4>, <V4>)
[ NsampAwbd,NstatAwbd, CstatAwbd, musimdatavecAwbd, varsimdatavecAwbd, timesimdatavecAwbd] = run_wkAllbdmodel(paramsAw, tsamp, Ninit, Nsim);

%% Plot the mean and variance for all three initial cell numbers
figure;
plot(timesimdatavecAwbd, musimdatavecAwbd,'r*', 'LineWidth',2)
hold on
for i = 1:length(Ninit)
    plot(tsamp, CstatAwbd(:,1,i), 'k.')
    text(tsamp(end-10), CstatAwbd(end-10, 1, i), ['N_{0}=', num2str(Ninit(i))])
end
xlabel('time (hours)')
ylabel('mean number of cells')
title('Expected and simulated mean for different N_{0} weak Allee')

figure;
plot(timesimdatavecAwbd, varsimdatavecAwbd,'g*', 'LineWidth',2)
hold on
for i = 1:length(Ninit)
    plot(tsamp, CstatAwbd(:,3,i), 'k.')
    text(tsamp(end-10), CstatAwbd(end-10, 3, i), ['N_{0}=', num2str(Ninit(i))])
end
xlabel('time (hours)')
ylabel('variance in number of cells')
title('Expected and simulated variance for different N_{0} weak Allee')
%% Weak Allee, Allee on birth probability
A = -1;
tau = 2;
paramsAw = [b,d, A, tau];
Ninit = [ 10; 15; 20; 30]; 
[NsampAwb,NstatAwb, CstatAwb, musimdatavecAwb, varsimdatavecAwb, timesimdatavecAwb] = run_wkAllbmodel(paramsAw, tsamp, Ninit, Nsim);

%% Plot the mean and variance for all three initial cell numbers for weak on b
figure;
plot(timesimdatavecAwb, musimdatavecAwb,'r*', 'LineWidth',2)
hold on
for i = 1:length(Ninit)
    plot(tsamp, CstatAwb(:,1,i), 'k.')
    text(tsamp(end-10), CstatAwb(end-10, 1, i), ['N_{0}=', num2str(Ninit(i))])
end
xlabel('time (hours)')
ylabel('mean number of cells')
title('Expected and simulated mean for different N_{0} weak Allee')

figure;
plot(timesimdatavecAwb, varsimdatavecAwb,'g*', 'LineWidth',2)
hold on
for i = 1:length(Ninit)
    plot(tsamp, CstatAwb(:,3,i), 'k.')
    text(tsamp(end-10), CstatAwb(end-10, 3, i), ['N_{0}=', num2str(Ninit(i))])
end
xlabel('time (hours)')
ylabel('variance in number of cells')
title('Expected and simulated variance for different N_{0} weak Allee')
%% Weak Allee, Allee on death probability
A = -1;
tau = 2;
paramsAw = [b,d, A, tau];
Ninit = [ 10; 15; 20; 30]; 
[NsampAwb,NstatAwd, CstatAwd, musimdatavecAwd, varsimdatavecAwd, timesimdatavecAwd] = run_wkAlldmodel(paramsAw, tsamp, Ninit, Nsim);

%% Plot the mean and variance for all three initial cell numbers for weak on b
figure;
plot(timesimdatavecAwd, musimdatavecAwd,'r*', 'LineWidth',2)
hold on
for i = 1:length(Ninit)
    plot(tsamp, CstatAwd(:,1,i), 'k.')
    text(tsamp(end-10), CstatAwd(end-10, 1, i), ['N_{0}=', num2str(Ninit(i))])
end
xlabel('time (hours)')
ylabel('mean number of cells')
title('Expected and simulated mean for different N_{0} weak Allee')

figure;
plot(timesimdatavecAwd, varsimdatavecAwd,'g*', 'LineWidth',2)
hold on
for i = 1:length(Ninit)
    plot(tsamp, CstatAwd(:,3,i), 'k.')
    text(tsamp(end-10), CstatAwd(end-10, 3, i), ['N_{0}=', num2str(Ninit(i))])
end
xlabel('time (hours)')
ylabel('variance in number of cells')
title('Expected and simulated variance for different N_{0} weak Allee')
%% Compare b-d, strong, and weak models mean and variance
Nsim = 1000;
tsamp = [ 0:4:336];
b = 0.0092;
d = 0.001;
A = 2;
tau = 3;
Ninit = [ 5; 10; 15]; 
paramsbd = [b,d];
paramsA= [b,d,A];
paramsAw = [b,d,-A, tau];

[ Nsampbd,Nstatbd, Cstatbd, musimdatavecbd, varsimdatavecbd, timesimdatavecbd] = run_bdmodel(paramsbd, tsamp, Ninit, Nsim);
[ NsampAbd,NstatAbd, CstatAbd, musimdatavecAbd, varsimdatavecAbd, timesimdatavecAbd] = run_strAllbdmodel(paramsA, tsamp, Ninit, Nsim);
[ NsampAwbd,NstatAwbd, CstatAwbd, musimdatavecAwbd, varsimdatavecAwbd, timesimdatavecAwbd] = run_wkAllbdmodel(paramsAw, tsamp, Ninit, Nsim);

%% PLot mean and variance of each model
figure;
% plot(timesimdatavecbd, musimdatavecbd,'r*', 'LineWidth',2)
% hold on
% plot(timesimdatavecAbd, musimdatavecAbd,'b*', 'LineWidth',2)
% plot(timesimdatavecAwbd, musimdatavecAwbd,'g*', 'LineWidth',2)
hold on
plot(tsamp, Cstatbd(:,1,1), 'r-','LineWidth',2)
text(tsamp(end-10), Cstatbd(end-10, 1, 1), ['N_{0}=', num2str(Ninit(1))])
plot(tsamp, CstatAbd(:,1,1), 'b-','LineWidth',2)
text(tsamp(end-10), CstatAbd(end-10, 1, 1), ['N_{0}=', num2str(Ninit(1))])
plot(tsamp, CstatAwbd(:,1,1), 'g-','LineWidth',2)
text(tsamp(end-10), CstatAwbd(end-10, 1, 1), ['N_{0}=', num2str(Ninit(1))])
for i = 2:length(Ninit)
    plot(tsamp, Cstatbd(:,1,i), 'r-','LineWidth',2)
    text(tsamp(end-10), Cstatbd(end-10, 1, i), ['N_{0}=', num2str(Ninit(i))])
end
hold on
for i = 2:length(Ninit)
    plot(tsamp, CstatAbd(:,1,i), 'b-','LineWidth',2)
    text(tsamp(end-10), CstatAbd(end-10, 1, i), ['N_{0}=', num2str(Ninit(i))])
end
hold on
for i = 2:length(Ninit)
    plot(tsamp, CstatAwbd(:,1,i), 'g-','LineWidth',2)
    text(tsamp(end-10), CstatAwbd(end-10, 1, i), ['N_{0}=', num2str(Ninit(i))])
end
xlim([0 tsamp(end)])
xlabel('time (hours)')
ylabel('mean number of cells')
title(['Mean for each model, constant b & d, A=+/-', num2str(A), ', \tau=', num2str(tau)])
legend ('simple birth-death model', 'strong Allee', 'weak Allee')
legend boxoff
%%
figure;
% plot(timesimdatavecbd, varsimdatavecbd,'r*', 'LineWidth',2)
% hold on
% plot(timesimdatavecAbd, varsimdatavecAbd,'b*', 'LineWidth',2)
% plot(timesimdatavecAwbd, varsimdatavecAwbd,'g*', 'LineWidth',2)
hold on
plot(tsamp, Cstatbd(:,3,1), 'r-','LineWidth',2)
text(tsamp(end-10), Cstatbd(end-10, 3, 1), ['N_{0}=', num2str(Ninit(1))])
plot(tsamp, CstatAbd(:,3,1), 'b-','LineWidth',2)
text(tsamp(end-10), CstatAbd(end-10, 3, 1), ['N_{0}=', num2str(Ninit(1))])
plot(tsamp, CstatAwbd(:,3,1), 'g-','LineWidth',2)
text(tsamp(end-10), CstatAwbd(end-10, 3, 1), ['N_{0}=', num2str(Ninit(1))])
for i = 2:length(Ninit)
    plot(tsamp, Cstatbd(:,3,i), 'r-', 'LineWidth',2)
    text(tsamp(end-10), Cstatbd(end-10, 3, i), ['N_{0}=', num2str(Ninit(i))])
end
hold on
for i = 2:length(Ninit)
    plot(tsamp, CstatAbd(:,3,i), 'b-','LineWidth',2)
    text(tsamp(end-10), CstatAbd(end-10, 3, i), ['N_{0}=', num2str(Ninit(i))])
end
hold on
for i = 2:length(Ninit)
    plot(tsamp, CstatAwbd(:,3,i), 'g-','LineWidth',2)
    text(tsamp(end-10), CstatAwbd(end-10, 3, i), ['N_{0}=', num2str(Ninit(i))])
end
xlim([0 tsamp(end)])
xlabel('time (hours)')
ylabel('Variance in number of cells')
title(['Variance for each model, constant b & d, A=+/-', num2str(A), ', \tau=', num2str(tau)])
legend ('simple birth-death model', 'strong Allee', 'weak Allee')
legend boxoff
%% Plot mean and variance for each model and each stochastic structure
% Already ran forward birth and death
% Now do all birth
[ NsampAb,NstatAb, CstatAb, musimdatavecAb, varsimdatavecAb, timesimdatavecAb] = run_strAllbmodel(paramsA, tsamp, Ninit, Nsim);
[ NsampAwb,NstatAwb, CstatAwb, musimdatavecAwb, varsimdatavecAwb, timesimdatavecAwb] = run_wkAllbmodel(paramsAw, tsamp, Ninit, Nsim);
% Now do all death
[ NsampAd,NstatAd, CstatAd, musimdatavecAd, varsimdatavecAd, timesimdatavecAd] = run_strAlldmodel(paramsA, tsamp, Ninit, Nsim);
[ NsampAwd,NstatAwd, CstatAwd, musimdatavecAwd, varsimdatavecAwd, timesimdatavecAwd] = run_wkAlldmodel(paramsAw, tsamp, Ninit, Nsim);

%% Strong Allee model mean and variance for each stochastic structure
figure;
hold on
plot(tsamp, CstatAb(:,1,1), 'c*','LineWidth',1)
text(tsamp(end-10), CstatAb(end-10, 1, 1), ['N_{0}=', num2str(Ninit(1))])
plot(tsamp, CstatAd(:,1,1), 'b--','LineWidth',4)
text(tsamp(end-10), CstatAd(end-10, 1, 1), ['N_{0}=', num2str(Ninit(1))])
plot(tsamp, CstatAbd(:,1,1), 'k.','LineWidth',8)
text(tsamp(end-10), CstatAbd(end-10,1, 1), ['N_{0}=', num2str(Ninit(1))])
for i = 2:length(Ninit)
    plot(tsamp, CstatAb(:,1,i), 'c*','LineWidth',1)
    text(tsamp(end-10), CstatAb(end-10, 1, i), ['N_{0}=', num2str(Ninit(i))])
end
hold on
for i = 2:length(Ninit)
    plot(tsamp, CstatAd(:,1,i), 'b--','LineWidth',4)
    text(tsamp(end-10), CstatAd(end-10, 1, i), ['N_{0}=', num2str(Ninit(i))])
end
hold on
for i = 2:length(Ninit)
    plot(tsamp, CstatAbd(:,1,i), 'k.','LineWidth',8)
    text(tsamp(end-10), CstatAbd(end-10, 1, i), ['N_{0}=', num2str(Ninit(i))])
end
xlim([0 tsamp(end)])
xlabel('time (hours)')
ylabel('mean number of cells')
title('Mean for each stochastic structure strong Allee model, constant b, d,& A')
legend ('Allee on P(birth)', 'Alee on P(death)', 'Allee equally on both')
legend boxoff

figure;
hold on
plot(tsamp, CstatAb(:,3,1), 'c*','LineWidth',1)
text(tsamp(end-10), CstatAb(end-10, 3, 1), ['N_{0}=', num2str(Ninit(1))])
plot(tsamp, CstatAd(:,3,1), 'b--','LineWidth',4)
text(tsamp(end-10), CstatAd(end-10, 1, 1), ['N_{0}=', num2str(Ninit(1))])
plot(tsamp, CstatAbd(:,3,1), 'k.','LineWidth',8)
text(tsamp(end-10), CstatAbd(end-10,3, 1), ['N_{0}=', num2str(Ninit(1))])
for i = 2:length(Ninit)
    plot(tsamp, CstatAb(:,3,i), 'c*','LineWidth',1)
    text(tsamp(end-10), CstatAb(end-10, 3, i), ['N_{0}=', num2str(Ninit(i))])
end
hold on
for i = 2:length(Ninit)
    plot(tsamp, CstatAd(:,3,i), 'b--','LineWidth',4)
    text(tsamp(end-10), CstatAd(end-10, 3, i), ['N_{0}=', num2str(Ninit(i))])
end
hold on
for i = 2:length(Ninit)
    plot(tsamp, CstatAbd(:,3,i), 'k.','LineWidth',8)
    text(tsamp(end-10), CstatAbd(end-10, 3, i), ['N_{0}=', num2str(Ninit(i))])
end
xlim([0 tsamp(end)])
xlabel('time (hours)')
ylabel('Variance in number of cells')
title('Variance for each stochastic structure strong Allee model, constant b, d,& A')
legend ('Allee on P(birth)', 'Alee on P(death)', 'Allee equally on both')
legend boxoff
%% Weak/Strong Allee model mean and variance for each stochastic structure
figure;
hold on
plot(tsamp, CstatAwb(:,1,1), 'y*','LineWidth',1)
text(tsamp(end-10), CstatAwb(end-10, 1, 1), ['N_{0}=', num2str(Ninit(1))])
plot(tsamp, CstatAwd(:,1,1), 'g--','LineWidth',4)
text(tsamp(end-10), CstatAwd(end-10, 1, 1), ['N_{0}=', num2str(Ninit(1))])
plot(tsamp, CstatAwbd(:,1,1), 'k.','LineWidth',8)
text(tsamp(end-10), CstatAwbd(end-10,1, 1), ['N_{0}=', num2str(Ninit(1))])
for i = 2:length(Ninit)
    plot(tsamp, CstatAwb(:,1,i), 'y*','LineWidth',1)
    text(tsamp(end-10), CstatAwb(end-10, 1, i), ['N_{0}=', num2str(Ninit(i))])
end
hold on
for i = 2:length(Ninit)
    plot(tsamp, CstatAwd(:,1,i), 'g--','LineWidth',4)
    text(tsamp(end-10), CstatAwd(end-10, 1, i), ['N_{0}=', num2str(Ninit(i))])
end
hold on
for i = 2:length(Ninit)
    plot(tsamp, CstatAwbd(:,1,i), 'k.','LineWidth',8)
    text(tsamp(end-10), CstatAwbd(end-10, 1, i), ['N_{0}=', num2str(Ninit(i))])
end
xlim([0 tsamp(end)])
xlabel('time (hours)')
ylabel('mean number of cells')
title('Mean for each stochastic structure weak Allee model, constant b, d, A &\tau')
legend ('Allee on P(birth)', 'Alee on P(death)', 'Allee equally on both')
legend boxoff

figure;
hold on
plot(tsamp, CstatAwb(:,3,1), 'y*','LineWidth',1)
text(tsamp(end-10), CstatAwb(end-10, 3, 1), ['N_{0}=', num2str(Ninit(1))])
plot(tsamp, CstatAwd(:,3,1), 'g--','LineWidth',4)
text(tsamp(end-10), CstatAwd(end-10, 1, 1), ['N_{0}=', num2str(Ninit(1))])
plot(tsamp, CstatAwbd(:,3,1), 'k.','LineWidth',8)
text(tsamp(end-10), CstatAwbd(end-10,3, 1), ['N_{0}=', num2str(Ninit(1))])
for i = 2:length(Ninit)
    plot(tsamp, CstatAwb(:,3,i), 'y*','LineWidth',1)
    text(tsamp(end-10), CstatAwb(end-10, 3, i), ['N_{0}=', num2str(Ninit(i))])
end
hold on
for i = 2:length(Ninit)
    plot(tsamp, CstatAwd(:,3,i), 'g--','LineWidth',4)
    text(tsamp(end-10), CstatAwd(end-10, 3, i), ['N_{0}=', num2str(Ninit(i))])
end
hold on
for i = 2:length(Ninit)
    plot(tsamp, CstatAwbd(:,3,i), 'k.','LineWidth',8)
    text(tsamp(end-10), CstatAwbd(end-10, 3, i), ['N_{0}=', num2str(Ninit(i))])
end
xlim([0 tsamp(end)])
xlabel('time (hours)')
ylabel('Variance in number of cells')
title('Variance for each stochastic structure weak Allee model, constant b, d, A & \tau')
legend ('Allee on P(birth)', 'Alee on P(death)', 'Allee equally on both')
legend boxoff