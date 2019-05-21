% This script find the distribution of cells actually seeded for an attempt
% at sorting N0=8 cells. This distribution will be scaled appropriately and
% used to draw from in the initial cell number used in the Gillespie
% Algorithm stochastic modeling simulations
close all; clear all; clc
[N, T] =xlsread('../data/N8variability.xls');

Ntot= N(:,1);
Ngreen = N(:,2);
Nred = N(:,3);

mu = mean(Ntot)
sigma = std(Ntot)

figure (2)
hold off
histogram(Ntot, 8)
hold on
text(3, 4, ['\mu= ', num2str(round(mu,2)), ', \sigma= ', num2str(round(sigma,2))])
xlabel ('Number of Cells Observed in Well')
xlim([0 8])
ylim([ 0 25])
ylabel('Frequency')
title( 'Distribution of Actual Cell Seeding when N_{0}=8')
%%
figure(1)
subplot(1,3,1)
hold off
histogram(Ntot, 8)
hold on
text(3, 4, ['\mu= ', num2str(round(mu,2)), ', \sigma= ', num2str(round(sigma,2))])
xlabel ('Number of Cells Observed in Well')
xlim([0 10])
ylabel('Frequency')
title( 'Distribution of Cell Seeding when N_{seed}=8')

distrib = normrnd(mu, sigma,[1000,1])

int_distrib = abs(round(distrib,0));


subplot(1,3,2)
hold off
histogram(abs(distrib))
hold on
text(3, 50, ['\mu= ', num2str(round(mu,2)), ', \sigma= ', num2str(round(sigma,2))])
xlabel ('Number of Cells Simulated')
xlim([0 10])
ylabel('Frequency')
title( 'Simulated N_{seed}=8')

subplot(1,3,3)
hold off
histogram(int_distrib, 8)
hold on
text(3, 50, ['\mu= ', num2str(round(mu,2)), ', \sigma= ', num2str(round(sigma,2))])
xlabel ('Number of Cells Simulated')
xlim([0 10])
ylabel('Frequency')
title( 'Simulated Rounded  N_{seed}=8')
distrib8 = [mu, sigma];

save('../out/distrib8.mat', 'distrib8')